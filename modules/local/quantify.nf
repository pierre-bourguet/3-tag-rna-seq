// Mapping and quantification: STAR on the genome, salmon on the STAR transcriptome BAM, salmon pseudo-alignment,
// bigwigs and sample-to-sample correlation

// The TranscriptomeSAM options are required for salmon alignment-mode counting downstream
process STAR_ALIGN {
    tag "${meta.id}"
    cpus 8
    memory { 12.GB + 4.GB * task.attempt }
    time '2h'
    publishDir "${params.outdir}/01_QC/STAR_logs_and_QC", mode: 'copy', pattern: '*_Log*.out'
    publishDir "${params.outdir}/03_STAR_bam", mode: 'copy', pattern: '*_Aligned.sortedByCoord.out.bam', enabled: params.save_bam

    input:
    tuple val(meta), path(fastq)
    path index

    output:
    tuple val(meta), path("${meta.id}_Signal*.out.bg"),                  emit: bedgraph
    tuple val(meta), path("${meta.id}_Log.final.out"),                   emit: log_final
    path "${meta.id}_Log.out",                                           emit: log
    tuple val(meta), path("${meta.id}_Aligned.toTranscriptome.out.bam"), emit: tx_bam
    path "${meta.id}_Aligned.sortedByCoord.out.bam",                     emit: bam

    script:
    """
    STAR --genomeDir ${index} --outFileNamePrefix ${meta.id}_ \\
        --readFilesIn ${fastq} --runThreadN ${task.cpus} \\
        --outSAMtype BAM SortedByCoordinate --outWigType bedGraph --outWigNorm RPM --outWigStrand Stranded \\
        --outFilterMismatchNoverLmax 0.04 --outFilterMultimapNmax 50 \\
        --quantMode TranscriptomeSAM --quantTranscriptomeBan Singleend --limitBAMsortRAM 8287822014
    """
}

// Sense (-l SF) and antisense (-l SR) counts, all reads and unique mappers (MAPQ 255) only
process SALMON_QUANT_ALIGNED {
    tag "${meta.id}"
    cpus 8
    memory { 24.GB + 24.GB * task.attempt }
    time '2h'
    errorStrategy 'retry'
    publishDir "${params.outdir}/02_counts/samples/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(bam)
    path transcripts

    output:
    tuple val(meta), path("*STAR_mapping_salmon_quant_${meta.id}*", type: 'dir'), emit: quant

    script:
    def id = meta.id
    """
    salmon quant -p ${task.cpus} -l SF -t ${transcripts} -a ${bam} -o STAR_mapping_salmon_quant_${id} --noLengthCorrection
    salmon quant -p ${task.cpus} -l SR -t ${transcripts} -a ${bam} -o STAR_mapping_salmon_quant_${id}_AS --noLengthCorrection

    samtools view --bam -q 255 ${bam} > unique_${bam}
    salmon quant -p ${task.cpus} -l SF -t ${transcripts} -a unique_${bam} -o unique_STAR_mapping_salmon_quant_${id} --noLengthCorrection
    salmon quant -p ${task.cpus} -l SR -t ${transcripts} -a unique_${bam} -o unique_STAR_mapping_salmon_quant_${id}_AS --noLengthCorrection
    rm unique_${bam}
    """
}

// Selective alignment against transcripts with the genome as decoy, independent of STAR
process SALMON_QUANT_PSEUDO {
    tag "${meta.id}"
    cpus 8
    memory { 4.GB + 8.GB * task.attempt }
    time '2h'
    errorStrategy 'retry'
    publishDir "${params.outdir}/02_counts/samples/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(fastq)
    path index

    output:
    tuple val(meta), path("quant_${meta.id}*", type: 'dir'), emit: quant

    script:
    """
    salmon quant -p ${task.cpus} -l SF -i ${index} -r ${fastq} -o quant_${meta.id} --noLengthCorrection
    salmon quant -p ${task.cpus} -l SR -i ${index} -r ${fastq} -o quant_${meta.id}_AS --noLengthCorrection
    """
}

// STAR RPM bedgraphs -> stranded and unstranded bigwigs, unique and unique+multi mappers
process BIGWIG {
    tag "${meta.id}"
    cpus 2
    memory '8 GB'
    time '1h'
    publishDir "${params.outdir}/03_STAR_bigwigs", mode: 'copy'

    input:
    tuple val(meta), path(bedgraphs)
    path chrom_sizes

    output:
    path "unique_multi_unstranded/${meta.id}_unique_multi.unstranded.bw", emit: unique_multi_unstranded
    path "unique_unstranded/${meta.id}_unique.unstranded.bw",             emit: unique_unstranded
    path "unique_multi_stranded/${meta.id}_unique_multi.str?.bw",         emit: unique_multi_stranded
    path "unique_stranded/${meta.id}_unique.str?.bw",                     emit: unique_stranded

    script:
    def id = meta.id
    """
    export PATH="${params.ucsc_tools}:\$PATH"
    for str in str1 str2; do
        bedGraphToBigWig ${id}_Signal.Unique.\${str}.out.bg ${chrom_sizes} ${id}_unique.\${str}.bw
        bedGraphToBigWig ${id}_Signal.UniqueMultiple.\${str}.out.bg ${chrom_sizes} ${id}_unique_multi.\${str}.bw
    done
    for type in unique unique_multi; do
        bigWigMerge ${id}_\${type}.str1.bw ${id}_\${type}.str2.bw ${id}_\${type}.unstranded.bg
        bedGraphToBigWig ${id}_\${type}.unstranded.bg ${chrom_sizes} ${id}_\${type}.unstranded.bw
        mkdir -p \${type}_stranded \${type}_unstranded
        mv ${id}_\${type}.str1.bw ${id}_\${type}.str2.bw \${type}_stranded/
        mv ${id}_\${type}.unstranded.bw \${type}_unstranded/
    done
    """
}

process DEEPTOOLS_CORRELATION {
    cpus { 6 + 12 * task.attempt }
    memory { 4.GB + 16.GB * task.attempt }
    time '2h'
    publishDir "${params.outdir}/01_QC/deeptools_plots", mode: 'copy', pattern: '*.{pdf,tab}'

    input:
    path bigwigs

    output:
    path "*.{pdf,tab}", emit: plots

    script:
    """
    multiBigwigSummary bins --bwfiles ${bigwigs} -p ${task.cpus} --outFileName summary.npz
    for method in spearman pearson; do
        plotCorrelation --corData summary.npz --corMethod \${method} --colorMap RdYlBu --skipZeros --removeOutliers \\
            -p heatmap -o Corr_\${method}.pdf --outFileCorMatrix Corr_\${method}.tab
    done
    plotPCA -in summary.npz -o PCA.pdf --outFileNameData PCA.tab
    """
}
