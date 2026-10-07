// Read preprocessing: adapter/polyA/quality trimming, UMI collapsing, downsampling, FastQC

// Trims 3' adapters, poly-X and low-quality tails; caps read length; drops short reads
process FASTP {
    tag "${meta.id}"
    cpus 8
    memory { 8.GB + 16.GB * task.attempt }
    time '1h'
    publishDir "${params.outdir}/01_QC/fastp_trimming", mode: 'copy', pattern: '*_fastp.{json,html,log}'

    input:
    tuple val(meta), path(fastq)
    path adapters

    output:
    tuple val(meta), path("${meta.id}_fastp_output.fastq"), emit: reads
    tuple val(meta), path("${meta.id}_fastp.json"),         emit: json
    path "${meta.id}_fastp.{html,log}",                     emit: report

    script:
    """
    fastp --thread ${task.cpus} -i ${fastq} -o ${meta.id}_fastp_output.fastq \\
        --trim_poly_x --adapter_fasta ${adapters} --cut_tail --max_len1 ${params.max_len} --length_required ${params.min_len} \\
        -h ${meta.id}_fastp.html -j ${meta.id}_fastp.json 2> ${meta.id}_fastp.log
    """
}

// Collapses PCR duplicates: reads with identical UMI (8 nt, end of read name) and near-identical sequence
// (clumpify: up to 2 mismatches). The UMI is written 3x in front of the sequence (with max-quality bases) so that
// any UMI mismatch exceeds clumpify's tolerance, then removed again.
process UMI_DEDUP {
    tag "${meta.id}"
    cpus 2
    memory { 10.GB + 20.GB * task.attempt }
    time '2h'
    errorStrategy 'retry'
    publishDir "${params.outdir}/01_QC/umi_dedup", mode: 'copy', pattern: '*_umi_dedup.{log,stats,length}'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}_umi_dedup_remove_umi.fastq"), emit: reads
    tuple val(meta), path("${meta.id}_umi_dedup.length"),           emit: length
    path "${meta.id}_umi_dedup.{log,stats}",                         emit: logs

    script:
    def id = meta.id
    """
    awk 'BEGIN {UMI = ""; q="IIIIIIII"} {r_n = (NR%4); if(r_n == 1) {UMI = substr(\$1,length(\$1)-7,length(\$1)); print;}; if(r_n==2) {printf "%s%s%s%s\\n", UMI, UMI, UMI, \$1}; if(r_n==3) {print}; if(r_n==0) {printf "%s%s%s%s\\n", q,q,q,\$1}}' ${fastq} > ${id}_triple_umi_in_R1.fastq

    clumpify.sh in=${id}_triple_umi_in_R1.fastq out=${id}_umi_dedup.fastq dedupe addcount tossjunk 2> ${id}_umi_dedup.log
    rm ${id}_triple_umi_in_R1.fastq

    # distribution of copies per collapsed read
    awk '(NR % 4) == 1' ${id}_umi_dedup.fastq | { grep copies || true; } | awk '{print substr(\$NF,8,100)}' | sort -nk1 | uniq -c | sort -nk2 > ${id}_umi_dedup.stats

    awk 'END {printf "%d\\n", NR / 4}' ${id}_umi_dedup.fastq > ${id}_umi_dedup.length

    awk '{if((NR%2)==1) {print} else {printf "%s\\n", substr(\$1,25,length(\$1))}}' ${id}_umi_dedup.fastq > ${id}_umi_dedup_remove_umi.fastq
    rm ${id}_umi_dedup.fastq
    """
}

// Random subsample to --max_n_read reads; samples at or below the limit pass through unchanged
process DOWNSAMPLE {
    tag "${meta.id}"
    cpus 1
    memory { 4.GB * task.attempt }
    time '1h'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}.fastq"), emit: reads

    script:
    """
    n=\$(awk 'END {printf "%d", NR / 4}' ${fastq})
    if [ "\$n" -le ${params.max_n_read} ]; then
        ln -s "\$(readlink -f ${fastq})" ${meta.id}.fastq
    else
        seqtk sample -2 -s ${params.seed} ${fastq} ${params.max_n_read} > ${meta.id}.fastq
    fi
    """
}

process FASTQC {
    tag "${meta.id}"
    cpus 1
    memory '2 GB'
    time '1h'
    publishDir "${params.outdir}/01_QC/fastqc", mode: 'copy'

    input:
    tuple val(meta), path(fastq)

    output:
    path "*_fastqc.{zip,html}", emit: reports

    script:
    """
    fastqc ${fastq}
    """
}
