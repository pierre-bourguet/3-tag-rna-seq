// Reference build: AtRTD3 transcriptome with TAIR10 ATTEs (TE genes removed) or with its TE genes, optional transgene,
// STAR and salmon indices, annotation tables for DESeq2. Output locations are set in PREPARE_REFERENCE.

// genome.fa (+ transgene), chromosome sizes, salmon decoy names
process REF_GENOME {
    cpus 1
    memory '4 GB'
    time '1h'
    publishDir "${params.reference_out}", mode: 'copy', pattern: '{genome.fa,genome.fa.fai,chrom_sizes.txt,decoys.txt}'

    input:
    path genome, stageAs: 'input/genome.fa'
    path transgene_genome, stageAs: 'input/transgene_genome.fa'

    output:
    path "genome.fa",               emit: genome
    path "genome_no_transgene.fa",  emit: genome_base
    path "genome.fa.fai",           emit: fai
    path "chrom_sizes.txt",         emit: chrom_sizes
    path "decoys.txt",              emit: decoys

    script:
    def tg = transgene_genome.name != 'NO_FILE' && params.transgene
    """
    cp ${genome} genome_no_transgene.fa
    cat genome_no_transgene.fa ${tg ? transgene_genome : ''} > genome.fa
    samtools faidx genome.fa
    cut -f 1,2 genome.fa.fai > chrom_sizes.txt
    grep ">" genome.fa | sed 's/ .*// ; s/^>//' > decoys.txt
    """
}

// AGAT input: AtRTD3 GTF (TE genes removed in ATTE mode) converted with gffread -T, plus the TAIR10 ATTE records
process REF_ANNOTATION {
    cpus 1
    memory '8 GB'
    time '1h'
    publishDir "${params.reference_out}/build", mode: 'copy', pattern: '{TEG_IDs.txt,TAIR10_GFF3_ATTEs.gff}'

    input:
    path gtf
    path tair10_gff

    output:
    path "agat_input.gff",        emit: gff
    path "TAIR10_GFF3_ATTEs.gff", emit: attes
    path "TEG_IDs.txt",           emit: teg_ids

    script:
    """
    awk '\$3=="transposable_element_gene" {print \$9}' ${tair10_gff} | sed 's/ID=// ; s/\\;.*//' > TEG_IDs.txt
    awk '\$3 == "transposable_element"' ${tair10_gff} > TAIR10_GFF3_ATTEs.gff

    if [ "${params.te_annotation}" = "ATTE" ]; then
        grep -v -f TEG_IDs.txt ${gtf} > transcripts.gtf
        gffread transcripts.gtf -T -o transcripts.gff
        cat transcripts.gff TAIR10_GFF3_ATTEs.gff > agat_input.gff
    else
        gffread ${gtf} -T -o agat_input.gff
    fi
    """
}

// Makes the GFF consistent (missing gene/transcript levels, IDs) so that STAR, gffread and salmon agree on transcripts
process AGAT_FIX {
    cpus 4
    memory { 16.GB * task.attempt }
    time '8h'
    publishDir "${params.reference_out}/build", mode: 'copy', pattern: '*.agat.log'

    input:
    path gff

    output:
    path "fixed.gff",  emit: gff
    path "*.agat.log", emit: log

    script:
    """
    agat_convert_sp_gxf2gxf.pl --gff ${gff} -o fixed.gff
    """
}

// annotation.gff (with transgene) and annotation_no_transgene.gff (for transcript extraction).
// ATTE mode: AGAT drops the ATTE records, so they are added back as gene/transcript/exon; sorted with en_US
// collation, as in the original reference.
process REF_FINAL_GFF {
    cpus 2
    memory '8 GB'
    time '1h'
    publishDir "${params.reference_out}", mode: 'copy', pattern: 'annotation.gff'

    input:
    path fixed
    path attes
    path transgene_gff, stageAs: 'input/transgene.gff'

    output:
    path "annotation.gff",               emit: gff
    path "annotation_no_transgene.gff",  emit: gff_base

    script:
    def tg = transgene_gff.name != 'NO_FILE' && params.transgene
    """
    if [ "${params.te_annotation}" = "ATTE" ]; then
        convert_atte_gff.awk ${attes} > converted_TAIR10_GFF3_ATTEs.gff
        cat ${fixed} converted_TAIR10_GFF3_ATTEs.gff | LC_ALL=en_US.UTF-8 sort -S 4G -T . -k1,1 -k4,4n > annotation_no_transgene.gff
    else
        cp ${fixed} annotation_no_transgene.gff
    fi
    cat annotation_no_transgene.gff ${tg ? transgene_gff : ''} > annotation.gff
    """
}

// Transcript sequences (spliced exons) from the genome, + transgene cDNAs
process AGAT_EXTRACT {
    cpus 2
    memory { 16.GB * task.attempt }
    time '8h'
    publishDir "${params.reference_out}", mode: 'copy', pattern: 'transcripts.fa'

    input:
    path gff
    path genome
    path transgene_cdna, stageAs: 'input/transgene_cdna.fa'

    output:
    path "transcripts.fa", emit: fasta

    script:
    def tg = transgene_cdna.name != 'NO_FILE' && params.transgene
    """
    agat_sp_extract_sequences.pl -g ${gff} -f ${genome} --mrna -o transcripts_no_transgene.fa
    cat transcripts_no_transgene.fa ${tg ? transgene_cdna : ''} > transcripts.fa
    """
}

process STAR_INDEX {
    cpus 16
    memory '48 GB'
    time '4h'
    publishDir "${params.reference_out}", mode: 'copy', pattern: 'annotation.gtf'
    publishDir "${params.index_out}", mode: 'copy', pattern: 'STAR'

    input:
    path genome
    path gff

    output:
    path "STAR",           emit: index
    path "annotation.gtf", emit: gtf

    script:
    """
    gffread ${gff} -T -o annotation.gtf
    STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir STAR \\
        --genomeFastaFiles ${genome} --sjdbGTFfile annotation.gtf --limitGenomeGenerateRAM 45292192010
    """
}

// Decoy-aware index (transcripts + genome) for salmon selective alignment
process SALMON_INDEX {
    cpus 16
    memory { 24.GB + 12.GB * task.attempt }
    time '4h'
    publishDir "${params.index_out}", mode: 'copy', pattern: 'salmon'

    input:
    path transcripts
    path genome
    path decoys

    output:
    path "salmon", emit: index

    script:
    """
    cat ${transcripts} ${genome} > gentrome.fa
    salmon index -t gentrome.fa -d ${decoys} -i salmon -p ${task.cpus}
    rm gentrome.fa
    """
}

// Tables read by deseq2/DESeq2_tagseq.R: protein-coding genes, the TE set of this mode (ATTE or TEG) with
// family/superfamily, TE x PCG overlaps in the same / opposite orientation, gene functional annotation
process DESEQ2_TABLES {
    cpus 1
    memory '4 GB'
    time '1h'
    publishDir "${params.reference_out}/deseq2", mode: 'copy'

    input:
    path tair10_gff
    path te_table
    path teg_table
    path gene_annotations

    output:
    path "*.tsv", emit: tables

    script:
    """
    awk -F'\\t' -v OFS='\\t' 'BEGIN {print "Chr","Start","End","Geneid","Type","Strand"}
        \$9 ~ /Note=protein_coding_gene/ || \$9 ~ /Note=pseudogene/ {
            split(\$9, a, ";"); id = a[1]; sub(/ID=/, "", id); type = a[2]; sub(/Note=/, "", type)
            if (type == "protein_coding_gene") print \$1, \$4, \$5, id, type, \$7 }' ${tair10_gff} > TAIR10_GFF_PCG.tsv

    if [ "${params.te_annotation}" = "ATTE" ]; then
        cp ${te_table} TE_table.tsv
        awk '\$3 == "transposable_element"' ${tair10_gff} > TE.gff
    else
        awk -F'\\t' -v OFS='\\t' 'NR == 1 {print "Transposon_Name","orientation_is_5prime","Transposon_min_Start","Transposon_max_End","Transposon_Family","Transposon_Super_Family"; next}
            {print \$4, (\$6 == "+" ? "true" : "false"), \$2 + 1, \$3, \$8, \$9}' ${teg_table} > TE_table.tsv
        awk '\$3 == "transposable_element_gene"' ${tair10_gff} > TE.gff
    fi

    tail -n+2 TAIR10_GFF_PCG.tsv > PCG.bed
    awk 'BEGIN {FS=";"} {print \$1}' TE.gff | sed 's/ID=//' | awk 'BEGIN {OFS="\\t"} {print \$1,\$4,\$5,\$9,\$3,\$7}' > TE.bed
    bedtools intersect -wao -a TE.bed -b PCG.bed | awk '\$13 != 0' > TE_PCG.tsv
    awk -F'\\t' -v OFS='\\t' '\$6 == \$12' TE_PCG.tsv > TE_PCG_intersect_same_orientation.tsv
    awk -F'\\t' -v OFS='\\t' '\$6 != \$12' TE_PCG.tsv > TE_PCG_intersect_opposite_orientation.tsv
    rm TE_PCG.tsv

    cp ${gene_annotations} gene_annotations.tsv
    """
}

// reference_manifest.tsv: absolute paths of everything a run needs, plus provenance
process REF_MANIFEST {
    cpus 1
    memory '1 GB'
    time '30m'
    publishDir "${params.reference_out}", mode: 'copy'

    input:
    path inputs, stageAs: 'inputs/*'
    path transgene_gff, stageAs: 'input/transgene.gff'
    val ready

    output:
    path "reference_manifest.tsv"

    script:
    def r = params.reference_out
    def i = params.index_out
    def tg = transgene_gff.name != 'NO_FILE' && params.transgene
    """
    ids=""
    if [ "${tg}" = "true" ]; then
        ids=\$(awk -F'\\t' '\$3 == "gene" {match(\$9, /ID=[^;]+/); print substr(\$9, RSTART + 3, RLENGTH - 3)}' ${transgene_gff} | paste -sd, -)
    fi
    {
        echo -e "# 3-tag-rna-seq reference, built \$(date +%F) with ${workflow.manifest.version} (\$(git -C ${projectDir} describe --tags --always --dirty 2>/dev/null || echo NA))"
        echo -e "te_annotation\\t${params.te_annotation}"
        echo -e "transgene\\t${params.transgene ?: 'none'}"
        echo -e "unique_count_ids\\t\$ids"
        echo -e "genome_fasta\\t${r}/genome.fa"
        echo -e "chrom_sizes\\t${r}/chrom_sizes.txt"
        echo -e "annotation_gff\\t${r}/annotation.gff"
        echo -e "annotation_gtf\\t${r}/annotation.gtf"
        echo -e "transcripts_fasta\\t${r}/transcripts.fa"
        echo -e "star_index\\t${i}/STAR"
        echo -e "salmon_index\\t${i}/salmon"
        echo -e "deseq2_pcg\\t${r}/deseq2/TAIR10_GFF_PCG.tsv"
        echo -e "deseq2_te\\t${r}/deseq2/TE_table.tsv"
        echo -e "deseq2_te_pcg_same\\t${r}/deseq2/TE_PCG_intersect_same_orientation.tsv"
        echo -e "deseq2_te_pcg_opposite\\t${r}/deseq2/TE_PCG_intersect_opposite_orientation.tsv"
        echo -e "deseq2_annotations\\t${r}/deseq2/gene_annotations.tsv"
        for f in inputs/*; do echo -e "md5:\$(readlink -f \$f)\\t\$(md5sum < \$f | cut -d' ' -f1)"; done
    } > reference_manifest.tsv
    """
}
