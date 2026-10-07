#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {downsample} from './downsample_pierre'
include {downsample as downsample2} from './downsample_pierre'
params.samples_list = "tagseq01_15978_try3.tsv"
params.outdir = "/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/tagseq_01_all_reads_new_barcodes"
params.genome_index = "/groups/berger/user/pierre.bourguet/genomics/indexes/bowtie2_atRTD3_29122021_and_TAIR10_TEs"
params.downsample = -1

process write_info {
    executor 'slurm'
    cpus  1
    memory '1 GB'
    time '1h'
    
    publishDir "${params.outdir}/${lib}/${sample_name}", mode: 'copy' 
    input:
        tuple file(filename), val(lib),  val(sample_name) 
    output:
        file("info.txt")
        file("raw_reads.length")

    script:
    """
        echo "${sample_name}" > info.txt
        echo "$filename" >> info.txt
        zcat $filename | wc -l | awk '{printf "%d\\n", \$1 / 4}' > raw_reads.length
    """
}

process trim_tagseq {
    executor 'slurm'
    cpus  2
    memory '2 GB'
    time '1h'
    module 'build-env/2020:fastp/0.20.1-gcc-8.2.0-2.31.1'   
    publishDir "${params.outdir}/${lib}/${sample_name}", mode: 'copy', pattern: '*.{log,length,json,html}'
    input:
        tuple file(filename), val(lib),  val(sample_name)
    output:
        val(lib)
        val(sample_name)
        file("fastp_output.fastq")
        file("fastp_output.html")
        file("fastp_output.json")
        file("fastp_output.log")
        file("fastp_output.length")
        file("fastp_fixed.fastq")
    script:
    """
        
        fastp -i $filename -o fastp_fixed.fastq --max_len1 50  2> fastp_fixed.log
        fastp -i fastp_fixed.fastq -o fastp_output.fastq --trim_poly_x -h fastp_output.html -j fastp_output.json --max_len1 70  2> fastp_output.log
        cat fastp_output.fastq | wc -l | awk '{printf "%d\\n", \$1 / 4}' > fastp_output.length
    """
}

process clean_umi_duplicates {
    executor 'slurm'
    cpus  2
    time '1h'
    memory { 10.GB + 20.GB * task.attempt }
    errorStrategy 'retry' 
    module 'build-env/2020:bbmap/38.26-foss-2018b'
    publishDir "${params.outdir}/${lib}/${sample_name}", mode: 'copy', pattern: '*.{log,stats,length}'
    input:
        val(lib)
        val(sample_name)
        file(filename)
    output:
        tuple val(lib), val(sample_name), file("umi_dedup_remove_umi.fastq")
        file("triple_umi_in_R1.fastq")
        file("umi_dedup.fastq")
        file("umi_dedup.log")
        file("umi_dedup.stats")
        file("umi_dedup_remove_umi.length")
    
    // This is done as preparation for the cleaning of the UMI using clumpify.sh
    // As clumpify allow up to two missmatch any difference in the UMI will indicate different instances of a read
    script:
    """
        cat $filename | awk 'BEGIN {UMI = ""; q="IIIIIIII"} {r_n = (NR%4); if(r_n == 1) {UMI = substr(\$1,length(\$1)-7,length(\$1)); print;}; if(r_n==2) {printf "%s%s%s%s\\n", UMI, UMI, UMI, \$1}; if(r_n==3) {print}; if(r_n==0) {printf "%s%s%s%s\\n", q,q,q,\$1}}' > triple_umi_in_R1.fastq
        
        clumpify.sh in=triple_umi_in_R1.fastq out=umi_dedup.fastq dedupe addcount 2> umi_dedup.log

        cat umi_dedup.fastq | awk '(NR % 4) == 1' | grep copies | awk '{print substr(\$NF,8,100)}' | sort -nk1 | uniq -c | sort -nk2 > umi_dedup.stats
        cat umi_dedup.fastq | wc -l | awk '{printf "%d\\n", \$1 / 4}' > umi_dedup.length

        cat umi_dedup.fastq | awk '{if((NR%2)==1) {print} else {printf "%s\\n", substr(\$1,25,length(\$1))}}' > umi_dedup_remove_umi.fastq
        cat umi_dedup_remove_umi.fastq | wc -l | awk '{printf "%d\\n", \$1 / 4}' > umi_dedup_remove_umi.length

    """ 
}

process map_reads {
    executor 'slurm'
    cpus  1
    memory { 10.GB + 20.GB * task.attempt }
    time '6h'
    module 'build-env/2020:bowtie2/2.3.5.1-foss-2018b'

    publishDir "${params.outdir}/${lib}/${sample_name}", mode: 'copy', pattern: '*.log' 
    input:
        tuple val(lib), val(sample_name), file(filename)
        val(bowtie_index)
    output:
        val(lib)
        val(sample_name)
        file("mapped_reads.map")
        file("mapped_reads.log")
    
    script:
    """
        bowtie2 -x ${bowtie_index} -U $filename -p 1 -a > mapped_reads.map 2> mapped_reads.log 
    """ 
}

process quantify_exp {
    module 'build-env/f2021:salmon/1.5.2-gompi-2020b'
    executor 'slurm'
    cpus  1
    memory '5 GB'
    time '1h'

    publishDir "${params.outdir}/${lib}/${sample_name}", mode: 'copy', pattern: 'quant*' 
    input:
        tuple val(lib), val(sample_name), file(filename), val(sample_size), val(rep)
        val(index)
    output:
        file("quant*")

    """
        salmon quant -l SF -i $index -r $filename -p 1 -o quant_${sample_size}_R${rep} --noLengthCorrection
    """
}

def indexed( items ) {
    return Channel.from( items.withIndex() ).map { item, idx -> tuple( idx, item ) }
}

workflow {
    samples = Channel.fromPath(params.samples_list)
        .splitText()
        .map { it.replaceFirst(/\n/,'') }
        .splitCsv(sep: '\t')
        .map {it -> [file(it[0]) , it[1], it[2]] }
    
    write_info(samples)
    trimmed_samples = trim_tagseq(samples)

    umi_cleaned = clean_umi_duplicates(trimmed_samples[0], trimmed_samples[1], trimmed_samples[2])
    ds_reads = downsample(umi_cleaned[0].combine(Channel.from(500000000)).combine(Channel.from(1)))
    
    map_reads(ds_reads.map {[it[0], it[1], it[2]]}, 
                        "/groups/berger/user/pierre.bourguet/genomics/indexes/bowtie2_atRTD3_29122021_and_TAIR10_TEs/atRTD3_29122021_and_TAIR10_TEs")

        quantify_exp(umi_cleaned[0],
            "/groups/berger/user/pierre.bourguet/genomics/scripts/salmon/salmon_index_atRTD3_29122021_and_TAIR10_TEs")
}
