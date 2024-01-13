#!/usr/bin/env nextflow

samples = Channel.fromPath(params.sample_list)
               .splitText()
               .map { it.replaceFirst(/\n/,'') }
               .splitCsv(sep:'\t')
               .map { it -> [file(it[0]), it[1]] }

//Pulling the AT bowtie2 index from singularity container
if(params.reference_genome == "tair10")
{
params.index="library://elin.axelsson/index/index_bowtie2_tair10:v2.4.1-release-47"
}

//Basic parameter check
if ( !params.index ) exit 1, "Error: no Bowtie2 index"

process pull_genome_index
{
     executor 'slurm'

     input:
     file params.index

     output:
     file params.reference_genome

     script:
     """
     singularity run ${params.index}
     """
}

process write_info {
     executor 'slurm'
     cpus 1
     memory '1 GB'
     time '1h'

     publishDir "${params.outdir}/${sample_name}", mode: 'copy'

     input:
     tuple file(filename), val(sample_name)

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
	cpus 2
	memory '2 GB'
	time '1h'
	module 'build-env/2020:fastp/0.20.1-gcc-8.2.0-2.31.1'

	publishDir "${params.outdir}/${sample_name}", mode: 'copy', pattern: '*.{log,length,json,html}'

	input:
	tuple file(filename), val(sample_name)

	output:
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
	cpus 2
	time '1h'
	memory { 10.GB + 20.GB * task.attempt }
	errorStrategy 'retry'
	module 'build-env/2020:bbmap/38.26-foss-2018b'

	publishDir "${params.outdir}/${sample_name}", mode: 'copy', pattern: '*.{log,stats,length}'

	input:
	val(sample_name)
	file(filename)

	output:
	tuple val(sample_name), file("umi_dedup_remove_umi.fastq")
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

process fastqc {
     executor 'slurm'
	cpus 8
	memory '64 GB'
	time '4h'
     module 'fastqc/0.11.9-java-11'

     publishDir "${params.outdir}/${sample_name}", mode: 'copy', pattern: '*.{zip,html,txt}'

     input:
     tuple val(sample_name), file(filename)

     output:
     file "*.{zip,html,txt}"

     script:
     """
     fastqc ${filename}
     """
}

process bowtie2_alignment {
	executor 'slurm'
	cpus 8
	memory '64 GB'
	time '4h'
	module 'build-env/f2022:bowtie2/2.5.1-gcc-12.2.0'

	publishDir "${params.outdir}/${sample_name}", mode: 'copy', pattern: '*.log'

	input:
	tuple val(sample_name), file(filename)
	file(genomeIndex)

	output:
	val(sample_name)
        file("mapped_reads.map")
        file("mapped_reads.log")

	script:
	"""
	bowtie2 -x ${genomeIndex}/${genomeIndex} -U $filename -p 8 -a > mapped_reads.map 2> mapped_reads.log
	"""
}

process create_salmon_index
{
     executor 'slurm'
     cpus 4
     memory '10 GB'
     time '1h'
     module 'build-env/f2021:salmon/1.5.2-gompi-2020b'

     publishDir "${params.outdir}/salmon_index/", mode: 'copy'

     input:
     val cdna_url
     val tes_url

     output:
     file("*")
     val("${params.outdir}/salmon_index")

     script:
     """
     wget $cdna_url -O cdna.fa
     wget $tes_url -O TEs.fa

     cat cdna.fa TEs.fa > concatenated_fasta.fa

     salmon index -t concatenated_fasta.fa -i ${params.outdir}/salmon_index
     """
}



process quantify_exp {
	executor 'slurm'
     cpus 4
     memory '10 GB'
     time '1h'
     module 'build-env/f2021:salmon/1.5.2-gompi-2020b'

     publishDir "${params.outdir}/${sample_name}", mode: 'copy', pattern: 'quant*'

     input:
     tuple val(sample_name), file(filename)
     val index

	output:
	file("quant*")

	script:
	"""
	salmon quant -l SF -i $index -r $filename -p 4 -o quant_${sample_name} --noLengthCorrection
	"""
}


workflow {
//Writing the raw read counts in the metadata
     write_info(samples)
//Pulling the bowtie2 AT genome index for bowtie2 alignment
     bw2index = pull_genome_index(params.index)
//Downloading cdna and te files from TAIR10 from these links
     params.cdna_url = "https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_blastsets/TAIR10_cdna_20110103_representative_gene_model_updated"
     params.tes_url = "https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_transposable_elements/TAIR10_TE.fas"
//Creating the salmon index from the downloaded (concatenated) files
     salmon_index = create_salmon_index(params.cdna_url, params.tes_url)
//Trimming poly
     trimmed_samples = trim_tagseq(samples)
	umi_cleaned = clean_umi_duplicates(trimmed_samples[0], trimmed_samples[1])
     fastqc(umi_cleaned[0])
	bowtie2_alignment(umi_cleaned[0], bw2index)


     quantify_exp(umi_cleaned[0], salmon_index[1])
}
