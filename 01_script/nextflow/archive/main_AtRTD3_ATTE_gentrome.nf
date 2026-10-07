#!/usr/bin/env nextflow

//Channel created from input file and mapping of file name and path
samples = Channel.fromPath(params.sample_list)
               .splitText()
               .map { it.replaceFirst(/\n/,'') }
               .splitCsv(sep:'\t')
               .map { it -> [file(it[0]), it[1]] }

//Channel for the R script to plot the read count summary in R
params.rscript = file("${baseDir}/Summarize_read_counts.R")
rscript = Channel.fromPath(params.rscript)

//This should remain 0 for anyone other than me using the pipeline
//levels=1 just adopts a x axis ordering custom difined only for my data
params.levels=0

//Default links of files to be imported
params.genome_url = "/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/02_genome_files/02_input/TAIR10_chr_all.fas"
params.gff_url = "/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/02_genome_files/03_output/atRTD3_TS_21Feb22_transfix_no_TEG_w_ATTE.gff"
params.gentrome_url = "/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/02_genome_files/03_output/salmon_decoys/gentrome_atRTD3_29122021_no_TEG_w_ATTE.fna"
params.decoys_url = "/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/02_genome_files/03_output/salmon_decoys/decoys.txt"


//This process counts the no. of reads in the raw files for each sample 
process write_info {
	executor 'slurm'
	cpus 1
	memory '1 GB'
	time '1h'

	publishDir "${params.outdir}/01_QC/raw_read_counts/", mode: 'copy'

	input:
	tuple file(filename), val(sample_name)

	output:
	file("${sample_name}_raw_read_counts.txt") 

	script:
	"""
	read_count=\$(zcat $filename | wc -l | awk '{printf "%d\\n", \$1 / 4}')
	echo -e "${sample_name}\t\$read_count" >> ${sample_name}_raw_read_counts.txt
	"""
}

//This process trims the read to 50 bp, removes the poly_x and counts the no. of reads after trimming
process trim_tagseq {
	executor 'slurm'
	cpus 8
        memory { 8.GB + 8.GB * task.attempt }
	time '1h'
	module 'build-env/2020:fastp/0.20.1-gcc-8.2.0-2.31.1'

	publishDir "${params.outdir}/01_QC/fastp_trimming/", mode: 'copy', pattern: '*.{log,length,json,html}'

	input:
	tuple file(filename), val(sample_name)

	output:
	val(sample_name)
        file("${sample_name}_fastp_output.fastq")
        file("${sample_name}_fastp_output.html")
        file("${sample_name}_fastp_output.json")
        file("${sample_name}_fastp_output.log")
        file("${sample_name}_trimmed_read_counts.txt")
        file("${sample_name}_fastp_fixed.fastq")

	script:
	"""
	fastp --thread 8 -i $filename -o "${sample_name}_fastp_fixed.fastq" --max_len1 50 2> "${sample_name}_fastp_fixed.log"
	fastp --thread 8 -i "${sample_name}_fastp_fixed.fastq" -o "${sample_name}_fastp_output.fastq" --trim_poly_x -h "${sample_name}_fastp_output.html" -j "${sample_name}_fastp_output.json" --max_len1 70  2> "${sample_name}_fastp_output.log"
	trimmed=\$(cat "${sample_name}_fastp_output.fastq" | wc -l | awk '{printf "%d\\n", \$1 / 4}')
	echo -e "${sample_name}\t\$trimmed" >> "${sample_name}_trimmed_read_counts.txt"
	"""
}

//After trimming, this process removes the UMIs and counts the read counts again
process clean_umi_duplicates {
	executor 'slurm'
	cpus 2
	time '2h'
	memory { 15.GB + 20.GB * task.attempt }
	errorStrategy 'retry'
	module 'build-env/2020:bbmap/38.26-foss-2018b'

	publishDir "${params.outdir}/01_QC/umi_dedup/", mode: 'copy', pattern: '*.{log,stats,length}'

	input:
	val(sample_name)
	file(filename)

	output:
	tuple val(sample_name), file("${sample_name}_umi_dedup_remove_umi.fastq")
	file("${sample_name}_triple_umi_in_R1.fastq")
	file("${sample_name}_umi_dedup.fastq")
	file("${sample_name}_umi_dedup.log") 
	file("${sample_name}_umi_dedup.stats")
	file("${sample_name}_umi_dedup_removed_read_counts.txt")

	// This is done as preparation for the cleaning of the UMI using clumpify.sh
	// As clumpify allow up to two missmatch any difference in the UMI will indicate different instances of a read

	script:
	"""
	cat $filename | awk 'BEGIN {UMI = ""; q="IIIIIIII"} {r_n = (NR%4); if(r_n == 1) {UMI = substr(\$1,length(\$1)-7,length(\$1)); print;}; if(r_n==2) {printf "%s%s%s%s\\n", UMI, UMI, UMI, \$1}; if(r_n==3) {print}; if(r_n==0) {printf "%s%s%s%s\\n", q,q,q,\$1}}' > ${sample_name}_triple_umi_in_R1.fastq

	clumpify.sh in=${sample_name}_triple_umi_in_R1.fastq out=${sample_name}_umi_dedup.fastq dedupe addcount tossjunk 2> ${sample_name}_umi_dedup.log

	cat ${sample_name}_umi_dedup.fastq | awk '(NR % 4) == 1' | grep copies | awk '{print substr(\$NF,8,100)}' | sort -nk1 | uniq -c | sort -nk2 > ${sample_name}_umi_dedup.stats
	cat ${sample_name}_umi_dedup.fastq | wc -l | awk '{printf "%d\\n", \$1 / 4}' > ${sample_name}_umi_dedup.length

	cat ${sample_name}_umi_dedup.fastq | awk '{if((NR%2)==1) {print} else {printf "%s\\n", substr(\$1,25,length(\$1))}}' > ${sample_name}_umi_dedup_remove_umi.fastq
	umi_removed=\$(cat ${sample_name}_umi_dedup_remove_umi.fastq | wc -l | awk '{printf "%d\\n", \$1 / 4}') 
	echo -e "${sample_name}\t\$umi_removed" >> "${sample_name}_umi_dedup_removed_read_counts.txt"
	"""
}

//This process collects all text files of read counts after every step
process combine_raw_counts {
	executor 'slurm'
        cpus 1
        memory '4 GB'
        time '1h'

	publishDir "${params.outdir}/01_QC", mode: 'copy'

	input:
	file(counts)
	file(trimmed_file)
	file(umi_removed)

	output:
	file ("all_raw_read_counts.txt")
	file ("all_trimmed_read_counts.txt")
	file ("all_umi_collapsed_read_counts.txt")

	script:
	'''
	cat *_raw_read_counts.txt >> "all_raw_read_counts.txt"
	cat *_trimmed_read_counts.txt >> "all_trimmed_read_counts.txt"
	cat *_umi_dedup_removed_read_counts.txt >> "all_umi_collapsed_read_counts.txt"
	'''
}

//This process summarizes the read counts by merging all files from process above and plotting bar plots by calling the Rscript 
process Summarize_read_counts {
	executor 'slurm'
	cpus 4
	memory '8 GB'
	time '1h'
	module 'build-env/f2022:r/4.2.0-foss-2021b'
	
	publishDir "${params.outdir}/01_QC/", mode: 'copy'

	input:
	file(rscript)
	file(raw)
	file(trimmed)
	file(umi_removed)

	output:
	file("all_read_counts_summarized.txt")
	file("all_read_counts_summarized.pdf")

	script:
	"""
	#Check if the Rscript exists
	if [ ! -f "${rscript}" ]; then
		echo "Error: Rscript- Summarize_read_counts.R '${rscript}' not found. Make sure that the file is in ${baseDir} or provide the path for Rscript using --rscript "
		exit 1
	fi

	Rscript ${rscript} $raw $trimmed $umi_removed ${params.levels}
	"""
}

//This process generates the fastqc reports for all files after UMI removal
process fastqc {
	executor 'slurm'
	cpus 8
	memory '4 GB'
	time '1h'
	module 'fastqc/0.11.9-java-11'

	publishDir "${params.outdir}/01_QC/fastqc/", mode: 'copy', pattern: '*.{zip,html,txt}'

	input:
	tuple val(sample_name), file(filename)

	output:
	file "*.{zip,html,txt}"

	script:
	"""
	fastqc ${filename}
	"""
}

//This process generates a salmon index using your input fasta as a background genome
process create_salmon_index
{
	executor 'slurm'
	cpus 32
	memory '12 GB'
	time '1h'
	module 'build-env/f2021:salmon/1.5.2-gompi-2020b'

	publishDir "${params.outdir}/04_indices/", mode: 'copy'

	input:
	val gentrome_url
	val decoys_url

	output:
	path("*")
	path("salmon_index")

	script:
	"""
	salmon index -t $gentrome_url -d $decoys_url -i salmon_index/ -p 32
	"""
}

//This process quantifies expression using the salmon index created above
process quantify_exp {
	executor 'slurm'
	cpus { 8 + 8 * task.attempt }
	time '1h'
	memory { 10.GB + 10.GB * task.attempt }
	errorStrategy 'retry'
	module 'build-env/f2021:salmon/1.5.2-gompi-2020b'

	publishDir "${params.outdir}/02_counts/samples/${sample_name}", mode: 'copy', pattern: 'quant*'

	input:
	tuple val(sample_name), file(filename)
	file index

	output:
	file("quant*")

	script:
	"""
	salmon quant -l SF -i $index -r $filename -o quant_${sample_name} --noLengthCorrection
	salmon quant -l SR -i $index -r $filename -o quant_${sample_name}_AS --noLengthCorrection
	"""
}

process download_genome {
	executor 'slurm'
        cpus 4
        memory '8 GB'
        time '1h'
	module 'build-env/2020:gffread/0.11.8-gcccore-8.3.0'
	module 'build-env/2020:samtools/1.10-gcc-8.3.0'

	input:
	val genome_url
        val gff_url

	output:
	tuple file("TAIR10_genome.fa"), file("TAIR10_annotations.gtf")
	file ("chrom_sizes.txt")

	script:
	"""
	cp $genome_url TAIR10_genome.fa
	samtools faidx TAIR10_genome.fa
	cut -f 1,2 TAIR10_genome.fa.fai > chrom_sizes.txt

	cp $gff_url TAIR10_annotations.gff

	gffread TAIR10_annotations.gff -T -o TAIR10_annotations.gtf
	"""
}

process create_star_index {
	executor 'slurm'
        cpus 16
        memory '64 GB'
        time '1h'
        module 'build-env/f2021:star/2.7.9a-gcc-10.2.0'

	publishDir "$params.outdir/04_indices/" , mode: 'copy'

	input:
	tuple file(genome), file(gff)

	output:
	file ("STAR_index/")

	script: 
	"""
	STAR --runThreadN 16 --runMode genomeGenerate \
	--genomeDir "STAR_index/" \
	--genomeFastaFiles $genome --sjdbGTFfile $gff \
	--limitGenomeGenerateRAM 45292192010
	"""
}
	
process star_alignment {
	executor 'slurm'
        cpus 8
        memory '16 GB'
        time '1h'
	module 'build-env/f2021:star/2.7.9a-gcc-10.2.0'

	publishDir "$params.outdir/01_QC/STAR_logs_and_QC/" , mode: 'copy', pattern: '*.{out}'

	input:
	file index_star
	tuple val(sample_name), file(fastq_file)

	output:
	path ("${sample_name}*.bam")
	tuple val (sample_name), file ("${sample_name}_Signal.Unique.str1.out.bg")
	path ("${sample_name}_Log.final.out")
	path("${sample_name}_Log.out")		

	script:
	"""
	STAR --genomeDir $index_star --outFileNamePrefix ${sample_name}_ \
	--readFilesIn  $fastq_file --runThreadN 8 \
	--outSAMtype BAM SortedByCoordinate --outWigType bedGraph --outWigNorm RPM --outWigStrand Unstranded \
	--outFilterMismatchNoverLmax 0.04 --outFilterMultimapNmax 20
	"""
}

process bedgraph2bigwig {
	executor 'slurm'
        cpus 2
        memory '2 GB'
        time '1h'
	module 'build-env/f2021:bedgraphtobigwig/385-linux-x86_64'

	publishDir "$params.outdir/03_STAR_bigwigs/", mode: 'copy', pattern: '{*.bw}'

	input:
	tuple val(sample_name), file(bg_file)
	file chr_sizes

	output:
	file ("${sample_name}.bw")

	script:
	"""
	bedGraphToBigWig $bg_file $chr_sizes "${sample_name}.bw"
	"""
}

process correlation_plots {
	executor 'slurm'
	cpus { 24 + 12 * task.attempt }
        memory { 16.GB + 16.GB * task.attempt }
        time '4h'
	module 'build-env/f2022:deeptools/3.5.4-foss-2022a'

	publishDir "$params.outdir/01_QC/deeptools_plots/", mode: 'copy', pattern: '*.{pdf,tab}'

	input:
	file(bigwig_files)

	output:
	file ("summary.npz")
	file ("*.pdf")
	file("*.tab")

	script:
	"""
	multiBigwigSummary bins --bwfiles $bigwig_files -p 24 --outFileName "summary.npz"

	plotCorrelation --corData summary.npz \
		--corMethod spearman \
		--colorMap RdYlBu \
		--skipZeros \
		--plotNumbers \
		--removeOutliers \
		-p heatmap \
		-o Corr_spearman.pdf \
		--outFileCorMatrix Corr_spearman.tab

	plotCorrelation --corData summary.npz \
		--corMethod pearson \
		--colorMap RdYlBu \
		--skipZeros \
		--plotNumbers \
		--removeOutliers \
		-p heatmap \
		-o Corr_pearson.pdf \
		--outFileCorMatrix Corr_pearson.tab

	plotPCA -in summary.npz \
		-o PCA.pdf \
		--outFileNameData PCA.tab 
	"""
}

workflow {

	//Writing the raw read counts in the metadata
	raw_counts = write_info(samples)
	
	//Downloading genome and annotation files from TAIR10
	downloaded_files = download_genome(params.genome_url, params.gff_url)

	//Creating STAR index from downloaded genome and annotation files
	star_index = create_star_index(downloaded_files[0])
	
	//Creating the salmon index from the downloaded (concatenated) files
	salmon_index = create_salmon_index(params.gentrome_url, params.decoys_url )
	
	//Trimming and removing poly_x
	trimmed_samples = trim_tagseq(samples)
	
	//UMI removal
	umi_cleaned = clean_umi_duplicates(trimmed_samples[0], trimmed_samples[1])

	//Generate fastqc reports
	fastqc(umi_cleaned[0])
	
	//Collecting all read counts
	read_counts = combine_raw_counts(raw_counts.collect(), trimmed_samples[5].collect(), umi_cleaned[5].collect())
	
	//Plotting the read count summary in R
	Summarize_read_counts(rscript, read_counts)
	
	//STAR alignment to get the sorted bam and bedGraph files
	aligned_bam = star_alignment(star_index, umi_cleaned[0])

	//Converting the bedGraph files coming out of STAR alignment to bigwig files
	bigwigs_to_plot = bedgraph2bigwig(aligned_bam[1], downloaded_files[1])

	//Plotting Pearson, Spearmann Correlation and PCA of the bw files
	correlation_plots(bigwigs_to_plot.toSortedList())

	//Quantify gene expression using salmon	
	quantify_exp(umi_cleaned[0], salmon_index[1])

}
