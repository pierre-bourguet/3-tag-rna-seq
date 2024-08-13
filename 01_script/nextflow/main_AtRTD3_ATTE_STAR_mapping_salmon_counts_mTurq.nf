#!/usr/bin/env nextflow

// Define a process parameter for the numerical argument
params.max_n_read = params.max_n_read ?: 10000000

// Access the argument in the code
def myArg = params.max_n_read
println "The maximum number of reads is: $myArg"

ch_max_n_read = Channel.value(params.max_n_read)

// Channel created from input file and mapping of file name and path
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
params.cdna = "/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/02_genome_files/03_output/agat_fix_gff/atRTD3_TS_21Feb22_transfix_no_TEG_w_ATTE_fixed_mTurq_transgene.fa"
params.genome = "/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/02_genome_files/03_output/TAIR10_chr_all_mTurq_transgene.fas"
params.gff = "/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/02_genome_files/03_output/agat_fix_gff/atRTD3_TS_21Feb22_transfix_no_TEG_w_ATTE_fixed_mTurq_transgene.gff"
params.gentrome = "/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/02_genome_files/03_output/salmon_decoys/gentrome_atRTD3_29122021_no_TEG_w_ATTE_fixed_mTurq_transgene.fna"
params.decoys = "/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/02_genome_files/03_output/salmon_decoys/decoys_transgene.txt"
params.adapters = "/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/02_genome_files/02_input/adapter_fasta.fa"

// import downsample process from downsample.nf
include {downsample} from './downsample'
params.downsample = -1

// This process counts the no. of reads in the raw files for each sample 
process write_info {
	executor 'slurm'
	cpus 1
	memory '1 GB'
	time '5m'

	publishDir "${params.outdir}/01_QC/raw_read_counts/", mode: 'copy'

	input:
	tuple path(filename), val(sample_name)

	output:
	path("${sample_name}_raw_read_counts.txt") 

	script:
	"""
	read_count=\$(zcat $filename | wc -l | awk '{printf "%d\\n", \$1 / 4}')
	echo -e "${sample_name}\t\$read_count" >> ${sample_name}_raw_read_counts.txt
	"""
}

// This process trims the reads to 50 bp, removes the poly_x and counts the no. of reads after trimming
process trim_tagseq {
	executor 'slurm'
	cpus 8
        memory { 8.GB + 16.GB * task.attempt }
	time '20m'
	module 'build-env/f2022:fastp/0.23.4-gcc-12.2.0'

	publishDir "${params.outdir}/01_QC/fastp_trimming/", mode: 'copy', pattern: '*.{log,length,json,html}'

	input:
	tuple path(filename), val(sample_name)
	val adapters

	output:
	val(sample_name)
        path("${sample_name}_fastp_output.fastq")
        path("${sample_name}_fastp.html")
        path("${sample_name}_fastp.json")
        path("${sample_name}_fastp.log")
        path("${sample_name}_trimmed_read_counts.txt")

	script:
	"""
	# trim low quality bases, poly X & adapters; remove low quality bases at the tail (3' end)
	fastp --thread 8 -i $filename -o "${sample_name}_fastp_output.fastq" \
	--trim_poly_x --adapter_fasta $adapters --cut_tail --max_len1 100 \
	 -h "${sample_name}_fastp.html" -j "${sample_name}_fastp.json" 2> "${sample_name}_fastp.log"

	trimmed=\$(cat "${sample_name}_fastp_output.fastq" | wc -l | awk '{printf "%d\\n", \$1 / 4}')
	echo -e "${sample_name}\t\$trimmed" >> "${sample_name}_trimmed_read_counts.txt"
	"""
}

// After trimming, this process removes the UMIs and counts the read counts again
process clean_umi_duplicates {
	executor 'slurm'
	cpus 2
	time '1h'
	memory { 10.GB + 20.GB * task.attempt }
	errorStrategy 'retry'
	module 'build-env/f2022:bbmap/39.08-gcc-12.2.0'

	publishDir "${params.outdir}/01_QC/umi_dedup/", mode: 'copy', pattern: '*.{log,stats,length}'

	input:
	val(sample_name)
	path(filename)

	output:
	tuple val(sample_name), file("${sample_name}_umi_dedup_remove_umi.fastq")
	path("${sample_name}_triple_umi_in_R1.fastq")
	path("${sample_name}_umi_dedup.fastq")
	path("${sample_name}_umi_dedup.log") 
	path("${sample_name}_umi_dedup.stats")
	path("${sample_name}_umi_dedup_removed_read_counts.txt")
	path("${sample_name}_umi_dedup.length")

	// This is done as preparation for the cleaning of the UMI using clumpify.sh
	// As clumpify allow up to two missmatch any difference in the UMI will indicate different instances of a read

	script:
	"""
	# Triple UMI in R1 Sequence
	cat $filename | awk 'BEGIN {UMI = ""; q="IIIIIIII"} {r_n = (NR%4); if(r_n == 1) {UMI = substr(\$1,length(\$1)-7,length(\$1)); print;}; if(r_n==2) {printf "%s%s%s%s\\n", UMI, UMI, UMI, \$1}; if(r_n==3) {print}; if(r_n==0) {printf "%s%s%s%s\\n", q,q,q,\$1}}' > ${sample_name}_triple_umi_in_R1.fastq

	# Deduplicate Reads Using UMI
	clumpify.sh in=${sample_name}_triple_umi_in_R1.fastq out=${sample_name}_umi_dedup.fastq dedupe addcount tossjunk 2> ${sample_name}_umi_dedup.log

	# Generate UMI Deduplication Statistics
	cat ${sample_name}_umi_dedup.fastq | awk '(NR % 4) == 1' | grep copies | awk '{print substr(\$NF,8,100)}' | sort -nk1 | uniq -c | sort -nk2 > ${sample_name}_umi_dedup.stats
	
	# Count Total Deduplicated Reads
	cat ${sample_name}_umi_dedup.fastq | wc -l | awk '{printf "%d\\n", \$1 / 4}' > ${sample_name}_umi_dedup.length

	# Remove UMI from Sequence
	cat ${sample_name}_umi_dedup.fastq | awk '{if((NR%2)==1) {print} else {printf "%s\\n", substr(\$1,25,length(\$1))}}' > ${sample_name}_umi_dedup_remove_umi.fastq
	
	# Count Reads After UMI Removal
	umi_removed=\$(cat ${sample_name}_umi_dedup_remove_umi.fastq | wc -l | awk '{printf "%d\\n", \$1 / 4}') 
	echo -e "${sample_name}\t\$umi_removed" >> "${sample_name}_umi_dedup_removed_read_counts.txt"
	"""
}

// This process collects all text files of read counts after every step
process combine_raw_counts {
	executor 'slurm'
        cpus 1
        memory '0.1 GB'
        time '5m'

	publishDir "${params.outdir}/01_QC", mode: 'copy'

	input:
	path(counts)
	path(trimmed_file)
	path(umi_removed)

	output:
	path ("all_raw_read_counts.txt")
	path ("all_trimmed_read_counts.txt")
	path ("all_umi_collapsed_read_counts.txt")

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
	cpus 1
	memory '1 GB'
	time '5m'
	module 'build-env/f2022:r/4.2.0-foss-2021b'
	
	publishDir "${params.outdir}/01_QC/", mode: 'copy'

	input:
	path(rscript)
	path(raw)
	path(trimmed)
	path(umi_removed)

	output:
	path("all_read_counts_summarized.txt")
	path("all_read_counts_summarized.pdf")

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

// This process generates the fastqc reports for all files after UMI removal
process fastqc {
	executor 'slurm'
	cpus 1
	memory '1 GB'
	time '20m'
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

process download_genome {
	executor 'slurm'
        cpus 4
        memory '8 GB'
        time '20m'
	module 'build-env/f2022:gffread/0.12.7-gcccore-12.2.0'
	module 'build-env/f2022:samtools/1.18-gcc-12.3.0'

	publishDir "${params.outdir}/04_indices/", mode: 'copy', pattern: '*.{gtf}'

	input:
	val genome
    val gff

	output:
	tuple file("TAIR10_genome.fa"), file("TAIR10_annotations.gtf")
	path ("chrom_sizes.txt")

	script:
	"""
	cp $genome TAIR10_genome.fa
	samtools faidx TAIR10_genome.fa
	cut -f 1,2 TAIR10_genome.fa.fai > chrom_sizes.txt

	cp $gff TAIR10_annotations.gff

	gffread TAIR10_annotations.gff -T -o TAIR10_annotations.gtf
	"""
}

process create_star_index {
	executor 'slurm'
        cpus 16
        memory '32 GB'
        time '1h'
	module 'build-env/f2022:star/2.7.11a-gcc-12.3.0'

	publishDir "$params.outdir/04_indices/" , mode: 'copy'

	input:
	tuple file(genome), file(gtf)

	output:
	file ("STAR_index/")

	script: 
	"""
	STAR --runThreadN 16 --runMode genomeGenerate \
	--genomeDir "STAR_index/" \
	--genomeFastaFiles $genome --sjdbGTFfile $gtf \
	--limitGenomeGenerateRAM 45292192010
	"""
}
	
process star_alignment {
	executor 'slurm'
        cpus 8
        memory { 12.GB + 4.GB * task.attempt }
        time '1h'
	module 'build-env/f2022:star/2.7.11a-gcc-12.3.0'

	publishDir "$params.outdir/01_QC/STAR_logs_and_QC/" , mode: 'copy', pattern: '*.{out}'

	input:
	file index_star
	tuple val(sample_name), file(fastq_file)

	output:
	tuple val (sample_name), path ("${sample_name}*.out.bg")
	path ("${sample_name}_Log.final.out")
	path ("${sample_name}_Log.out")
    tuple val(sample_name), path ("${sample_name}_Aligned.toTranscriptome.out.bam")

	script:
	"""
	STAR --genomeDir $index_star --outFileNamePrefix ${sample_name}_ \
	--readFilesIn $fastq_file --runThreadN 8 \
	--outSAMtype BAM SortedByCoordinate --outWigType bedGraph --outWigNorm RPM --outWigStrand Stranded \
	--outFilterMismatchNoverLmax 0.04 --outFilterMultimapNmax 50 \
    --quantMode TranscriptomeSAM --quantTranscriptomeBan Singleend --limitBAMsortRAM 8287822014
	"""
}

// This process generates a salmon index using your input fasta as a background genome
process create_salmon_index
{
	executor 'slurm'
	cpus 16
	memory { 24.GB + 12.GB * task.attempt }
	time '10m'
	module 'build-env/f2022:salmon/1.10.1-gcc-12.2.0'

	publishDir "${params.outdir}/04_indices/", mode: 'copy'

	input:
	val gentrome
	val decoys

	output:
	path("*")
	path("salmon_index")

	script:
	"""
	salmon index -t $gentrome -d $decoys -i salmon_index/ -p 16
	"""
}

// This process quantifies expression using the salmon index created above
process quantify_exp {
	executor 'slurm'
	cpus 8
	time '1h'
	memory { 4.GB + 8.GB * task.attempt }
	errorStrategy 'retry'
	module 'build-env/f2022:salmon/1.10.1-gcc-12.2.0'

	publishDir "${params.outdir}/02_counts/samples/${sample_name}", mode: 'copy', pattern: '{quant}*'

	input:
	tuple val(sample_name), path(filename)
	path index

	output:
	path("quant*")

	script:
	"""
	salmon quant -p 8 -l SF -i $index -r $filename -o quant_${sample_name} --noLengthCorrection
	salmon quant -p 8 -l SR -i $index -r $filename -o quant_${sample_name}_AS --noLengthCorrection
	"""
}

// This process quantifies expression using the STAR-aligned BAM file. It requires no more than 8 threads, memory is limiting
process quantify_exp_STAR_mapping {
	executor 'slurm'
	cpus 8
	time '1h'
	memory { 24.GB + 24.GB * task.attempt }
	errorStrategy 'retry'
	module 'build-env/f2022:salmon/1.10.1-gcc-12.2.0'

	publishDir "${params.outdir}/02_counts/samples/${sample_name}", mode: 'copy', pattern: '{STAR_mapping_salmon_quant}*'

	input:
    tuple val(sample_name), path (bamfile)
    val cdna

	output:
	path("STAR_mapping_salmon_quant*")

	script:
	"""
    salmon quant -p 8 -l SF -t $cdna -a $bamfile -o STAR_mapping_salmon_quant_${sample_name} --noLengthCorrection
	salmon quant -p 8 -l SR -t $cdna -a $bamfile -o STAR_mapping_salmon_quant_${sample_name}_AS --noLengthCorrection
	"""
}

process bedgraph2bigwig {
	executor 'slurm'
		cpus 2
		memory '8 GB'
		time '1h'
	module 'build-env/f2021:bedgraphtobigwig/385-linux-x86_64'


	publishDir "$params.outdir/03_STAR_bigwigs/", mode: 'copy', pattern: '**/*bw'

	input:
	tuple val(sample_name), path (bg_files)
	path chr_sizes

	output:
	path ("unique_multi_unstranded/${sample_name}*multi.unstranded.bw")
	path ("unique_unstranded/${sample_name}_unique.unstranded.bw")
	path ("unique_multi_stranded/${sample_name}*multi.str?.bw")
	path ("unique_stranded/${sample_name}_unique.str?.bw")

	script:
	"""
	# convert STAR bedgraphs to bigwig
	bedGraphToBigWig "${sample_name}_Signal.Unique.str1.out.bg" $chr_sizes "${sample_name}_unique.str1.bw"
	bedGraphToBigWig "${sample_name}_Signal.UniqueMultiple.str1.out.bg" $chr_sizes "${sample_name}_unique_multi.str1.bw"
	bedGraphToBigWig "${sample_name}_Signal.Unique.str2.out.bg" $chr_sizes "${sample_name}_unique.str2.bw"
	bedGraphToBigWig "${sample_name}_Signal.UniqueMultiple.str2.out.bg" $chr_sizes "${sample_name}_unique_multi.str2.bw"
	
	# add stranded bigwigs from both strands to compute unstranded bedgraphs
	wget --quiet -O bigWigMerge http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigMerge
	bigWigMerge "${sample_name}_unique.str1.bw" "${sample_name}_unique.str2.bw" "${sample_name}_unique.unstranded.bg"
	bigWigMerge "${sample_name}_unique_multi.str1.bw" "${sample_name}_unique_multi.str2.bw" "${sample_name}_unique_multi.unstranded.bg"

	# convert unstranded bedgraphs to bigwigs
	bedGraphToBigWig "${sample_name}_unique.unstranded.bg" $chr_sizes "${sample_name}_unique.unstranded.bw"
	bedGraphToBigWig "${sample_name}_unique_multi.unstranded.bg" $chr_sizes "${sample_name}_unique_multi.unstranded.bw"

	# organize outputs
	mkdir -p unique_stranded unique_multi_stranded unique_unstranded unique_multi_unstranded
	mv "${sample_name}_unique.str1.bw" "${sample_name}_unique.str2.bw" unique_stranded
	mv "${sample_name}_unique_multi.str1.bw" "${sample_name}_unique_multi.str2.bw" unique_multi_stranded
	mv "${sample_name}_unique.unstranded.bw" unique_unstranded
	mv "${sample_name}_unique_multi.unstranded.bw" unique_multi_unstranded
	"""
}

process correlation_plots {
	executor 'slurm'
	cpus { 6 + 12 * task.attempt }
    memory { 4.GB + 16.GB * task.attempt }
    time '1h'
	module 'build-env/f2022:deeptools/3.5.4-foss-2022a'

	publishDir "$params.outdir/01_QC/deeptools_plots/", mode: 'copy', pattern: '*.{pdf,tab}'

	input:
	path (bigwig_files)

	output:
	path ("summary.npz")
	path ("*.pdf")
	path ("*.tab")

	script:
	"""
	multiBigwigSummary bins --bwfiles $bigwig_files -p 24 --outFileName "summary.npz"

	plotCorrelation --corData summary.npz \
		--corMethod spearman \
		--colorMap RdYlBu \
		--skipZeros \
		--removeOutliers \
		-p heatmap \
		-o Corr_spearman.pdf \
		--outFileCorMatrix Corr_spearman.tab

	plotCorrelation --corData summary.npz \
		--corMethod pearson \
		--colorMap RdYlBu \
		--skipZeros \
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

	// Writing the raw read counts in the metadata
	raw_counts = write_info(samples)
	
	// Trimming and removing poly_x
	trimmed_samples = trim_tagseq(samples, params.adapters)
	
	// UMI removal
	umi_cleaned = clean_umi_duplicates(trimmed_samples[0], trimmed_samples[1])

	// downsampling reads
	ds_reads = downsample(umi_cleaned[0].combine(ch_max_n_read).combine(Channel.from(1)))

	// Generate fastqc reports
	fastqc(ds_reads[0])
	
	// Collecting all read counts
	read_counts = combine_raw_counts(raw_counts.collect(), trimmed_samples[5].collect(), umi_cleaned[5].collect())
		
	//Plotting the read count summary in R
	Summarize_read_counts(rscript, read_counts)

	// Downloading genome and annotation files from TAIR10
	downloaded_files = download_genome(params.genome, params.gff)

	// Creating STAR index from downloaded genome and annotation files
	star_index = create_star_index(downloaded_files[0])
	
	// Creating the salmon index from the downloaded (concatenated) files
	salmon_index = create_salmon_index(params.gentrome, params.decoys )
	
	// STAR alignment to get the sorted bam and bedGraph files
	aligned_bam = star_alignment(star_index, ds_reads[0])

	// Quantify gene expression using salmon	
	quantify_exp(ds_reads[0], salmon_index[1])

	// Quantify gene expression using salmon	
	quantify_exp_STAR_mapping(aligned_bam[3], params.cdna)

	// Converting the bedGraph files coming out of STAR alignment to bigwig files
	bigwigs_to_plot = bedgraph2bigwig(aligned_bam[0], downloaded_files[1])

	// Plotting Pearson, Spearmann Correlation and PCA of the bw files
	correlation_plots(bigwigs_to_plot[0].toSortedList())

}
