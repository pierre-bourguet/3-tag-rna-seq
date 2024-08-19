# Overview

This pipeline is meant to analyze short-read transcriptome libraries built with a custom protocol established by Yoav Voichek and Pieter Clauw in the Nordborg lab.
In short, messenger RNAs are reverse-transcribed using a oligo d(T) primer, RNA:cDNA duplexes are cleaved by Tn5, PCR amplifies the 3' ends of transcripts and adds Illumina adapters for sequencing. Further details of the library construction protocol can be found here:

http://UPDATELINKLATER

The pipeline does the following:
- trim reads in 3' to remove adapters, polyA, low quality base calls, reads shorter than 50bp
- collapse duplicates based on 8-bp UMIs (mismatch allowed: 0 for UMIs, 2 in the transcript)
- randomly subsample to a given number of reads (default: 10 millions)
- create indexes based on provided genome files
- map with STAR
- output bigwig files, stranded or unstranded, including or excluding multi-mapping reads
- quantify sense and antisense reads with Salmon in alignment-mode
- generate plots to inspect quality metrics: multiplexing, trimming, duplicates, multiQC, sample to sample correlation
- analyze differential gene expression with DESeq2, producing tables of differentially expressed genes, heatmaps, barplots, PCA

# Usage
## Prepare genome files
Execute the code in 
```shell
cd 02_genome_files/01_script/
./prepare_genome_files.sh
./gene_annotations_for_deseq2.sh
```
The scripts are not perfectly wrapped, you have to execute each block in an interactive session and make sure everything goes well, but you have to do this only once. If TAIR blocks your wget queries, you can download files from a web browser.

## Mapping and counting reads
- Prepare a sample_list file. Example:
```shell
head 03_sample_lists/sample_list_tagseq_03_cdca7.tsv
```
- Run the main script. Example:
```shell
sample=test13
nextflow run 01_script/nextflow/main_AtRTD3_ATTE_STAR_mapping_salmon_counts.nf --sample_list 03_sample_lists/sample_list_tagseq_03_sandbox.tsv --outdir 04_output/"$sample" -profile cbe -w "$SCRATCHDIR"/nf_tmp_"$sample" --max_n_read 5000000
```

- Run post processing steps to aggregate read counts, plot quality metrics and run MultiQC. Example:
```shell
cd 01_script/post_processing/
sbatch 01.0_post_processing.sbatch ../../04_output/test13
```

## DESeq2
An analysis of differentially expressed genes (DEGs). You can remove outliers (e.g. : samples with low read counts), or exclude all samples matching a string pattern, to test their influence on the outcome. By default, the analysis focuses on protein-coding genes and transposable elements, but you can control this by modifying the DESeq2 environment file.

- (optional) Modify the DESeq2 environment file, to change annotations of interest. Default annotations are protein-coding genes and TAIR10 transposable elements.
```shell
vi 01_script/post_processing/02.2_DESeq2_environment.R
```
- Run the DESeq2 script. It does the following:
	- run DESeq2, normalize with estimated size factors (ESF) and transform counts using a regularized log (rlog).
	- find differentially expressed genes (DEGs) for all treatments relative to a control condition (absolute log2FC >= 1, P < 0.1)	
	- produce plots: PCA, heatmap of sample-to-sample correlations (using variance-stabilized counts), number of DEGs per sample, heatmaps of normalized & transformed counts at DEGs
	- write tables: normalized (ESF) & transformed (rlog) counts, pairwise comparisons, DEGs. RPM tables are also provided based on Salmon outputs.
	- saves the environment and creates a script in "07_analysis" where you can start you own analysis using this environment.

The order of arguments are: 
```shell
sbatch 02.0_DESeq2.sbatch output_directory \
    outliers \
    outlier_pattern \
    sample_list_file \
    control_condition
```
Example:
```shell
sbatch 02.0_DESeq2.sbatch "../../04_output/tagseq_03_cdca7_complementation_AtRTD3_ATTE/02_counts/" \
    "WT_R6,cdca7_ab_R1,cdca7_ab_R6,cdca7_ab_dCter_R4" \
    "none" \
    "../../03_sample_lists/sample_list_tagseq_03_cdca7.tsv" \
    "WT"
```

# Contributions
Yoav Voichek developed the original nextflow pipeline. Vikas Shukla improved it and put it on github. I forked from Vikas's pipeline and did the following modifications:
- reference transcriptome changed from TAIR10 to AtRTD3
- transposable_element_gene annotations removed and replaced by TAIR10 AT.TE annotations. This improves the counting of 3' fragments and avoids ambiguity issues with counting reads at overlapping annotations
- to improve mapping rates, 3' adapters and polyA are trimmed explicitly and reads shorter than 50bp are discarded
- incorporated read downsampling from Yoav's original code
- map with STAR and count reads with salmon in alignment mode, instead of salmon pseudoalignment
- count reads in sense and antisense
- output stranded and unstranded bigwig files
- reorganized output folder architecture
- created an alternative version of the genome files to include a transgenic construct

I added a post processing step which generates:
- files with aggregated counts
- multiQC
- tables and plots to summarize statistics for trimming, collapsing, mapping

I incorporated the DESeq2 downstream analysis.