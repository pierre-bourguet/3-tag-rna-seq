# Overview

This pipeline is meant to analyze short-read transcriptome libraries built with a custom protocol established by Yoav Voichek and Pieter Clauw in the Nordborg lab.
In short, messenger RNAs are reverse-transcribed using a oligo d(T) primer, RNA:cDNA duplexes are cleaved by Tn5, PCR amplifies the 3' ends of transcripts and adds Illumina adapters for sequencing. Further details of the library construction protocol can be found here:

http://UPDATELINKLATER

The pipeline does the following:
- trim reads in 3' to remove adapters, polyA, low quality base calls
- collapse duplicates based on 8-bp UMIs (no UMI mismatch allowed, two mismatches allowed in the transcript)
- randomly subsample to a given number of reads (default: 10 millions)
- create indexes based on provided genome files
- map with STAR
- output stranded and unstranded bigwig files
- quantify sense and antisense reads with Salmon in alignment-mode
- generate plots to inspect quality metrics: multiplexing, trimming, duplicates, multiQC, sample to sample correlation
- analyze relative to a control condition with DESeq2, producing tables of differentially expressed genes, heatmaps, barplots, PCA

# Usage
## Prepare genome files
Execute the code in 
```
cd 02_genome_files/01_script/
./prepare_genome_files.sh
./gene_annotations_for_deseq2.sh
```
The scripts are not perfectly wrapped, you have to execute each block in an interactive session and make sure everything goes well, but you have to do this only once. If TAIR blocks your wget queries, you can download files from a web browser.

## Mapping and counting reads
- Prepare a sample_list file. Example:
```
head 03_sample_lists/sample_list_tagseq_03_cdca7.tsv
```
- Run the main script. Example:
```
sample=test13
nextflow run 01_script/nextflow/main_AtRTD3_ATTE_STAR_mapping_salmon_counts.nf --sample_list 03_sample_lists/sample_list_tagseq_03_sandbox.tsv --outdir 04_output/"$sample" -profile cbe -w "$SCRATCHDIR"/nf_tmp_"$sample" --max_n_read 5000000
```

- Run post processing steps to aggregate read counts, plot quality metrics and run MultiQC. Example:
```
cd 01_script/post_processing/
sbatch 01.0_post_processing.sbatch ../../04_output/test13
```

## DESeq2
A standard analysis of differentially expressed genes (DEGs). You can remove outliers (e.g. : samples with low read counts), or even exclude all samples matching a string pattern, to test the influence of some samples on the DESeq2 outcome. By default, the analysis focuses on protein-coding genes and transposable elements, but you can control this by modifying the DESeq2 environment file.

- Prepare a DESeq2 sample_list file to describe your experiment design. Examples in:
```
cat 02_sample_lists/DESeq2_tagseq_03_cdca7.tsv
```
- (optional) Modify the DESeq2 environment file, to change annotations of interest. Default annotations are protein-coding genes and TAIR10 transposable elements.
```
vi 01_script/post_processing/02.2_DESeq2_environment.R
```
- Run a DESeq2 analysis. It runs DESeq2, normalizes and transforms counts, plots a PCA and a sample correlation heatmap, writes count tables, identifies differentially expressed genes (DEGs) for all treatment versus control comparisons and write tables, plots heatmaps of normalized & transformed counts at DEGs, write barplots of the number of DEGs per sample. It finally saves the environment and creates a script in "07_analysis" where you can start you own analysis.
The order of arguments are: 
```
sbatch 02.0_DESeq2.sbatch output_directory \
    outliers \
    outlier_pattern \
    sample_list_file \
    sample_DESeq2_file \
    control_condition
```
Example:
```
sbatch 02.0_DESeq2.sbatch "../../04_output/tagseq_03_cdca7_complementation_AtRTD3_ATTE/02_counts/" \
    "WT_R6,cdca7_ab_R1,cdca7_ab_R6,cdca7_ab_dCter_R4" \
    "none" \
    "../../03_sample_lists/sample_list_tagseq_03_cdca7.tsv" \
    "../../03_sample_lists/DESeq2_tagseq_03_cdca7.tsv" \
    "WT"
```

# Contributions
Yoav Voichek developed the original pipeline. Vikas Shukla further improved it and put it on github. I created a fork from Vikas's pipeline. The main modifications I did are the following:
- reference transcriptome changed from TAIR10 to AtRTD3
- transposable_element_gene annotations removed and replaced by TAIR10 AT.TE annotations. This improves the counting of 3' fragments and avoids ambiguity issues with counting reads at overlapping annotations
- 3' adapters and polyA are trimmed more efficiently, increasing mapping rates by ~100%
- incorporated read downsampling from Yoav's original code
- map with STAR and count reads with salmon in alignment mode, instead of salmon pseudoalignment
- count reads in sense and antisense
- output stranded and unstranded bigwig files
- created an alternative version of the genome files to include a transgenic construct
- incorporated the DESeq2 downstream analysis

I also added a post processing step which generates:
- files with aggregated counts
- multiQC
- tables and plots to summarize statistics for trimming, collapsing, mapping
