# Overview

This pipeline is meant to analyze short-read transcriptome libraries built with a custom protocol established by Yoav Voichek and Pieter Clauw in the Nordborg lab.
In short, messenger RNAs are reverse-transcribed using a oligo d(T) primer, RNA:cDNA duplexes are cleaved by Tn5, PCR amplifies the 3' ends of transcripts and adds Illumina adapters for sequencing. Further details of the library construction protocol can be found here:

# Usage
## Mapping
- Prepare a sample_list file. Examples in:
02_sample_lists/

- Run the main script from the 01_script/nextflow/ folder. Commands are found in 01_script/running_nextflow.sh. Example:
sample=test13
nextflow run 01_script/nextflow/main_AtRTD3_ATTE_STAR_mapping_salmon_counts.nf --sample_list 03_sample_lists/sample_list_tagseq_03_sandbox.tsv --outdir 04_output/"$sample" -profile cbe -w "$SCRATCHDIR"/nf_tmp_"$sample" --reference_genome tair10 -resume

- Run post processing steps to aggregate read counts, plot quality metrics and run MultiQC. Commands are found in 01_script/post_processing/00_post_processing.sh. Example:
cd 01_script/post_processing/
sbatch 01.0_post_processing.sbatch ../../04_output/test13

## DESeq2
- Prepare a DESeq2 sample_list file to describe your experiment design. Examples in:
02_sample_lists/DESeq2_tagseq_03_cdca7.tsv

- (optional) Modify the DESeq2 environment file, to change annotations of interest. Default annotations are protein-coding genes and TAIR10 transposable elements.
vi 01_script/post_processing/02.2_DESeq2_environment.R

- Run a DESeq2 analysis. It runs DESeq2, normalizes and transforms counts, plots a PCA and a sample correlation heatmap, writes count tables, identifies differentially expressed genes (DEGs) for all treatment versus control comparisons and write tables, plots heatmaps of normalized & transformed counts at DEGs, write barplots of the number of DEGs per sample. It finally saves the environment and creates a script in "07_analysis" where you can start you own analysis.
The order of arguments are: 
sbatch 02.0_DESeq2.sbatch output_directory outliers sample_list_file sample_DESeq2_file control_condition
Example:
sbatch 02.0_DESeq2.sbatch "../../04_output/tagseq_03_cdca7_complementation_AtRTD3_ATTE/02_counts/" "WT_R6,cdca7_ab_R1,cdca7_ab_R6,cdca7_ab_dCter_R4" "none" "../../03_sample_lists/sample_list_tagseq_03_cdca7.tsv" "../../03_sample_lists/DESeq2_tagseq_03_cdca7.tsv" "WT"

# Contributions
Yoav Voichek developed the original pipeline. Vikas Shukla further improved it and put it on github. I created a fork from Vikas's pipeline. The main modifications I did are the following:
- 3' adapters and polyA are trimmed more efficiently, increasing mapping rates by ~100%
- reference transcriptome changed from TAIR10 to AtRTD3
- transposable_element_gene annotations removed and replaced by TAIR10 AT.TE annotations. This improves the counting of 3' fragments and avoids ambiguity issues with counting reads at overlapping annotations
- incorporated read downsampling from Yoav's original code
- map with STAR and count reads with salmon in alignment mode, instead of salmon pseudoalignment
- count reads in sense and antisense
- created an alternative version of the genome files to include a transgenic construct

I also added a post processing step which generates:
- files with aggregated counts
- multiQC
- tables and plots to summarize statistics for trimming, collapsing, mapping

I incorporated an automated DESeq2 analysis of the data, where one can control outliers to remove after a first analysis (e.g. : samples with low read counts). One can also test the influence of different samples on the DESeq2 outcome by excluding samples based on string patterns. One can control the type of annotations analyzed by modifying the DESeq2 environment file, where default includes protein-coding genes and transposable elements.