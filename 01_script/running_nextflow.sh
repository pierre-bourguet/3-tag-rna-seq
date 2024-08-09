#!/bin/bash

# test fastq with 100k reads
sample=test13
nextflow run 01_script/nextflow/main_AtRTD3_ATTE_STAR_mapping_salmon_counts.nf --sample_list 03_sample_lists/sample_list_tagseq_03_sandbox.tsv --outdir 04_output/"$sample" -profile cbe -w "$SCRATCHDIR"/nf_tmp_"$sample" --reference_genome tair10 -resume

# tag-seq01
sample=tagseq_01_cdca7_mutants_AtRTD3_ATTE
nextflow run 01_script/nextflow/main_AtRTD3_ATTE_STAR_mapping_salmon_counts.nf --sample_list 03_sample_lists/sample_list_tagseq_01_cdca7.tsv --outdir 04_output/"$sample" -profile cbe -w "$SCRATCHDIR"/nf_tmp_"$sample" --reference_genome tair10 -resume --max_n_read 5000000

# tag-seq03
sample=tagseq_03_cdca7_complementation_AtRTD3_ATTE
nextflow run 01_script/nextflow/main_AtRTD3_ATTE_STAR_mapping_salmon_counts_mTurq.nf --sample_list 03_sample_lists/sample_list_tagseq_03_cdca7.tsv --outdir 04_output/"$sample" -profile cbe -w "$SCRATCHDIR"/nf_tmp_"$sample" --reference_genome tair10 -resume

# tag-seq04
sample=tagseq_04_cdca7_complementation_AtRTD3_ATTE
nextflow run 01_script/nextflow/main_AtRTD3_ATTE_STAR_mapping_salmon_counts_mTurq.nf --sample_list 03_sample_lists/sample_list_tagseq_04_cdca7.tsv --outdir 04_output/"$sample" -profile cbe -w "$SCRATCHDIR"/nf_tmp_"$sample" --reference_genome tair10 -resume

