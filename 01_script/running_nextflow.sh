#!/bin/bash

# test files with with 100k reads
sample=test13
nextflow run 01_script/nextflow/main_AtRTD3_ATTE_STAR_mapping_salmon_counts.nf --sample_list 03_sample_lists/sample_list_tagseq_03_sandbox.tsv --outdir 04_output/"$sample" -profile cbe -w "$SCRATCHDIR"/nf_tmp_"$sample" -resume

# tag-seq01
sample=tagseq_01_cdca7_mutants_AtRTD3_ATTE
nextflow run 01_script/nextflow/main_AtRTD3_ATTE_STAR_mapping_salmon_counts.nf --sample_list 03_sample_lists/sample_list_tagseq_01_cdca7.tsv --outdir 04_output/"$sample" -profile cbe -w "$SCRATCHDIR"/nf_tmp_"$sample" --max_n_read 5000000 -resume

# tag-seq03
sample=tagseq_03_cdca7_complementation_AtRTD3_ATTE_150bp_3M_min_50bp
nextflow run 01_script/nextflow/main_AtRTD3_ATTE_STAR_mapping_salmon_counts_mTurq.nf --sample_list 03_sample_lists/sample_list_tagseq_03_cdca7.tsv --outdir 04_output/"$sample" -profile cbe -w "$SCRATCHDIR"/nf_tmp_"$sample" --max_n_read 300000000 -resume

# tag-seq04
sample=tagseq_04_cdca7_complementation_AtRTD3_ATTE
nextflow run 01_script/nextflow/main_AtRTD3_ATTE_STAR_mapping_salmon_counts_mTurq.nf --sample_list 03_sample_lists/sample_list_tagseq_04_cdca7.tsv --outdir 04_output/"$sample" -profile cbe -w "$SCRATCHDIR"/nf_tmp_"$sample" --max_n_read 300000000 -resume

# tag-seq05
sample=tagseq_05_ddm1_EMS_mutants_10M
nextflow run 01_script/nextflow/main_AtRTD3_ATTE_STAR_mapping_salmon_counts.nf --sample_list 03_sample_lists/sample_list_tagseq_05_kanno_ddm1_EMS.tsv --outdir 04_output/"$sample" -profile cbe -w "$SCRATCHDIR"/nf_tmp_"$sample" --max_n_read 10000000

# tag-seq06
sample=tagseq_06_ddm1_alleles
nextflow run 01_script/nextflow/main_AtRTD3_ATTE_STAR_mapping_salmon_counts.nf --sample_list 03_sample_lists/tagseq_06_ddm1_alleles.tsv --outdir 04_output/"$sample" -profile cbe -w "$SCRATCHDIR"/nf_tmp_"$sample" --max_n_read 100000000

# tag-seq07
sample=tagseq_07_remodelers_0.5M
nextflow run 01_script/nextflow/main_AtRTD3_ATTE_STAR_mapping_salmon_counts.nf --sample_list 03_sample_lists/tagseq_07_remodelers.tsv --outdir 04_output/"$sample" -profile cbe -w "$SCRATCHDIR"/nf_tmp_"$sample" --max_n_read 500000 -resume

# tag-seq05 (misnumbered, but is tag-seq 05 according to the ones I generated myself)
sample=tagseq_05_CMT3_CD_and_ddm1_alleles
nextflow run 01_script/nextflow/main_AtRTD3_ATTE_STAR_mapping_salmon_counts.nf --sample_list 03_sample_lists/sample_list_tagseq_05_standard.tsv --outdir 04_output/"$sample" -profile cbe -w "$SCRATCHDIR"/nf_tmp_"$sample" --max_n_read 10000000

# tag-seq05 (misnumbered, but is tag-seq 05 according to the ones I generated myself) with more reads (demultiplexed only on the .2 index)
sample=tagseq_05_CMT3_CD_and_ddm1_demultiplexed_w_index_2_only
nextflow run 01_script/nextflow/main_AtRTD3_ATTE_STAR_mapping_salmon_counts.nf --sample_list 03_sample_lists/sample_list_tagseq_05_index_2_only.tsv --outdir 04_output/"$sample" -profile cbe -w "$SCRATCHDIR"/nf_tmp_"$sample" --max_n_read 10000000 -resume
