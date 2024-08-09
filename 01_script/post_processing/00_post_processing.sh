#!/bin/bash

cd 01_script/post_processing/

# tag-seq 01
# aggregating counts, multiQC
sbatch -p m 01.0_post_processing.sbatch ../../04_output/tagseq_01_cdca7_mutants_AtRTD3_ATTE
# deseq2
sbatch -p m 02.0_DESeq2.sbatch "../../04_output/tagseq_01_cdca7_mutants_AtRTD3_ATTE_100bp/02_counts/" "empty_R1,ddm1_a_long_b_2_R1,F2_ddm1_R1,F2_ddm1_a_2_R1,ddm1_ab_2_R1,ddm1_a_1_R3,a_long_b_R1,b_2_R1,ab_2_R2,ab_1_R1,a_long_2_R2" "mom1,ddm1_mom1,F2_WT,F2_a_2" "../../03_sample_lists/sample_list_tagseq_01_cdca7.tsv" "../../03_sample_lists/DESeq2_tagseq_01_cdca7.tsv" "Col_0"
# all samples with less than 50% uniquely mapping reads were considered outliers (10 samples), but I kept at least 2 replicates per genotype, meaning sometimes samples with low % of unique reads were kept (F2_ddm1_a_2_R3,ab_2_R1). If two samples out of 3 replicates had low unique %, I kept the one with the highest absolute read number. Also removed 2 samples with not enough reads (ddm1_a_long_b_2_R1,ab_2_R2) and a_long_2_R2 (looks contaminated).

# tag-seq 03
sbatch -p m 01.0_post_processing.sbatch ../../04_output/tagseq_03_cdca7_complementation_AtRTD3_ATTE/
sbatch -p m 02.0_DESeq2.sbatch "../../04_output/tagseq_03_cdca7_complementation_AtRTD3_ATTE/02_counts/" "WT_R6,cdca7_ab_R1,cdca7_ab_R6,cdca7_ab_dCter_R4" "none" "../../03_sample_lists/sample_list_tagseq_03_cdca7.tsv" "../../03_sample_lists/DESeq2_tagseq_03_cdca7.tsv" "WT"
# outliers have low number of input reads

# tag-seq 04
sbatch -p m 01.0_post_processing.sbatch ../../04_output/tagseq_04_cdca7_complementation_AtRTD3_ATTE/
sbatch -p m 02.0_DESeq2.sbatch "../../04_output/tagseq_04_cdca7_complementation_AtRTD3_ATTE/02_counts/" "NC_R1,WT_R5,WT_R4" "none" "../../03_sample_lists/sample_list_tagseq_04_cdca7.tsv" "../../03_sample_lists/DESeq2_tagseq_04_cdca7.tsv" "WT"
# outliers have low number of input reads