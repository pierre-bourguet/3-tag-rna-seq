#!/bin/bash

cd 01_script/post_processing/

# tag-seq 01
# aggregating counts, multiQC
sbatch -p m 01.0_post_processing.sbatch ../../04_output/tagseq_01_cdca7_mutants_AtRTD3_ATTE_150bp_5M_min_50bp
# deseq2

# outliers: low counts, low unique mappers, contaminated samples
sbatch -p c 02.0_DESeq2.sbatch \
	"../../04_output/tagseq_01_cdca7_mutants_AtRTD3_ATTE_150bp_5M_min_50bp/02_counts/" \
	"empty_R1,ab_1_R1,ab_2_R2,ddm1_a_long_b_2_R1,F2_ddm1_R1,F2_ddm1_a_2_R1,ddm1_ab_2_R1,ddm1_a_1_R3,a_long_b_R1,b_2_R1,a_2_R3,a_long_1_R3,a_long_2_R1,a_long_2_R2,ddm1_2_G2_R2" \
	"mom1,F2_WT,F2_a_2" \
	"../../03_sample_lists/sample_list_tagseq_01_cdca7.tsv" \
	"Col_0"

# all samples with less than 50% uniquely mapping reads were considered outliers (10 samples), but I kept at least 2 replicates per genotype, meaning sometimes samples with low % of unique reads were kept (F2_ddm1_a_2_R3,ab_2_R1). If two samples out of 3 replicates had low unique %, I kept the one with the highest absolute read number. Also removed 2 samples with not enough reads (ddm1_a_long_b_2_R1,ab_2_R2) and samples that look contaminated (F2_ddm1_R1, a_2_R3, a_long_1_R3, a_long_2_R1, a_long_2_R2, b_2_R1, ddm1_2_G2_R2).


# tag-seq 03
sbatch -p m 01.0_post_processing.sbatch ../../04_output/tagseq_03_cdca7_complementation_AtRTD3_ATTE_150bp_3M_min_50bp/
sbatch -p c 02.0_DESeq2.sbatch "../../04_output/tagseq_03_cdca7_complementation_AtRTD3_ATTE_150bp_3M_min_50bp_unique_mappers/02_counts/" "WT_R6,cdca7_ab_R1,cdca7_ab_R6,cdca7_ab_dCter_R4" "none" "../../03_sample_lists/sample_list_tagseq_03_cdca7.tsv" "WT"
# outliers have low number of input reads

# tag-seq 03 individualized
sbatch -p c 02.0_DESeq2_individualized.sbatch "../../04_output/tagseq_03_cdca7_complementation_AtRTD3_ATTE_150bp_3M_min_50bp_unique_mappers_DESeq2_rerun/02_counts/" "WT_R6,cdca7_ab_R1,cdca7_ab_R6,cdca7_ab_dCter_R4" "none" "../../03_sample_lists/sample_list_tagseq_03_cdca7_individualized.tsv" "WT"

# tag-seq 04
sbatch -p c 01.0_post_processing.sbatch ../../04_output/tagseq_04_cdca7_complementation_AtRTD3_ATTE_150bp_3M_min_50bp_unique_mappers/
sbatch -p c 02.0_DESeq2.sbatch "../../04_output/tagseq_04_cdca7_complementation_AtRTD3_ATTE_150bp_3M_min_50bp_unique_mappers/02_counts/" "NC_R1,WT_R5,WT_R4" "none" "../../03_sample_lists/sample_list_tagseq_04_cdca7.tsv" "WT"
# outliers have low number of input reads

# tagseq 05 ddm1 EMS mutants
sbatch -p m 01.0_post_processing.sbatch ../../04_output/tagseq_05_ddm1_EMS_mutants_10M/
sbatch -p m 02.0_DESeq2.sbatch "../../04_output/tagseq_05_ddm1_EMS_mutants_10M/02_counts/" \
	"P624L_R1" \
	"none" \
	"../../03_sample_lists/sample_list_tagseq_05_kanno_ddm1_EMS.tsv" \
	"WT"

# tagseq 06 ddm1 alleles
sbatch -p c 01.0_post_processing.sbatch ../../04_output/tagseq_06_ddm1_alleles/
sbatch -p c 02.0_DESeq2.sbatch "../../04_output/tagseq_06_ddm1_alleles/02_counts/" \
	"none" \
	"none" \
	"../../03_sample_lists/tagseq_06_ddm1_alleles.tsv" \
	"WT"

# tagseq 07 remodelers
sbatch -p c 01.0_post_processing.sbatch ../../04_output/tagseq_07_remodelers/
sbatch -p c 02.0_DESeq2.sbatch "../../04_output/tagseq_07_remodelers/02_counts/" \
	"none" \
	"none" \
	"../../03_sample_lists/tagseq_07_remodelers.tsv" \
	"Col"

# tagseq 07 remodelers 0.5M
sbatch -p c 01.0_post_processing.sbatch ../../04_output/tagseq_07_remodelers_0.5M/
sbatch -p c 02.0_DESeq2.sbatch "../../04_output/tagseq_07_remodelers_0.5M/02_counts/" \
	"ddm1_chr11_17_R1,mom1_chr11_17_R1,mom1_R3" \
	"none" \
	"../../03_sample_lists/tagseq_07_remodelers.tsv" \
	"Col"

# tagseq 05 CMT3 CD and ddm1 alleles
target_dir=tagseq_05_CMT3_CD_and_ddm1_alleles
target_dir=tagseq_05_CMT3_CD_and_ddm1_demultiplexed_w_index_2_only

sbatch -p c 01.0_post_processing.sbatch ../../04_output/${target_dir}/

sample_list="../../03_sample_lists/sample_list_tagseq_05_standard.tsv"
sample_list="../../03_sample_lists/sample_list_tagseq_05_index_2_only.tsv"

sbatch -p c 02.0_DESeq2.sbatch "../../04_output/${target_dir}/02_counts/" \
	"none" \
	"none" \
	"${sample_list}" \
	"Col_0_A"

# subset only CMT3 samples
sbatch -p c 02.0_DESeq2.sbatch "../../04_output/${target_dir}/02_counts/" \
	"none" \
	"bmi,Col,ddm1,h2az,ring,ros" \
	"${sample_list}" \
	"WT"

# subset only ddm1 prc1 samples
sbatch -p c 02.0_DESeq2.sbatch "../../04_output/${target_dir}/02_counts/" \
	"ddm1_2_R1,ddm1_2_R2,ddm1_2_R3" \
	"WT,cmt3,Col_0_B,ros1,h2az" \
	"${sample_list}" \
	"Col_0_A"

# subset only ddm1 h2az ros1 samples
sbatch -p c 02.0_DESeq2.sbatch "../../04_output/${target_dir}/02_counts/" \
	"none" \
	"WT,cmt3,Col_0_A,G1,G2,ring,bmi1a,ros1_3" \
	"${sample_list}" \
	"Col_0_B"