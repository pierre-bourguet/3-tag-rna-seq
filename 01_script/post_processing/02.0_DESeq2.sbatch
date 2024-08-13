#!/bin/bash
#SBATCH --output=%j_DESeq2.out
#SBATCH --error=%j_DESeq2.out
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --nodes=1

echo -e "\n$(date) . . . Starting . . ."

# importing arguments
count_dir=$1            # e.g. ../../04_output/tagseq_01_cdca7_mutants_AtRTD3_ATTE_STAR_12/02_counts/
outliers=$2             # comma separated values of samples to exclude
outlier_pattern=$3      # comma separated values of samples to exclude based on a pattern
sample_list=$4          # e.g. ../../03_sample_lists/DESeq2_tagseq_01_cdca7.tsv
reference_condition=$5  # e.g. Col_0

# Print imported arguments as a control with a newline between each
echo -e "count directory: $count_dir\n"
echo -e "outliers: $outliers\n"
echo -e "outlier patterns: $outlier_pattern\n"
echo -e "sample info file: $sample_list\n"
echo -e "reference condition: $reference_condition\n"

ml build-env/f2022 r-bundle-bioconductor/3.18-foss-2023a-r-4.3.2

echo -e "\n$(date) . . . running 02.1_DESeq2_tagseq.R . . ."

Rscript 02.1_DESeq2_tagseq.R $count_dir $outliers $outlier_pattern $sample_list $reference_condition

echo -e "$(date) . . . done . . ."
