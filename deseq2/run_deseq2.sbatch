#!/bin/bash
#SBATCH --output=%j_DESeq2.out
#SBATCH --error=%j_DESeq2.out
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --nodes=1

# Usage: sbatch run_deseq2.sbatch <count_dir> <outliers> <outlier_patterns> <sample_sheet> <reference_condition> [--key=value ...]
#   e.g. sbatch run_deseq2.sbatch <run>/02_counts "WT_R6,mut_R1" "none" <run>/00_demultiplex/samplesheet.csv WT --outdir=<run>/06_DESeq2_subset
# Arguments and options: see the header of DESeq2_tagseq.R

# sbatch runs a copy of this script from the spool directory, so its own location cannot be derived from $0
deseq2_dir=/groups/berger/user/pierre.bourguet/shared/pipelines/3-tag-rna-seq/deseq2

echo -e "$(date) . . . Starting . . ."
echo -e "count directory: $1\noutliers: $2\noutlier patterns: $3\nsample sheet: $4\nreference condition: $5\noptions: ${*:6}\n"

ml build-env/f2022 r-bundle-bioconductor/3.18-foss-2023a-r-4.3.2

Rscript "${deseq2_dir}/DESeq2_tagseq.R" "$@"

echo -e "$(date) . . . done . . ."
