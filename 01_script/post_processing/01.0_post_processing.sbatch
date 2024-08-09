#!/bin/bash
#SBATCH --output=%j_post_processing.out
#SBATCH --error=%j_post_processing.out
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=12G
#SBATCH --ntasks=2
#SBATCH --nodes=1

# $1 should be path to output files as in ../../04_output/tagseq_01_cdca7_mutants_AtRTD3_ATTE/
input_dir=${1%/}

# Check if the folder exists
if [ ! -d "$input_dir" ]; then
  echo "Error: Directory $input_dir does not exist."
  exit 1
fi

echo -e "$(date) . . . Starting. This script should take 45-60 minutes to run with 96 samples . . . "

echo -e "$(date) . . . Processing 3' tag-seq output from $1 . . . "

# aggregate count data & collect alignment rates
wd=$PWD
./01.1_aggregate_data.sh "$input_dir"
cd $wd

# plot statistics for QC metrics
echo -e "\n$(date) . . . plot statistics for QC metrics with plot_statistics.R . . ."
ml build-env/f2022 r-bundle-bioconductor/3.19-foss-2023b-r-4.4.1
Rscript 01.2_plot_statistics.R "$input_dir"

# run multiQC
echo -e "\n$(date) . . . run multiQC . . ."
ml build-env/f2022 multiqc/1.14-foss-2021b
multiqc "$1" -o "$input_dir"/05_multiqc


echo -e "$(date) . . . done . . ."
