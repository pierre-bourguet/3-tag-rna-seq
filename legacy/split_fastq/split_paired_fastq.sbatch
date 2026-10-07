#!/bin/bash
#SBATCH --job-name=split_paired_fastq
#SBATCH --output=%j_split_paired_fastq.out
#SBATCH --error=%j_split_paired_fastq.out
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --ntasks=2

# Load the seqtk module
ml build-env/f2021 seqtk/1.3-gcc-10.2.0

# Set input files and output directory
input_file_r1=$1
input_file_r2=$2
output_dir=$3
mkdir -p ${output_dir}

# Number of lines per file (100 million reads x 4 lines per read)
lines_per_file=$((100000000 * 4))

# Process Read1
zcat ${input_file_r1} | split -d -l ${lines_per_file} - ${output_dir}/fastq_part_r1_ &

# Process Read2
zcat ${input_file_r2} | split -d -l ${lines_per_file} - ${output_dir}/fastq_part_r2_

wait


# Compress each output part for both Read1 and Read2
ml build-env/f2022 pigz/2.7-gcccore-12.2.0 
for file in ${output_dir}/fastq_part_r1_*
do
    pigz $file
done

for file in ${output_dir}/fastq_part_r2_*
do
    pigz $file
done

echo -e "$date
Splitting completed."

