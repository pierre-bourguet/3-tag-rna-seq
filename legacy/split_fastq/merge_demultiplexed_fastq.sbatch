#!/bin/bash
#SBATCH --job-name=merge_fastq_files
#SBATCH --output=%j_merge_fastq_files.out
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --ntasks=1
#SBATCH --nodes=1

echo -e "$(date) Starting..."

# Define the base directory where all the part directories are located
base_dir=$1

# directory prefix as input
dir_prefix=$2

# Directory to store the merged files
output_dir=$3
mkdir -p "${output_dir}"

# Find all unique file names across all directories that start with the specified prefix
find "${base_dir}" -type d -name "${dir_prefix}*" -print0 |
xargs -0 -I {} find {} -type f -name '*.fastq.gz' | sed 's#.*/##' | sort -u |
while read filename; do
    echo "Merging files named ${filename}"

    # Correctly using find to gather all files and cat to merge them into a single file
    find "${base_dir}" -type d -name "${dir_prefix}*" -exec find {} -type f -name "${filename}" \; |
    xargs cat > "${output_dir}/${filename}"

    echo "Merged ${filename} and saved in ${output_dir}"
done

echo -e "$(date) All files have been merged successfully."
