#!/bin/bash

cd /groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/split_fastq

########## split the fastq into files of 100 M reads

sbatch split_paired_fastq.sbatch \
    /scratch-cbe/users/pierre.bourguet/tagseq_05_standard/235HWMLT3_1_R19186_20250723/demultiplexed/354741/354741_S1_R1_001.fastq.gz \
    /scratch-cbe/users/pierre.bourguet/tagseq_05_standard/235HWMLT3_1_R19186_20250723/demultiplexed/354741/354741_S1_R2_001.fastq.gz \
    /scratch-cbe/users/pierre.bourguet/tagseq_05_standard/split_fastq

sbatch split_paired_fastq.sbatch \
    /scratch-cbe/users/pierre.bourguet/tagseq_05_index_2_only/235HWMLT3_1_R19186_20250724/demultiplexed/354741/354741_S1_R1_001.fastq.gz \
    /scratch-cbe/users/pierre.bourguet/tagseq_05_index_2_only/235HWMLT3_1_R19186_20250724/demultiplexed/354741/354741_S1_R2_001.fastq.gz \
    /scratch-cbe/users/pierre.bourguet/tagseq_05_index_2_only/split_fastq

########## rename the splitted fastq files in order to have them compatible with the demultiplex script

for i in $(seq -w 00 14); do
    # Create the directory
    mkdir -p "part_$i"
    # Move the corresponding files into the directory
    mv "fastq_part_r1_$i.gz" "part_$i/"
    mv "fastq_part_r2_$i.gz" "part_$i/"
done

# Base directory containing all part_XX directories
cd /scratch-cbe/users/pierre.bourguet/tagseq_05_index_2_only
base_dir="split_fastq"

# Loop over each directory in the base directory
# Check if base_dir is set and is a directory  
if [ -z "$base_dir" ]; then  
    echo "Error: Base directory variable 'base_dir' is not set." >&2  
    exit 1  
elif [ ! -d "$base_dir" ]; then  
    echo "Error: Base directory '$base_dir' not found or is not a directory." >&2  
    exit 1  
fi  
  
for dir in ${base_dir}/part_*; do  
    # Check if the item is a directory  
    # The glob can return the pattern itself if no matches are found  
    if [ "$dir" == "${base_dir}/part_*" ]; then  
        echo "Warning: No directories matching '${base_dir}/part_*' found."  
        break # Exit the loop if no matching directories  
    fi  
  
    if [ -d "$dir" ]; then  
        echo "Processing directory: $dir"  
  
        # Rename Read1 files  
        for r1_file in ${dir}/fastq_part_r1_*.gz; do  
            # Check if the file exists (glob might return pattern if no match)  
            if [ -e "$r1_file" ]; then  
                new_r1_name="${dir}/R1.fastq_part${r1_file##*_r1_}"  
                echo "Renaming $r1_file to $new_r1_name"  
                mv "$r1_file" "$new_r1_name"  
            else  
                echo "No R1 files found in $dir matching fastq_part_r1_*.gz"  
            fi  
        done  
  
        # Rename Read2 files  
        for r2_file in ${dir}/fastq_part_r2_*.gz; do  
            # Check if the file exists (glob might return pattern if no match)  
            if [ -e "$r2_file" ]; then  
                new_r2_name="${dir}/R2.fastq_part${r2_file##*_r2_}"  
                echo "Renaming $r2_file to $new_r2_name"  
                mv "$r2_file" "$new_r2_name"  
            else  
                echo "No R2 files found in $dir matching fastq_part_r2_*.gz"  
            fi  
        done  
    fi  
done


########## run the demultiplexing on each splitted file

# Base directory where your part_XX folders are located
BASE_DIR=$base_dir
output_dir=/scratch-cbe/users/pierre.bourguet/tagseq_05_index_2_only/demultiplexed

# Loop through each part directory
for part_dir in $BASE_DIR/part_*; do
    # Extract part number to use as experiment name
    part_name=$(basename $part_dir)
    echo "Processing $part_name ..."

    # Run your demultiplexing script
    sbatch -p c /groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/000_demultiplex_script/demultiplex_libraries.sbatch $part_dir $part_name $output_dir
done

########## merge back everything together
sample=tagseq_05_standard
sample=tagseq_05_index_2_only
sbatch -p c merge_demultiplexed_fastq.sbatch /scratch-cbe/users/pierre.bourguet/${sample}/demultiplexed part /scratch-cbe/users/pierre.bourguet/${sample}/demultiplex_merged

# sandbox to prepare the sample files
/scratch-cbe/users/pierre.bourguet/tagseq_06_ddm1_alleles/demultiplex_merged/
/scratch-cbe/users/pierre.bourguet/tagseq_07_remodelers/demultiplex_merged/
/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/split_fastq
/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/03_sample_lists/samples_Kanno

for i in $(seq 1 36) ; do sample=`cut -f1 /groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/03_sample_lists/samples_Kanno/tagseq_06_ddm1_alleles | head -${i} | tail -n 1` && ls /scratch-cbe/users/pierre.bourguet/tagseq_06_ddm1_alleles/demultiplex_merged/${sample}_* >> paths ; done

for i in $(seq 1 30) ; do sample=`cut -f1 /groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/03_sample_lists/samples_Kanno/tagseq_07_remodelers | head -${i} | tail -n 1` && ls /scratch-cbe/users/pierre.bourguet/tagseq_07_remodelers/demultiplex_merged/${sample}_* >> paths ; done