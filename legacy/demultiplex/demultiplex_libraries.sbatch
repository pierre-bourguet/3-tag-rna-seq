#!/usr/bin/env bash
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --output=%j_demultiplex.out

# $1 should be the input folder where the fastq.gz files are located
# $2 should be an experiment name
# $3 should be an output directory

outputdir="${3%/}"/"$2"

echo -e "$(date) . . . Starting . . .
input folder: $1
output: $outputdir\n"

python3 /groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/000_demultiplex_script/demultiplex_tagseq.py --input $1 --output $outputdir

wait

# compress files
cd $outputdir
ml build-env/f2022 pigz/2.7-gcccore-12.2.0
pigz *.fastq

echo -e "$(date) . . . done . . ."
