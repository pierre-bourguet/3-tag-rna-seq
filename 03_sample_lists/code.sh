cd /groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/sample_lists

####### this file has code to create a sample_list from the demultiplexed folder ######

### You need one table as follow (sample   well):

w_s456_H2A.W_KDPK       A1
h2aw_s456       A10
w_s456_pHTA6_HTA9       A2
w_s456_pHTA9_HTA9_HTA6Cter      A3
w_s456_H2A.W_dCter      A4
WT      A5
s456    A6
w_s456_H2A.W_dCter      A7
WT      A8
w_s456_H2A.W    A9

# saved in tagseq_02_sample_to_well for example

### python script that will look for corresponding files in the demultiplexed folder, and create a table for the tagseq pipeline

# tagseq01 cdca7
python create_sample_table.py /scratch-cbe/users/pierre.bourguet/3_tag_seq_demultiplexed_tagseq_01_merged/ tagseq_01_sample_to_well tagseq_01_cdca7

# tagseq_02 h2a.w suvh456 complementation
python create_sample_table.py /scratch-cbe/users/pierre.bourguet/3_tag_seq_demultiplexed_tagseq_02_w_s456/ tagseq_02_sample_to_well tagseq_02_w_s456
# then we rename the mystery sample to mystery_R1
# also remove the unknown sample manually

# tag_seq 03
python create_sample_table.py /scratch-cbe/users/pierre.bourguet/3_tag_seq_demultiplexed_tagseq_03_cdca7_merged/ tagseq_03_sample_to_well tagseq_03_cdca7

# tagseq 04
python create_sample_table.py /scratch-cbe/users/pierre.bourguet/3_tag_seq_demultiplexed_tagseq_04_merged/ tagseq_04_sample_to_well tagseq_04_cdca7

# to create a table for DESeq2
echo -e "run_accession,condition,sample" > DESeq2_tagseq_03_cdca7.tsv
cut -f2 sample_list_tagseq_03_cdca7.tsv | while IFS= read -r line; do     echo "$line" | sed -E 's/^(.*)_(R[0-9]+)$/\1,\2/ ; s/^/,/'; done >> DESeq2_tagseq_03_cdca7.tsv

# to create a table for DESeq2
echo -e "run_accession,condition,sample" > DESeq2_tagseq_01_cdca7.tsv
cut -f2 sample_list_tagseq_01_cdca7.tsv | while IFS= read -r line; do     echo "$line" | sed -E 's/^(.*)_(R[0-9]+)$/\1,\2/ ; s/^/,/'; done >> DESeq2_tagseq_01_cdca7.tsv

######### match tag-seq samples with individual plants: the goal here is to map back transgene expression values to the T1 plant

# tagseq 03
mkdir -p tagseq_03_replicate_to_individual_plant && cd tagseq_03_replicate_to_individual_plant
vi tagseq_03_well_to_individual_plant
sort -k3,3 tagseq_03_well_to_individual_plant > tmp && mv tmp tagseq_03_well_to_individual_plant
sed 's/.*\/// ; s/_.*fastq.gz//' ../sample_list_tagseq_03_cdca7.tsv | sort -k1,1 > well_to_tagseq_replicate
join -1 3 -2 1 tagseq_03_well_to_individual_plant well_to_tagseq_replicate | tr " " "\t" > tagseq_03_replicate_to_individual_plant.tsv

# tagseq 04
mkdir -p tagseq_04_replicate_to_individual_plant && cd tagseq_04_replicate_to_individual_plant
vi tagseq_04_well_to_individual_plant
sort -k3,3 tagseq_04_well_to_individual_plant > tmp && mv tmp tagseq_04_well_to_individual_plant
sed 's/.*\/// ; s/_.*fastq.gz//' ../sample_list_tagseq_04_cdca7.tsv | sort -k1,1 > well_to_tagseq_replicate
join -1 3 -2 1 tagseq_04_well_to_individual_plant well_to_tagseq_replicate | tr " " "\t" > tagseq_04_replicate_to_individual_plant.tsv