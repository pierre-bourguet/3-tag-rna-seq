# Sort the first file on the second column and the second file on the first column
sort -k2,2 tagseq03_plant_ID_to_well > sorted_plant_ID_to_well
sort -k1,1 tagseq03_well_to_replicate > sorted_well_to_replicate
sort -k1,1 tagseq03_plant_ID_to_WGBS_replicate > sorted_tagseq03_plant_ID_to_WGBS_replicate
head sorted_tagseq03_plant_ID_to_WGBS_replicate

# Join the files
join -1 2 -2 1 -o 1.1,0,2.2 sorted_plant_ID_to_well sorted_well_to_replicate | sort -k1,1 > plant_ID-well-tagseq_sample
head plant_ID-well-tagseq_sample
join -1 1 -2 1 plant_ID-well-tagseq_sample sorted_tagseq03_plant_ID_to_WGBS_replicate > tagseq03_plant_ID_well_tagseq_WGBS.tsv