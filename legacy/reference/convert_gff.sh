#!/bin/bash

# Input and output file paths
input_file="../TAIR10_GFF3_ATTEs.gff"
output_file="converted_TAIR10_GFF3_ATTEs.gff"

# Create or clear the output file and add GFF version header
echo "##gff-version 3" > $output_file

# Read each line from the input file
while IFS=$'\t' read -r chr source type start end score strand phase attributes; do
    # Extract the ID from attributes for the Parent field
    id=$(echo "$attributes" | grep -oP 'ID=[^;]+' | cut -d'=' -f2)

    # Remove Name and Alias attributes
    attributes=$(echo "$attributes" | sed -E 's/;?Name=[^;]+//g' | sed -E 's/;?Alias=[^;]+//g')

    # Generate gene attributes with transcript_id
    gene_attributes="$attributes;gene_id=$id;transcript_id=${id}.1"
    
    # Generate transcript attributes with modified ID
    transcript_id="${id}.1"
    transcript_attributes=$(echo "$attributes" | sed -E "s/ID=[^;]+/ID=${transcript_id}/")
    transcript_attributes="$transcript_attributes;Parent=$id;gene_id=$id;transcript_id=${transcript_id}"
    
    # Generate exon attributes with modified ID
    exon_id="exon-${id}"
    exon_attributes=$(echo "$attributes" | sed -E "s/ID=[^;]+/ID=${exon_id}/")
    exon_attributes="$exon_attributes;Parent=${transcript_id};gene_id=$id;transcript_id=${transcript_id}"

    # Create the new rows
    gene_row="$chr\t$source\tgene\t$start\t$end\t$score\t$strand\t$phase\t$gene_attributes"
    transcript_row="$chr\t$source\ttranscript\t$start\t$end\t$score\t$strand\t$phase\t$transcript_attributes"
    exon_row="$chr\t$source\texon\t$start\t$end\t$score\t$strand\t$phase\t$exon_attributes"

    # Append the new rows to the output file
    echo -e "$gene_row" >> $output_file
    echo -e "$transcript_row" >> $output_file
    echo -e "$exon_row" >> $output_file

done < "$input_file"

echo "Conversion complete. Output written to $output_file"

