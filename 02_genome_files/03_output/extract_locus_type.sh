#!/bin/bash

# Input and output file paths
input_file="Araport11_gene_ATTE_annotations.tsv"
output_file="output_with_locus_type_column.tsv"

# Create or clear the output file
> $output_file

# Read each line from the input file
while IFS=$'\t' read -r chr start end id strand attributes loc; do
    # Default locus_type to NA
    locus_type="NA"

    # Check if the row contains a transposable element
    if [[ "$id" =~ AT[0-9]TE ]]; then
        locus_type="transposable_element"
    else
        # Extract locus_type if present in the attributes
        if [[ "$attributes" =~ locus_type=([^;]+) ]]; then
            locus_type="${BASH_REMATCH[1]}"
        fi
    fi

    # Append the new row to the output file
    echo -e "$chr\t$start\t$end\t$id\t$strand\t$attributes\t$loc\t$locus_type" >> $output_file
done < "$input_file"

echo "Extraction complete. Output written to $output_file"

