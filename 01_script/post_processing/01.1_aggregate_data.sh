#!/bin/bash

mkdir -p "$1"/02_counts/normalized_counts/no_filter_no_transcript_merge
cd "$1"/02_counts/samples # $1 is the nextflow output folder

############################## merge counts and RPM files into salmon_counts.tsv and salmon_RPM.tsv
echo -e "\n$(date) . . . merge counts & RPM files . . . "

# Initialize the header files with the first column name
echo -n "Geneid" > header_sense
echo -n "Geneid" > header_antisense
echo -n "Geneid" > header_star_sense
echo -n "Geneid" > header_star_antisense
echo -n "Geneid" > header_unique_star_sense
echo -n "Geneid" > header_unique_star_antisense

# Find one template sample file to extract the Geneid column
salmon_template=$(find . -path '*quant*' -not -path '*STAR_mapping*' -name 'quant.sf' | head -1)
star_template=$(find . -path '*quant*' -not -path '*/quant_*/*' -name 'quant.sf' -o -path '*/STAR_mapping*/*' -name 'quant.sf' | head -1)
unique_star_template=$(find . -path '*quant*' -not -path '*/quant_*/*' -name 'quant.sf' -o -path '*/unique_STAR_mapping*/*' -name 'quant.sf' | head -1)

# Extract the Geneid column from the salmon & star template file
tail -n+2 "$salmon_template" | cut -f1 | sed 's/\.1$//' > salmon_counts.tsv
tail -n+2 "$star_template" | cut -f1 | sed 's/\.1$//' > star_counts.tsv
tail -n+2 "$unique_star_template" | cut -f1 | sed 's/\.1$//' > unique_star_counts.tsv

# Copy Geneid column
cp salmon_counts.tsv salmon_counts_AS.tsv
cp salmon_counts.tsv salmon_RPM.tsv
cp salmon_counts.tsv salmon_RPM_AS.tsv

cp star_counts.tsv star_counts_AS.tsv
cp star_counts.tsv star_RPM.tsv
cp star_counts.tsv star_RPM_AS.tsv

cp unique_star_counts.tsv unique_star_counts_AS.tsv
cp unique_star_counts.tsv unique_star_RPM.tsv
cp unique_star_counts.tsv unique_star_RPM_AS.tsv

# Prepare associative arrays to hold sample names and file paths for sense and antisense
declare -A sample_paths_sense=()
declare -A sample_paths_antisense=()
declare -A star_sample_paths_sense=()
declare -A star_sample_paths_antisense=()
declare -A unique_star_sample_paths_sense=()
declare -A unique_star_sample_paths_antisense=()

# Collect all quant.sf files and corresponding sample names, separating those in folders starting with "STAR"
while read -r quant_file; do
    dir_name=$(dirname "$quant_file")
    sample_name=$(basename "$dir_name")

    if [[ "$sample_name" == STAR* ]]; then
        if [[ "$sample_name" == *_AS ]]; then
            star_sample_paths_antisense["$sample_name"]="$quant_file"
        else
            star_sample_paths_sense["$sample_name"]="$quant_file"
        fi
    elif [[ "$sample_name" == unique* ]]; then
        if [[ "$sample_name" == *_AS ]]; then
            unique_star_sample_paths_antisense["$sample_name"]="$quant_file"
        else
            unique_star_sample_paths_sense["$sample_name"]="$quant_file"
        fi
    else
        if [[ "$sample_name" == *_AS ]]; then
            sample_paths_antisense["$sample_name"]="$quant_file"
        else
            sample_paths_sense["$sample_name"]="$quant_file"
        fi
    fi
done < <(find . -path '*quant*' -name 'quant.sf')

# Sort sample names and process files in sorted order for sense
echo -e "\n$(date) . . . Sort sample names and process files in sorted order for sense . . . "

for sample_name in $(echo "${!sample_paths_sense[@]}" | tr ' ' '\n' | sort); do
    quant_file=${sample_paths_sense[$sample_name]}

    # Add the sample name to the sense header with a tab separator
    echo -ne "\t$sample_name" >> header_sense

    # Add the corresponding counts and RPM values to the sense summary files
    awk 'NR > 1 {print $5}' "$quant_file" | paste salmon_counts.tsv - > tmp && mv tmp salmon_counts.tsv
    awk 'NR > 1 {print $4}' "$quant_file" | paste salmon_RPM.tsv - > tmp && mv tmp salmon_RPM.tsv

done

# Sort sample names and process files in sorted order for antisense
echo -e "\n$(date) . . . Sort sample names and process files in sorted order for antisense . . . "

for sample_name in $(echo "${!sample_paths_antisense[@]}" | tr ' ' '\n' | sort); do
    quant_file=${sample_paths_antisense[$sample_name]}

    # Add the sample name to the antisense header with a tab separator
    echo -ne "\t$sample_name" >> header_antisense

    # Add the corresponding counts and RPM values to the antisense summary files
    awk 'NR > 1 {print $5}' "$quant_file" | paste salmon_counts_AS.tsv - > tmp && mv tmp salmon_counts_AS.tsv
    awk 'NR > 1 {print $4}' "$quant_file" | paste salmon_RPM_AS.tsv - > tmp && mv tmp salmon_RPM_AS.tsv

done

# Sort sample names and process files in sorted order for STAR sense
echo -e "\n$(date) . . . Sort sample names and process files in sorted order for STAR sense . . . "

for sample_name in $(echo "${!star_sample_paths_sense[@]}" | tr ' ' '\n' | sort); do
    quant_file=${star_sample_paths_sense[$sample_name]}

    # Add the sample name to the STAR sense header with a tab separator
    echo -ne "\t$sample_name" >> header_star_sense

    # Add the corresponding counts and RPM values to the STAR sense summary files
    awk 'NR > 1 {print $5}' "$quant_file" | paste star_counts.tsv - > tmp && mv tmp star_counts.tsv
    awk 'NR > 1 {print $4}' "$quant_file" | paste star_RPM.tsv - > tmp && mv tmp star_RPM.tsv

done

# Sort sample names and process files in sorted order for STAR antisense
echo -e "\n$(date) . . . Sort sample names and process files in sorted order for STAR antisense . . . "

for sample_name in $(echo "${!star_sample_paths_antisense[@]}" | tr ' ' '\n' | sort); do
    quant_file=${star_sample_paths_antisense[$sample_name]}

    # Add the sample name to the STAR antisense header with a tab separator
    echo -ne "\t$sample_name" >> header_star_antisense

    # Add the corresponding counts and RPM values to the STAR antisense summary files
    awk 'NR > 1 {print $5}' "$quant_file" | paste star_counts_AS.tsv - > tmp && mv tmp star_counts_AS.tsv
    awk 'NR > 1 {print $4}' "$quant_file" | paste star_RPM_AS.tsv - > tmp && mv tmp star_RPM_AS.tsv

done


# Sort sample names and process files in sorted order for unique STAR sense
echo -e "\n$(date) . . . Sort sample names and process files in sorted order for unique STAR sense . . . "

for sample_name in $(echo "${!unique_star_sample_paths_sense[@]}" | tr ' ' '\n' | sort); do
    quant_file=${unique_star_sample_paths_sense[$sample_name]}

    # Add the sample name to the unique STAR sense header with a tab separator
    echo -ne "\t$sample_name" >> header_unique_star_sense

    # Add the corresponding counts and RPM values to the unique STAR sense summary files
    awk 'NR > 1 {print $5}' "$quant_file" | paste unique_star_counts.tsv - > tmp && mv tmp unique_star_counts.tsv
    awk 'NR > 1 {print $4}' "$quant_file" | paste unique_star_RPM.tsv - > tmp && mv tmp unique_star_RPM.tsv

done

# Sort sample names and process files in sorted order for unique STAR antisense
echo -e "\n$(date) . . . Sort sample names and process files in sorted order for unique STAR antisense . . . "

for sample_name in $(echo "${!unique_star_sample_paths_antisense[@]}" | tr ' ' '\n' | sort); do
    quant_file=${unique_star_sample_paths_antisense[$sample_name]}

    # Add the sample name to the unique STAR antisense header with a tab separator
    echo -ne "\t$sample_name" >> header_unique_star_antisense

    # Add the corresponding counts and RPM values to the unique STAR antisense summary files
    awk 'NR > 1 {print $5}' "$quant_file" | paste unique_star_counts_AS.tsv - > tmp && mv tmp unique_star_counts_AS.tsv
    awk 'NR > 1 {print $4}' "$quant_file" | paste unique_star_RPM_AS.tsv - > tmp && mv tmp unique_star_RPM_AS.tsv

done


echo -e "\n$(date) . . . format counts and RPM files . . . "

# Add a new line at the end of the headers
echo "" >> header_sense
echo "" >> header_antisense
echo "" >> header_star_sense
echo "" >> header_star_antisense
echo "" >> header_unique_star_sense
echo "" >> header_unique_star_antisense

# remove the quant_ prefix from headers
sed -i 's/quant_//g' header_sense
sed -i 's/quant_//g' header_antisense
sed -i 's/STAR_mapping_salmon_quant_//g' header_star_sense
sed -i 's/STAR_mapping_salmon_quant_//g' header_star_antisense
sed -i 's/unique_STAR_mapping_salmon_quant_//g' header_unique_star_sense
sed -i 's/unique_STAR_mapping_salmon_quant_//g' header_unique_star_antisense

# Combine headers and data files
cat header_sense salmon_counts.tsv > tmp && mv tmp ../salmon_counts.tsv
cat header_sense salmon_RPM.tsv > tmp && mv tmp ../normalized_counts/no_filter_no_transcript_merge/salmon_RPM.tsv

cat header_antisense salmon_counts_AS.tsv > tmp && mv tmp ../salmon_counts_AS.tsv
cat header_antisense salmon_RPM_AS.tsv > tmp && mv tmp ../normalized_counts/no_filter_no_transcript_merge/salmon_RPM_AS.tsv

cat header_star_sense star_counts.tsv > tmp && mv tmp ../star_counts.tsv
cat header_star_sense star_RPM.tsv > tmp && mv tmp ../normalized_counts/no_filter_no_transcript_merge/star_RPM.tsv

cat header_star_antisense star_counts_AS.tsv > tmp && mv tmp ../star_counts_AS.tsv
cat header_star_antisense star_RPM_AS.tsv > tmp && mv tmp ../normalized_counts/no_filter_no_transcript_merge/star_RPM_AS.tsv

cat header_unique_star_sense unique_star_counts.tsv > tmp && mv tmp ../unique_star_counts.tsv
cat header_unique_star_sense unique_star_RPM.tsv > tmp && mv tmp ../normalized_counts/no_filter_no_transcript_merge/unique_star_RPM.tsv

cat header_unique_star_antisense unique_star_counts_AS.tsv > tmp && mv tmp ../unique_star_counts_AS.tsv
cat header_unique_star_antisense unique_star_RPM_AS.tsv > tmp && mv tmp ../normalized_counts/no_filter_no_transcript_merge/unique_star_RPM_AS.tsv

# Cleanup
rm header_sense header_antisense header_star_sense header_star_antisense header_unique_star_sense header_unique_star_antisense
rm *tsv

# Replace the last two genes (the two transgene cDNAs) with unique counts
# in sense
head -n -2 star_counts.tsv > temp_star_counts.tsv
tail -n 2 unique_star_counts.tsv >> temp_star_counts.tsv
mv temp_star_counts.tsv star_counts.tsv
# in antisense
head -n -2 star_counts_AS.tsv > temp_star_counts_AS.tsv
tail -n 2 unique_star_counts_AS.tsv >> temp_star_counts_AS.tsv
mv temp_star_counts_AS.tsv star_counts_AS.tsv

############################## collect salmon & STAR alignment rates

echo -e "\n$(date) . . . creating STAR mapping summary file . . .\n"

cd ../../01_QC
touch tmp_STAR tmp_header
for i in STAR_logs_and_QC/*Log.final.out; do
    basename ${i%Log.final.out} >> tmp_header
    # % mapped reads
    sed -n '10p;25p;32p' $i | cut -f2 -d "|" | sed 's/\t// ; s/%//' > tmp
    # raw number of mapped reads
    sed -n '6p;9p;24p;29p;31p' $i | cut -f2 -d "|" | sed 's/\t//' >> tmp
    paste tmp_STAR tmp > tmp2
    mv tmp2 tmp_STAR
done
sed 's/_$//' tmp_header > tmp_header2 && mv tmp_header2 tmp_header
tr "\n" "\t" < tmp_header > tmp2 && mv tmp2 tmp_header
echo "" >> tmp_header

echo -e "\n$(date) . . . creating salmon mapping summary file . . .\n"

touch tmp_salmon tmp_salmon_AS
for i in ../02_counts/samples/*/quant_*/logs/salmon_quant.log; do
    if [[ "$i" =~ quant_.*_AS/ ]]; then
        grep "Mapping rate" "$i" | cut -d "=" -f2 | sed 's/\%//' > tmp
        paste tmp_salmon_AS tmp > tmp2
        mv tmp2 tmp_salmon_AS
    else
        grep "Mapping rate" "$i" | cut -d "=" -f2 | sed 's/\%//' > tmp
        paste tmp_salmon tmp > tmp2
        mv tmp2 tmp_salmon
    fi
done
rm tmp

echo -e "\n$(date) . . . merging STAR and salmon alignment rates . . .\n"
echo -e "parameter\nSTAR:%_mapped_unique\nSTAR:%_mapped_multi\nSTAR:%_of_reads_unmapped:_too_short\nSTAR:input_reads\nSTAR:mapped_unique\nSTAR:mapped_multi\nSTAR:unmapped:_too_many_mismatches\nSTAR:unmapped:_too_short\nsalmon:%_of_reads_mapped_in_sense\nsalmon:%_of_reads_mapped_in_antisense" > tmp_row_names
cat tmp_STAR tmp_salmon tmp_salmon_AS | sed 's/\t//' | cat tmp_header - | paste tmp_row_names - | sed 's/\t$//' > tmp && mv tmp pipeline_statistics.tsv
rm tmp_*

echo -e "\n$(date) . . . aggregate_data.sh done . . .\n"