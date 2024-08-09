#!/bin/bash

#### create a bed file with protein coding genes and pseudogenes from TAIR10, for DESeq2  #########################################################

echo -e "Chr\tStart\tEnd\tStrand\tType\tGeneid" > 03_output/TAIR10_GFF_PCG_pseudogene.tsv

awk 'BEGIN {FS="\t"; OFS="\t"} ; $9 ~ /Note=protein_coding_gene/ || $9 ~ /Note=pseudogene/ {print $9}' 02_input/TAIR10_GFF3_genes_transposons.gff | awk 'BEGIN {FS=";"; OFS="\t"} ; {print $2, $1}' | sed 's/ID=// ; s/Note=//' > part2

awk 'BEGIN {FS="\t"; OFS="\t"} ; $9 ~ /Note=protein_coding_gene/ || $9 ~ /Note=pseudogene/ {print $1,$4,$5,$7}' 02_input/TAIR10_GFF3_genes_transposons.gff | paste - part2 >> 03_output/TAIR10_GFF_PCG_pseudogene.tsv

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, $6, $5, $4}' 03_output/TAIR10_GFF_PCG_pseudogene.tsv > tmp && mv tmp 03_output/TAIR10_GFF_PCG_pseudogene.tsv

rm part2

# with protein coding genes only

awk 'BEGIN {FS="\t"; OFS="\t"} ; $5=="Type" || $5=="protein_coding_gene"' 03_output/TAIR10_GFF_PCG_pseudogene.tsv > 03_output/TAIR10_GFF_PCG.tsv

#### import Araport11 gene annotations and subcellular predictions and format #########################################################

# download the data
wget https://www.arabidopsis.org/download/file?path=Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.current.gff.gz https://www.arabidopsis.org/download/file?path=Genes/Araport11_genome_release/Araport11-Subcellular_Predictions
gzip -d Araport11_GFF3_genes_transposons.current.gff.gz

## filter genes and format the gff files
# extract the Geneid
awk 'BEGIN {FS="\t"; OFS="\t"} ; $3 ~ /gene/ {print $9}' Araport11_GFF3_genes_transposons.current.gff | awk 'BEGIN {FS=";"; OFS="\t"} ; {print $1}' | sed 's/ID=//' > tmp
# extract the functional annotations and merge with geneids
awk 'BEGIN {FS="\t"; OFS="\t"} ; $3 ~ /gene/ {print $1,$4,$5,$7,$3,$9}' Araport11_GFF3_genes_transposons.current.gff | paste tmp - | sort -k1,1 > tmp2
# extract subcellular predictions
grep "\.1" Araport11-Subcellular_Predictions | sed 's/\.1//' | sort -k1,1 > tmp3

# left join on the functional annotations (there are more than subcellular predictions)
awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {a[$1]=$2; next} {print $2,$3,$4,$1,$5,$7,$8 a[$1] ? a[$1] : "NA"}' tmp3 tmp2 > Araport11_gene_annotations.tsv
rm tmp tmp2 tmp3

# add ATTEs to the file
tail -n+2 TAIR10_Transposable_Elements.txt | sort -k6,6 -k5,5 | cut -f1 | sed 's/AT/Chr/ ; s/.......$//' > tmp
tail -n+2 TAIR10_Transposable_Elements.txt | sort -k6,6 -k5,5 | sed 's/false/\-/; s/true/\+/' | awk 'BEGIN{OFS="\t"} {print $3,$4,$1,$2,$5,$6}' | paste tmp - | cat Araport11_gene_annotations.tsv - | sort -k1,1 -k2,2n > Araport11_gene_ATTE_annotations.tsv
mv Araport11_gene_annotations.tsv Araport11_gene_ATTE_annotations.tsv ../03_output/
rm tmp

# add a locus type column and reorganize the order of columns
./extract_locus_type.sh
awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$4,$5,$8,$6,$7}' output_with_locus_type_column.tsv > tmp && mv tmp Araport11_gene_ATTE_annotations.tsv
rm output_with_locus_type_column.tsv


#### find TEs intersecting PCGs #########################################################

ml build-env/f2021 bedtools/2.30.0-gcc-10.2.0

# convert PCG file to bed
tail -n+2 03_output/TAIR10_GFF_PCG.tsv > PCG_bed
# convert ATTE file to bed using awk
awk 'BEGIN {FS=";"; OFS="\t"} ; {print $1}' 03_output/TAIR10_GFF3_ATTEs.gff | sed 's/ID=//' | awk 'BEGIN{OFS="\t"} {print $1,$4,$5,$9,$3,$7}' > ATTE_bed

# bedtools intersect, filter overlaps with TE and PCG annotations in the same orientation
bedtools intersect -wao -a ATTE_bed -b PCG_bed | awk '$13!=0' | awk 'BEGIN{FS=OFS="\t"} $6==$12 {print $0}' > 03_output/TAIR10_ATTE_PCG_intersect_same_orientation.tsv
bedtools intersect -wao -a ATTE_bed -b PCG_bed | awk '$13!=0' | awk 'BEGIN{FS=OFS="\t"} $6!=$12 {print $0}' > 03_output/TAIR10_ATTE_PCG_intersect_opposite_orientation.tsv
rm PCG_bed ATTE_bed