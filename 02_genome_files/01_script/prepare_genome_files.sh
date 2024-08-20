#!/bin/bash

cd /groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/AtRTD3_no_TEG_w_ATTE
mkdir -p 01_script 02_input 03_output


#### Download data #########################################################

# get the AtRTD3 data
wget -P 02_input https://ics.hutton.ac.uk/atRTD/RTD3/atRTD3_29122021.fa https://ics.hutton.ac.uk/atRTD/RTD3/atRTD3_07082020.bed https://ics.hutton.ac.uk/atRTD/RTD3/atRTD3_TS_21Feb22_transfix.gtf https://ics.hutton.ac.uk/atRTD/RTD3/AtRTD3_gene_transcript.csv

# rename html files
for i in 02_input/*html ; do mv $i ${i%.html} ; done

# get the TAIR10 data
wget -P 02_input https://www.arabidopsis.org/download/file?path=Genes/TAIR10_genome_release/TAIR10_transposable_elements/TAIR10_TE.fas https://www.arabidopsis.org/download/file?path=Genes/TAIR10_genome_release/TAIR10_transposable_elements/TAIR10_Transposable_Elements.txt https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes_transposons.gff https://www.arabidopsis.org/download/file?path=Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz



#### remove TEGs from AtRTD3 #########################################################

# isolate TEGs
awk 'BEGIN{OFS="\t"} $3=="transposable_element_gene" {print $9}' 02_input/TAIR10_GFF3_genes_transposons.gff | sed 's/ID=// ; s/\;.*//' > 02_input/TEG_IDs.txt

# convert AtRTD3 fasta to tabular mode
./01_script/FastaToTbl.sh 02_input/atRTD3_29122021.fa | sort -k1,1 > atRTD3_29122021.fa.table

# remove TEGs from the table and convert back to fasta
grep -v -f 02_input/TEG_IDs.txt atRTD3_29122021.fa.table | ./01_script/TblToFasta.sh > 03_output/atRTD3_29122021_no_TEG.fa
rm atRTD3_29122021.fa.table

# remove TEGs from AtRTD3 GTF file
grep -v -f 02_input/TEG_IDs.txt 02_input/atRTD3_TS_21Feb22_transfix.gtf > 03_output/atRTD3_TS_21Feb22_transfix_no_TEG.gtf


####  merge with ATTE #########################################################

# reformat ATTE fasta
sed -i 's/|/ | /g' 02_input/TAIR10_TE.fas

# merge fastas
cat 03_output/atRTD3_29122021_no_TEG.fa 02_input/TAIR10_TE.fas > 03_output/atRTD3_29122021_no_TEG_w_ATTE.fa

# convert gtf to gff
ml build-env/2020 gffread/0.11.8-gcccore-8.3.0
gffread 03_output/atRTD3_TS_21Feb22_transfix_no_TEG.gtf -T -o 03_output/atRTD3_TS_21Feb22_transfix_no_TEG.gff

# filter TAIR10 GFF for ATTE elements
awk '$3 == "transposable_element"' 02_input/TAIR10_GFF3_genes_transposons.gff > 03_output/TAIR10_GFF3_ATTEs.gff

# merge gff files
cat 03_output/atRTD3_TS_21Feb22_transfix_no_TEG.gff 03_output/TAIR10_GFF3_ATTEs.gff > 03_output/atRTD3_TS_21Feb22_transfix_no_TEG_w_ATTE.gff


#### merge with mTurq-3xcMyc transgene #########################################################

# import manually created gff and fasta files
#cp /groups/berger/user/pierre.bourguet/genomics/Araport11/mTurq-3xcMyc_transgene/mTurq_3xcMyc_UBQ10_ter_pAlli_Venus.gff /groups/berger/user/pierre.bourguet/genomics/Araport11/mTurq-3xcMyc_transgene/mTurq_3xcMyc_UBQ10_ter_pAlli_Venus.fas /groups/berger/user/pierre.bourguet/genomics/Araport11/mTurq-3xcMyc_transgene/transgene_cDNAs.fas 02_input
# i remade them instead to improve counting of transgene reads. See tagseq 03 PPT summary.

# split the transgene fasta file into lines of 60 characters maximum, in case that is a problem later on
cd 02_input
awk 'length($0) > 60 && $0 !~ /^>/ {
    while (length($0) > 60) {
        print substr($0, 1, 60);
        $0 = substr($0, 61);
    }
    print $0;
    next
}
{ print }' transgene_cDNAs.fas > transgene_cDNAs.fas.tmp
mv transgene_cDNAs.fas.tmp transgene_cDNAs.fas 
cd ..

# genome fasta file (for bowtie)
gzip -d 02_input/TAIR10_chr_all.fas.gz
cat 02_input/TAIR10_chr_all.fas 02_input/mTurq_3xcMyc_UBQ10_ter_pAlli_Venus.fas > 03_output/TAIR10_chr_all_mTurq_transgene.fas 

# transcript fasta
cat 03_output/atRTD3_29122021_no_TEG_w_ATTE.fa 02_input/transgene_cDNAs.fas > 03_output/atRTD3_29122021_no_TEG_w_ATTE_mTurq_transgene.fa

## GFF

# create the GFF entry for the transgene: use vim to paste this into a file
mTurq_3xcMyc_UBQ10_ter_and_pAlli_Venus	TAIR10	chromosome	1	3147	.	.	.	ID=mTurq_3xcMyc_UBQ10_ter_and_pAlli_Venus;Name=mTurq_3xcMyc_UBQ10_ter_and_pAlli_Venus
mTurq_3xcMyc_UBQ10_ter_and_pAlli_Venus	TAIR10	gene	1	1505	.	+	.	ID=mTurq_3xcMyc;Note=protein_coding_gene;Name=mTurq_3xcMyc
mTurq_3xcMyc_UBQ10_ter_and_pAlli_Venus	TAIR10	mRNA	1	1505	.	+	.	ID=mTurq_3xcMyc.1;Parent=mTurq_3xcMyc;Name=mTurq_3xcMyc.1;Index=1
mTurq_3xcMyc_UBQ10_ter_and_pAlli_Venus	TAIR10	protein	1	873	.	+	.	ID=mTurq_3xcMyc.1-Protein;Name=mTurq_3xcMyc.1;Derives_from=mTurq_3xcMyc.1
mTurq_3xcMyc_UBQ10_ter_and_pAlli_Venus	TAIR10	exon	1	1505	.	+	.	Parent=mTurq_3xcMyc.1
mTurq_3xcMyc_UBQ10_ter_and_pAlli_Venus	TAIR10	CDS	1	873	.	+	0	Parent=mTurq_3xcMyc.1,mTurq_3xcMyc.1-Protein;
mTurq_3xcMyc_UBQ10_ter_and_pAlli_Venus	TAIR10	three_prime_UTR	874	1505	.	+	.	Parent=mTurq_3xcMyc.1
mTurq_3xcMyc_UBQ10_ter_and_pAlli_Venus	TAIR10	gene	1510	3147	.	+	.	ID=pAlli_Venus;Note=protein_coding_gene;Name=pAlli_Venus
mTurq_3xcMyc_UBQ10_ter_and_pAlli_Venus	TAIR10	mRNA	1510	3147	.	+	.	ID=pAlli_Venus.1;Parent=pAlli_Venus;Name=pAlli_Venus.1;Index=1
mTurq_3xcMyc_UBQ10_ter_and_pAlli_Venus	TAIR10	protein	2001	2720	.	+	.	ID=pAlli_Venus.1-Protein;Name=pAlli_Venus.1;Derives_from=pAlli_Venus.1
mTurq_3xcMyc_UBQ10_ter_and_pAlli_Venus	TAIR10	exon	1510	3147	.	+	.	Parent=pAlli_Venus.1
mTurq_3xcMyc_UBQ10_ter_and_pAlli_Venus	TAIR10	five_prime_UTR	1510	2000	.	+	.	Parent=pAlli_Venus.1
mTurq_3xcMyc_UBQ10_ter_and_pAlli_Venus	TAIR10	CDS	2001	2720	.	+	0	Parent=pAlli_Venus.1,pAlli_Venus.1-Protein;
mTurq_3xcMyc_UBQ10_ter_and_pAlli_Venus	TAIR10	three_prime_UTR	2721	3147	.	+	.	Parent=pAlli_Venus.1

# merge with AtRTD3
cat 03_output/atRTD3_TS_21Feb22_transfix_no_TEG_w_ATTE.gff 02_input/mTurq_3xcMyc_UBQ10_ter_pAlli_Venus.gff > 03_output/atRTD3_TS_21Feb22_transfix_no_TEG_w_ATTE_mTurq_transgene.gff 


#### salmon decoys #########################################################

##### concatenate AtRTD3 transcripts with chromosome and mTurq transgene fasta files for Salmon indexing
mkdir 03_output/salmon_decoys

## decoy files
# without transgene
cat 02_input/TAIR10_chr_all.fas | grep ">" | sed 's/\ .*// ; s/^>//' > 03_output/salmon_decoys/decoys.txt
# with transgene
cat 03_output/TAIR10_chr_all_mTurq_transgene.fas | grep ">" | sed 's/\ .*// ; s/^>//' > 03_output/salmon_decoys/decoys_transgene.txt

## gentrome files
# without transgene
cat 03_output/atRTD3_29122021_no_TEG_w_ATTE.fa 02_input/TAIR10_chr_all.fas > 03_output/salmon_decoys/gentrome_atRTD3_29122021_no_TEG_w_ATTE.fna
# with transgene
cat 03_output/atRTD3_29122021_no_TEG_w_ATTE_mTurq_transgene.fa 03_output/TAIR10_chr_all_mTurq_transgene.fas > 03_output/salmon_decoys/gentrome_atRTD3_29122021_no_TEG_w_ATTE_mTurq_transgene.fna


#### To have STAR mapping compatible with salmon counting: repair the GFF file using AGAT #########################################################
# this is because the AtRTD3 GFF file doesn't always match the fasta file: there is sometimes a few bp difference between two transcripts, causing salmon to crash

mkdir -p 03_output/agat_fix_gff
cd 03_output/agat_fix_gff

# use singularity to load agat: you need an interactive session with more resources than usual for this (eg 16 CPUs, 16 Gb RAM)
singularity pull docker://quay.io/biocontainers/agat:1.0.0--pl5321hdfd78af_0
singularity run agat_1.0.0--pl5321hdfd78af_0.sif

## GFF
# fix the gff file
agat_convert_sp_gxf2gxf.pl --gff ../atRTD3_TS_21Feb22_transfix_no_TEG_w_ATTE.gff -o atRTD3_TS_21Feb22_transfix_no_TEG_fixed.gff
# ATTEs get lost in the process, so I reintroduce them now

## create a proper GFF file for ATTEs for compatibility with STAR index (ATTEs are ignored by STAR otherwise)

# convert ATTE GFF into "gene", "transcript" and "exon" entries
./convert_gff.sh ../TAIR10_GFF3_ATTEs.gff

# merge fixed GFF with ATTE GFF and sort based on position
cat atRTD3_TS_21Feb22_transfix_no_TEG_fixed.gff converted_TAIR10_GFF3_ATTEs.gff | sort -k1,1 -k4,4n > atRTD3_TS_21Feb22_transfix_no_TEG_w_ATTE_fixed.gff

# add the transgene
cat atRTD3_TS_21Feb22_transfix_no_TEG_w_ATTE_fixed.gff ../../02_input/mTurq_3xcMyc_UBQ10_ter_pAlli_Venus.gff > atRTD3_TS_21Feb22_transfix_no_TEG_w_ATTE_fixed_mTurq_transgene.gff

## fasta
## regenerate a fasta file from the fixed GFF file
agat_sp_extract_sequences.pl -g atRTD3_TS_21Feb22_transfix_no_TEG_w_ATTE_fixed.gff -f ../../02_input/TAIR10_chr_all.fas --mrna -o atRTD3_TS_21Feb22_transfix_no_TEG_w_ATTE_fixed.fa

# merge with the transgene fasta
cat atRTD3_TS_21Feb22_transfix_no_TEG_w_ATTE_fixed.fa ../../02_input/transgene_cDNAs.fas > atRTD3_TS_21Feb22_transfix_no_TEG_w_ATTE_fixed_mTurq_transgene.fa

## regenerate gentrome files
# without transgene
cat atRTD3_TS_21Feb22_transfix_no_TEG_w_ATTE_fixed.fa ../../02_input/TAIR10_chr_all.fas > ../salmon_decoys/gentrome_atRTD3_29122021_no_TEG_w_ATTE_fixed.fna
# with transgene
cat atRTD3_TS_21Feb22_transfix_no_TEG_w_ATTE_fixed_mTurq_transgene.fa ../../03_output/TAIR10_chr_all_mTurq_transgene.fas > ../salmon_decoys/gentrome_atRTD3_29122021_no_TEG_w_ATTE_fixed_mTurq_transgene.fna
