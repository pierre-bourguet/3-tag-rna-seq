#!/usr/bin/env bash

echo "Geneid" > header
file1=$(ls */quant_null_Rnull/quant.sf | head -1)
tail -n+2 $file1 | cut -f1 > tagseq_01_TPM.tsv
cp tagseq_01_TPM.tsv tagseq_01_counts.tsv
echo -e "Reads in\nReads out" > tagseq_01_umi_dedup.tsv

for j in */
	do

	echo -e "doing $j files in j loop"

	# header
	echo $j | sed 's/_quant.*$// ; s/\/quant_null_Rnull\///' | paste header - > tmp
	mv tmp header

	# salmon TPMs
	tail -n+2 ${j}quant_null_Rnull/quant.sf | cut -f 4 | paste tagseq_01_TPM.tsv - > tmp
	mv tmp tagseq_01_TPM.tsv

	# salmon counts
	tail -n+2 ${j}quant_null_Rnull/quant.sf | cut -f 5 | paste tagseq_01_counts.tsv - > tmp
	mv tmp tagseq_01_counts.tsv

	# UMI deduplication
	grep "Reads In\|Reads Out" ${j}umi_dedup.log | sed 's/.*   //' | paste tagseq_01_umi_dedup.tsv - > tmp
	mv tmp tagseq_01_umi_dedup.tsv

done

cat header tagseq_01_TPM.tsv > tmp && mv tmp tagseq_01_TPM.tsv
cat header tagseq_01_counts.tsv > tmp && mv tmp tagseq_01_counts.tsv
cat header tagseq_01_umi_dedup.tsv > tmp && mv tmp tagseq_01_umi_dedup.tsv
rm header tmp2