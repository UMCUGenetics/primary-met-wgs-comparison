#!/bin/bash

base_dir=/hpc/cuppen/projects/P0025_PCAWG_HMF/passengers/
manifest=$base_dir/processed/manifest/04_20211203/manifest_HMF_PCAWG.txt.gz
out_txt=$base_dir/analysis/sv_comparison/features/linx.vis_sv_data.merged.txt.gz

counter=0
colnames=$(zcat $manifest | head -n1)
zcat $manifest | tail -n+2 | while read $colnames; do
#gunzip -c $manifest | tac | head -n -1 | while read $colnames; do
	counter=$((counter+1))

	in_txt=$dir/$linx_vis_sv_data_tsv

	if [[ ! -f $in_txt ]]; then 
		echo -ne "Skipping [$counter]: $sample\r"
		continue
	fi

	echo -ne "Processing [$counter]: $sample\r"

	if [[ $counter -eq 1 ]]; then
		head -n1 $in_txt | 
		cut -d$'\t' -f2- |
		awk '{print "cohort""\t""sample""\t"$0}' | 
		gzip -c > $out_txt
	fi

	tail -n+2 $in_txt | 
	cut -f2- |		
	awk -v cohort="$cohort" -v sample="$sample" '{print cohort"\t"sample"\t"$0}' | 
	gzip -c >> $out_txt
done
echo -e '\n'
