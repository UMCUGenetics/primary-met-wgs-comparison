#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --partition=cpu
#SBATCH --output=slurm.out

base_dir=/hpc/cuppen/projects/P0025_PCAWG_HMF/passengers/processed/

manifest=$base_dir/manifest/05_20220216/manifest_HMF_PCAWG.txt.gz
out_txt=$base_dir/metadata/purity/05_20220216/purple.purity.txt.gz

counter=0
colnames=$(zcat $manifest | head -n1)
zcat $manifest | tail -n+2 | while read $colnames; do
	counter=$((counter+1))

	in_txt=$dir/$purity_tsv

	if [[ ! -f $in_txt ]]; then continue; fi

	echo -ne "Processing [$counter]: $sample\r"

	if [[ $counter -eq 1 ]]; then
		head -n1 $in_txt | 
		awk '{print "sample""\t"$0}' | 
		gzip -c > $out_txt
	fi

	tail -n+2 $in_txt | 
	awk -v sample="$sample" '{print sample"\t"$0}' | 
	gzip -c >> $out_txt
done
echo -e '\n'
