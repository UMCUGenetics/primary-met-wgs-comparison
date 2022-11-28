#!/bin/bash

base_dir=/hpc/cuppen/projects/P0025_PCAWG_HMF/

wd=$base_dir/passengers/processed/mut_contexts/04_fixed_smnvs/
job_dir=$wd/jobs/; mkdir -p $job_dir
out_dir=$wd/matrices/; mkdir -p $out_dir
RSCRIPT=$wd/scripts/01_extract_sigs.R

manifest=$base_dir/passengers/processed/manifest/05_20220216/manifest_HMF_PCAWG.txt.gz

counter=0
colnames=$(zcat $manifest | head -n1)
zcat $manifest | tail -n+2 | while read $colnames; do
	counter=$((counter+1))

	job_file=$job_dir/xs_${sample}.job
	som_vcf=$dir/$som_vcf
	sv_vcf=$dir/$sv_vcf
	
	if [[ ! -f $som_vcf || ! -f $sv_vcf ]]; then
		echo "[$counter] Somatic SMNV SV or does not exist for: $sample"
		continue
	fi

	if [[ ! -f ${job_file}.done ]]; then

echo "#!/bin/bash
#SBATCH --time=00:45:00
#SBATCH --mem=10G
#SBATCH --output=${job_file}.o

guixr load-profile ~/.guix-profile --<<EOF
Rscript $RSCRIPT $sample $out_dir $som_vcf $sv_vcf && touch ${job_file}.done
EOF
" > $job_file

		echo -n "[$counter] Submitting: $sample; "
		sbatch $job_file
	else
    	#continue
      	echo "[$counter] Skipping: $sample"
  	fi

  	#if [[ $counter -eq 1 ]]; then break; fi
done


