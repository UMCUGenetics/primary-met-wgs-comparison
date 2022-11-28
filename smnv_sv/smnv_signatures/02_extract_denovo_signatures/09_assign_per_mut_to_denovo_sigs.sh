#!/bin/bash

parent_out_dir=/hpc/cuppen/projects/P0025_PCAWG_HMF/passengers/processed/sigs_denovo/extractions/12_fixed_seeds/
fit_metadata=$parent_out_dir/metadata/fit_metadata.txt.gz
job_dir=$parent_out_dir/sig_contrib/muts_assigned.jobs/; mkdir -p $job_dir
rscript=$parent_out_dir/scripts/09_assign_per_mut_to_denovo_sigs.R

counter=0
colnames=$(zcat $fit_metadata | head -n1)
zcat $fit_metadata | tail -n+2 | while read $colnames; do
	counter=$((counter+1))

	job_file=${job_dir}/${sample}.job
	done_file=${job_file}.done

	echo "#!/bin/bash
#SBATCH --output=${job_file}.o
#SBATCH --mem=16G
#SBATCH --time=1:00:00
guixr load-profile ~/.guix-profile/ --<<EOF
Rscript ${rscript} \
$sample \
${sample_dir}/${som_vcf_path} \
${sig_profile_dir}/${sig_path_SBS} \
${sig_profile_dir}/${sig_path_DBS} \
${sig_profile_dir}/${sig_path_ID} && touch ${done_file}
EOF
" > $job_file
	if [[ ! -f $done_file ]]; then
		sbatch $job_file
   	else
      	echo Done file exists: $(basename $done_file)
      	continue
   	fi

   	#if [[ $counter -eq 20 ]]; then break; fi
done


