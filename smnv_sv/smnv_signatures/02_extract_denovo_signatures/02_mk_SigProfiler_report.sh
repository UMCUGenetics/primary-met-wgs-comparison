#!/bin/bash

parent_out_dir=/hpc/cuppen/projects/P0025_PCAWG_HMF/passengers/processed/sigs_denovo/extractions/12_fixed_seeds/
rscript=$parent_out_dir/scripts/02_mk_SigProfiler_report.R

for i in $parent_out_dir/output/*/; do
	
	job_file=$i/mk_report.job
	done_file=$i/mk_report.done

	echo "#!/bin/bash
#SBATCH --output=${job_file}.o
#SBATCH --mem=4G
#SBATCH --time=0:15:00
guixr load-profile ~/.guix-profile/ --<<EOF
Rscript ${rscript} ${i} && touch ${done_file}
EOF
" > $job_file

	if [[ ! -f $done_file ]]; then
		sbatch $job_file
   	else
      	echo Done file exists: $done_file
   	fi
done

