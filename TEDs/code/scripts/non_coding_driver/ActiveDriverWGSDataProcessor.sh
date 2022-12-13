#!/bin/bash


wd="/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis"
dt=$(date +%y%m%d_%T | sed -e "s/:/-/g")
RSCRIPT=${wd}/dna-rep-ann/final-update/ActiveDriverWGSDataProcessor.R



treatment_nr=(0 3 16 16 9 6 3 3 6 6 1 1 4 2 2 2 1 2 1 5 3 7 4 1 4 4 3 1 3)

for i in {1..28}; do #28 cancer types selected for non-coding driver analysis

for j in $(seq 1 ${treatment_nr[${i}]}); do #corresponding to the number treatment groups available for each cancer type

for k in {1..6}; do #6 genomic elements we included for non-coding driver analysis


job_dir="${wd}/slurm_out/jobs/${dt}.genestat/"; mkdir -p $job_dir
err_dir="${wd}/slurm_out/errs/${dt}.genestat/"; mkdir -p $err_dir
outjob_dir="${wd}/slurm_out/outs/${dt}.genestat/"; mkdir -p $outjob_dir


job_script=${job_dir}/"${i}.${j}.${k}".jobscript.txt

        if [[ ! -p ${job_script} ]]; then
        touch ${job_script};
        fi

echo "#!/bin/bash

#SBATCH --job-name="${i}.${j}.${k}"
#SBATCH --output=${outjob_dir}/"${i}.${j}.${k}".output.txt
#SBATCH --error=${err_dir}/"${i}.${j}.${k}".error.txt
#SBATCH --time=72:00:00
#SBATCH --mem=20G
#SBATCH --mail-user=a.movasati@uu.nl


guixr load-profile /gnu/store/l78z66k3wxhjr42cmvhyj9750ds89n07-profile --<<EOF


Rscript $RSCRIPT ${i} ${j} ${k}

EOF
" > ${job_script}

sbatch ${job_script}

done

done

done

