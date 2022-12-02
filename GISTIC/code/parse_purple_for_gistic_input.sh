#!/bin/bash

###adjust me###
sample_list=$(cat "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/json_files_Nov.txt")
job_dir="/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/jobs/"

Rscript_link=/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/parse_purple_for_gistic_input.R


for json in $sample_list
do

                BASEDIR=$(dirname "$json")
                echo "$BASEDIR"
                CANCERTYPE=$(echo $BASEDIR | awk -F "json_files/" '{print $2}')
                echo "$CANCERTYPE"
                DRUGTYPE=$(basename "$json" .json)
                echo "$DRUGTYPE"
                mkdir -p "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/${CANCERTYPE}"
                outdir="/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/${CANCERTYPE}/"
                echo "$tempout"

                job_name=${job_dir}/GIST_${CANCERTYPE}_${DRUGTYPE}.sh
                if [[ ! -p $job_name ]]; then
                    touch $job_name
                fi
                echo "#!/bin/sh" > $job_name
                echo guixr load-profile /home/cog/avanhoeck/.guix-profile '--<<EOF' >> $job_name
                echo Rscript $Rscript_link $json $outdir/${CANCERTYPE}_${DRUGTYPE}_gistic_input.csv >> $job_name
                echo "EOF" >> $job_name

                sbatch -o $job_dir/output_${CANCERTYPE}_${DRUGTYPE}.csv -e $job_dir/error_${CANCERTYPE}_${DRUGTYPE}.txt --mail-user=a.vanhoeck@umcutrecht.nl --time=01:00:00 --mem=5G $job_name


done
