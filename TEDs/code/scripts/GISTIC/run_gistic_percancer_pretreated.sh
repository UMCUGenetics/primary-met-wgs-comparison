#!/bin/bash

###adjust me###
SAMPLE_PATH="/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/"
job_dir="/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/jobs/"
#mkdir -p $job_dir
#mkdir -p $job_dir/jobs

sample_list=$(find "${SAMPLE_PATH}" -maxdepth 4 -iname "*_gistic_input.csv")

for csv in $sample_list
do

                BASEDIR=$(dirname "$csv")
                echo "$BASEDIR"
                FILENAME=$(basename "$csv" _gistic_input.csv)
                echo $FILENAME
                FILE_ABS_PATH="$(cd "$(dirname "$csv")"; pwd)/$(basename "$csv")"
                echo $FILE_ABS_PATH
                #FOLDERNAME=$(echo $FILENAME | cut -d"_" -f1)
                #basedir=${job_dir}/$FOLDERNAME/${FILENAME}
                #segfile= $FILE_ABS_PATH
                #refgenefile=/hpc/cuppen/projects/P0009_Resistance_Invivo/HMF_PCAWG_analysis/analysis/copyNumber/GISTIC2/ftp.broadinstitute.org/pub/GISTIC2.0/refgenefiles//hg19.UCSC.add_miR.140312.refgene.mat

                job_name=${job_dir}/run_gistic_${FILENAME}.sh
                if [[ ! -p $job_name ]]; then
                    touch $job_name
                fi
                echo "#!/bin/sh" > $job_name
                echo "cd  /hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/GISTIC/GISTIC_tool/" >> $job_name

                echo basedir=${BASEDIR}/GISTIC_${FILENAME} >> $job_name
                echo "mkdir -p \$basedir" >> $job_name

                echo segfile=$FILE_ABS_PATH >> $job_name
                echo refgenefile=/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/GISTIC/GISTIC_tool/refgenefiles/hg19.UCSC.add_miR.140312.refgene.mat >> $job_name

                echo "./gistic2 -b $basedir -seg $csv -refgene $refgenefile -genegistic 1 -gcm extreme -maxseg 4000 -broad 1 -brlen 0.98 -conf 0.95 -rx 0 -cap 3 -saveseg 0 -armpeel 1 -smallmem 0 -res 0.01 -ta 0.1 -td 0.1 -savedata 0 -savegene" 1 -qvt 0.1 >> $job_name

                sbatch -o $job_dir/output_gistic_${FILENAME}.csv -e $job_dir/error_gistic_${FILENAME}.txt --mail-user=a.vanhoeck@umcutrecht.nl --time=24:00:00 --mem=16G $job_name


done

