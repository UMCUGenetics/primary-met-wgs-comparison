#!/bin/bash


job_dir=/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/jobs/
Rscript_link=/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/parse_gistic_output_percancer.R


SAMPLE_PATH=/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/
sample_list=$(find "${SAMPLE_PATH}" -maxdepth 4 -iname "scores.gistic")

for gisticrun in $sample_list
do
	echo $gisticrun
	Gisticpath=$(dirname $gisticrun)
	echo $Gisticpath
	GisticID=$(basename $gisticrun)
	echo $GisticID
	
	drugcancertype=$(echo $Gisticpath | awk -F "GISTIC" '{print $3}')
	drugcancertype="${drugcancertype:1}"
	echo $drugcancertype
	
	job_name=${job_dir}/GIS_${drugcancertype}.job
	if [[ ! -p $job_name ]]; then
		touch $job_name
	fi
	
	echo "#!/bin/sh" > $job_name
	echo guixr load-profile ~/.guix-profile '--<<EOF' >> $job_name
	echo Rscript $Rscript_link $Gisticpath $Gisticpath "Yes" >> $job_name
	echo 'EOF' >> $job_name
	
	sbatch -o $job_dir/output_${drugcancertype}.txt -e $job_dir/error_${drugcancertype}.txt --mail-user=a.vanhoeck@umcutrecht.nl --time=02:00:00 --mem=40G $job_name
	sleep .01
done
