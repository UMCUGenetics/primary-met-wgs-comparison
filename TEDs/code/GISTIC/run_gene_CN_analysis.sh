#!/bin/bash


job_dir=/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/jobs/
check_dir=/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/CN_gene_database/checklist/
mkdir -p $job_dir
cd $job_dir

sample_list=/hpc/cuppen/projects/P0009_Resistance_Invivo/HMF_PCAWG_analysis/analysis/resistance/resistance_project/scripts/CN_gene/todolist.txt
#sample_list=/hpc/cuppen/projects/P0009_Resistance_Invivo/HMF_PCAWG_analysis/analysis/resistance/resistance_project/scripts/CN_gene/list_rest.txt
out_dir=/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/CN_gene_database//output/

Rscript_link=/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/CN_gene.R



cat $sample_list | while read sample_ID
do
	echo $sample_ID
	if [ -f $check_dir/${sample_ID}_annotation.done ]; then
   		echo "$check_dir${sample_ID}_annotation.done exists - exit pipeline"
		exit 0
	else 
		echo "remove previous job file"
		rm -f ${job_dir}/mim_${sample_ID}.job
		
	fi
	
	job_name=${job_dir}/mim_${sample_ID}.job
	if [[ ! -p $job_name ]]; then
		touch $job_name
	fi
	
	echo "#!/bin/sh" > $job_name
	echo guixr load-profile ~/.guix-profile '--<<EOF' >> $job_name
	echo Rscript $Rscript_link $sample_ID $out_dir >> $job_name
	echo 'EOF' >> $job_name
	
	sbatch -o $job_dir/output/output_${sample_ID}.txt -e $job_dir/error/error_${sample_ID}.txt --mail-user=a.vanhoeck@umcutrecht.nl --time=02:00:00 --mem=5G $job_name
	sleep .01
done
