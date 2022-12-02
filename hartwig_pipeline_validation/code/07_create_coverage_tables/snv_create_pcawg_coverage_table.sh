#!/bin/bash

#SBATCH --time=01:00:00

### Usage
# sbatch snv_create_pcawg_coverage_table.sh.sh

# get current date as variable
current_date=$(date +"%Y_%m_%d")

# define dirs
base_dir=$(dirname $(dirname $PWD))

manifest_path=$base_dir/data/processed/metadata/manifest/bam_metrics_paths.txt.gz
job_dir=$base_dir/results/07_create_coverage_tables/${current_date}/scripts/jobs/; mkdir -p $job_dir; cd $job_dir

output_dir=$base_dir/results/07_create_coverage_tables/${current_date}/results; mkdir -p $output_dir

# path to Rscript
R_SCRIPT=$base_dir/code/07_create_coverage_tables/snv_create_pcawg_coverage_table.R

counter=0
zcat $manifest_path | tail -n +2 | while read sample snv_pcawg_wgsmetrics; do
	let counter=counter+1

	job_file=$job_dir/cov_${sample}.job

	echo "--------------------------------------------------------------"

	echo "now analysing $snv_pcawg_wgsmetrics ..."

	# extract PCAWG mean coverage from bam metrics
	grep -A 1 "MEAN_COVERAGE" $snv_pcawg_wgsmetrics | \
	cut -f2 > $output_dir/${sample}_pcawg_coverage.txt

	snv_pcawg_mean_coverage=$output_dir/${sample}_pcawg_coverage.txt

	echo "#!/bin/bash
	#SBATCH --job-name=$(basename $job_file)
	#SBATCH --output=$(basename $job_file).o
	#SBATCH --time=00:20:00
	#SBATCH --mem=10G

	Rscript $R_SCRIPT $sample $snv_pcawg_mean_coverage $output_dir && touch ${job_file}.done

	rm $snv_pcawg_mean_coverage
	" > $job_file

	if [[ ! -f ${job_file}.done ]]; then
		sbatch --job-name=create_pcawg_coverage_table $job_file
	fi

	sleep 0.5s

	# if [[ $counter -eq 2 ]]; then break; fi
done

### combine all the separate table into one big table
cd $output_dir

# create the file + header for the combined SV length cutoff
echo -e "patient_id\tcohort\tmean_coverage\tcoverage_label" \
> $output_dir/pcawg_coverage_info.tsv

# create combine_tables.sh, a Bash script that combines all the short *.temp tables into a big table
echo '#!/bin/bash
#SBATCH --job-name=pcawg_coverage_combine_tables
#SBATCH --output=pcawg_coverage_combine_tables.o
#SBATCH --time=02:00:00
#SBATCH --mem=8G

echo "combining tables..."

for file in *.temp; do
	cat $file >> pcawg_coverage_info.tsv
	rm $file
done

echo "...done"
' > pcawg_coverage_combine_tables.sh

sbatch --dependency=singleton --job-name=create_pcawg_coverage_table pcawg_coverage_combine_tables.sh
