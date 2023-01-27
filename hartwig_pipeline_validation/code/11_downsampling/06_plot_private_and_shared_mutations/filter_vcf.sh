#!/usr/bin/env bash

#SBATCH --output=filter_vcf_%j.out
#SBATCH --time=02:00:00
#SBATCH --mem=2G
#SBATCH --mail-user=<add_mail?>
#SBATCH --mail-type=END

### Description
# This script spawns a bunch of daughter scripts that are all submitted to the HPC scheduler.
# It reads and parses a manifest file, and passes the parsed values to the underlying job script (usually a Rscript or Python script).
# The number of jobs that are spawned is determined by number of lines in the manifest file - 1.

### Input
# This script needs a (filtered) manifest file

### Output
# A bunch of job scripts that are submitted to the HPC scheduler.

### Usage
# sbatch filter_vcf.sh -m <manifest_date>
# e.g. sbatch filter_vcf.sh -s true

# get the commandline input
while getopts ":s:" opt; do
  case ${opt} in
	s )
	sage_filtered=$OPTARG
	 ;;
    \? ) echo "Usage: cmd [-s]"
     ;;
    : ) echo "Invalid option: $OPTARG requires an argument" 1>&2
	 ;;
  esac
done

# get today's date
current_date=$(date +"%Y_%m_%d")

# define dirs
base_dir=$(dirname $(dirname $(dirname $PWD)))

if [[ $sage_filtered == true ]]; then
	job_dir=$base_dir/results/11_downsampling/06_plot_private_and_shared_mutations/${current_date}/sage_filtered/scripts/jobs

	output_dir=$base_dir/results/11_downsampling/06_plot_private_and_shared_mutations/${current_date}/sage_filtered/results
else
	job_dir=$base_dir/results/11_downsampling/06_plot_private_and_shared_mutations/${current_date}/original/scripts/jobs

	output_dir=$base_dir/results/11_downsampling/06_plot_private_and_shared_mutations/${current_date}/original/results
fi

# create output and job dirs
mkdir -p $job_dir; cd $job_dir
mkdir -p $output_dir

# define input files
manifest_path=$base_dir/data/processed/metadata/manifest/downsampled_vcf_paths.txt.gz

R_SCRIPT=${base_dir}/code/11_downsampling/06_plot_private_and_shared_mutations/filter_vcf.R

# get header line of manifest file
manifest_column_names=$(zcat $manifest_path | head -n 1)

echo "
Job started on $(date)

The column names are: $manifest_column_names

Results can be found at $output_dir
"

# init counter
counter=0

# read the manifest line by line, create one job script for one line and submit that job if it is not .done yet
# manifest colnames: sample	som_vcf	purity	driver_catalog	fusion_catalog	breakpoints
zcat $manifest_path | tail -n +2 | while read $manifest_column_names; do
	let counter=counter+1

	job_file=$job_dir/filtvcf_${sample}.job

	# choose correct input vcf and output directory
	if [[ $sage_filtered == true ]]; then
		input_vcf=$som_vcf_sage_filtered
	else
		input_vcf=$som_vcf
	fi

	cat > $job_file <<- EOM
	#!/usr/bin/env bash

	#SBATCH --output=$(basename $job_file).out
	#SBATCH --job-name=$(basename $job_file)
	#SBATCH --time=00:10:00
	#SBATCH --mem=2G
	#SBATCH --mail-user=<add_mail>
	#SBATCH --mail-type=FAIL

	zcat $input_vcf | \\
	SnpSift filter "FILTER == 'PASS'" - | \\
	SnpSift extractFields -s "," -e "NA" - CHROM POS REF ALT PURPLE_AF SUBCL | \\
	Rscript ${R_SCRIPT} ${current_date} ${sample} ${sage_filtered} ${output_dir} | \\
	ifne tee $output_dir/${sample}_var.txt

	gzip $output_dir/${sample}_var.txt && touch ${job_file}.done
	EOM

	if [[ ! -f ${job_file}.done ]]; then
		sbatch \
			--output=$(basename $job_file).out \
			--job-name=$(basename $job_file) \
			--time=00:10:00 \
			--mem=2G \
			--mail-user=<add_mail> \
			--mail-type=FAIL \
			$job_file
	else
		echo "Warning: There is a .done file in output folder for sample $sample"
	fi

	sleep 0.5s

	# if [[ $counter -eq 2 ]]; then break; fi
done
