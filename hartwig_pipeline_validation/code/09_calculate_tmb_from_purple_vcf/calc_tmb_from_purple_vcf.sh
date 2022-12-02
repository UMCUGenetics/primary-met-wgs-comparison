#!/bin/bash

#SBATCH --time=04:00:00

### Usage
# sbatch calc_tmb_from_purple_vcf.sh -m <manifest_file>

# get commandline input
while getopts ":m:" opt; do
  case ${opt} in
    m )
    manifest_file=$OPTARG
      ;;
    \? ) echo "Usage: cmd [-m]"
      ;;
		: ) echo "Invalid option: $OPTARG requires an argument" 1>&2
	    ;;
  esac
done

# get current date as variable
current_date=$(date +"%Y_%m_%d")

# define dirs
base_dir=$(dirname $(dirname $PWD))

manifest_path=$base_dir/data/processed/metadata/manifest/$manifest_file
job_dir=$base_dir/results/09_calculate_tmb_from_purple_vcf/${current_date}/scripts/jobs/; mkdir -p $job_dir; cd $job_dir

output_dir=$base_dir/results/09_calculate_tmb_from_purple_vcf/${current_date}/results/; mkdir -p $output_dir

# path to the Rscript
R_SCRIPT=$base_dir/code/09_calculate_tmb_from_purple_vcf/calc_tmb_from_purple_vcf.R

# get header line of manifest file
manifest_column_names=$(zcat $manifest_path | head -n 1)

echo """
the column names are: $manifest_column_names\n
using this manifest file: $(basename $manifest_path)
"""

# init counter
counter=0

# read the manifest line by line, create one job script for one line and submit that job if it is not .done yet
# sample purple_vcf clonality
zcat $manifest_path | tail -n +2 | while read $manifest_column_names; do
	let counter=counter+1

	job_file=$job_dir/tmb_${sample}.job

	echo "#!/bin/bash
	#SBATCH --job-name=$(basename $job_file)
	#SBATCH --time=00:20:00
	#SBATCH --mem=10G

	Rscript $R_SCRIPT $sample $purple_vcf $output_dir && touch ${job_file}.done
	" > $job_file

	if [[ ! -f ${job_file}.done ]]; then
		sbatch -o $(basename $job_file).out --job-name=extract_tmb_from_vcf $job_file
	fi

	sleep 0.5s

	# if [[ $counter -eq 2 ]]; then break; fi
done

### combine all the separate table into one big table
cd $output_dir

# create the file + header for the combined SV length cutoff
echo -e "sample_id\tmut_type\tclonality\tn_per_mut_type_and_clonality" \
> $output_dir/hartwig_variant_calls_clonality_annotated.tsv

# create combine_tables.sh, a Bash script that combines all the short *.temp tables into a big table
echo '#!/bin/bash
#SBATCH --job-name=combine_tables
#SBATCH --output=combine_tables.o
#SBATCH --time=02:00:00
#SBATCH --mem=8G

echo "combining tables..."

for file in *.temp; do
	cat $file >> hartwig_variant_calls_clonality_annotated.tsv
	rm $file
done

echo "...done"
' > combine_tables.sh

sbatch --dependency=singleton --job-name=extract_tmb_from_vcf combine_tables.sh
