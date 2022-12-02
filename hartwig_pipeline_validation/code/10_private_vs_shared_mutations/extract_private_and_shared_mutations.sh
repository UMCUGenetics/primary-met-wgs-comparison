#!/bin/bash

#SBATCH --time=04:00:00

### Usage
# sbatch extract_private_and_shared_mutations.sh -m <manifest_file>

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
job_dir=$base_dir/results/10_private_vs_shared_mutations/${current_date}/scripts/jobs/; mkdir -p $job_dir; cd $job_dir

output_dir=$base_dir/results/10_private_vs_shared_mutations/${current_date}/results/; mkdir -p $output_dir

# path to the Rscript
R_SCRIPT=$base_dir/code/10_private_vs_shared_mutations/extract_private_and_shared_mutations.R

# get header line of manifest file
manifest_column_names=$(zcat $manifest_path | head -n 1)

echo "the column names are: $manifest_column_names"
echo "using this manifest file: $(basename $manifest_path)"

# init counter
counter=0

# read the manifest line by line, create one job script for one line and submit that job if it is not .done yet
zcat $manifest_path | tail -n +2 | while read $manifest_column_names; do
	let counter=counter+1

	job_file=$job_dir/privacy_${sample}.job

	echo "#!/bin/bash
	#SBATCH --job-name=$(basename $job_file)
	#SBATCH --time=00:20:00
	#SBATCH --mem=8G

	Rscript $R_SCRIPT $sample $som_hartwig    $som_pcawg_snv $som_pcawg_indel $output_dir && touch ${job_file}.done
	" > $job_file

	if [[ ! -f ${job_file}.done ]]; then
		sbatch -o $(basename $job_file).out --mem=8G --job-name=extract_private_and_shared_mutations $job_file
	fi

	sleep 0.5s

	# if [[ $counter -eq 2 ]]; then break; fi
done

### combine all the separate table into one big table
cd $output_dir

# create the file + header for the combined SV length cutoff
echo -e "sample_id\tprivacy\tmut_type\thartwig_clonality\tn" \
> $output_dir/variant_privacy.tsv

# create combine_tables.sh, a Bash script that combines all the short *.temp tables into a big table
echo '#!/bin/bash
#SBATCH --job-name=combine_tables
#SBATCH --output=combine_tables.o
#SBATCH --time=02:00:00
#SBATCH --mem=8G

echo "combining tables..."

for file in *.temp; do
	cat $file >> variant_privacy.tsv
	rm $file
done

echo "...done"
' > combine_tables.sh

sbatch --dependency=singleton --job-name=extract_private_and_shared_mutations combine_tables.sh
