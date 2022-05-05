#!/bin/bash

#SBATCH --time=02:00:00

### Usage
# sbatch combine_driver_catalogs.sh -m 2022_04_12_linx_output_paths.txt.gz

# get the linx version as commandline input to pass down to Rscript
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
base_dir=/base/path

manifest_path=$base_dir/path/to/$manifest_file
job_dir=$base_dir/path/to/job/scripts/${current_date}_scripts/jobs/; mkdir -p $job_dir; cd $job_dir

output_dir=$base_dir/path/to/output/${current_date}/; mkdir -p $output_dir

# path to Rscript
R_SCRIPT=$base_dir/path/to/Rscript/combine_driver_catalogs.R

# get header line of manifest file
manifest_column_names=$(zcat $manifest_path | head -n 1)

echo "the column names are: $manifest_column_names"

# initialize counter
counter=0

# colnames: sample linx_driver_catalog linx_fusion_catalog linx_breakpoints_catalog
zcat $manifest_path | tail -n +2 | while read $manifest_column_names; do
	let counter=counter+1

	echo "parsing $(basename $linx_driver_catalog) driver catalog..."

	job_file=$job_dir/combdriver_${sample}.job

	echo "#!/bin/bash
	#SBATCH --job-name=$(basename $job_file)
	#SBATCH --time=00:20:00
	#SBATCH --mem=10G


	Rscript $R_SCRIPT $sample $linx_driver_catalog $output_dir && touch ${job_file}.done
	" > $job_file

	if [[ ! -f ${job_file}.done ]]; then
		sbatch -o $(basename $job_file).out --job-name=combine_driver_catalogs $job_file
	fi

	echo "...done"

	sleep 0.5s

	# if [[ $counter -eq 2 ]]; then break; fi
done

# combine all the separate table into one big table
cd $output_dir

echo -e "sample_id\tchromosome\tchromosomeBand\tgene\tdriver\tcategory\tlikelihoodMethod\tdriverLikelihood\tdndsLikelihood\tmissense\tnonsense\tsplice\tinframe\tframeshift\tbiallelic\tminCopyNumber\tmaxCopyNumber" \
> linx_drivers_prefilter.tsv

echo '#!/bin/bash
#SBATCH --job-name=combine_driver_catalog_tables
#SBATCH --output=combine_driver_catalog_tables.o
#SBATCH --time=02:00:00
#SBATCH --mem=8G

echo "combining tables..."

for file in *.temp; do
	cat $file >> linx_drivers_prefilter.tsv
	rm $file
done

echo "...done"
' > combine_driver_catalog_tables.sh

sbatch --dependency=singleton --job-name=combine_driver_catalogs combine_driver_catalog_tables.sh
