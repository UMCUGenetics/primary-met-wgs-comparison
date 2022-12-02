############### Create driver catalog table ###############
# author: remy (sascha)
# date: 19/07/2021
# last updated: 22/07/2021

### Description
# This script reads the individual driver catalogs, add the sample name a separate column, 
# then writes this table to processed/linx_driver/2_combine_driver_catalogs/results

### Input
# LINX driver catalog file (LINX output)

### Output
# the same catalog but with a sample_id column added

# libraries
source(paste0(here::here(), '/code/r_objects/libs.R'))

#========= Path prefixes =========#
base_dir <- list(
  path=paste0(here::here(), '/data')
)

for(i in base_dir){
  if(dir.exists(i)){
    base_dir <- i
    break
  }
}

# pass input arguments to vector
args = commandArgs(trailingOnly = TRUE)

sample <- args[1]
driver_path_input <- args[2]
output_dir <- args[3]

# read linx driver catalog into R
driver_catalog <- read_tsv(driver_path_input,
                           skip = 1, # skip faulty column names
                           col_names = c('chromosome', 'chromosomeBand',
                                         'gene', 'driver', 'category', 
                                         'likelihoodMethod', 'driverLikelihood',
                                         'dndsLikelihood', 'missense',
                                         'nonsense', 'splice', 'inframe',
                                         'frameshift', 'biallelic',
                                         'minCopyNumber', 'maxCopyNumber'))

# add sample_id as a separate column for distinction between samples
driver_catalog <- driver_catalog %>%
  mutate(sample_id = sample, .before = chromosome)

# write this sample table to disk
write.table(driver_catalog,
            file = paste0(output_dir, sample, "_driver_catalog.temp"),
            sep = "\t", quote = FALSE, append = FALSE, row.names = FALSE,
            col.names = FALSE)