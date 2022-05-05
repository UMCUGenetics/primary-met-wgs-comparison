############### Create individual fusion catalog table ###############
# author: remy (sascha)
# date: 03/08/2021
# last updated: 03/08/2021

### Description
# This script reads the individual fusion catalogs, adds the sample name as separate column, 
# then writes this table to processed/linx_driver/2_combine_driver_catalogs/results

# libraries
library(tidyverse)

#========= Path prefixes =========#
base_dir <- list(
  hpc='/base/path',
  mnt='/base/path',
  umc='/base/path'
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
fusion_path_input <- args[2]
output_dir <- args[3]

# read linx driver catalog into R
fusion_catalog <- read_tsv(fusion_path_input,
                           skip = 1, # skip faulty column names
                           col_names = c('fivePrimeBreakendId', 'threePrimeBreakendId',
                                         'name', 'reported', 'reportedType', 
                                         'phased', 'likelihood',
                                         'chainLength', 'chainLinks',
                                         'chainTerminated', 'domainsKept', 'domainsLost',
                                         'skippedExonsUp', 'skippedExonsDown',
                                         'fusedExonUp', 'fusedExonDown',
                                         'geneStart', 'geneContextStart',
                                         'transcriptStart', 'geneEnd',
                                         'geneContextEnd', 'transcriptEnd', 'junctionCopyNumber'))

# add sample_id as a separate column for distinction between samples
fusion_catalog <- fusion_catalog %>%
  mutate(sample_id = sample, .before = fivePrimeBreakendId)

# write this sample table to disk
write.table(fusion_catalog,
            file = paste0(output_dir, sample, "_fusion_catalog.temp"),
            sep = "\t", quote = FALSE, append = FALSE, row.names = FALSE,
            col.names = FALSE)