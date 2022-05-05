############### Filter the driver catalogs ###############
# author: remy (sascha)
# date: 15/10/2021
# last updated: 21/12/2021

### Description
# This script filters the combined LINX drivers catalog for MUTATION drivers that were also identified by the dndscv
# method, dropping all the other MUTATION drivers that were identified by the HMF method

# libraries
library(tidyverse) # data manipulation and plotting

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

# input variables
driver_results_date <- '2021_12_03'
filter_resistance_drivers <- FALSE

# ------------------------------ METADATA
metadata <- read_tsv(paste0(base_dir, '/path/to/metadata.tsv')) %>%
  # select relevant columns
  select(sample_id, cohort, cancer_type, cancer_type_code)

# ------------------------------ LINX DRIVERS
linx_drivers <- read_tsv(paste0(base_dir, '/path/to/', driver_results_date ,'/linx_drivers_prefilter.tsv'),
                         col_names = TRUE)

# ----------------------------- TREATMENT RESISTANCE DRIVERS
resistance_drivers <- read_tsv(file = paste0(base_dir, '/path/to/resistance_drivers_list.tsv'))

# little bit of string cleaning to match the driver column in 'linx_drivers_high' table
resistance_drivers <- resistance_drivers %>%
  mutate(driver = if_else(type_alt == 'coding mutations', 'mutation', type_alt)) %>%
  mutate(driver = str_to_upper(driver))

# filter out drivers with a likelihood > 0.5
linx_drivers_high <- linx_drivers %>%
  filter(driverLikelihood > 0.5)

# right join the metadata table to keep only the relevant samples that are not blacklisted
linx_drivers_high <- linx_drivers_high %>%
  right_join(metadata, by = 'sample_id')

# collapse factors: HOM_DISRUPTION is classified as DEL, PARTIAL_AMP is classified as AMP.
linx_drivers_high <- linx_drivers_high %>%
  mutate(driver = fct_collapse(driver,
                               DEL = c('DEL', 'HOM_DISRUPTION'),
                               AMP = c('AMP', 'PARTIAL_AMP')))

# remove the found resistance drivers from the linx_drivers table
if (filter_resistance_drivers) {
  linx_drivers_high <- linx_drivers_high %>%
    anti_join(resistance_drivers, by = c('gene' = 'gene_driver', 'cancer_type_code', 'driver')) 
}

linx_drivers_high_final <- linx_drivers_high

# write final table to disk
write_tsv(linx_drivers_high_final, 
          file = paste0(base_dir, '/output/path/', 
                        driver_results_date ,'/linx_drivers_nores.tsv'))
