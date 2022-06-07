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
  hpc='/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers',
  mnt='/home/sascha/hpc_mnt/projects/P0025_PCAWG_HMF/drivers',
  umc='/home/cog/sbrunner/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers'
)

for(i in base_dir){
  if(dir.exists(i)){
    base_dir <- i
    break
  }
}

# input variables
driver_results_date <- '2021_12_03' # old: 2021_07_22, new PANEL == 60: 2021_10_25, new PANEL == 70: 2021_12_03, new PANEL == 90: 2021_12_07
filter_resistance_drivers <- FALSE

# ------------------------------ METADATA
metadata <- read_tsv(paste0(base_dir, '/processed/linx_driver/metadata/metadata.tsv')) %>%
  # select relevant columns
  select(sample_id, cohort, cancer_type, cancer_type_code)

# ------------------------------ LINX DRIVERS
linx_drivers <- read_tsv(paste0(base_dir, '/processed/linx_driver/2_combine_driver_catalogs/results/', driver_results_date ,'/linx_drivers_prefilter.tsv'),
                         col_names = TRUE)

# ------------------------------ DNDSCV DRIVERS (not necessary with PANEL == 70 anymore as of 03/12/2021)
# dndscv_drivers <- read_tsv(file = paste0(base_dir, '/processed/linx_driver/full_dataset_drivers_dndscv_intogen.tsv')) %>%
#   # remove the pancancer drivers (not informative)
#   filter(cancer_type != 'pancancer')

# ----------------------------- TREATMENT RESISTANCE DRIVERS
resistance_drivers <- read_tsv(file = paste0(base_dir, '/processed/linx_driver/6_plot_drivers_per_sample/2021_12_21_resistance_drivers_list.tsv'))

# little bit of string cleaning to match the driver column in 'linx_drivers_high' table
resistance_drivers <- resistance_drivers %>%
  mutate(driver = if_else(type_alt == 'coding mutations', 'mutation', type_alt)) %>%
  mutate(driver = str_to_upper(driver))

# filter out drivers with a likelihood > 0.5
# filter out SMAD3 from anaylsis, because this is liekly an artefact from the HOTSPOT analysis
linx_drivers_high <- linx_drivers %>%
  filter(driverLikelihood > 0.5) %>%
  filter(gene != 'SMAD3')

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

### dndscv filter (not necessary anymore as of 03/12/2021, see above)
# # split driver table into DNDS MUTATION drivers and rest
# linx_drivers_high_mut <- linx_drivers_high %>%
#   filter(driver == 'MUTATION' & likelihoodMethod == 'DNDS')
# linx_drivers_high_nonmut <- linx_drivers_high %>%
#   filter(!(driver == 'MUTATION' & likelihoodMethod == 'DNDS'))
# 
# # only keep those DNDS MUTATION drivers that are in the dndscv driver set
# linx_drivers_high_mut <- linx_drivers_high_mut %>%
#   semi_join(dndscv_drivers, by = c('cancer_type', 'gene' = 'gene_name'))
# 
# # drop duplicate rows
# linx_drivers_high_mut <- linx_drivers_high_mut[!duplicated(linx_drivers_high_mut),]
# 
# # combine all drivers into one final catalog, only keep first 17 columns (only linx driver columns)
# linx_drivers_high_final <- linx_drivers_high_mut %>%
#   union(., linx_drivers_high_nonmut) %>%
#   select(1:17)

linx_drivers_high_final <- linx_drivers_high

# write final table to disk
write_tsv(linx_drivers_high_final, 
          file = paste0(base_dir, '/processed/linx_driver/2_combine_driver_catalogs/results/', 
                        driver_results_date ,'/linx_drivers.tsv'))
