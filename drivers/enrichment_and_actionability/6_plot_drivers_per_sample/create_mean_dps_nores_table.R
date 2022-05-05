############### Create summarized DPS tables ###############
# author: remy (sascha)
# date: 21/12/2021
# last updated: 03/01/2022

### Description
# This script creates a summary table of the drivers per sample (DPS) either including
# or excluding identified treatment resistance drivers in the LINX catalog.

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

# define dates of the input files
driver_results_date <- '2021_12_03' # old: 2021_07_22, new: 2021_10_25, PANEL == 70: 2021_12_03, PANEL == 90: 2021_12_07
fusion_results_date <- '2021_12_03' # old: 2021_08_03, new: 2021_10_25, PANEL == 70: 2021_12_03, PANEL == 90: 2021_12_07

# read in metadata
metadata <- read_tsv(paste0(base_dir, '/path/to/metadata.tsv')) %>%
  # select relevant columns
  select(sample_id, cohort, cancer_type)

# read in the completed table
linx_drivers_high <- read_tsv(file = paste0(base_dir, '/path/to/', 
                                            driver_results_date, '/full_drivers_per_sample_table.tsv'))
linx_drivers_high_nores <- read_tsv(file = paste0(base_dir, '/path/to/', 
                                                  driver_results_date, '/full_drivers_per_sample_table_nores.tsv'))

# add new column 'group' that shows if the table has or has no resistance drivers
linx_drivers_high <- linx_drivers_high %>%
  mutate(group = '+res')
linx_drivers_high_nores <- linx_drivers_high_nores %>%
  mutate(group = '-res')

# combine both tables
linx_drivers_high <- linx_drivers_high %>%
  union(., linx_drivers_high_nores)

# join the metadata
linx_drivers_high <- linx_drivers_high %>%
  left_join(metadata, by = 'sample_id')

# ------------------ create total drivers per sample with/without resistance drivers

total_dps <- linx_drivers_high %>%
  group_by(sample_id, cancer_type, group) %>%
  summarise(total_dps = sum(driver_per_clonality), .groups = 'drop') %>%
  mutate(group = if_else(group == '+res', 'with_res', 'no_res')) %>%
  pivot_wider(., names_from = 'group', values_from = 'total_dps')

write_tsv(total_dps, file = paste0(base_dir, '/path/to/total_dps_res.tsv'))
  
# ------------------ create mean dps table

# summarise table, first sum up all drivers then calculate the mean DPS and standard deviation for each cancer type in each cohort, 
# with and without resistance drivers
linx_drivers_high_summed <- linx_drivers_high %>%
  group_by(sample_id, cancer_type, cohort, group) %>%
  summarise(total_drivers = sum(driver_per_clonality)) %>%
  group_by(cancer_type, cohort, group) %>%
  summarise(mean_dps = mean(total_drivers),
            sd_dps = sd(total_drivers))

linx_drivers_high_summed_pancan <- linx_drivers_high %>%
  group_by(sample_id, cancer_type, cohort, group) %>%
  summarise(total_drivers = sum(driver_per_clonality)) %>%
  group_by(cohort, group) %>%
  summarise(mean_dps = mean(total_drivers),
            sd_dps = sd(total_drivers)) %>%
  mutate(cancer_type = 'pancancer')

linx_drivers_high_summed <- linx_drivers_high_summed %>%
  union(., linx_drivers_high_summed_pancan)

# write table to disk
write_tsv(linx_drivers_high_summed, file = paste0(base_dir, '/output/path/to/resistance_mean_dps.tsv'))
