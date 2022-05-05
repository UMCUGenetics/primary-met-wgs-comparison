############### Complete drivers per sample table ###############
# author: remy (sascha)
# date: 16/11/2021
# last updated: 15/12/2021

### Description
# This script creates the completed drivers per sample table that includes all driver types (amplification, deletion, fusion, point mutation)
# and their counts per sample. For the point mutation drivers this table also includes the amount of clonal/subclonal drivers.
# this table is used by the plot_drivers_per_sample.R script.

# libraries
library(tidyverse) # data manipulation and plotting
library(naniar) # replace with NA function
library(googlesheets4)

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
driver_results_date <- '2021_12_03'
fusion_results_date <- '2021_12_03' 

exclude_hypermutators <- FALSE
filter_resistance <- FALSE

# read in metadata
metadata <- read_tsv(paste0(base_dir, '/path/to/metadata.tsv')) %>%
  # select relevant columns
  select(sample_id, cohort, is_hypermutated, cancer_type, cancer_type_code)

# exclude hypermutators?
if (exclude_hypermutators) {
  metadata <-  metadata %>%
    filter(is_hypermutated == FALSE)
}

# ----------------------------- RESISTANCE DRIVERS

resistance_drivers <- read_tsv(paste0(base_dir, '/path/to/2021_12_21_resistance_drivers_list.tsv')) %>%
  mutate(type_alt = recode(type_alt, 'coding mutations' = 'MUTATION', 'Amp' = 'AMP', 'Del' = 'DEL'))

# ----------------------------- AMP/DEL/MUT DRIVERS

# read in linx drivers
linx_drivers <- read_tsv(paste0(base_dir, '/path/to/',
                                driver_results_date ,'/linx_drivers.tsv'),
                         col_names = TRUE)

# filter out drivers with a likelihood > 0.5
linx_drivers_high <- linx_drivers %>%
  { if (filter_resistance) filter(., clonality2 == clonality | clonality2 == 'none') else . } %>% # optionally filter resistance drivers
  filter(driverLikelihood > 0.5)

# ----------------------------- FUSIONS

# read in linx fusions
linx_fusions_high <- read_tsv(paste0(base_dir, '/path/to/',
                                     fusion_results_date ,'/linx_fusions.tsv'),
                         col_names = TRUE) %>%
  filter(reported == TRUE)

# left join the metadata table, add driver column, rename 'name' column to 'gene' (to match the 'gene' column)
linx_fusions_high <- linx_fusions_high %>%
  left_join(metadata, by = 'sample_id') %>%
  mutate(driver = 'FUSION') %>%
  rename(gene = 'name')

# add the fusions to the drivers table
# get the names of columns that are the same in linx drivers and fusions table
common_column_names <- names(linx_drivers_high)[names(linx_drivers_high) %in% names(linx_fusions_high)]

# filter linx drivers table based on the character vector above
linx_drivers_high <- linx_drivers_high %>%
  select(common_column_names)
linx_fusions_high <- linx_fusions_high %>%
  select(common_column_names)

# combine fusions and the other drivers
linx_drivers_high_plot <- linx_drivers_high %>%
  union(., linx_fusions_high)

# join the linx_drivers_hig table to the metadata table to include all samples
linx_drivers_high_plot <- metadata %>%
  left_join(linx_drivers_high_plot, by = 'sample_id')

# --------------------------------------------- read the driver clonality table in

driver_clonality_date <- '2021_12_03' # softfilter: 2021_10_26 strictfilter: 2021_10_29, panel == 70: 2021_12_03

# read in driver clonality
driver_clonality <- read_tsv(paste0(base_dir, '/path/to/',
                                    driver_clonality_date, '/driver_clonality_clean.tsv')) %>%
  distinct()

# join the clonality to the linx_drivers_high_plot table
linx_drivers_high_plot <- linx_drivers_high_plot %>%
  left_join(driver_clonality, by = c('sample_id', 'gene', 'driver'))

# filter out MUTATION drivers (because they need to be completed for clonality as well)
# then complete the table to include all the leftover driver types for every sample
# mark driver per driver type with 0 for those sample that dont have the specific driver
complete_table <- metadata %>% select(sample_id)

linx_drivers_high_plottt <- linx_drivers_high_plot %>%
  # filter(driver != "MUTATION") %>%
  right_join(complete_table) %>%
  group_by(sample_id, driver) %>%
  summarise(driver_per_clonality = n(), .groups = 'drop') %>%
  complete(., sample_id, driver, fill = list(driver_per_clonality = 0)) %>%
  drop_na(driver) %>%
  mutate(clonality2 = 'none') %>%
  filter(driver != 'MUTATION') # for some reason need to exclude MUTATION drivers at the end here (not in the beginning)

# filter for MUTATION drivers only, then calculate number of clonal and subclonal drivers per sample
linx_drivers_high_plottt2 <- linx_drivers_high_plot %>%
  filter(driver == 'MUTATION') %>%
  right_join(complete_table) %>%
  group_by(sample_id, driver, clonality2) %>%
  summarise(driver_per_clonality = n(), .groups = 'drop') %>%
  complete(., sample_id, driver, clonality2, fill = list(driver_per_clonality = 0)) %>%
  drop_na(driver, clonality2) %>%
  filter(!driver %in% c('AMP', 'DEL', 'FUSION'))

# combine both plots, now the table is COMPLETE
linx_drivers_high_plot <- union(linx_drivers_high_plottt, linx_drivers_high_plottt2)

# write this table to disk
if (filter_resistance) {
  write_tsv(linx_drivers_high_plot, file = paste0(base_dir, '/path/to/', 
                                                  driver_results_date, '/full_drivers_per_sample_table_nores.tsv'))
} else {
  write_tsv(linx_drivers_high_plot, file = paste0(base_dir, '/path/to/', 
                                                  driver_results_date, '/full_drivers_per_sample_table.tsv'))
}

################################### SUPPLEMENTARY TABLES
# refresh token for gsheets
gs4_auth(email = 'your_email@gmail.com')

# get the correct IDs for supp tables (HMF IDs for HMF samples)
hmf_ids <- read_tsv(paste0(base_dir, '/path/to/metadata.tsv')) %>%
  # select relevant columns
  select(sample_id, supp_table_id, cohort, cancer_type)

# drivers per sample and type of alteration
dps_supp <- linx_drivers_high_plot %>%
  group_by(sample_id, driver) %>%
  summarise(n = sum(driver_per_clonality), .groups = 'drop') %>%
  pivot_wider(., names_from = 'driver', values_from = 'n', names_prefix = 'n_') %>%
  inner_join(hmf_ids, by = 'sample_id') %>%
  select(supp_table_id, cancer_type, n_AMP, n_DEL, n_FUSION, n_MUTATION) %>%
  mutate(n_TOTAL = n_AMP + n_DEL + n_FUSION + n_MUTATION) %>%
  rename(sample_id = supp_table_id)
  

# push to gsheets
googlesheets4::write_sheet(dps_supp,
                           ss = 'https://docs.google.com/spreadsheets/...',
                           sheet = 'drivers_per_sample')
