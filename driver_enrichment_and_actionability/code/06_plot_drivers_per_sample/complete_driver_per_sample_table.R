############### Complete drivers per sample table ###############
# author: remy (sascha)
# date: 16/11/2021
# last updated: 03/10/2022

### Description
# This script creates the completed drivers per sample table that includes all driver types (amplification, deletion, fusion, point mutation)
# and their counts per sample. For the point mutation drivers this table also includes the amount of clonal/subclonal drivers.
# this table is used by the plot_drivers_per_sample.R script.

# libraries
source(paste0(here::here(), '/code/r_objects/libs.R'))

#========= Path prefixes =========#
base_dir <- list(
  path=paste0(here::here())
)

for(i in base_dir){
  if(dir.exists(i)){
    base_dir <- i
    break
  }
}

# define dates of the input files
driver_results_date <- '2021_12_03' # prefix of your generated LINX driver table
fusion_results_date <- '2021_12_03' # prefix of your generated fusion driver table
driver_clonality_date <- '2021_12_03' # prefix of your generated driver clonality file
filter_resistance <- FALSE # TRUE: filter resistance drivers, FALSE: dont

# ----------------------------- METADATA

# read in metadata
metadata <- read_tsv(paste0(base_dir, '/data/processed/metadata/metadata.tsv'),
                     col_types = 'iccccccccccccccllllcccdiiildclccllllllliiiiiiiiiii') %>%
  # select relevant columns
  select(sample_id, cohort, cancer_type, cancer_type_code)

# ----------------------------- AMP/DEL/MUT DRIVERS

# read in linx drivers
if (filter_resistance) {
  linx_drivers <- read_tsv(paste0(base_dir, '/results/02_combine_driver_catalogs/', driver_results_date ,'/driver/results/linx_drivers_nores.tsv'),
                                col_names = TRUE)
} else {
  linx_drivers <- read_tsv(paste0(base_dir, '/results/02_combine_driver_catalogs/', driver_results_date ,'/driver/results/linx_drivers.tsv'),
                                col_names = TRUE)
}

# ----------------------------- FUSIONS

# read in linx fusions
linx_fusions_high <- read_tsv(paste0(base_dir, '/results/02_combine_driver_catalogs/', fusion_results_date ,'/fusion/results/linx_fusions.tsv'),
                         col_names = TRUE) %>%
  filter(reported == TRUE)

# left join the metadata table, add driver column, rename 'name' column to 'gene' (to match the 'gene' column)
linx_fusions_high <- linx_fusions_high %>%
  left_join(metadata, by = 'sample_id') %>%
  mutate(driver = 'FUSION') %>%
  rename(gene = 'name')

# add the fusions to the drivers table
# get the names of columns that are the same in linx drivers and fusions table
common_column_names <- names(linx_drivers)[names(linx_drivers) %in% names(linx_fusions_high)]

# filter linx drivers table based on the character vector above
linx_drivers <- linx_drivers %>%
  select(common_column_names)
linx_fusions_high <- linx_fusions_high %>%
  select(common_column_names)

# combine fusions and the other drivers
linx_drivers_plot <- linx_drivers %>%
  union(., linx_fusions_high)

# join the linx_drivers_hig table to the metadata table to include all samples
linx_drivers_plot <- metadata %>%
  left_join(linx_drivers_plot, by = 'sample_id')

# ----------------------------- DRIVER CLONALITY

# read in driver clonality
driver_clonality <- read_tsv(paste0(base_dir, '/results/03_parse_driver_clonality/', driver_clonality_date, '/results/driver_clonality_clean.tsv')) %>%
  distinct()

# join the clonality to the linx_drivers_plot table
linx_drivers_plot <- linx_drivers_plot %>%
  left_join(driver_clonality, by = c('sample_id', 'gene', 'driver'))

# filter out MUTATION drivers (because they need to be completed for clonality as well)
# then complete the table to include all the leftover driver types for every sample
# mark driver per driver type with 0 for those sample that dont have the specific driver
complete_table <- metadata %>% select(sample_id)

linx_drivers_plottt <- linx_drivers_plot %>%
  right_join(complete_table) %>%
  group_by(sample_id, driver) %>%
  summarise(driver_per_clonality = n(), .groups = 'drop') %>%
  complete(., sample_id, driver, fill = list(driver_per_clonality = 0)) %>%
  drop_na(driver) %>%
  mutate(clonality2 = 'none') %>%
  filter(driver != 'MUTATION') # for some reason need to exclude MUTATION drivers at the end here (not in the beginning)

# filter for MUTATION drivers only, then calculate number of clonal and subclonal drivers per sample
linx_drivers_plottt2 <- linx_drivers_plot %>%
  filter(driver == 'MUTATION') %>%
  right_join(complete_table) %>%
  group_by(sample_id, driver, clonality2) %>%
  summarise(driver_per_clonality = n(), .groups = 'drop') %>%
  complete(., sample_id, driver, clonality2, fill = list(driver_per_clonality = 0)) %>%
  drop_na(driver, clonality2) %>%
  filter(!driver %in% c('AMP', 'DEL', 'FUSION'))

# combine both plots, now the table is COMPLETE
linx_drivers_plot <- union(linx_drivers_plottt, linx_drivers_plottt2)

# write this table to disk
if (filter_resistance) {
  write_tsv(linx_drivers_plot, file = paste0(base_dir, '/results/02_combine_driver_catalogs/', driver_results_date ,'/driver/results/full_drivers_per_sample_table_nores.tsv'))
} else {
  write_tsv(linx_drivers_plot, file = paste0(base_dir, '/results/02_combine_driver_catalogs/', driver_results_date ,'/driver/results/full_drivers_per_sample_table.tsv'))
}