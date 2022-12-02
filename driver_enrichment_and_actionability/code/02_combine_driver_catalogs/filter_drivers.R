############### Filter the driver catalogs ###############
# author: remy (sascha)
# date: 15/10/2021
# last updated: 21/12/2021

### Description
# This script filters the combined LINX drivers catalog.
# Applied filters are:
# 1. Minimum driver likelihood of 0.5
# 2. exclude the gene SMAD3, because it is likely an artefact from the HOTSPOT method (see LINX driver catalog on Github)
# 3. The user can also specify wether to exclude treatment resistance drivers that are specified in a separate list.

### Input
# metadata.tsv
# combined linx drivers prefilter table from step 2.
# optional: a list of resistance drivers with three columns: gene_name, cancer_type_code, driver_type (e.g. coding mutations, Amp, Del)

### Output
# A filtered combine LINX driver catalog file of the whole cohort.

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

# input variables
driver_results_date <- '2021_12_03'# the date of your combined LINX driver catalog (date name of your output folder)
filter_resistance_drivers <- FALSE # exclude treatment resistance drivers? (specified in an external list with cancer type column and gene column that match the respective columns in the driver catalog file)

# ------------------------------ METADATA

metadata <- read_tsv(paste0(base_dir, '/data/processed/metadata/metadata.tsv')) %>%
  # select relevant columns
  select(sample_id, cohort, cancer_type, cancer_type_code)

# ------------------------------ LINX DRIVERS

linx_drivers <- read_tsv(paste0(base_dir, '/results/02_combine_driver_catalogs/', driver_results_date ,'/driver/results/linx_drivers_prefilter.tsv'),
                         col_names = TRUE)

# ----------------------------- TREATMENT RESISTANCE DRIVERS

resistance_drivers <- read_tsv(file = paste0(base_dir, '/data/processed/metadata/resistance_drivers_list.tsv'))

# little bit of string cleaning to match the driver column in 'linx_drivers_high' table
resistance_drivers <- resistance_drivers %>%
  mutate(driver = if_else(type_alt == 'coding mutations', 'mutation', type_alt)) %>%
  mutate(driver = str_to_upper(driver))

# filter out drivers with a likelihood > 0.5
# filter out SMAD3 from analysis, because this is likely an artefact from the HOTSPOT analysis
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

linx_drivers_high_final <- linx_drivers_high

# write final table to disk
if (filter_resistance_drivers) {
  write_tsv(linx_drivers_high_final, 
            file = paste0(base_dir, '/results/02_combine_driver_catalogs/', 
                          driver_results_date ,'/driver/results/linx_drivers_nores.tsv'))
} else {
  write_tsv(linx_drivers_high_final, 
            file = paste0(base_dir, '/results/02_combine_driver_catalogs/', 
                          driver_results_date ,'/driver/results/linx_drivers.tsv'))
}