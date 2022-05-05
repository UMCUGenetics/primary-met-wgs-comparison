############### Investigate LINX driver output ###############
# author: remy (sascha)
# date: 21/07/2021
# last updated: 09/12/2021

### Description
# This script creates all the plots...

# libraries
library(tidyverse) # data manipulation and plotting
library(ggpubr) # wilcoxon test for plots
library(naniar) # replace with NA function

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

## ggplot custom theme
theme_bigfont <- theme(plot.title = element_text(size=22),
                       plot.subtitle = element_text(size = 16),
                       axis.text.x = element_text(size=15),
                       axis.text.y = element_text(size=15), 
                       axis.title = element_text(size=18),
                       legend.text = element_text(size = 14),
                       strip.text.x = element_text(size = 15))

# define dates of the input files
driver_results_date <- '2021_12_03'
fusion_results_date <- '2021_12_03'
exclude_hypermutators <- FALSE

# read in metadata
metadata <- read_tsv(paste0(base_dir, '/path/to/metadata.tsv')) %>%
  # select relevant columns
  select(sample_id, cohort, is_hypermutated, cancer_type, msi_status)

# ----------------------------- FUSIONS
linx_drivers_high <- read_tsv(paste0(base_dir, '/processed/linx_driver/2_combine_driver_catalogs/results/', driver_results_date ,'/linx_drivers.tsv'),
                         col_names = TRUE)

# ----------------------------- FUSIONS
linx_fusions_high <- read_tsv(paste0(base_dir, '/processed/linx_driver/2_combine_driver_catalogs/results/', fusion_results_date ,'/linx_fusions.tsv'),
                         col_names = TRUE) %>%
  filter(reported == TRUE)

# left join the metadata table, add driver column, rename 'name' column to 'gene'
linx_fusions_high <- linx_fusions_high %>%
  left_join(metadata, by = 'sample_id') %>%
  mutate(driver = 'FUSION') %>%
  rename(gene = 'name')

linx_fusions_high %>%
  group_by(cohort, gene) %>%
  count() %>%
  ungroup() %>%
  filter(cohort == 'PCAWG') %>%
  mutate(prop = n / sum(n) * 100) %>%
  arrange(desc(n))

# --------------------------------------------- DRIVER CLONALITY

driver_clonality_date <- '2021_12_03' 

# read in driver clonality
driver_clonality <- read_tsv(paste0(base_dir, '/path/to/', driver_clonality_date, '/driver_clonality_clean.tsv')) %>%
  distinct()

# left join the linx_drivers_hig table to the metadata table to include all samples
linx_drivers_high <- metadata %>%
  left_join(linx_drivers_high, by = 'sample_id')

# add the fusions to the drivers
# get the names of columns that are the same in linx drivers and fusions table
common_column_names <- names(linx_drivers_high)[names(linx_drivers_high) %in% names(linx_fusions_high)]

# filter linx drivers table based on the character vector above
linx_drivers_high_plot <- linx_drivers_high %>%
  select(common_column_names)
linx_fusions_high_plot <- linx_fusions_high %>%
  select(common_column_names)

# combine fusions and the other drivers
linx_drivers_high_plot <- linx_drivers_high_plot %>%
  union(., linx_fusions_high_plot)

# filter hypermutators?
if (exclude_hypermutators) {
  linx_drivers_high_plot <-  linx_drivers_high_plot %>%
    filter(is_hypermutated == FALSE)
}

# join the clonality to the linx_drivers_high_plot table
linx_drivers_high_plot <- linx_drivers_high_plot %>%
  left_join(driver_clonality, by = c('sample_id', 'gene', 'driver'))

### how many samples don't have any driver?
tt <- linx_drivers_high_plot %>%
  group_by(sample_id) %>%
  count(driver) %>%
  filter(is.na(driver))

nrow(tt) / n_distinct(linx_drivers_high$sample_id)

tt %>%
  inner_join(metadata, by = 'sample_id') %>%
  group_by(cancer_type) %>%
  count()

tt <- metadata %>%
  inner_join(linx_drivers_high, by = c('sample_id', 'cohort', 'cancer_type'))

# check how many ERBB@ & RSF1 mUT and AMP we have
linx_drivers_high %>%
  filter(gene == 'RSF1') %>%
  filter(cancer_type_code == 'BRCA') %>%
  group_by(cancer_type, gene, driver, likelihoodMethod) %>%
  count()
