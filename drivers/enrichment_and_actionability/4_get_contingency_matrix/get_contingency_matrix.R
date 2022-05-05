############### Get Fishers exact test results for driver genes ###############
# author: remy (sascha)
# date: 30/07/2021
# last updated: 08/12/2021

### Description
# This script creates a table of 2x2 Fishers exact test results that was conducted for each cancer type 
# and gene between cohort sample groups.
# Note: Only the 22 cancer types with at least 10 samples in either PCAWG of HMF cohort are analysed.

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

# parse commandline input
args <- commandArgs(trailingOnly = TRUE)
driver_results_date <- args[1]
driver_clonality_date <- args[2]
exclude_hypermutators <- args[3]
filter_clonality <- args[4]
clonality <- args[5]
output_dir <- args[6]

# set todays date as the output_date
output_date <- format(Sys.Date(), '%Y_%m_%d')

# ------------------------------- Custom functions

# function to transform the linx_driver_high table into a 1x4 matrix that is used a input for the 
# cramers V and Fishers exact test function
get_contingency_matrix <- function(cancer_type_input, driver_input, gene_input) {
  
  # filter for certain driver only
  linx_drivers_high_mut <- linx_drivers_high %>%
    filter(cancer_type == cancer_type_input) %>%
    filter(driver == driver_input) %>%
    filter(gene == gene_input) %>%
    select(-is_metastatic, -is_hypermutated, -cancer_type, -cohort, -strata)
  
  # print(linx_drivers_high_mut)
  
  # join the filtered linx driver table to the metadata table to stratify the dataset and then count the amount
  # of samples that fall into the 4 contingency table categories:
  # mut_false_group | wt_false_group
  # mut_true_group  | wt_true_group
  linx_drivers_high_mut_metadata <- metadata %>%
    filter(cancer_type == cancer_type_input) %>%
    left_join(linx_drivers_high_mut, by = 'sample_id') %>%
    mutate(allele = if_else(gene == gene_input, 'mut', 'wt', missing = 'wt'), .after = cancer_type) %>%
    select(allele, strata) %>%
    group_by(strata, allele) %>%
    mutate(n_mutated = n()) %>%
    ungroup() %>%
    distinct() %>%
    tidyr::complete(allele, strata, fill = list(n_mutated = 0)) %>%
    mutate(strata = as.character(strata)) %>%
    mutate(strata = replace(strata, strata == 'TRUE', 'true_group'),
           strata = replace(strata, strata == 'FALSE', 'false_group')) %>%
    unite('allele', allele, strata, sep = '_')
  
  # print(linx_drivers_high_mut_metadata)
  
  # handle implicitly missing rows by creating a dummy table which contains all the possible allele combination
  full_colnames <- c('mut_false_group', 'wt_false_group', 'mut_true_group', 'wt_true_group')
  full_tbl <- tibble(allele = full_colnames, dummy_col = 0)
  
  # full join the dummy table to complete the sumarized table
  linx_drivers_high_mut_metadata <- linx_drivers_high_mut_metadata %>%
    full_join(full_tbl, by = 'allele') %>%
    replace_na(., replace = list(n_mutated = 0)) %>%
    select(-dummy_col) %>%
    pivot_wider(names_from = allele, values_from = n_mutated)

  # select the right order of the columns for the fisher.matrix()
  linx_drivers_high_mut_metadata <- linx_drivers_high_mut_metadata %>%
    select(mut_false_group, wt_false_group, mut_true_group, wt_true_group)

  linx_drivers_high_mut_metadata <- as.data.frame(linx_drivers_high_mut_metadata)

  rownames(linx_drivers_high_mut_metadata) <- paste0(str_replace(cancer_type_input, ' ', '_'), '_', driver_input,
                                                     '_', gene_input)

  return(linx_drivers_high_mut_metadata)

}

# ------------------------------ METADATA
metadata <- read_tsv(paste0(base_dir, '/path/to/metadata.tsv')) %>%
  # select relevant columns
  select(sample_id, cohort, is_metastatic, is_hypermutated, cancer_type, cancer_type_code)

# exclude hypermutators?
if (exclude_hypermutators == 'yes') {
  metadata <- metadata %>%
    filter(is_hypermutated == FALSE)
}

# ----------------- Splitting the dataset by specified strata

# transform strata column into a logical if needed
metadata <- metadata %>%
  mutate(strata = if_else(cohort == 'HMF', TRUE, FALSE))

# ------------------------------ GENE PANEL
# read in driver gene panel (important to know total number of genes that are investigated by LINX)
driver_gene_panel <- read_tsv(file = paste0(base_dir, '/path/to/DriverGenePanel.37.tsv'))

# ------------------------------ FILTERED LINX DRIVERS
linx_drivers <- read_tsv(paste0(base_dir, '/path/to/', 
                                driver_results_date ,'/linx_drivers.tsv'),
                         col_names = TRUE)

# right join the metadata table to only keep the cancer types with enough samples
linx_drivers_high <- linx_drivers %>%
  right_join(metadata, by = 'sample_id')

# collapse AMP & PARTIAL_AMP into AMP, DEL & HOM_DISRUPTION into DEL before generating the contingency table
linx_drivers_high <- linx_drivers_high %>%
  mutate(driver = fct_collapse(driver, 
                               AMP = c('AMP', 'PARTIAL_AMP'),
                               DEL = c('DEL', 'HOM_DISRUPTION')))

# ------------------------------ DRIVER CLONALITY
driver_clonality <- read_tsv(file = paste0(base_dir, '/path/to/', 
                                           driver_clonality_date, '/driver_clonality_clean.tsv'))

# join driver clonality to the linx_driver catalog
linx_drivers_high <- linx_drivers_high %>%
  left_join(driver_clonality, by = c('sample_id', 'driver', 'gene')) %>%
  mutate(clonality2 = if_else(is.na(clonality2), 'none', clonality2))

# filter clonal / subclonal?
if (filter_clonality == 'yes') {
  linx_drivers_high <- linx_drivers_high %>%
    filter(clonality2 == clonality | clonality2 == 'none')
}

# define iterators
all_genes <- unique(driver_gene_panel$gene)
all_cancers <- unique(metadata$cancer_type)
all_drivers <- unique(linx_drivers_high$driver)[complete.cases(unique(linx_drivers_high$driver))] # drop NAs

# initiate list
driver_contingency_list <- list()

# calculate Fishers exact test p-value, FDR adjusted p-value (per cancer type, per driver mutation type),
# cramers V and the Odds ratio for each gene in the filtered linx_driver_high table per cancer type per driver mutation type
# stratified by the column you specified under the section 'Splitting the dataset by specified strata'
# then combine it all in one big table
for (cancer in all_cancers) {
  for (driver in all_drivers) {
    for (gene in all_genes) {
      print(paste('Running...', cancer, driver, gene))
      driver_contingency_list[[cancer]][[driver]][[gene]] <- get_contingency_matrix(cancer_type_input = cancer,
                                                                                    driver_input = driver,
                                                                                    gene_input = gene)
    }
    driver_contingency_list[[cancer]][[driver]] <- do.call(rbind, driver_contingency_list[[cancer]][[driver]])

  }

  driver_contingency_list[[cancer]] <- do.call(rbind, driver_contingency_list[[cancer]])
}

# combine all the contingency dfs into one
driver_contingency_df <- do.call(rbind, driver_contingency_list)

# add the cancer_type, driver and gene column again, parsed from the row names
driver_contingency_df <- driver_contingency_df %>%
  mutate(driver = rownames(driver_contingency_df), .before = mut_false_group) %>%
  separate(col = 'driver', into = c('cancer_type', 'driver', 'gene'), sep = '\\.')

# changed output file name based on whether the table includes of excludes hypermutators
if (filter_clonality == 'yes') {
  output_file <- paste0(output_dir, '/', output_date, '_', clonality, '_driver_contingency_matrix_with_hyp.tsv')
} else {
  output_file <- paste0(output_dir, '/', output_date, '_driver_contingency_matrix_with_hyp.tsv')
}

# write table to disk
write_tsv(driver_contingency_df, file = output_file,
          col_names = TRUE, append = FALSE)
