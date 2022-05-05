############### Get Fishers exact test results for driver genes ###############
# author: remy (sascha)
# date: 03/08/2021
# last updated: 29/09/2021

### Description
# This script creates a table of 2x2 Fishers exact test results that was conducted for each cancer type 
# and gene between primary and metastatic cancer sample groups.
# Note: Only the 22 cancer types with at least 10 samples in either PCAWg of HMF cohort are analysed.

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
results_date <- args[1]
exclude_hypermutators <- args[2]
output_dir <- args[3]

# set todays date as the output_date
output_date <- format(Sys.Date(), '%Y_%m_%d')

# ------------------------------- Custom functions

# function to transform the linx_fusions_high table into a 1x4 matrix that is used a input for the 
# cramers V and Fishers exact test function
get_contingency_matrix <- function(cancer_type_input, driver_input, gene_input) {
  
  # filter for certain driver only
  linx_fusions_high_mut <- linx_fusions_high %>%
    filter(cancer_type == cancer_type_input) %>%
    filter(driver == driver_input) %>%
    filter(gene == gene_input) %>%
    select(-is_metastatic, -is_hypermutated, -cancer_type, -cohort, -strata)
  
  # print(linx_fusions_high_mut)
  
  # join the filtered linx driver table to the metadata table to stratify the dataset and then count the amount
  # of smaples that fall into the 4 contingency table categories:
  # mut_false_group | wt_false_group
  # mut_true_group  | wt_true_group
  linx_fusions_high_mut_metadata <- metadata %>%
    filter(cancer_type == cancer_type_input) %>%
    left_join(linx_fusions_high_mut, by = 'sample_id') %>%
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
  
  # print(linx_fusions_high_mut_metadata)
  
  # handle implicitly missing rows by creating a dummy table which contains all the possible allele combination
  full_colnames <- c('mut_false_group', 'wt_false_group', 'mut_true_group', 'wt_true_group')
  full_tbl <- tibble(allele = full_colnames, dummy_col = 0)
  
  # full join the dummy table to complete the sumarized table
  linx_fusions_high_mut_metadata <- linx_fusions_high_mut_metadata %>%
    full_join(full_tbl, by = 'allele') %>%
    replace_na(., replace = list(n_mutated = 0)) %>%
    select(-dummy_col) %>%
    pivot_wider(names_from = allele, values_from = n_mutated)
  
  # print(linx_fusions_high_mut_metadata)
  
  # select the right order of the columns for the fisher.matrix()
  linx_fusions_high_mut_metadata <- linx_fusions_high_mut_metadata %>%
    select(mut_false_group, wt_false_group, mut_true_group, wt_true_group)
  
  linx_fusions_high_mut_metadata <- as.data.frame(linx_fusions_high_mut_metadata)
  
  rownames(linx_fusions_high_mut_metadata) <- paste0(str_replace(cancer_type_input, ' ', '_'), '_', driver_input,
                                                     '_', gene_input)
  
  return(linx_fusions_high_mut_metadata)
  
}

# read in metadata
metadata <- read_tsv(paste0(base_dir, '/path/to/metadata.tsv')) %>%
  # select relevant columns
  select(sample_id, cohort, is_metastatic, is_hypermutated, cancer_type)

# ----------------- Splitting the dataset by specified strata

# transform strata column into a logical if needed
metadata <- metadata %>%
  mutate(strata = if_else(cohort == 'HMF', TRUE, FALSE))

# ----------------

# read in linx drivers
linx_fusions <- read_tsv(paste0(base_dir, '/path/to/', results_date ,'/linx_fusions.tsv'),
                         col_names = TRUE) %>%
  filter(reported == TRUE)

linx_fusions_high <- linx_fusions

# left join the metadata table, add driver column, rename 'name' column to 'gene'
linx_fusions_high <- linx_fusions %>%
  left_join(metadata, by = 'sample_id') %>%
  mutate(driver = 'FUSION') %>%
  rename(gene = 'name')

## Rename fusions (stolen from luan)
## Rename one sided promiscuous fusions, add new column 'fusion_name'
linx_fusions_high$fusion_name <- (function(){
  v <- linx_fusions_high$gene
  
  is_promiscuous_5 <- linx_fusions_high$reportedType=='PROMISCUOUS_5'
  v[is_promiscuous_5] <- paste0(linx_fusions_high$geneStart[is_promiscuous_5],'_*')
  
  is_promiscuous_3 <- linx_fusions_high$reportedType=='PROMISCUOUS_3'
  v[is_promiscuous_3] <- paste0('*_',linx_fusions_high$geneEnd[is_promiscuous_3])
  
  return(v)
})()

## Split fusions with 2 promiscuous genes into 2 entries
if(sum(linx_fusions_high$reportedType=='PROMISCUOUS_BOTH')>0){
  linx_fusions_high <- (function(){
    df_promiscuous_both <- subset(linx_fusions_high, reportedType=='PROMISCUOUS_BOTH')
    
    rbind(
      subset(linx_fusions_high, reportedType!='PROMISCUOUS_BOTH'),
      within(df_promiscuous_both, fusion_name <- paste0(geneStart,'_*')),
      within(df_promiscuous_both, fusion_name <- paste0('*_',geneEnd))
    )
  })()
}

# drop "gene" column, rename 'fusion_name' to 'gene'
linx_fusions_high <- linx_fusions_high %>%
  select(-gene) %>%
  rename(gene = 'fusion_name')

# exclude hypermutators
if (exclude_hypermutators == 'yes') {
  linx_fusions_high <- linx_fusions_high %>%
    filter(is_hypermutated == FALSE)
}

# define iterators
all_genes <- unique(linx_fusions_high$gene)
all_cancers <- unique(metadata$cancer_type) # cancer types from metadata, because it is smaller than the number in linx_drivers_high 
all_drivers <- unique(linx_fusions_high$driver)

# initiate list
fusion_contingency_list <- list()

# calculate Fishers exact test p-value, FDR adjusted p-value (per cancer type, per driver mutation type),
# cramers V and the Odds ratio for each gene in the filtered linx_driver_high table per cancer type per driver mutation type
# then combine it all in one big table
for (cancer in all_cancers) {
  for (driver in all_drivers) {
    for (gene in all_genes) {
      print(paste('Running...', cancer, driver, gene))
      fusion_contingency_list[[cancer]][[driver]][[gene]] <- get_contingency_matrix(cancer_type_input = cancer,
                                                                               driver_input = driver, 
                                                                               gene_input = gene)
    }
    fusion_contingency_list[[cancer]][[driver]] <- do.call(rbind, fusion_contingency_list[[cancer]][[driver]])

  }

  fusion_contingency_list[[cancer]] <- do.call(rbind, fusion_contingency_list[[cancer]])
}

# combine all the contingency dfs into one
fusion_contingency_df <- do.call(rbind, fusion_contingency_list)

# add the cancer_type, driver and gene column again, parsed from the row names
fusion_contingency_df <- fusion_contingency_df %>%
  mutate(driver = rownames(fusion_contingency_df), .before = mut_false_group) %>%
  separate(col = driver, into = c('cancer_type', 'driver', 'gene'), sep = '\\.')

# changed output file name based on whether the table includes of excludes hypermutators
if (exclude_hypermutators == 'yes') {
  output_file <- paste0(output_dir, '/', output_date, '_fusion_contingency_matrix_no_hyp.tsv')
} else {
  output_file <- paste0(output_dir, '/', output_date, '_fusion_contingency_matrix_with_hyp.tsv')
}

# write table to disk
write_tsv(fusion_contingency_df, file = output_file,
          col_names = TRUE, append = FALSE)
