############### Get Fishers exact test results for driver genes ###############
# author: remy (sascha)
# date: 03/08/2021
# last updated: 21/10/2022

### Description
# This script creates a matrix of 2x2 contingency tables for each fusion gene in each cancer type
# It counts how many times a fusion driver is reported in the PCAWG and Hartwig cohorts
# and stores these values in 4 columns of a matrix. Each row represents a fusion driver in a cancer type.
# The user can specify whether to do the contingecy table generation by cancer subtype, metastatic location, or primary progression.

### Input
# metadata.tsv
# LINX fusions catalog from step 2.

### Output
# A contingency table where each row is a fusion driver in a specific cancer type.
# The first three columns are the cancer type, the gene name(s) and the driver (in this case FUSION).
# The next four column denote counts of how many PCAWG samples do have the driver or not, and how many Hartwig samples do have the driver or not.

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

# parse commandline input
args <- commandArgs(trailingOnly = TRUE)
results_date <- args[1]
by_subtype <- args[2]
by_subtype_met_location <- args[3]
by_progression <- args[4]
by_met_location <- args[5]
output_dir <- args[6]

# checks
if (sum(str_detect(c(by_subtype, by_progression, by_met_location), pattern = 'yes')) > 1) {
  stop('Cannot perform contingency table generation by multiple strata at the same time.\n
       Choose either -s, -p, -m or none of the above to be "yes" and run script again.')
}

# set todays date as the current_date
current_date <- format(Sys.Date(), '%Y_%m_%d')

# ------------------------------- Custom functions

# function to transform the linx_fusions table into a 1x4 matrix that is used a input for the 
# cramers V and Fishers exact test function
get_contingency_matrix <- function(cancer_type_input, driver_input, gene_input) {
  
  # filter for certain driver only
  linx_fusions_mut <- linx_fusions %>%
    filter(cancer_type == cancer_type_input) %>%
    filter(driver == driver_input) %>%
    filter(gene == gene_input) %>%
    select(-cancer_type, -cohort, -strata)
  
  # join the filtered linx driver table to the metadata table to stratify the dataset and then count the amount
  # of smaples that fall into the 4 contingency table categories:
  # mut_false_group | wt_false_group
  # mut_true_group  | wt_true_group
  linx_fusions_mut_metadata <- metadata %>%
    filter(cancer_type == cancer_type_input) %>%
    left_join(linx_fusions_mut, by = 'sample_id') %>%
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
  
  # handle implicitly missing rows by creating a dummy table which contains all the possible allele combination
  full_colnames <- c('mut_false_group', 'wt_false_group', 'mut_true_group', 'wt_true_group')
  full_tbl <- tibble(allele = full_colnames, dummy_col = 0)
  
  # full join the dummy table to complete the sumarized table
  linx_fusions_mut_metadata <- linx_fusions_mut_metadata %>%
    full_join(full_tbl, by = 'allele') %>%
    replace_na(., replace = list(n_mutated = 0)) %>%
    select(-dummy_col) %>%
    pivot_wider(names_from = allele, values_from = n_mutated)
  
  # select the right order of the columns for the fisher.matrix()
  linx_fusions_mut_metadata <- linx_fusions_mut_metadata %>%
    select(mut_false_group, wt_false_group, mut_true_group, wt_true_group)
  
  linx_fusions_mut_metadata <- as.data.frame(linx_fusions_mut_metadata)
  
  rownames(linx_fusions_mut_metadata) <- paste0(str_replace(cancer_type_input, ' ', '_'), ';', driver_input,
                                                     ';', gene_input)
  
  return(linx_fusions_mut_metadata)
  
}

# ----------------- METADATA

metadata <- read_tsv(paste0(base_dir, '/data/processed/metadata/metadata.tsv'),
                     col_types = 'icccccccccccccccllllcccdiiildclccllllllliiiiiiiiiii') %>%
  # select relevant columns
  select(patient_id, sample_id, cohort, cancer_type, cancer_type_code, 
         cancer_subtype, progression_status_code, is_blacklisted_subtype, metastatic_location) %>%
  # add has_subtype column for easier filtering
  mutate(has_subtype = if_else(cancer_type_code %in% c('BRCA', 'COREAD', 'UCEC'), TRUE, FALSE))

# ----------------- SUBTYPES

# do analysis by subtype
if (by_subtype == 'yes') {
  
  # duplicate the rows of samples with subtype annotation, fuse together cancer_type and cancer_subtype column
  # do analysis by subtype and metastatic location
  if (by_subtype_met_location == 'yes') {
    
    subtypes_met_loc <- metadata %>%
      filter(has_subtype & !is_blacklisted_subtype & !(metastatic_location %in% c('Unknown', 'Primary', NA_character_))) %>%
      mutate(cancer_type = paste0(cancer_type, '_', cancer_subtype, '_', metastatic_location))
    
    # duplicate PCAWG samples rows to give them a pseudo cancer_type + metastatic_location column for easier data handling
    pcawg <- metadata %>%
      filter(cohort == 'PCAWG' & has_subtype & !is_blacklisted_subtype)
    
    pcawg_local <- pcawg %>%
      mutate(cancer_type = paste0(cancer_type, '_', cancer_subtype, '_Local'))
    
    pcawg_distant <- pcawg %>%
      mutate(cancer_type = paste0(cancer_type, '_', cancer_subtype, '_Distant'))
    
    pcawg_lymph <- pcawg %>%
      mutate(cancer_type = paste0(cancer_type, '_', cancer_subtype, '_Lymph'))
    
    # combine artificial pcawg duplicate rows
    subtypes <- subtypes_met_loc %>%
      union(., pcawg_local) %>%
      union(., pcawg_distant) %>%
      union(., pcawg_lymph)
    
    # filter out those cancer subtypes with low hartwig samples in a certain metastatic location category (<5 samples)
    low_hartwig_samples <- subtypes %>%
      group_by(cancer_type, cohort) %>%
      count() %>%
      ungroup() %>%
      pivot_wider(names_from = 'cohort', values_from = 'n', values_fill = 0) %>%
      filter(pmax(Hartwig, PCAWG) < 5) %>%
      pull(cancer_type)
    
    subtypes <- subtypes %>%
      filter(!cancer_type %in% low_hartwig_samples)
    
  } else {
    subtypes <- metadata %>%
      filter(has_subtype & !is_blacklisted_subtype) %>%
      mutate(cancer_type = paste0(cancer_type, '_', cancer_subtype))
  }
  
  # add subtype rows to the metadata df
  metadata <- bind_rows(metadata, subtypes)
  
  # filter for cancer types that only have subtypes
  metadata <- metadata %>%
    filter(has_subtype)
}

# ---------------------- PRIMARY PROGRESSION

# do analysis by comparing different primary progression phenotypes to metastatic samples
if (by_progression == 'yes') {
  
  progression_cancer_types <- metadata %>%
    filter(cancer_type_code %in% c('PRAD', 'PANET')) %>%
    mutate(cancer_type = cancer_type_code)
  
  # create new contrast groups (primary progression stable vs. relapse, relapse vs. met, stable vs. met) for all cancer types with annotation
  progression <- metadata %>%
    filter(progression_status_code %in% c('PROG/RL', 'ST/RM')) %>%
    mutate(cancer_type = paste0(cancer_type_code, '_', progression_status_code))
  
  # duplicate metastatic rows and change the cancer type accordingly to create artificial cancer types contrasts
  
  ## Prostate
  prad_met <- metadata %>%
    filter(cancer_type_code == 'PRAD',
           cohort == 'Hartwig')
  
  prad_met_stable <- prad_met %>%
    mutate(cancer_type = 'PRAD_ST/RM')
  
  prad_met_relapse <- prad_met %>%
    mutate(cancer_type = 'PRAD_PROG/RL')
  
  ## Pancreas neuroendocrine
  panet_met <- metadata %>%
    filter(cancer_type_code == 'PANET',
           cohort == 'Hartwig')
  
  panet_met_stable <- panet_met %>%
    mutate(cancer_type = 'PANET_ST/RM')
  
  panet_met_relapse <- panet_met %>%
    mutate(cancer_type = 'PANET_PROG/RL')
  
  # combine all of the tables
  metadata <- progression_cancer_types %>%
    union(., progression) %>%
    union(., prad_met_relapse) %>%
    union(., prad_met_stable) %>%
    union(., panet_met_relapse) %>%
    union(., panet_met_stable)
  
}

# ---------------------- METASTATIC LOCATION

# do analysis by metastatic location
if (by_met_location == 'yes') {
  
  met_loc <- metadata %>%
    filter(!(metastatic_location %in% c('Unknown', 'Primary', NA_character_)))
  
  met_loc_cancer_types <- unique(met_loc$cancer_type)
  
  met_loc <- met_loc %>%
    mutate(cancer_type = paste0(cancer_type, '_', metastatic_location)) %>%
    mutate(cancer_type_code = paste0(cancer_type_code, '_', metastatic_location))
  
  # duplicate PCAWG samples rows to create artificial metastatic location contrast groups
  pcawg_local <- metadata %>%
    filter(cohort == 'PCAWG') %>%
    mutate(cancer_type = paste0(cancer_type, '_Local')) %>%
    mutate(cancer_type_code = paste0(cancer_type_code, '_Local'))
  
  pcawg_distant <- metadata %>%
    filter(cohort == 'PCAWG') %>%
    mutate(cancer_type = paste0(cancer_type, '_Distant')) %>%
    mutate(cancer_type_code = paste0(cancer_type_code, '_Distant'))
  
  pcawg_lymph <- metadata %>%
    filter(cohort == 'PCAWG') %>%
    mutate(cancer_type = paste0(cancer_type, '_Lymph')) %>%
    mutate(cancer_type_code = paste0(cancer_type_code, '_Lymph'))
  
  # filter the dataset for cancer types that actually have metastatic location annotation
  metadata <- metadata %>%
    filter(cancer_type %in% met_loc_cancer_types)
  
  # add artificial pcawg duplicate rows to metadata
  metadata <- metadata %>%
    union(., met_loc) %>%
    union(., pcawg_local) %>%
    union(., pcawg_distant) %>%
    union(., pcawg_lymph)
  
  # filter out those cancer_types with low hartwig samples in a certain metastatic location category
  low_hartwig_samples <- metadata %>%
    group_by(cancer_type, cohort) %>%
    count() %>%
    ungroup() %>%
    pivot_wider(names_from = 'cohort', values_from = 'n', values_fill = 0) %>%
    filter(pmin(Hartwig, PCAWG) < 5) %>%
    pull(cancer_type)
  
  metadata <- metadata %>%
    filter(!cancer_type %in% low_hartwig_samples)
}

# ----------------- Splitting the dataset by specified binary strata

# transform strata column into a logical if needed
metadata <- metadata %>%
  mutate(strata = if_else(cohort == 'Hartwig', TRUE, FALSE))

# ---------------- LINX FUSIONS

linx_fusions <- read_tsv(paste0(base_dir, '/results/02_combine_driver_catalogs/', results_date ,'/fusions/results/linx_fusions.tsv'),
                         col_names = TRUE) %>%
  filter(reported == TRUE)

# left join the metadata table, add driver column, rename 'name' column to 'gene'
linx_fusions <- linx_fusions %>%
  left_join(metadata, by = 'sample_id') %>%
  mutate(driver = 'FUSION') %>%
  rename(gene = 'name')

## Rename fusions (stolen from luan)
## Rename one sided promiscuous fusions, add new column 'fusion_name'
linx_fusions$fusion_name <- (function(){
  v <- linx_fusions$gene
  
  is_promiscuous_5 <- linx_fusions$reportedType == 'PROMISCUOUS_5'
  v[is_promiscuous_5] <- paste0(linx_fusions$geneStart[is_promiscuous_5], '_*')
  
  is_promiscuous_3 <- linx_fusions$reportedType == 'PROMISCUOUS_3'
  v[is_promiscuous_3] <- paste0('*_', linx_fusions$geneEnd[is_promiscuous_3])
  
  return(v)
})()

## Split fusions with 2 promiscuous genes into 2 entries
if(sum(linx_fusions$reportedType == 'PROMISCUOUS_BOTH')>0){
  linx_fusions <- (function(){
    df_promiscuous_both <- subset(linx_fusions, reportedType == 'PROMISCUOUS_BOTH')
    
    rbind(
      subset(linx_fusions, reportedType != 'PROMISCUOUS_BOTH'),
      within(df_promiscuous_both, fusion_name <- paste0(geneStart, '_*')),
      within(df_promiscuous_both, fusion_name <- paste0('*_', geneEnd))
    )
  })()
}

# drop "gene" column, rename 'fusion_name' to 'gene'
linx_fusions <- linx_fusions %>%
  select(-gene) %>%
  rename(gene = 'fusion_name')

# define iterators
all_genes <- unique(linx_fusions$gene)
all_cancers <- unique(metadata$cancer_type)
all_drivers <- unique(linx_fusions$driver)

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
    fusion_contingency_list[[cancer]][[driver]] <- do.call(rbind, unname(fusion_contingency_list[[cancer]][[driver]]))

  }

  fusion_contingency_list[[cancer]] <- do.call(rbind, unname(fusion_contingency_list[[cancer]]))
}

# combine all the contingency dfs into one
fusion_contingency_df <- do.call(rbind, unname(fusion_contingency_list))

# add the cancer_type, driver and gene column again, parsed from the row names
fusion_contingency_df <- fusion_contingency_df %>%
  mutate(driver = rownames(fusion_contingency_df), .before = mut_false_group) %>%
  separate(col = driver, into = c('cancer_type', 'driver', 'gene'), sep = '\\;')

# changed output file name based on input
output_file <- paste0(output_dir, '/', current_date, '_fusion_contingency_matrix')

if (by_subtype == 'yes') {
  if (by_subtype_met_location == 'yes') {
    output_file <- paste0(output_file, '_by_subtype_and_met_location.tsv')
  } else {
    output_file <- paste0(output_file, '_by_subtype.tsv')
  }
} else if (by_progression == 'yes') {
  output_file <- paste0(output_file, '_by_progression.tsv')
} else if (by_met_location == 'yes') {
  output_file <- paste0(output_file, '_by_met_location.tsv')
} else {
  output_file <- paste0(output_file, '.tsv')
}
# write table to disk
write_tsv(fusion_contingency_df, file = output_file,
          col_names = TRUE, append = FALSE)
