############### LINX driver enrichment analysis ###############
# author: Remy (sascha)
# date: 21/07/2021
# last updated: 31/10/2022

### Description
# This script creates the dotplot heatmap and volcano plot for the driver prevalence analysis.
# The user can choose whether to perform the analysis:
# 1. At cancer type level (by_subtype == FALSE; by_progression == FALSE, by_met_location == FALSE)
# 2. by cancer subtypes (by_subtype == TRUE; by_progression == FALSE, by_met_location == FALSE)
# 3. by primary progression (by_subtype == FALSE; by_progression == TRUE, by_met_location == FALSE)
# 4. by metastatic location (by_subtype == FALSE; by_progression == FALSE, by_met_location == TRUE).
# If the user chooses to perform the analysis at cancer type level, there is a subtype filter built in that
# excludes all significant hits at cancer type level that are not also significant at subtype level in 
# cancer types that have subtype annotation.

### Input
# metadata.tsv
# contingency matrices generated in step 4.

### Output
# Figure 5b
# Supplementary Figure 5a

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

output_dir <- list(
  path=paste0(here::here(), '/results')
)

for(i in output_dir){
  if(dir.exists(i)){
    output_dir <- i
    break
  }
}

# create output dir
dir.create(path = paste0(output_dir, '/05_plot_enrichment'), recursive = TRUE)
output_dir <- paste0(output_dir, '/05_plot_enrichment')

# get todays date
current_date <- format(Sys.Date(), '%Y_%m_%d')

## custom fonts
font_add_google(
  name = 'Inter',
  family = 'Helvetica'
)

## ggplot custom themes
theme_bigfont_heatmap <- theme(plot.title = element_text(size = 22),
                               plot.subtitle = element_text(size = 16),
                               axis.text.x = element_text(size = 15),
                               axis.text.y = element_text(size = 10), 
                               axis.title = element_text(size = 18),
                               legend.text = element_text(size = 12),
                               strip.text.x = element_text(size = 15))

theme_bigfont_volcano <- theme(plot.title = element_text(size = 22),
                               plot.subtitle = element_text(size = 16),
                               axis.text.x = element_text(size = 15),
                               axis.text.y = element_text(size = 15), 
                               axis.title = element_text(size = 18),
                               legend.text = element_text(size = 12),
                               strip.text.x = element_text(size = 15))

# ---------------------- Custom functions & objects

# source the functions to calculate the fishers exact test and effect sizes (cramer's V and odds Ratio)
# from test matrices
source(paste0(base_dir, '/code/00_func/fisher_wilcoxon_matrix_functions.R'))

# define cancer type order
source(paste0(base_dir, '/code/r_objects/cancer_type_order.R'))

# import cohort color
source(paste0(base_dir, '/code/r_objects/color_codes.R'))

# ---------------------- Define vars

contingency_matrix_results_date <- '2022_11_15' # old: 2021_07_30, new (with intogen filter): 2021_11_18, PANEL == 70: 2022_09_29, PANEL == 90: 2021_12_07, with subtypes/met loc: 2022_11_15
fusion_contingency_matrix_results_date <- '2022_11_15' # old: 2021_08_04, new: 2021_10_26, PANEL == 70: 2022_09_29, PANEL == 90: 2021_12_07, with subtypes/met loc: 2022_11_15
q_val_cutoff <- 0.01 # any value between 0 - 1
plot_clonality <- FALSE # TRUE: only show 'clonal' or 'subclonal' mutations in plot (specify below)
clonality <- 'clonal' # clonal/subclonal
strata <- 'cohort' # basically any binary strata column you have in the dataset
by_subtype <- FALSE # TRUE: perform analysis by cancer subtype, FALSE: dont
by_subtype_met_location <- FALSE # TRUE: perform analysis by cancer subtype and metastatic location, FALSE: dont, by_subtype must be TRUE!
by_progression <- FALSE # TRUE: perform analysis by primary progression type, FALSE: dont
by_met_location <- FALSE # TRUE: perform analysis by metastatic location, FALSE: dont

# checks
if (sum(c(plot_clonality, by_subtype, by_progression, by_met_location)) > 1) {
  stop('Multiple conditionals are TRUE, choose only one while leaving the others FALSE.')
}

# ---------------------- METADATA

metadata <- read_tsv(paste0(base_dir, '/data/processed/metadata/metadata.tsv'),
                     col_types = 'icccccccccccccccllllcccdiiildclccllllllliiiiiiiiiii') %>%
  # select relevant columns
  select(patient_id, sample_id, cohort, cancer_type, cancer_type_code, 
         cancer_subtype, progression_status_code, is_blacklisted_subtype, metastatic_location) %>%
  # add has_subtype column for easier filtering
  mutate(has_subtype = if_else(cancer_type_code %in% c('BRCA', 'COREAD', 'UCEC'), TRUE, FALSE)) %>%
  mutate(cancer_maintype = cancer_type, .before = cancer_type) %>%
  mutate(cancer_maintype_code = cancer_type_code)

# ---------------------- SUBTYPES

# do you want to do the analysis by cancer subtypes?
if (by_subtype) {
  
  # duplicate the rows of samples with subtype annotation, fuse together cancer_type and cancer_subtype column
  if (by_subtype_met_location) {
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
    
    # add artificial pcawg duplicate rows
    subtypes <- subtypes_met_loc %>%
      union(., pcawg_local) %>%
      union(., pcawg_distant) %>%
      union(., pcawg_lymph)
    
    # filter out those cancer subtypes with low hartwig samples in a certain biopsy loc category
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
  
  # add subtypes as separate column
  subtypes_plot <- metadata %>%
    distinct(cancer_maintype, cancer_type, cancer_type_code)
}

# ---------------------- PRIMARY PROGRESSION

if (by_progression) {
  
  progression_cancer_types <- metadata %>%
    filter(cancer_type_code %in% c('PRAD', 'PANET')) %>%
    mutate(cancer_type = cancer_type_code)
  
  # create new contrast groups (primary progression stable vs. relapse, relapse vs. met, stable vs. met) for all cancer types with annotation
  progression <- metadata %>%
    filter(progression_status_code %in% c('PROG/RL', 'ST/RM')) %>%
    mutate(cancer_type = paste0(cancer_type_code, '_', progression_status_code))
  
  # duplicate met rows and change the cancer type accordingly
  
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
  
  # add subtypes as separate column
  subtypes_plot <- metadata %>%
    distinct(cancer_maintype, cancer_type, cancer_type_code)
}

# ---------------------- METASTATIC LOCATION

# do analysis by metastatic location
if (by_met_location) {
  met_loc <- metadata %>%
    filter(!(metastatic_location %in% c('Unknown', 'Primary', NA_character_)))
  
  met_loc_cancer_types <- unique(met_loc$cancer_type)
  
  met_loc <- met_loc %>%
    mutate(cancer_type = paste0(cancer_type, '_', metastatic_location))
  
  # duplicate PCAWG samples rows to give them a pseudo cancer_type + metastatic_location column for easier data handling
  pcawg_local <- metadata %>%
    filter(cohort == 'PCAWG') %>%
    mutate(cancer_type = paste0(cancer_type, '_Local'))
  
  pcawg_distant <- metadata %>%
    filter(cohort == 'PCAWG') %>%
    mutate(cancer_type = paste0(cancer_type, '_Distant'))
  
  pcawg_lymph <- metadata %>%
    filter(cohort == 'PCAWG') %>%
    mutate(cancer_type = paste0(cancer_type, '_Lymph'))
  
  # filter the dataset for cancer types that actually have met loc annotation
  metadata <- metadata %>%
    filter(cancer_type %in% met_loc_cancer_types)
  
  # add artificial pcawg duplicate rows
  metadata <- metadata %>%
    union(., met_loc) %>%
    union(., pcawg_local) %>%
    union(., pcawg_distant) %>%
    union(., pcawg_lymph)
  
  # filter out those cancer_types with low hartwig samples in a certain biopsy loc category
  low_hartwig_samples <- metadata %>%
    group_by(cancer_type, cohort) %>%
    count() %>%
    ungroup() %>%
    pivot_wider(names_from = 'cohort', values_from = 'n', values_fill = 0) %>%
    filter(pmin(Hartwig, PCAWG) < 5) %>%
    pull(cancer_type)
  
  metadata <- metadata %>%
    filter(!cancer_type %in% low_hartwig_samples)
  
  # create subtypes df (important later when joining non-significant subtypes results to significant hits)
  subtypes_plot <- metadata %>%
    distinct(cancer_maintype, cancer_type, cancer_type_code)
}

# ---------------------- calculate sample size from metadata

sample_size <- metadata %>%
  group_by(cancer_type, cancer_type_code, across(all_of(strata))) %>%
  summarise(sample_size_tot = n(), .groups = 'drop') %>%
  pivot_wider(names_from = c('cohort'), values_from = 'sample_size_tot', names_prefix = 'n_') %>%
  replace_na(., replace = list(n_Hartwig = 0, n_PCAWG = 0)) %>%
  mutate(sample_size = n_Hartwig + n_PCAWG) %>%
  mutate(cancer_type_label = str_c(cancer_type, ' (', n_PCAWG, ' vs. ', n_Hartwig, ')')) %>%
  mutate(cancer_type_code_label = str_c(cancer_type_code, ' (', n_PCAWG, ' vs. ', n_Hartwig, ')')) %>%
  ungroup()

# ---------------------- CONTINGENCY MATRIX

# base input path for matrix
driver_matrix_path <- paste0(base_dir, '/results/04_get_contingency_matrix/results/', contingency_matrix_results_date ,'/', 
                             contingency_matrix_results_date, '_driver_contingency_matrix')
fusion_matrix_path <- paste0(base_dir, '/results/04_get_contingency_matrix/results/', contingency_matrix_results_date ,'/', 
                             contingency_matrix_results_date, '_fusion_contingency_matrix')

if (by_subtype) {
  if (by_subtype_met_location) {
    
    drivers_contingency_matrix_df <- read_tsv(paste0(driver_matrix_path, '_by_subtype_and_met_location.tsv'),
                                              col_names = TRUE)
    
    fusion_contingency_matrix_df <- read_tsv(paste0(fusion_matrix_path, '_by_subtype_and_met_location.tsv'),
                                             col_names = TRUE)
  } else {
    
    drivers_contingency_matrix_df <- read_tsv(paste0(driver_matrix_path, '_by_subtype.tsv'),
                                              col_names = TRUE) 
    
    fusion_contingency_matrix_df <- read_tsv(paste0(fusion_matrix_path, '_by_subtype.tsv'),
                                             col_names = TRUE) 
  }
} else if (by_progression) {
  
  drivers_contingency_matrix_df <- read_tsv(paste0(driver_matrix_path, '_by_progression.tsv'),
                                            col_names = TRUE)
  
  fusion_contingency_matrix_df <- read_tsv(paste0(fusion_matrix_path, '_by_progression.tsv'),
                                           col_names = TRUE)
} else if (by_met_location) {
  
  drivers_contingency_matrix_df <- read_tsv(paste0(driver_matrix_path, '_by_met_location.tsv'),
                                            col_names = TRUE) 
  
  fusion_contingency_matrix_df <- read_tsv(paste0(fusion_matrix_path, '_by_met_location.tsv'),
                                           col_names = TRUE)
} else {
  
  drivers_contingency_matrix_df <- read_tsv(paste0(driver_matrix_path, '.tsv'),
                                            col_names = TRUE) 
  
  fusion_contingency_matrix_df <- read_tsv(paste0(fusion_matrix_path, '.tsv'),
                                           col_names = TRUE)
  
  # ------------------ SUBTYPES DATA
  
  # read in the subtypes contingency matrices as well to exclude sig. hits at cancer type level 
  # that are not found in at least one subtype in cancer types that have subtype annotation
  
  # duplicate the rows of samples with subtype annotation, fuse together cancer_type and cancer_subtype column
  subtypes <- metadata %>%
    filter(has_subtype & !is_blacklisted_subtype) %>%
    mutate(cancer_type = paste0(cancer_type, '_', cancer_subtype))
  
  # add subtype rows to the metadata df
  metadata <- bind_rows(metadata, subtypes)
  
  # read in subtypes contingency matrix results to apply subtypes filter
  drivers_contingency_matrix_df_subtypes <- read_tsv(paste0(driver_matrix_path, '_by_subtype.tsv'),
                                                     col_names = TRUE) 
  
  fusion_contingency_matrix_df_subtypes <- read_tsv(paste0(fusion_matrix_path, '_by_subtype.tsv'),
                                                    col_names = TRUE)
  
  # combine both matrices
  contingency_matrix_df_subtypes <- rbind(drivers_contingency_matrix_df_subtypes, fusion_contingency_matrix_df_subtypes)
  
  # remove irrelevant subtypes from contingency table
  contingency_matrix_df_subtypes <- contingency_matrix_df_subtypes %>%
    mutate(cancer_type = str_replace(cancer_type, pattern = regex('_'), replacement = ' ')) %>%
    semi_join(subtypes, by = 'cancer_type')
  
  # create subtypes df (important later when joining non-significant subtypes results to significant hits)
  subtypes_plot <- metadata %>%
    distinct(cancer_maintype, cancer_type, cancer_type_code)
}

# combine both matrices
contingency_matrix_df <- rbind(drivers_contingency_matrix_df, fusion_contingency_matrix_df)

# add subtypes for the subtype filters
if (exists('contingency_matrix_df_subtypes')) {
  contingency_matrix_df <- contingency_matrix_df %>%
    mutate(cancer_type = str_replace(cancer_type, pattern = regex('_'), replacement = ' '))
  contingency_matrix_df <- rbind(contingency_matrix_df, contingency_matrix_df_subtypes)
}

# replace the first underscore in cancer_type
if (by_subtype) {
  contingency_matrix_df <- contingency_matrix_df %>%
    mutate(cancer_type = str_replace(cancer_type, pattern = regex('_'), replacement = ' ')) %>%
    right_join(subtypes_plot, by = 'cancer_type') %>%
    mutate(cancer_type_label = paste0(cancer_type, ' (', mut_false_group, '/', (mut_false_group + wt_false_group), ' vs. ', mut_true_group, '/', (mut_true_group + wt_true_group), ')'), .after = cancer_type) %>%
    mutate(sample_size = mut_false_group + mut_true_group + wt_false_group + wt_true_group)
} else if (by_met_location) {
  contingency_matrix_df <- contingency_matrix_df %>%
    mutate(cancer_type = str_replace(cancer_type, pattern = regex('_'), replacement = ' ')) %>%
    mutate(cancer_type = str_replace(cancer_type, pattern = regex(' (Distant|Local|Lymph)'), replacement = '_\\1')) %>%
    left_join(subtypes_plot, by = 'cancer_type') %>%
    mutate(cancer_type_label = paste0(cancer_type, ' (', mut_false_group, '/', (mut_false_group + wt_false_group), ' vs. ', mut_true_group, '/', (mut_true_group + wt_true_group), ')'), .after = cancer_type) %>%
    mutate(sample_size = mut_false_group + mut_true_group + wt_false_group + wt_true_group)
} else if (by_progression) {
  contingency_matrix_df <- contingency_matrix_df %>%
    left_join(subtypes_plot, by = 'cancer_type') %>%
    mutate(cancer_type_label = paste0(cancer_type, ' (', mut_false_group, '/', (mut_false_group + wt_false_group), ' vs. ', mut_true_group, '/', (mut_true_group + wt_true_group), ')'), .after = cancer_type) %>%
    mutate(sample_size = mut_false_group + mut_true_group + wt_false_group + wt_true_group)
} else {
  contingency_matrix_df <- contingency_matrix_df %>%
    right_join(subtypes_plot, by = 'cancer_type') %>%
    mutate(cancer_type_label = paste0(cancer_type, ' (', mut_false_group, '/', (mut_false_group + wt_false_group), ' vs. ', mut_true_group, '/', (mut_true_group + wt_true_group), ')'), .after = cancer_type) %>%
    mutate(sample_size = mut_false_group + mut_true_group + wt_false_group + wt_true_group)
}

# create the contingency matrix with only the four columns needed in the correct order
contingency_matrix <- contingency_matrix_df %>%
  select(mut_false_group, wt_false_group, mut_true_group, wt_true_group) %>%
  as.matrix(.)

# calculate the fisher p-value, the odds ratio and cramers V from that matrix, row by row
contingency_matrix_df$fisher_pval <- fisherTest.matrix(contingency_matrix)
contingency_matrix_df$odds_ratio <- oddsRatio.matrix(contingency_matrix)
contingency_matrix_df$cramers_v <- cramerV.matrix(contingency_matrix)

# multiply cramers V by -1 so that the metastatic cancers get the positive numbers
contingency_matrix_df <- contingency_matrix_df %>%
  mutate(cramers_v = -1 * cramers_v)

# filter out cases where there are less than 5 mutated samples in both primary or metastatic group
min_freq_idx <- pmax(contingency_matrix_df$mut_false_group, contingency_matrix_df$mut_true_group) > 4
contingency_matrix_df_minfreq <- contingency_matrix_df %>%
  filter(min_freq_idx)

# adjust p-value for multiple testing across all the rows using the FDR method
contingency_matrix_df_minfreq <- contingency_matrix_df_minfreq %>%
  group_by(cancer_type, driver) %>%
  mutate(fisher_pval_adj = p.adjust(fisher_pval, method = 'fdr'), .after = fisher_pval) %>%
  ungroup()

# put cancer types in correct order
contingency_matrix_df_minfreq <- contingency_matrix_df_minfreq %>%
  { if (sum(c(by_subtype, by_progression, by_met_location)) == 0) mutate(., cancer_maintype = factor(cancer_maintype, levels = cancer_type_order)) else . }

# filter for significantly enriched drivers (either in primary or metastatic group)
contingency_matrix_df_sig <- contingency_matrix_df_minfreq %>%
  filter(fisher_pval_adj < q_val_cutoff)

# prepare the df for plotting
contingency_matrix_df_sig <- contingency_matrix_df_sig %>%
  # add exclusivity column to show which drivers are exclusive to either prim or met
  mutate(exclusivity = case_when(mut_false_group == 0 ~ 'Hartwig only',
                                 mut_true_group == 0 ~ 'PCAWG only',
                                 TRUE ~ 'Both'))

# subtypes hits filter: only keep significant hits at cancer type level if the same hit is also significant in at least one subtype
if (sum(c(by_subtype, by_progression, by_met_location)) == 0) {
  contingency_matrix_df_sig_subtype_hits <- contingency_matrix_df_sig %>%
    filter(str_detect(cancer_type, pattern = regex('\\_')))
  
  contingency_matrix_df_sig_list <- contingency_matrix_df_sig %>%
    mutate(split_col = if_else(cancer_type_code %in% c('BRCA'), 'split_group1', 'split_group2')) %>%
    group_by(split_col) %>%
    group_split(.keep = FALSE)
  
  contingency_matrix_df_sig_list[[1]] <- contingency_matrix_df_sig_list[[1]] %>%
    semi_join(contingency_matrix_df_sig_subtype_hits, by = c('cancer_maintype', 'gene', 'driver'))
  
  contingency_matrix_df_sig <- do.call(rbind, contingency_matrix_df_sig_list)
  
  # remove subtypes hits
  contingency_matrix_df_sig <- contingency_matrix_df_sig %>%
    filter(!str_detect(cancer_type, pattern = '_'))
}

# shortened facet labels
if (plot_clonality) {
  facet_labels <- c(
    `AMP` = 'AMP',
    `DEL` = 'DEL',
    `FUSION` = 'FUS',
    `MUTATION` = paste0('MUT (', clonality, ')')
  )
  plot_subtitle <- paste0('Showing ', clonality, ' point mutations')
} else {
  facet_labels <- c(
    `AMP` = 'AMP',
    `DEL` = 'DEL',
    `FUSION` = 'FUS',
    `MUTATION` = paste0('MUT (total)')
  )
  plot_subtitle <- 'Showing all point mutations'
}

########### calculate basic statistics about the results
# freq_driver_muts <- contingency_matrix_df_sig %>%
#   mutate(freq_pcawg = mut_false_group / wt_false_group * 100,
#          freq_hartwig = mut_true_group / wt_true_group * 100, .after = wt_true_group)
# 
# enrichment_stats <- contingency_matrix_df_sig %>%
#   mutate(enriched_in = if_else(cramers_v > 0, 'enriched in hartwig', 'enriched in PCAWG'))
# 
# total_enrichment <- nrow(enrichment_stats)
# unique_enriched_driver_genes <- n_distinct(enrichment_stats$gene)
# unique_cancer_types <- n_distinct(enrichment_stats$cancer_type)
# 
# recurrent_drivers <- enrichment_stats %>%
#   group_by(gene) %>%
#   summarise(n = n_distinct(cancer_type), .groups = 'drop') %>%
#   mutate(percent_cancer_type = n / n_distinct(enrichment_stats$cancer_type) * 100)
# 
# percent_enriched_in <- enrichment_stats %>%
#   group_by(enriched_in) %>%
#   count() %>%
#   mutate(percent = n / total_enrichment * 100)
# 
# percent_exclusive <- enrichment_stats %>%
#   group_by(exclusivity) %>%
#   count() %>%
#   mutate(percent = n / total_enrichment * 100)
# 
# percent_mut_type_cohort <- enrichment_stats %>%
#   group_by(driver, enriched_in) %>%
#   count() %>%
#   group_by(driver) %>%
#   mutate(n_per_driver = sum(n)) %>%
#   mutate(percent = n / n_per_driver * 100)
###########

# order subtypes
cancer_type_code_order_subtype <- str_c(cancer_type_code_order, '_subtype')

# give correct names to colors
names(cancer_type_palette) <- cancer_type_code_order_subtype

# get the fisher pval adj data and combine with the complete contingency table
if (by_subtype | by_met_location | by_progression) {
  contingency_matrix_df_table <- contingency_matrix_df_minfreq %>%
    select(cancer_type_label, gene, driver, fisher_pval_adj) %>%
    right_join(contingency_matrix_df, by = c('cancer_type_label', 'gene', 'driver'))
  
  eff_table <- contingency_matrix_df_sig %>%
    distinct(cancer_maintype, cancer_type_code, gene, driver) %>%
    inner_join(., contingency_matrix_df_table, by = c('cancer_maintype', 'cancer_type_code', 'gene', 'driver'))
  
  # some prep for display
  eff_table <- eff_table %>%
    # make the label prettier for subtype & met location plots
    { if (by_subtype) {
      if (by_subtype_met_location) {
        mutate(., cancer_type_label = str_replace(cancer_type_label, pattern = regex('\\_(.+)\\_(.+)$'), replacement = ' \\1: Primary vs. \\2'))
      } else {
        mutate(., cancer_type_label = str_replace(cancer_type_label, pattern = regex('\\_(.+)$'), replacement = ': \\1'))
      } 
    } else if (by_met_location) {
      mutate(., cancer_type_label = str_replace(cancer_type_label, pattern = regex('\\_(.+)$'), replacement = ': Primary vs. \\1'))
    } else { . }} %>%
    mutate(cramers_v = round(cramers_v, 2)) %>%
    mutate(fisher_pval = case_when(
      is.na(fisher_pval) ~ 'Not calc (min. freq < 5)',
      fisher_pval < 0.01 ~ '< 0.01',
      fisher_pval < 0.05 ~ '< 0.05',
      TRUE ~ 'n.s.'
    )) %>%
    mutate(fisher_pval_adj = case_when(
      is.na(fisher_pval_adj) ~ 'Not calc (min. freq < 5)',
      fisher_pval_adj < 0.01 ~ '< 0.01',
      fisher_pval_adj < 0.05 ~ '< 0.05',
      TRUE ~ 'n.s.'
    )) %>%
    arrange(driver, gene, cancer_type_label)
  
  eff_table_plot <- eff_table %>%
    select(cancer_maintype, cancer_type_label, gene, driver, cramers_v, fisher_pval, fisher_pval_adj) %>%
    group_by(driver) %>%
    gt(.) %>%
    data_color(
      .,
      columns = cancer_maintype,
      colors = cancer_type_palette,
      alpha = 0.3,
      autocolor_text = FALSE
    ) %>%
    gt::tab_style(
      style = list(
        cell_borders(
          sides = 'top',
          color = 'black',
          style = 'solid',weight = 10
        )
      ),
      locations = list(
        cells_row_groups()
      )
    ) %>%
    gt::cols_hide(columns = 'cancer_maintype') %>%
    gt::cols_label(
      cancer_type_label = 'Cancer type', 
      gene = 'Gene', 
      driver = 'Driver', 
      cramers_v = "Cramer's V", 
      fisher_pval = 'P-value', 
      fisher_pval_adj = 'P-value (FDR adjusted)'
    )
  
  # save table as HTML
  eff_table_plot %>%
    gtsave(.,
           path = output_dir,
           filename = paste0(current_date, '_eff_size_minfreq5_q', q_val_cutoff,
                             if (by_subtype) { 
                               if (by_subtype_met_location) { '_by_subtype_and_met_location.html' }
                               else { '_by_subtype.html' } }
                             else if (by_progression) { '_by_progression.html' } 
                             else { '_by_met_location.html' }))
  
}

### dotplot heatmap

dotplot_heatmap <- contingency_matrix_df_sig %>%
  # make the label prettier for subtype & met location plots
  { if (by_subtype) {
    if (by_subtype_met_location) {
      mutate(., cancer_type = str_replace(cancer_type, pattern = regex('\\_(.+)\\_(.+)$'), replacement = ' \\1 Primary vs. \\2'))
    } else {
      mutate(., cancer_type = str_replace(cancer_type, pattern = regex('\\_(.+)$'), replacement = ' \\1'))
    } 
  } else if (by_met_location) {
    mutate(., cancer_type = str_replace(cancer_type, pattern = regex('\\_(.+)$'), replacement = ' Primary vs. \\1'))
  } else { . }} %>%
  ggplot(., aes(x = cancer_type, 
                y = fct_reorder(gene, desc(gene)))) +
  geom_point(aes(size = sample_size,
                 shape = exclusivity,
                 fill = cramers_v), color = 'lightgrey') +
  facet_grid(~driver, scales = 'free_x', space = 'free_x', labeller = as_labeller(facet_labels)) +
  scale_fill_gradient2(low = cohort_colors[1], mid = '#FFFFFF', high = cohort_colors[2], breaks = c(0.6,0.5, 0.4, 0.3, 0.2, 0.1, 0, -0.1, -0.2, -0.3, -0.4, -0.5, -0.6),
                       guide = guide_colorbar(title = "Eff size (Cramer's V)", order = 1,
                                              show.limits = TRUE, label.position = 'left', barheight = 10,
                                              frame.colour = 'black', ticks.colour = 'black')) +
  scale_size_continuous(range = c(2, 6), breaks = c(50, 250, 500, 1000), limits = c(10, 1000),
                        guide = guide_legend(title = 'Sample size', order = 2, override.aes = list(color = 'black'))) + # for q < 0.1
  scale_shape_manual(values = c(21, 24, 25), guide = guide_legend(title = '', order = 3, override.aes = list(size = 4, fill = 'black'))) +
  labs(
    x = '',
    y = '') +
  theme_bw() +
  theme_bigfont_heatmap +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = unit(c(1,1,1,3), 'cm'))

### volcano plot
  
# prepare df for volcano plot
contingency_matrix_df_volcano <- contingency_matrix_df_minfreq %>%
 { if (q_val_cutoff == 0.01 & !by_met_location)  filter(., driver != 'FUSION') else . } %>% # optional filter for prettier display of low q-val plot
  mutate(significance = if_else(fisher_pval_adj < q_val_cutoff, 'yes', 'no'))

if (by_subtype) {
  if (by_subtype_met_location) {
    contingency_matrix_df_sig_volcano <- contingency_matrix_df_volcano %>%
      filter(significance == 'yes') %>%
      separate(., col = 'cancer_type', into = c('cancer_type', 'cancer_subtype', 'metastatic_location'), sep = '_') %>%
      mutate(cancer_type_code = if_else(!is.na(cancer_subtype), paste0(cancer_type_code, '_', cancer_subtype, '_', metastatic_location), cancer_type_code)) %>%
      mutate(., gene_label = paste0(gene, ' (', cancer_type_code, ')'))
  } else {
    contingency_matrix_df_sig_volcano <- contingency_matrix_df_volcano %>%
      filter(significance == 'yes') %>%
      select(cancer_type, driver, gene, significance) %>%
      inner_join(contingency_matrix_df_sig, by = c('cancer_type', 'driver', 'gene')) %>%
      mutate(., gene_label = if_else(str_detect(cancer_type, pattern = regex('\\_')), paste0(gene, ' (', cancer_type_code, str_replace(cancer_type, pattern = regex('^.+(\\_.+$)'), replacement = '\\1'), ')'), paste0(gene, ' (', cancer_type_code, ')')))
  }
} else if (by_progression) {
  contingency_matrix_df_sig_volcano <- contingency_matrix_df_volcano %>%
    filter(significance == 'yes') %>%
    select(cancer_type, driver, gene, significance) %>%
    inner_join(contingency_matrix_df_sig, by = c('cancer_type', 'driver', 'gene')) %>%
    mutate(., gene_label = paste0(gene, ' (', cancer_type, ')'))
} else if (by_met_location) {
  contingency_matrix_df_sig_volcano <- contingency_matrix_df_volcano %>%
    filter(significance == 'yes') %>%
    select(cancer_type, driver, gene, significance) %>%
    inner_join(contingency_matrix_df_sig, by = c('cancer_type', 'driver', 'gene')) %>%
    mutate(., gene_label = if_else(str_detect(cancer_type, pattern = regex('\\_')), paste0(gene, ' (', cancer_type_code, str_replace(cancer_type, pattern = regex('^.+(\\_.+$)'), replacement = '\\1'), ')'), paste0(gene, ' (', cancer_type_code, ')')))
} else {
  contingency_matrix_df_sig_volcano <- contingency_matrix_df_volcano %>%
    filter(significance == 'yes') %>%
    select(cancer_type, driver, gene, significance) %>%
    inner_join(contingency_matrix_df_sig, by = c('cancer_type', 'driver', 'gene')) %>%
    mutate(., gene_label = if_else(str_detect(cancer_type, pattern = regex('\\_')), paste0(gene, ' (', cancer_type_code, str_replace(cancer_type, pattern = regex('^.+(\\_.+$)'), replacement = '\\1'), ')'), paste0(gene, ' (', cancer_type_code, ')')))
  
  # remove subtype results from df
  contingency_matrix_df_volcano <- contingency_matrix_df_volcano %>%
    filter(!str_detect(cancer_type, pattern = '_'))
}

# give correct names to colors
names(cancer_type_palette) <- cancer_type_order

# volcano plot; cramer's v
volcano_plot_cramersv <- ggplot(contingency_matrix_df_volcano, aes(x = cramers_v, y = -log10(fisher_pval_adj))) +
  geom_point(aes(x = cramers_v, color = significance, size = abs(cramers_v)),
  ) +
  geom_vline(xintercept = 0, color = '#B3B3B3', linetype = 'dashed') +
  geom_hline(yintercept = -log10(q_val_cutoff), color = '#B3B3B3', linetype = 'dashed') +
  labs(color = 'Significance') +
  geom_point(
    data = contingency_matrix_df_sig_volcano,
    aes(size = abs(cramers_v)),
    color = 'black',
    shape = 21,
    stroke = 1.1
  ) +
  geom_text_repel(data = contingency_matrix_df_sig_volcano, 
                  aes(label = gene_label, color = significance), key_glyph = 'rect', size = 5,
                  xlim = c(-1.1, 1.1), 
                  max.overlaps = Inf,
                  force = 9, force_pull = 0, 
                  max.time = 10, direction = 'both', 
                  min.segment.length = 0, point.padding = 0.5, box.padding = 0.5) +
  facet_wrap(~driver, labeller = as_labeller(facet_labels)) +
  scale_x_continuous(breaks = seq(-1, 1, by = 0.2), limits = c(-1, 1)) +
  scale_y_continuous(breaks = c(2, 5, 10, 15, 20), limits = c(0, 20)) +
  scale_color_manual(values = c('#E5E5E5', '#8B1A1A'), guide = 'none') +
  labs(
    x = "Effect size (Cramer's V)\nEnriched in PCAWG <=> Enriched in Hartwig" ,
    y = '-Log10(q-value)',
    color = '') +
  theme_bw() +
  theme_bigfont_volcano +
  guides(color = guide_legend(nrow = 11, byrow = TRUE)) +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'none')

# save plots
# as PDF
pdf(
  file = paste0(output_dir, '/', current_date, '_driver_fusion_dotplot_heatmap_minfreq5_q', q_val_cutoff,
                if (by_subtype) { 
                  if (by_subtype_met_location) { '_by_subtype_and_met_location.pdf' }
                  else { '_by_subtype.pdf' } }
                else if (by_progression) { '_by_progression.pdf' } 
                else if (by_met_location) { '_by_met_location.pdf'}
                else { '.pdf' }),
  width = if (q_val_cutoff == 0.01) { if (by_met_location) { 20 } else { 12 } } else { 16 },
  height = if (by_subtype_met_location) { 12 } else if (by_met_location) { 10 } else { 8 },
  useDingbats = FALSE,
  compress = FALSE
)
print(dotplot_heatmap)
dev.off()

pdf(
  file = paste0(output_dir, '/', current_date, '_driver_volcanoplot_minfreq5_q', q_val_cutoff,
                if (by_subtype) { 
                  if (by_subtype_met_location) { '_by_subtype_and_met_location.pdf' }
                  else { '_by_subtype.pdf' } }
                else if (by_progression) { '_by_progression.pdf' } 
                else if (by_met_location) { '_by_met_location.pdf'}
                else { '.pdf' }),
  width = 16, 
  height = if (q_val_cutoff == 0.01) { 8 } else { 16 },
  useDingbats = FALSE,
  compress = FALSE
  )
print(volcano_plot_cramersv)
dev.off()