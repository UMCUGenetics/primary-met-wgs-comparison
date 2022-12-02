############### Plot actionable variants ###############
# author: remy (sascha)
# date: 27/01/2022
# last updated: 25/10/2022

### Description
# input: actionable variants table that was parsed from VCF files by parse_actionable_variants.R script.
# This script plots the barplots including on and off label actionable variants. On-label variants are approved or have the potential to be approved
# for a certain cancer type. Off-label variants were approved or have the potential to be approved targets in other cancer types/

### Input
# metadata.tsv
# Clean Hartwig actionability table from step 8.

### Output
# Supplementary Figure 6

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

# get todays date
current_date <- format(Sys.Date(), '%Y_%m_%d')

## custom fonts
font_add_google(
  name = 'Inter',
  family = 'Helvetica'
)

## ggplot custom themes
resize_factor <- 1.2
theme_bigfont_barplot <- theme(plot.title = element_text(size = 16 * resize_factor),
                               plot.subtitle = element_text(size = 14 * resize_factor),
                               axis.text.x = element_text(size = 12 * resize_factor),
                               axis.text.y = element_text(size = 11 * resize_factor), 
                               axis.title = element_text(size = 18 * resize_factor),
                               legend.title = element_text(size = 16 * resize_factor),
                               legend.text = element_text(size = 14 * resize_factor),
                               strip.text.x = element_text(size = 11 * resize_factor),
                               strip.text.y = element_text(size = 10 * resize_factor))

theme_smallfont <- theme(plot.title = element_text(size = 14),
                         plot.subtitle = element_text(size = 11),
                         axis.text.x = element_text(size = 9),
                         axis.text.y = element_text(size = 9), 
                         axis.title = element_text(size = 11),
                         legend.title = element_text(size = 5),
                         legend.text = element_text(size = 10),
                         strip.text.y = element_text(size = 8))

# define dates of the input files and other variables
actionable_variants_results_date <- '2022_01_26' # old: 2022_01_20 new: 2022_01_26
response_filter <- 'sensitive' # can be 'sensitive' or 'resistant', restricts anaylsis to either sensitive or resistant actionable variants
q_val_cutoff <- 0.05 # significance threshold for all analysis
filter_sig <- F # TRUE: show only significant effect size differences as labels in plot, FALSE: show all effect size labels
by_subtype <- F # TRUE: perform analysis by cancer subtype, FALSE: perform analysis at cancer type level

# ------------------------------- Custom functions & objects

# source the functions to calculate the fishers exact test and wilcoxon test p-values and effect sizes
# from test matrices
source(paste0(base_dir, '/code/00_func/fisher_wilcoxon_matrix_functions.R'))

# define cancer type order
source(paste0(base_dir, '/code/r_objects/cancer_type_order.R'))
cancer_type_order <- c('Pancancer', 'Breast carcinoma', 'Breast carcinoma_ER+/HER2+', 'Breast carcinoma_ER+/HER2-', 'Breast carcinoma_ER-/HER2+', 'Breast carcinoma_TNBC', cancer_type_order[2:length(cancer_type_order)])
cancer_type_code_order <- c('BRCA', 'BRCA_ER+/HER2+', 'BRCA_ER+/HER2-', 'BRCA_ER-/HER2+', 'BRCA_TNBC', cancer_type_code_order[2:length(cancer_type_code_order)])

# import cohort color
source(paste0(base_dir, '/code/r_objects/color_codes.R'))

# define ASCO evidence tier level colors
tier_label_levels <- c('A_On_label', 'A_Off_label', 'B_On_label', 'B_Off_label')
label_colors <- c('#1e6220', '#2F913F', '#35E358', '#C4FFD0')
names(label_colors) <- tier_label_levels

# get contingency matrix functions
get_contingency_matrix <- function(cancer_type_input, label_input) {
  
  counts <- actionable_variants_clean_plot %>%
    select(cohort, cancer_type, label, n) %>%
    complete(., cohort, cancer_type, label, fill = list(n = 0)) %>%
    filter(cancer_type == cancer_type_input) %>%
    filter(label == label_input) %>%
    pivot_wider(., names_from = 'cohort', values_from = 'n', names_prefix = 'n_') %>%
    inner_join(sample_size, by = 'cancer_type') %>%
    pivot_wider(., names_from = 'cohort', values_from = 'sample_size') %>%
    mutate(n_non_Hartwig = Hartwig - n_Hartwig,
           n_non_PCAWG = PCAWG - n_PCAWG) %>%
    select(n_PCAWG, n_non_PCAWG, n_Hartwig, n_non_Hartwig)
  
  counts <- as.data.frame(counts)
  
  rownames(counts) <- paste0(str_replace(cancer_type_input, ' ', '_'), '_', label_input)
  
  return(counts)
  
}

get_global_contingency_matrix <- function(cancer_type_input) {
  
  counts <- actionable_variants_clean_plot_global %>%
    right_join(sample_size, by = c('cohort', 'cancer_type', 'cancer_type_code', 'cancer_type_label', 'sample_size')) %>%
    select(cohort, cancer_type, n) %>%
    complete(., cohort, cancer_type, fill = list(n = 0)) %>%
    filter(cancer_type == cancer_type_input) %>%
    pivot_wider(., names_from = 'cohort', values_from = 'n', names_prefix = 'n_') %>%
    inner_join(sample_size, by = 'cancer_type') %>%
    pivot_wider(., names_from = 'cohort', values_from = 'sample_size') %>%
    mutate(n_non_Hartwig = Hartwig - n_Hartwig,
           n_non_PCAWG = PCAWG - n_PCAWG) %>%
    select(n_PCAWG, n_non_PCAWG, n_Hartwig, n_non_Hartwig)
  
  counts <- as.data.frame(counts)
  
  rownames(counts) <- paste0(str_replace(cancer_type_input, ' ', '_'))
  
  return(counts)
  
}

# -------------------------------  OUTPUT DIR

# create output dir
dir.create(path = paste0(output_dir, '/08_actionability/', actionable_variants_results_date, '/results/plots'), recursive = TRUE)
output_dir <- paste0(output_dir, '/08_actionability/', actionable_variants_results_date, '/results/plots')

# ------------------------------- METADATA

metadata <- read_tsv(paste0(base_dir, '/data/processed/metadata/metadata.tsv'),
                     col_types = 'icccccccccccccccllllcccdiiildclccllllllliiiiiiiiiii') %>%
  # select relevant columns
  select(sample_id, cohort, cancer_type, cancer_type_code, cancer_subtype, is_blacklisted_subtype) %>%
  # add has_subtype column for easier filtering
  mutate(has_subtype = if_else(cancer_type_code %in% c('BRCA', 'COREAD', 'UCEC'), TRUE, FALSE)) %>%
  mutate(cancer_maintype = cancer_type, .before = cancer_type) %>%
  mutate(cancer_maintype_code = cancer_type_code)

# split BRCA tumors into subtypes
if (!by_subtype) {
  metadata <- metadata %>%
    mutate(cancer_type_code = if_else(cancer_type_code == 'BRCA' & !is.na(cancer_subtype), paste0(cancer_type_code, '_', cancer_subtype), cancer_type_code)) %>%
    mutate(cancer_type = if_else(cancer_type == 'Breast carcinoma' & !is.na(cancer_subtype), paste0(cancer_type, '_', cancer_subtype), cancer_type)) %>%
    filter(cancer_type_code != 'BRCA')
}

# ------------------------------- SUBTYPES

# do you want to do the analysis by cancer subtypes?
if (by_subtype) {
  
  # duplicate the rows of samples with subtype annotation, fuse together cancer_type and cancer_subtype column
  subtypes <- metadata %>%
    filter(has_subtype & !is_blacklisted_subtype) %>%
    mutate(cancer_type = paste0(cancer_type, '_', cancer_subtype)) %>%
    mutate(cancer_type_code = paste0(cancer_type_code, '_', cancer_subtype))
  
  # add subtype rows to the metadata df
  metadata <- bind_rows(metadata, subtypes)
  
  # filter for cancer types that only have subtypes
  metadata <- metadata %>%
    filter(has_subtype)
  
  # add subtypes as separate column
  subtypes_plot <- metadata %>%
    distinct(cancer_maintype, cancer_type, cancer_type_code)
}

# ------------------------------- ACTIONABILITY

if (response_filter == 'sensitive') {
  response_type_filter <- c('sensitive', 'predicted - sensitive')
} else {
  response_type_filter <- c('resistant', 'predicted - resistant')
}

actionability <- read_tsv(file = paste0(base_dir, '/results/08_actionability/', actionable_variants_results_date, '/results/hartwig_actionability.tsv')) %>%
  # remove liquid tumors because they are weird
  filter(!str_detect(cancerType, pattern = 'leukemia|CLL')) %>%
  # remove variants for which we dont have data for (expression and hypermethylation)
  filter(!mutated_variant %in% c('over exp', 'dec exp', 'hypermethylation')) %>%
  # filter only for specified response type (either 'sensitive' or 'resistant')
  filter(responseType %in% response_type_filter) %>%
  select(efficacyEvidence, everything())

# only keep relevant rows
actionability <- actionability %>%
  select(ckbEntryId, tumorProfile, gene, mutated_variant, cancerType, ampCapAscoEvidenceLevel) %>%
  distinct() %>%
  # only keep the top ampCapAscoEvidenceLevel for each actionable variant in each cancer type
  # (e.g. If there is a level A and B entry for EGFR3 mut in LUAD, then only level A is kept)
  group_by(ckbEntryId, cancerType, tumorProfile) %>%
  slice_min(ampCapAscoEvidenceLevel, n = 1) %>%
  ungroup()

# calculate sample size by strata
sample_size <- metadata %>%
  group_by(cancer_type, cancer_type_code, cohort) %>%
  summarise(sample_size = n(), .groups = 'drop') %>%
  pivot_wider(names_from = 'cohort', values_from = 'sample_size') %>%
  mutate(cancer_type_label = paste0(cancer_type, '\n(', PCAWG, ' vs. ', Hartwig, ')')) %>%
  pivot_longer(., cols = c(Hartwig, PCAWG), names_to = 'cohort', values_to = 'sample_size')

# calc cohort size
cohort_size <- metadata %>%
  group_by(cohort) %>%
  summarise(cohort_size = n(), .groups = 'drop') %>%
  pivot_wider(names_from = 'cohort', values_from = 'cohort_size') %>%
  mutate(cancer_type_label = paste0('Pancancer\n(', PCAWG, ' vs. ', Hartwig, ')')) %>%
  pivot_longer(., cols = c(Hartwig, PCAWG), names_to = 'cohort', values_to = 'sample_size') %>%
  mutate(cancer_type = 'Pancancer', cancer_type_code = 'Pancancer')

# combine sample_size and cohort_size
sample_size <- sample_size %>%
  union(., cohort_size)

# ------------------------------- ACTIONABLE VARIANTS

actionable_variants_clean <- read_tsv(file = paste0(base_dir, '/results/08_actionability/', actionable_variants_results_date, '/results/actionable_variants.tsv')) %>%
  # exclude wild-type AVs
  filter(mutated_variant != 'wild_type') %>%
  # exclude sample level actionable variants (MSI & TMB variants)
  filter(!(str_detect(mutated_variant, pattern = regex('MSI|TMB')))) %>%
  # join the actionability table to get the cancer type labels of treatments that are approved or in early/late clinical trial
  inner_join(actionability, by = c('gene', 'mutated_variant')) %>%
  # join metadata to get the actual cancer type labels of the sample and cohort label
  inner_join(metadata, by = 'sample_id') %>%
  distinct()

# add a new column 'On/off-label' indicating whether the approved cancer type overlaps with the actual cancer type label
actionable_variants_cleantt <- actionable_variants_clean %>%
  mutate(on_off_label = if_else(cancer_maintype == cancerType, 'On_label', 'Off_label')) %>%
  # join the sample_size of each cancer type
  inner_join(sample_size, by = c('cohort', 'cancer_type', 'cancer_type_code'))

# for each sample get the highest actionability level, ranked by ASCO evidence level and On-/Off-label
actionable_variants_cleantt2 <- actionable_variants_cleantt %>%
  group_by(sample_id, on_off_label) %>%
  slice_min(ampCapAscoEvidenceLevel, n = 1) %>%
  group_by(sample_id) %>%
  slice_max(on_off_label, n = 1) %>%
  distinct(sample_id, cancer_type, ampCapAscoEvidenceLevel, on_off_label)

################### count number of samples per cohort, ASCO level and label
### global
actionable_variants_clean_plot_pancan_global <- actionable_variants_cleantt2 %>%
  inner_join(metadata, by = c('sample_id', 'cancer_type')) %>%
  group_by(cohort) %>%
  count() %>%
  ungroup() %>%
  mutate(cancer_type = 'Pancancer', cancer_type_code = 'Pancancer') %>%
  inner_join(sample_size, by = c('cohort', 'cancer_type', 'cancer_type_code')) %>%
  mutate(percent = n / sample_size * 100)

actionable_variants_clean_plot_global <- actionable_variants_cleantt2 %>%
  inner_join(metadata, by = c('sample_id', 'cancer_type')) %>%
  group_by(cohort, cancer_type) %>%
  count() %>%
  ungroup() %>%
  inner_join(sample_size, by = c('cohort', 'cancer_type')) %>%
  mutate(percent = n / sample_size * 100)

# add pancancer % to the per cancer % table
actionable_variants_clean_plot_global <- actionable_variants_clean_plot_pancan_global %>%
  union(., actionable_variants_clean_plot_global)

### per label level
actionable_variants_clean_plot_pancan <- actionable_variants_cleantt2 %>%
  right_join(metadata, by = c('sample_id', 'cancer_type')) %>%
  group_by(cohort, ampCapAscoEvidenceLevel, on_off_label) %>%
  count() %>%
  ungroup() %>%
  mutate(cancer_type = 'Pancancer', cancer_type_code = 'Pancancer') %>%
  inner_join(sample_size, by = c('cohort', 'cancer_type', 'cancer_type_code')) %>%
  mutate(percent = n / sample_size * 100) %>%
  unite(., col = 'label', c(ampCapAscoEvidenceLevel, on_off_label), sep = '_')

# count number of samples per cohort, cancer_type. ASCO level and label
actionable_variants_clean_plot <- actionable_variants_cleantt2 %>%
  right_join(metadata, by = c('sample_id', 'cancer_type')) %>%
  group_by(cohort, cancer_type, ampCapAscoEvidenceLevel, on_off_label) %>%
  count() %>%
  ungroup() %>%
  inner_join(sample_size, by = c('cohort', 'cancer_type')) %>%
  mutate(percent = n / sample_size * 100) %>%
  unite(., col = 'label', c(ampCapAscoEvidenceLevel, on_off_label), sep = '_')

# add pancancer % to the per cancer % table
actionable_variants_clean_plot <- actionable_variants_clean_plot_pancan %>%
  union(., actionable_variants_clean_plot)

################### complete the AV table to include all sample - label combinations
complete_cohort <- metadata %>% select(sample_id, cancer_type)

complete_tier_label_table <- actionable_variants_cleantt %>%
  distinct(sample_id, gene, mutated_variant, cancer_type, ampCapAscoEvidenceLevel, on_off_label) %>%
  group_by(sample_id, cancer_type, ampCapAscoEvidenceLevel, on_off_label) %>%
  count() %>%
  ungroup() %>%
  right_join(complete_cohort, by = c('sample_id', 'cancer_type')) %>%
  tidyr::replace_na(., replace = list(ampCapAscoEvidenceLevel = 'A', on_off_label = 'Off_label', n = 0)) %>%
  select(-cancer_type) %>%
  distinct() %>%
  complete(sample_id, ampCapAscoEvidenceLevel, on_off_label, fill = list(n = 0)) %>%
  inner_join(metadata, by = 'sample_id') %>%
  relocate(., n, .after = last_col()) %>%
  unite(., col = 'label', ampCapAscoEvidenceLevel, on_off_label)

################### calculate fisher p-value per cancer type and label between cohorts

# get all cancer type names and labels
all_labels <- unique(complete_tier_label_table$label)
all_cancers <- unique(sample_size$cancer_type)

### global

# initiate list
global_contingency_list <- list()

# calculate Fishers exact test p-value, FDR adjusted p-value
# cramers V and the Odds ratio for each label in the filtered complete_tier_label_table per cancer type
# stratified by cohort, then combine it all in one big table
for (cancer in all_cancers) {
  print(paste('Running...', cancer))
  global_contingency_list[[cancer]] <- get_global_contingency_matrix(cancer_type_input = cancer) 
}

# combine all the contingency dfs into one
global_contingency_df <- do.call(rbind, global_contingency_list)

# add the cancer_type, driver and gene column again, parsed from the row names
global_contingency_df <- global_contingency_df %>%
  rownames_to_column(., var = 'cancer_type')

# rename for stats calculation
global_stats <- global_contingency_df %>%
  mutate(label = 'Total', .after = cancer_type)

# filter out cases where there are less than 5 mutated samples in both primary or metastatic group
min_freq_idx <- pmax(global_contingency_df$n_PCAWG, global_contingency_df$n_Hartwig) > 4
global_contingency_df <- global_contingency_df %>%
  filter(min_freq_idx)

# create the contingency matrix with only the four columns needed in the correct order
contingency_matrix <- global_contingency_df %>%
  select(n_PCAWG, n_non_PCAWG, n_Hartwig, n_non_Hartwig) %>%
  as.matrix(.)

# calculate the fisher p-value, the odds ratio and cramers V from that matrix, row by row
global_contingency_df$fisher_pval <- fisherTest.matrix(contingency_matrix)
global_contingency_df$odds_ratio <- oddsRatio.matrix(contingency_matrix)
global_contingency_df$cramers_v <- cramerV.matrix(contingency_matrix)

# adjust p-value for multiple testing across all the rows using the FDR method
global_contingency_df <- global_contingency_df %>%
  group_by(cancer_type) %>%
  mutate(fisher_pval_adj = p.adjust(fisher_pval, method = 'fdr'), .after = fisher_pval) %>%
  ungroup()

# assign global comparison contingency table to a new variable for supplementary table later
global_contingency_df_supp <- global_contingency_df

# mark significantly different cases and calculate the percentage increase/decrease of AVs for those cases (PCAWG -> Hartwig)
global_contingency_df <- global_contingency_df %>%
  mutate(global_significance = if_else(fisher_pval_adj < q_val_cutoff, '*', 'ns'),
         cohort = 'Hartwig') %>%
  mutate(prop_Hartwig = n_Hartwig / (n_Hartwig + n_non_Hartwig),
         prop_PCAWG = n_PCAWG / (n_PCAWG + n_non_PCAWG)) %>%
  mutate(fold_change = round(prop_Hartwig / prop_PCAWG, 1),
         fold_change_label = str_c(fold_change, 'x')) %>%
  mutate(abs_diff = round((prop_Hartwig - prop_PCAWG), 1),
         abs_diff_label = if_else(abs_diff > 0, str_c('+', abs_diff), as.character(abs_diff))) %>%
  select(cohort, cancer_type, global_significance, fold_change, fold_change_label, abs_diff, abs_diff_label)

### per cancer type & label level

# initiate list
contingency_list <- list()

# calculate Fishers exact test p-value, FDR adjusted p-value (per cancer type, per driver mutation type),
# cramers V and the Odds ratio for each label in the filtered complete_tier_label_table per cancer type
# stratified by cohort, then combine it all in one big table
for (cancer in all_cancers) {
  for (label in all_labels) {
    print(paste('Running...', cancer, label))
    contingency_list[[cancer]][[label]] <- get_contingency_matrix(cancer_type_input = cancer,
                                                                  label_input = label)
  }
  contingency_list[[cancer]] <- do.call(rbind, contingency_list[[cancer]])
}

# combine all the contingency dfs into one
contingency_df <- do.call(rbind, contingency_list)

# add the cancer_type, driver and gene column again, parsed from the row names
contingency_df <- contingency_df %>%
  rownames_to_column(., var = 'label') %>%
  separate(col = 'label', into = c('cancer_type', 'label'), sep = '\\.')

# rename for stats calculation
per_label_stats <- contingency_df

# filter out cases where there are less than 5 mutated samples in both primary or metastatic group
min_freq_idx <- pmax(contingency_df$n_PCAWG, contingency_df$n_Hartwig) > 4
contingency_df <- contingency_df %>%
  filter(min_freq_idx)

# create the contingency matrix with only the four columns needed in the correct order
contingency_matrix <- contingency_df %>%
  select(n_PCAWG, n_non_PCAWG, n_Hartwig, n_non_Hartwig) %>%
  as.matrix(.)

# calculate the fisher p-value, the odds ratio and cramers V from that matrix, row by row
contingency_df$fisher_pval <- fisherTest.matrix(contingency_matrix)
contingency_df$odds_ratio <- oddsRatio.matrix(contingency_matrix)
contingency_df$cramers_v <- cramerV.matrix(contingency_matrix)

# adjust p-value for multiple testing across all the rows using the FDR method
contingency_df <- contingency_df %>%
  group_by(cancer_type) %>%
  mutate(fisher_pval_adj = p.adjust(fisher_pval, method = 'fdr'), .after = fisher_pval) %>%
  ungroup()

# add columns for barplot borders annotation (only sig. A_On_label are annotated)
contingency_df <- contingency_df %>%
  mutate(significance = if_else(fisher_pval_adj < q_val_cutoff & label == 'A_On_label', 'yes', 'no')) %>%
  mutate(effect_size = case_when(significance == 'yes' & odds_ratio < 1 ~ 'Enriched in Hartwig',
                                 significance == 'yes' & odds_ratio > 1 ~ 'Enriched in PCAWG',
                                 TRUE ~ 'none')) %>%
  select(cancer_type, label, significance, effect_size) %>%
  complete(., cancer_type, label, fill = list(significance = 'no', effect_size = 'none'))

# get the correct cancer type label order
cancer_type_label_order <- sample_size %>%
  { if (by_subtype) filter(.,  cancer_type != 'Pancancer') %>% arrange(., factor(cancer_type)) else arrange(., factor(cancer_type, levels = cancer_type_order)) } %>%
  pull(cancer_type_label) %>%
  unique()

# impose order on cancer type label & tier level, filter NA level
actionable_variants_clean_plottt <- actionable_variants_clean_plot %>%
  { if (by_subtype) { filter(., cancer_type != 'Pancancer') } else { . }} %>% # filter out meaningless 'Pancancer' group when performing analysis by subtype
  left_join(contingency_df, by = c('cancer_type', 'label')) %>%
  left_join(global_contingency_df, by = c('cohort', 'cancer_type')) %>%
  mutate(cancer_type_label = factor(cancer_type_label, levels = cancer_type_label_order)) %>%
  filter(label != 'NA_NA') %>%
  mutate(label = factor(label, levels = rev(names(label_colors))))

# calculate mean of mean % increase in AVs across cancer types
mean_diff <- global_stats %>%
  union(., per_label_stats) %>%
  filter(cancer_type != 'Pancancer') %>%
  mutate(prop_Hartwig = n_Hartwig / (n_Hartwig + n_non_Hartwig),
         prop_PCAWG = n_PCAWG / (n_PCAWG + n_non_PCAWG)) %>%
  mutate(fold_change = round(prop_Hartwig / prop_PCAWG, 1)) %>%
  mutate(abs_diff = round((prop_Hartwig - prop_PCAWG), 1)) %>%
  filter(fold_change != Inf) %>%
  group_by(label) %>%
  summarise(mean_fold_change = mean(fold_change),
            sd_fold_change = sd(fold_change),
            mean_abs_diff = mean(abs_diff),
            sd_abs_diff = sd(abs_diff))

# get the alternating index number for the grey background rectangles
if (by_subtype) {
  rect_idx_vec <- sort(unique(actionable_variants_clean_plottt$cancer_type)[order(match(unique(actionable_variants_clean_plottt$cancer_type), cancer_type_order))])
  rect_idx <- seq(2, length(rect_idx_vec), by = 2)
  rect_idx <- rect_idx_vec[rect_idx]
} else {
  rect_idx <- seq(2, length(cancer_type_label_order), by = 2)
  rect_idx_vec <- unique(actionable_variants_clean_plottt$cancer_type)[order(match(unique(actionable_variants_clean_plottt$cancer_type), cancer_type_order))]
  rect_idx <- rect_idx_vec[rect_idx]
}

# plot stacked barplot: Proportion of A & B On/Off-label AVs, facetted by cancer type
global_av_plot <- actionable_variants_clean_plottt %>%
  ggplot(., aes(x = cohort, y = percent)) +
  geom_rect(
    data = . %>% 
      filter(cancer_type %in% rect_idx) %>%
      distinct(cancer_type, cancer_type_label), 
    aes(fill = cancer_type),
    inherit.aes = FALSE,
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = Inf, 
    alpha = 0.3, fill = '#a0a0a0'
    ) +
  geom_col(
    aes(fill = label, color = effect_size), 
    position = 'stack', alpha = 0.9, size = 1.1
    ) +
  geom_text(
    data = . %>%
      { if (filter_sig) filter(., global_significance != 'ns') else . } %>%
      distinct(cohort, fold_change_label, cancer_type_label),
    aes(label = fold_change_label), y = 104, nudge_x = 0.5, size = 4
  ) +
  facet_wrap(
    cancer_type_label ~ ., ncol = 1, 
    strip.position = 'left'
    ) +
  scale_x_discrete(
    position = 'top',
    expand = c(0.5,0.5)
    ) +
  scale_y_continuous(
    limits = c(0, 110),
    breaks = seq(0, 100, by = 25)
    ) +
  scale_fill_manual(
    values = label_colors,
    labels = c('A On-label', 'A Off-label', 'B On-label', 'B Off-label')
    ) +
  scale_color_manual(
    values = c('Enriched in Hartwig' = '#9966CC', 'Enriched in PCAWG' = '#F58134', 'none' = NA), 
    na.value = 'transparent', 
    labels = c('Hartwig', 'PCAWG', '')
    ) +
  coord_flip() +
  labs(
    title = paste0('Showing ', response_filter, ' actionable variants'),
    subtitle = paste0('global q < ', q_val_cutoff),
    x = '',
    y = '% cohort',
    fill = ''
    ) +
  guides(color = guide_legend(title = 'A On-label enrichment in',
                              override.aes = list(fill = NA))) +
  theme_classic() +
  theme_bigfont_barplot +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    panel.spacing = unit(0, 'cm'))

# as PDF
pdf(
  file = paste0(output_dir, '/', current_date, '_', response_filter, '_av_barplot_noleuk_fc',
                if (by_subtype) { if (filter_sig) { '_by_subtype_sig.pdf' } else { '_by_subtype_all.pdf' } } 
                else { if (filter_sig) { '_sig.pdf' } else { '_all.pdf' } }),
  width = 10,
  height = 12,
  useDingbats = FALSE,
  compress = FALSE
)
print(global_av_plot)
dev.off()

################### Tier level AV counts per sample (Total, A & B On/Off-label)

################ avg number of actionable variants per sample by tier level and On/Off-label
# calculate total number of AVs per sample + add this total to the original table
av_per_sample_tier_label <- complete_tier_label_table %>%
  group_by(., sample_id, cohort, cancer_type, cancer_type_code, cancer_subtype, cancer_maintype, cancer_maintype_code, is_blacklisted_subtype, has_subtype) %>%
  summarise(n = sum(n), .groups = 'drop') %>%
  mutate(label = 'Total')

av_per_sample_tier_label <- complete_tier_label_table %>%
  union(., av_per_sample_tier_label)

############# calculate mean of mean AVs per sample to get an idea on pancancer level
av_per_sample_tier_label_mean_of_mean <- av_per_sample_tier_label %>%
  group_by(cancer_type, cohort, label) %>%
  summarise(mean_av = mean(n), .groups = 'drop') %>%
  group_by(cohort, label) %>%
  summarise(mean_av = round(mean(mean_av), 1), .groups = 'drop') %>%
  pivot_wider(., names_from = 'cohort', values_from = 'mean_av') %>%
  mutate(mean_diff = round(Hartwig - PCAWG, 2)) %>%
  mutate(mean_diff_label = if_else(mean_diff > 0, paste0('+ ', mean_diff), paste0('- ', abs(mean_diff)))) %>%
  pivot_longer(., cols = c('Hartwig', 'PCAWG'), names_to = 'cohort', values_to = 'mean_av')
############

# in non-subtype mode, duplicate the data, but add 'Pancancer' as the cancer type, then add this duplicated table to the original table
if (!by_subtype) {
  pancan_av_per_sample_tier_label <- av_per_sample_tier_label %>%
    mutate(cancer_type = 'Pancancer', cancer_type_code = 'Pancancer')
  
  av_per_sample_tier_label <- av_per_sample_tier_label %>%
    union(., pancan_av_per_sample_tier_label) 
}

# join sample_size table to the av_per_sample_tier_label for the plotting
av_per_sample_tier_label <- av_per_sample_tier_label %>%
  inner_join(sample_size, by = c('cohort', 'cancer_type', 'cancer_type_code')) %>%
  mutate(cancer_type_label = factor(cancer_type_label, levels = cancer_type_label_order))

# calculate mean av pancan & per cancer type + mean av difference PCAWG -> Hartwig
mean_av_sample_tier_label <- av_per_sample_tier_label %>%
  group_by(cohort, cancer_type_label, label) %>%
  summarise(mean_av = mean(n), .groups = 'drop') %>%
  pivot_wider(., names_from = 'cohort', values_from = 'mean_av') %>%
  mutate(mean_diff = round(Hartwig - PCAWG, 2)) %>%
  mutate(mean_diff_label = if_else(mean_diff > 0, paste0('+ ', mean_diff), paste0('- ', abs(mean_diff)))) %>%
  pivot_longer(., cols = c('Hartwig', 'PCAWG'), names_to = 'cohort', values_to = 'mean_av')

if (response_filter == 'sensitive') {
  total_y_limit <- 7
  a_on_label_limit <- 5
  a_off_label_limit <- 6
  b_on_label_limit <- 5
  b_off_label_limit <- 6
} else {
  total_y_limit <- 5
  a_on_label_limit <- 5
  a_off_label_limit <- 5
  b_on_label_limit <- 5
  b_off_label_limit <- 4
}

legend_x_pos <- 1.7

# av per sample violin plots
total_avps_plot <- av_per_sample_tier_label %>%
  filter(label == 'Total') %>%
  ggplot(., aes(x = cohort, y = log2(n + 1))) +
  geom_rect(data = . %>% 
              filter(cohort == 'Hartwig') %>%
              filter(cancer_type %in% rect_idx) %>%
              distinct(cancer_type_label),
            inherit.aes = FALSE,
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, 
            alpha = 0.3, fill = '#a0a0a0') +
  geom_jitter(color = 'gray40', size = 0.1, alpha = 0.25) +
  geom_violin(aes(color = cohort), fill = NA, position = position_dodge(1), size = 1.05, scale = 'width', width = 0.5, adjust = 2.5) +
  geom_text(data = . %>% 
              compare_means(
                data = ., formula = n ~ cohort,
                method = 'wilcox.test', group.by = 'cancer_type_label', p.adjust.method = 'fdr',
                symnum.args = list(cutpoints = c(0, q_val_cutoff, 1), symbols = c('*', 'ns'))) %>% 
              left_join(mean_av_sample_tier_label, by = 'cancer_type_label', 'group1' = 'cohort') %>%
              filter(label == 'Total', cohort == 'Hartwig') %>%
              { if (filter_sig) filter(., p.signif != 'ns') else . },
            aes(label = mean_diff_label, y = (total_y_limit - 2) * 1.3), nudge_x = 0.5) +
  geom_point(data = mean_av_sample_tier_label %>% filter(label == 'Total'), aes(y = log2(mean_av + 1))) +
  facet_wrap(cancer_type_label ~ ., ncol = 1, strip.position = 'left') +
  coord_flip(
    ylim = c(0, total_y_limit - 1.9),
    clip = 'off'
  ) +
  scale_y_continuous(
    breaks = c(0, log2(2 ^ seq(0, total_y_limit - 2, by = 1) + 1)), 
    labels = c(0, 2 ^ seq(0, total_y_limit - 2, by = 1))
  ) +
  scale_color_manual(
    values = rev(cohort_colors),
    guide = guide_legend(title = '', reverse = TRUE, override.aes = list(fill = cohort_colors))
    ) +
  labs(
    title = 'Total',
    subtitle = paste0(str_to_title(response_filter), 
                      ' actionable variants | q < ', q_val_cutoff),
    x = '',
    y = 'Actionable variants per sample',
    color = ''
  ) +
  theme_classic() +
  theme_smallfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = c(legend_x_pos, 0.5),
    panel.spacing = unit(0, 'cm'),
    plot.margin = margin(1,5,1,0, unit = 'cm')
    )

# A On-label plot
a_on_label_avps_plot <- av_per_sample_tier_label %>%
  filter(label == 'A_On_label') %>%
  ggplot(., aes(x = cohort, y = log2(n + 1))) +
  geom_rect(data = . %>% 
              filter(cohort == 'Hartwig') %>%
              filter(cancer_type %in% rect_idx) %>%
              distinct(cancer_type_label), 
            inherit.aes = FALSE,
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, 
            alpha = 0.3, fill = '#a0a0a0') +
  geom_jitter(color = 'gray40', size = 0.1, alpha = 0.25) +
  geom_violin(aes(color = cohort), fill = NA, position = position_dodge(1), size = 1.05, scale = 'width', width = 0.5, adjust = 2.5) +
  geom_text(data = . %>% 
              compare_means(
                data = ., formula = n ~ cohort,
                method = 'wilcox.test', group.by = 'cancer_type_label', p.adjust.method = 'fdr',
                symnum.args = list(cutpoints = c(0, q_val_cutoff, 1), symbols = c('*', 'ns'))) %>% 
              left_join(mean_av_sample_tier_label, by = 'cancer_type_label', 'group1' = 'cohort') %>%
              filter(label == 'A_On_label', cohort == 'Hartwig') %>%
              { if (filter_sig) filter(., p.signif != 'ns') else . },
            aes(label = mean_diff_label, y = (a_on_label_limit - 2) * 1.3), nudge_x = 0.5) +
  stat_summary(aes(group = cohort), fun=mean, geom="point", shape=16, size=2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  facet_wrap(cancer_type_label ~ ., ncol = 1, strip.position = 'left') +
  coord_flip(
    ylim = c(0, a_on_label_limit - 1.9),
    clip = 'off'
  ) +
  scale_y_continuous(
    breaks = c(0, log2(2 ^ seq(0, a_on_label_limit - 2, by = 1) + 1)), 
    labels = c(0, 2 ^ seq(0, a_on_label_limit - 2, by = 1))
  ) +
  scale_color_manual(
    values = rev(cohort_colors),
    guide = guide_legend(title = '', reverse = TRUE, override.aes = list(fill = cohort_colors))
    ) +
  labs(
    title = 'A On-label',
    subtitle = paste0(str_to_title(response_filter), 
                      ' actionable variants | q < ', q_val_cutoff),
    x = '',
    y = 'Actionable variants per sample',
    color = ''
  ) +
  theme_classic() +
  theme_smallfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = c(legend_x_pos, 0.5),
    panel.spacing = unit(0, 'cm'),
    plot.margin = margin(1,5,1,0, unit = 'cm')
    )

# A Off-label plot
a_off_label_avps_plot <- av_per_sample_tier_label %>%
  filter(label == 'A_Off_label') %>%
  ggplot(., aes(x = cohort, y = log2(n + 1))) +
  geom_rect(data = . %>% 
              filter(cohort == 'Hartwig') %>%
              filter(cancer_type %in% rect_idx) %>%
              distinct(cancer_type_label),
            inherit.aes = FALSE,
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, 
            alpha = 0.3, fill = '#a0a0a0') +
  geom_jitter(color = 'gray40', size = 0.1, alpha = 0.25) +
  geom_violin(aes(color = cohort), fill = NA, position = position_dodge(1), size = 1.05, scale = 'width', width = 0.5, adjust = 2.5) +
  geom_text(data = . %>% 
              compare_means(
                data = ., formula = n ~ cohort,
                method = 'wilcox.test', group.by = 'cancer_type_label', p.adjust.method = 'fdr',
                symnum.args = list(cutpoints = c(0, q_val_cutoff, 1), symbols = c('*', 'ns'))) %>% 
              left_join(mean_av_sample_tier_label, by = 'cancer_type_label', 'group1' = 'cohort') %>%
              filter(label == 'A_Off_label', cohort == 'Hartwig') %>%
              { if (filter_sig) filter(., p.signif != 'ns') else . },
            aes(label = mean_diff_label, y = (a_off_label_limit - 2) * 1.3), nudge_x = 0.5) +
  stat_summary(aes(group = cohort), fun=mean, geom="point", shape=16, size=2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  facet_wrap(cancer_type_label ~ ., ncol = 1, strip.position = 'left') +
  coord_flip(
    ylim = c(0, a_off_label_limit - 1.9),
    clip = 'off'
  ) +
  scale_y_continuous(
    breaks = c(0, log2(2 ^ seq(0, a_off_label_limit - 2, by = 1) + 1)), 
    labels = c(0, 2 ^ seq(0, a_off_label_limit - 2, by = 1))
  ) +
  scale_color_manual(
    values = rev(cohort_colors),
    guide = guide_legend(title = '', reverse = TRUE, override.aes = list(fill = cohort_colors))
    ) +
  labs(
    title = 'A Off-label',
    subtitle = paste0(str_to_title(response_filter), 
                      ' actionable variants | q < ', q_val_cutoff),
    x = '',
    y = 'Actionable variants per sample',
    color = ''
  ) +
  theme_classic() +
  theme_smallfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = c(legend_x_pos, 0.5),
    panel.spacing = unit(0, 'cm'),
    plot.margin = margin(1,5,1,0, unit = 'cm')
    )

# B On-label plot
b_on_label_avps_plot <- av_per_sample_tier_label %>%
  filter(label == 'B_On_label') %>%
  ggplot(., aes(x = cohort, y = log2(n + 1))) +
  geom_rect(data = . %>% 
              filter(cohort == 'Hartwig') %>%
              filter(cancer_type %in% rect_idx) %>%
              distinct(cancer_type_label),
            inherit.aes = FALSE,
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, 
            alpha = 0.3, fill = '#a0a0a0') +
  geom_jitter(color = 'gray40', size = 0.1, alpha = 0.25) +
  geom_violin(aes(color = cohort), fill = NA, position = position_dodge(1), size = 1.05, scale = 'width', width = 0.5, adjust = 2.5) +
  geom_text(data = . %>% 
              compare_means(
                data = ., formula = n ~ cohort,
                method = 'wilcox.test', group.by = 'cancer_type_label', p.adjust.method = 'fdr',
                symnum.args = list(cutpoints = c(0, q_val_cutoff, 1), symbols = c('*', 'ns'))) %>% 
              left_join(mean_av_sample_tier_label, by = 'cancer_type_label', 'group1' = 'cohort') %>%
              filter(label == 'B_On_label', cohort == 'Hartwig') %>%
              { if (filter_sig) filter(., p.signif != 'ns') else . },
            aes(label = mean_diff_label, y = (b_on_label_limit - 2) * 1.3), nudge_x = 0.5) +
  stat_summary(aes(group = cohort), fun=mean, geom="point", shape=16, size=2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  facet_wrap(cancer_type_label ~ ., ncol = 1, strip.position = 'left') +
  coord_flip(
    ylim = c(0, b_on_label_limit - 1.9),
    clip = 'off'
  ) +
  scale_y_continuous(
    breaks = c(0, log2(2 ^ seq(0, b_on_label_limit - 2, by = 1) + 1)), 
    labels = c(0, 2 ^ seq(0, b_on_label_limit - 2, by = 1))
  ) +
  scale_color_manual(
    values = rev(cohort_colors),
    guide = guide_legend(title = '', reverse = TRUE, override.aes = list(fill = cohort_colors))
    ) +
  labs(
    title = 'B On-label',
    subtitle = paste0(str_to_title(response_filter), 
                      ' actionable variants | q < ', q_val_cutoff),
    x = '',
    y = 'Actionable variants per sample',
    color = ''
  ) +
  theme_classic() +
  theme_smallfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = c(legend_x_pos, 0.5),
    panel.spacing = unit(0, 'cm'),
    plot.margin = margin(1,5,1,0, unit = 'cm')
    )

# B Off-label plot
b_off_label_avps_plot <- av_per_sample_tier_label %>%
  filter(label == 'B_Off_label') %>%
  ggplot(., aes(x = cohort, y = log2(n + 1))) +
  geom_rect(data = . %>% 
              filter(cohort == 'Hartwig') %>%
              filter(cancer_type %in% rect_idx) %>%
              distinct(cancer_type_label),
            inherit.aes = FALSE,
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, 
            alpha = 0.3, fill = '#a0a0a0') +
  geom_jitter(color = 'gray40', size = 0.1, alpha = 0.25) +
  geom_violin(aes(color = cohort), fill = NA, position = position_dodge(1), size = 1.05, scale = 'width', width = 0.5, adjust = 2.5) +
  geom_text(data = . %>% 
              compare_means(
                data = ., formula = n ~ cohort,
                method = 'wilcox.test', group.by = 'cancer_type_label', p.adjust.method = 'fdr',
                symnum.args = list(cutpoints = c(0, q_val_cutoff, 1), symbols = c('*', 'ns'))) %>% 
              left_join(mean_av_sample_tier_label, by = 'cancer_type_label', 'group1' = 'cohort') %>%
              filter(label == 'B_Off_label', cohort == 'Hartwig') %>%
              { if (filter_sig) filter(., p.signif != 'ns') else . },
            aes(label = mean_diff_label, y = (b_off_label_limit - 2) * 1.3), nudge_x = 0.5) +
  stat_summary(aes(group = cohort), fun=mean, geom="point", shape=16, size=2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  facet_wrap(cancer_type_label ~ ., ncol = 1, strip.position = 'left') +
  coord_flip(
    ylim = c(0, b_off_label_limit - 1.9),
    clip = 'off'
  ) + # coord_flip() needs to be before scale_y_continuous()!
  scale_y_continuous(
    breaks = c(0, log2(2 ^ seq(0, b_off_label_limit - 2, by = 1) + 1)), 
    labels = c(0, 2 ^ seq(0, b_off_label_limit - 2, by = 1))) +
  scale_color_manual(
    values = rev(cohort_colors),
    guide = guide_legend(title = '', reverse = TRUE, override.aes = list(fill = cohort_colors))
    ) +
  labs(
    title = 'B Off-label',
    subtitle = paste0(str_to_title(response_filter), 
                      ' actionable variants | q < ', q_val_cutoff),
    x = '',
    y = 'Actionable variants per sample',
    color = ''
  ) +
  theme_classic() +
  theme_smallfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = c(legend_x_pos, 0.5),
    panel.spacing = unit(0, 'cm'),
    plot.margin = margin(1,5,1,0, unit = 'cm')
    )

# save plots
# as PDF
pdf(
  file = paste0(output_dir, '/', current_date, '_', response_filter, '_total_avps_noleuk',
                if (by_subtype) { if (filter_sig) { '_by_subtype_sig.pdf' } else { '_by_subtype_all.pdf' } } 
                else { if (filter_sig) { '_sig.pdf' } else { '_all.pdf' } }),
  width = 6,
  height = 9,
  useDingbats = FALSE,
  compress = FALSE
)
print(total_avps_plot)
dev.off()

pdf(
  file = paste0(output_dir, '/', current_date, '_', response_filter, '_a_on_label_avps_noleuk',
                if (by_subtype) { if (filter_sig) { '_by_subtype_sig.pdf' } else { '_by_subtype_all.pdf' } } 
                else { if (filter_sig) { '_sig.pdf' } else { '_all.pdf' } }),
  width = 6,
  height = 9,
  useDingbats = FALSE,
  compress = FALSE
)
print(a_on_label_avps_plot)
dev.off()

pdf(
  file = paste0(output_dir, '/', current_date, '_', response_filter, '_a_off_label_avps_noleuk',
                if (by_subtype) { if (filter_sig) { '_by_subtype_sig.pdf' } else { '_by_subtype_all.pdf' } } 
                else { if (filter_sig) { '_sig.pdf' } else { '_all.pdf' } }),
  width = 6,
  height = 9,
  useDingbats = FALSE,
  compress = FALSE
)
print(a_off_label_avps_plot)
dev.off()

pdf(
  file = paste0(output_dir, '/', current_date, '_', response_filter, '_b_on_label_avps_noleuk',
                if (by_subtype) { if (filter_sig) { '_by_subtype_sig.pdf' } else { '_by_subtype_all.pdf' } } 
                else { if (filter_sig) { '_sig.pdf' } else { '_all.pdf' } }),
  width = 6,
  height = 9,
  useDingbats = FALSE,
  compress = FALSE
)
print(b_on_label_avps_plot)
dev.off()

pdf(
  file = paste0(output_dir, '/', current_date, '_', response_filter, '_b_off_label_avps_noleuk',
                if (by_subtype) { if (filter_sig) { '_by_subtype_sig.pdf' } else { '_by_subtype_all.pdf' } } 
                else { if (filter_sig) { '_sig.pdf' } else { '_all.pdf' } }),
  width = 6,
  height = 9,
  useDingbats = FALSE,
  compress = FALSE
)
print(b_off_label_avps_plot)
dev.off()

##################### Which actionable variants contribute to the significant increase in the different tier levels?

# read in treatment information to include as a label in the figure
treatment <- read_tsv(file = paste0(base_dir, '/results/08_actionability/', actionable_variants_results_date, '/results/hartwig_actionability.tsv')) %>%
  # remove liquid tumors because they are weird
  filter(!str_detect(cancerType, pattern = 'leukemia|CLL')) %>%
  # remove variants for which we dont have data for (expression and hypermethylation)
  filter(!mutated_variant %in% c('over exp', 'dec exp', 'hypermethylation')) %>%
  # filter only for specified response type (either 'sensitive' or 'resistant')
  filter(responseType %in% response_type_filter) %>%
  # select relevant columns
  select(ckbEntryId, gene, mutated_variant, cancerType, ampCapAscoEvidenceLevel, treatment) %>%
  distinct()

# combine treatment information with our parse variants from our cohort
treatment <- actionable_variants_cleantt %>%
  inner_join(treatment, by = c('ckbEntryId', 'gene', 'mutated_variant', 'cancerType', 'ampCapAscoEvidenceLevel')) %>%
  # create the tier level and gene variant columns again by uniting two cols
  mutate(tier_level = paste0(ampCapAscoEvidenceLevel, '_', on_off_label)) %>%
  mutate(gene_variant = paste0(gene, '_', mutated_variant)) %>%
  distinct(cancer_type_label, cancer_type_code, tier_level, gene_variant, treatment) %>%
  # collapse the different treatments that are available for a particular variant into one row separate by a line break
  group_by(cancer_type_label, tier_level, gene_variant) %>%
  summarise(treatment_label = str_c(treatment, collapse = '\n'), .groups = 'drop')

# custom functions to get the frequencies of different tier level AVs per cancer type
get_av_frequencies <- function(cancer_type_input, tier_level_input, on_off_label_input = 'On_label') {
  
  output_table <- actionable_variants_cleantt %>%
    filter(cancer_type == cancer_type_input,
           ampCapAscoEvidenceLevel == tier_level_input,
           on_off_label == on_off_label_input) %>%
    distinct(sample_id, cancer_type, cohort, gene, mutated_variant) %>%
    group_by(cancer_type, cohort, gene, mutated_variant) %>%
    count() %>%
    ungroup() %>%
    unite(., col ='gene_variant', c(gene, mutated_variant)) %>%
    complete(cancer_type, cohort, gene_variant, fill = list(n = 0)) %>%
    arrange(gene_variant) %>%
    inner_join(sample_size, by = c('cancer_type', 'cohort')) %>%
    mutate(freq = n / sample_size * 100) %>%
    mutate(tier_level = paste0(tier_level_input, '_', on_off_label_input))

  return(output_table)
}

# define cancer type iterator
a_on_label_cancer_types <- all_cancers
a_off_label_cancer_types <- all_cancers
b_on_label_cancer_types <- all_cancers
b_off_label_cancer_types <- all_cancers

# get the AVs for the specified 
a_on_label_variants <- map(a_on_label_cancer_types, ~get_av_frequencies(.x, tier_level_input = 'A', on_off_label_input = 'On_label'))
a_off_label_variants <- map(a_off_label_cancer_types, ~get_av_frequencies(.x, tier_level_input = 'A', on_off_label_input = 'Off_label'))
b_on_label_variants <- map(b_on_label_cancer_types, ~get_av_frequencies(.x, tier_level_input = 'B', on_off_label_input = 'On_label'))
b_off_label_variants <- map(b_off_label_cancer_types, ~get_av_frequencies(.x, tier_level_input = 'B', on_off_label_input = 'Off_label'))

# concatenate sublists into one
variant_freq_list <- c(
  a_on_label_variants,
  a_off_label_variants,
  b_on_label_variants,
  b_off_label_variants
)

# transform list of dfs to one large df
variant_freq_df <- do.call(rbind, variant_freq_list)

# make table wider and calculate the frequency difference between cohorts
variant_freq_df2 <- variant_freq_df %>%
  pivot_wider(., names_from = 'cohort', values_from = c('n', 'sample_size', 'freq')) %>%
  mutate(freq_diff = freq_Hartwig - freq_PCAWG) %>%
  mutate(n_non_PCAWG = sample_size_PCAWG - n_PCAWG,
         n_non_Hartwig = sample_size_Hartwig - n_Hartwig) %>%
  arrange(desc(freq_diff))

# filter out cases where there are less than 4 mutated samples in both primary or metastatic group
min_freq_idx <- pmax(variant_freq_df2$n_PCAWG, variant_freq_df2$n_Hartwig) > 4
variant_minfreq_df <- variant_freq_df2 %>%
  filter(min_freq_idx)

### calculate fisher's exact test p-value
# create fisher matrix
variant_freq_mat <- variant_minfreq_df %>%
  mutate(cancer_type_gene_variant = paste0(cancer_type_code, '_', tier_level, '_', gene_variant)) %>%
  select(cancer_type_gene_variant, n_PCAWG, n_non_PCAWG, n_Hartwig, n_non_Hartwig) %>%
  column_to_rownames(., var = 'cancer_type_gene_variant') %>%
  as.matrix()

# calculate exact test p-value
variant_minfreq_df$fisher_pval <- fisherTest.matrix(variant_freq_mat)

# adjust p-value
variant_minfreq_df <- variant_minfreq_df %>%
  group_by(cancer_type) %>%
  mutate(fisher_pval_adj = p.adjust(fisher_pval, method = 'fdr')) %>%
  ungroup() %>%
  mutate(significance_level = if_else(fisher_pval_adj < q_val_cutoff, '*', NA_character_))

# prepare ylim for plot
pcawg_freq_limit <- -100
hartwig_freq_limit <- 100

# get maintype code
cancer_maintype_code <- metadata %>%
  distinct(cancer_type_code, cancer_maintype_code)

# arrange data for plotting
variant_freq_df_final <- variant_minfreq_df %>%
  filter(pmax(freq_Hartwig, freq_PCAWG) > 5) %>%
  filter(abs(freq_diff) > 5) %>%
  semi_join(global_contingency_df %>% filter(global_significance != 'ns'), by = 'cancer_type') %>%
  select(-n_Hartwig, -n_PCAWG, -sample_size_Hartwig, -sample_size_PCAWG) %>%
  pivot_longer(., cols = c(freq_Hartwig, freq_PCAWG), names_to = 'cohort', names_prefix = 'freq_', values_to = 'freq') %>%
  mutate(freq = if_else(cohort == 'PCAWG', freq * -1, freq)) %>%
  { if (by_subtype) { inner_join(., cancer_maintype_code, by = 'cancer_type_code') %>%
      mutate(., cancer_maintype_code = factor(cancer_type_code, levels = cancer_type_code_order)) } 
    else { mutate(., cancer_type_code = factor(cancer_type_code, levels = cancer_type_code_order)) }} %>%
  mutate(cancer_type_gene_variant = paste0(cancer_type_code, ': ', str_replace(gene_variant, pattern = '_', replacement = ' '))) %>%
  mutate(tier_level = factor(tier_level, levels = tier_label_levels)) %>%
  arrange(factor(cancer_type_label, levels = cancer_type_label_order), cancer_type_gene_variant)

if (by_subtype) {
  
  tmp1 <- subtypes_plot %>%
    inner_join(variant_freq_df, by = c('cancer_type', 'cancer_type_code'))
  
  tmp2 <- variant_freq_df_final %>%
    filter(fisher_pval_adj < q_val_cutoff) %>%
    inner_join(subtypes_plot, by = c('cancer_type', 'cancer_type_code')) %>%
    select(cancer_maintype, cancer_type, gene_variant, tier_level, significance_level)
  
  variant_freq_df_finaltt <- tmp1 %>%
    semi_join(tmp2, by = c('cancer_maintype', 'gene_variant', 'tier_level'))
  
  subtype_gene_variant_sample_size <- variant_freq_df_finaltt %>%
    select(cancer_type_code, cohort, gene_variant, n, sample_size) %>%
    distinct() %>%
    pivot_wider(., names_from = 'cohort', values_from = c('n', 'sample_size')) %>%
    mutate(cancer_type_gene_variant = paste0(str_replace(cancer_type_code, pattern = regex('\\_'), replacement = ' '), ': ', str_replace(gene_variant, pattern = regex('\\_'), replacement = ' '), ' (', n_PCAWG, '/', sample_size_PCAWG, ' vs. ', n_Hartwig, '/', sample_size_Hartwig, ')')) %>%
    select(cancer_type_code, gene_variant, cancer_type_gene_variant)
  
  variant_freq_df_final <- variant_freq_df_finaltt %>%
    left_join(tmp2 %>% select(-cancer_maintype), by = c('cancer_type', 'gene_variant', 'tier_level')) %>%
    inner_join(subtype_gene_variant_sample_size, by = c('cancer_type_code', 'gene_variant')) %>%
    distinct() %>%
    mutate(freq = if_else(cohort == 'PCAWG', freq * -1, freq))
  
}

if (by_subtype) {
  cancer_type_gene_variant_order <- variant_freq_df_final %>%
    arrange(cancer_maintype, gene_variant) %>%
    pull(cancer_type_gene_variant) %>%
    unique()
  
} else {
  # get the order of the cancer type gene variants to match general cancer type order
  cancer_type_gene_variant_order <- unique(variant_freq_df_final$cancer_type_gene_variant)
  
}

# gene variants slice the highest 
variant_freq_df_final <- variant_freq_df_final %>%
  mutate(tier_level = factor(tier_level, levels = tier_label_levels)) %>%
  { if (by_subtype) { group_by(., cancer_maintype, gene_variant) } else { group_by(., cancer_type, gene_variant) }} %>%
  slice_min(., order_by = tier_level, n = 1) %>%
  ungroup()

variant_freq_df_final <- variant_freq_df_final %>%
  # impose order on the x axis variable
  mutate(cancer_type_gene_variant = factor(cancer_type_gene_variant, levels = rev(cancer_type_gene_variant_order)))

# prepare treatment labels for plotting
treatment <- treatment %>%
  inner_join(variant_freq_df_final, by = c('cancer_type_label', 'tier_level', 'gene_variant')) %>%
  filter(cohort == 'Hartwig') %>%
  mutate(tier_level = factor(tier_level, levels = tier_label_levels))

# plot AV frequencies
av_frequency_plot <- variant_freq_df_final %>%
  ggplot(., aes(x = cancer_type_gene_variant, y = freq)) +
  geom_rect(data = . %>%
              distinct(tier_level),
            aes(fill = tier_level),
            inherit.aes = FALSE,
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, 
            alpha = 0.7) +
  scale_fill_manual(
    values = label_colors,
    labels = c('A On-label', 'A Off-label', 'B On-label', 'B Off-label'),
    guide = guide_legend(
      title = '',
      override.aes = list(alpha = 1),
      order = 1)
    ) +
  new_scale_fill() +
  geom_bar(aes(fill = cohort), stat = 'identity', position = 'identity', color = 'black') +
  scale_fill_manual(
    values = rev(cohort_colors),
    guide = guide_legend(title = '',
                         order = 2)
    ) +
  geom_text(aes(y = if_else(cohort == 'PCAWG', freq - 8, freq + 8), label = paste0(abs(round(freq, 1)), '%')), size = 4) +
  geom_text(data = . %>% filter(cohort == 'Hartwig'), aes(y = if_else(cohort == 'PCAWG', freq - 20, freq + 20), label = significance_level),
            nudge_x = -0.25, size = 7.5) +
  geom_label(data = treatment %>% filter(!is.na(significance_level)), aes(y = 90, label = treatment_label),
             size = 4, 
             ) +
  facet_grid(tier_level ~ .,
             scales = 'free_y', space = 'free_y') +
  coord_flip(
    ylim = c(pcawg_freq_limit - 5, hartwig_freq_limit + 5),
    clip = 'off'
    ) +
  scale_y_continuous(
    breaks = seq(pcawg_freq_limit, hartwig_freq_limit, by = 10),
    labels = c(seq(hartwig_freq_limit, 0, by = -10), seq(10, hartwig_freq_limit, by = 10))
  ) +
  labs(
    x= '',
    y = 'Frequency [%]'
  ) +
  theme_bw() +
  theme_bigfont_barplot +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000', size = 0.5),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    plot.margin = margin(0,25,0,0, unit = 'points')
    )

# save the plot
# as PDF
pdf(
  file = paste0(output_dir, '/', current_date, '_', response_filter, '_av_frequency_q', q_val_cutoff, '_noleuk',
                if (by_subtype) { '_subtype.pdf' } else { '.pdf' }),
  width = if (by_subtype) { 14 } else { 14 },
  height = if (by_subtype) { 8 } else { 10 },
  useDingbats = FALSE,
  compress = FALSE
)
print(av_frequency_plot)
dev.off()