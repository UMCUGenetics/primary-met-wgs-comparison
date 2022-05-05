############### Plot actionable variants ###############
# author: remy (sascha)
# date: 27/01/2022
# last updated: 06/04/2022

### Description
# input: actionable variants table that was parsed from VCF files by parse_actionable_variants.R script.
# This script plots the barplots including on and off label actionable variants. On-label variants are approved or have the potential to be approved
# for a certain cancer type. Off-label variants were approved or have the potential to be approved targets in other cancer types/

# libraries
library(tidyverse) # data manipulation and plotting
library(ggpubr) # compare means  while plotting
library(ggrepel) # use repelling labels in plots
library(ggnewscale) # second color scales
library(naniar) # replace with NA function
library(svglite) # handle inkscape SVG issue
library(googlesheets4) # push supp tables to gsheets

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

output_dir <- list(
  local='/output/path',
  umc='/output/path'
)

for(i in output_dir){
  if(dir.exists(i)){
    output_dir <- i
    break
  }
}

current_date <- format(Sys.Date(), '%Y_%m_%d')

## ggplot custom themes
theme_bigfont <- theme(plot.title = element_text(size=22),
                       plot.subtitle = element_text(size = 16),
                       axis.text.x = element_text(size=15),
                       axis.text.y = element_text(size=15), 
                       axis.title = element_text(size=18),
                       legend.text = element_text(size = 14),
                       strip.text.x = element_text(size = 15))

theme_smallfont <- theme(plot.title = element_text(size=17),
                         plot.subtitle = element_text(size = 11),
                         axis.text.x= element_text(size=10),
                         axis.text.y= element_text(size=10), 
                         axis.title=element_text(size=13),
                         legend.text = element_text(size = 12),
                         strip.text.x = element_text(size = 10))

# define dates of the input files and other variables
actionable_variants_results_date <- '2022_01_26' # old: 2022_01_20 new: 2022_01_26
linx_drivers_resutls_date <- '2021_12_03'
fusion_results_date <- '2021_12_03'
purity_results_date <- '2021_12_03'
response_filter <- 'sensitive' # can be 'sensitive' or 'resistant'
q_val_cutoff <- 0.05
p_val_cutoff <- 0.01


# source the functions to calculate the fishers exact test and wilcoxon test p-values and effect sizes
# from test matrices
source(file = paste0(base_dir, '/path/to/0_func/fisher_wilcoxon_matrix_functions.R'))

# ------------------------------- Custom functions

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
    # filter(cancer_type == 'Kidney clear cell carcinoma') %>%
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

# --------------------------------------------- ACTIONABILITY
if (response_filter == 'sensitive') {
  response_type_filter <- c('sensitive', 'predicted - sensitive')
} else {
  response_type_filter <- c('resistant', 'predicted - resistant')
}

actionability <- read_tsv(file = paste0(base_dir, '/path/to/hartwig_actionability.tsv')) %>%
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
  # (e.g. If there is a level A and B entry for EGFR3 mut in NSCLC, then only level A is kept)
  group_by(ckbEntryId, cancerType, tumorProfile) %>%
  slice_min(ampCapAscoEvidenceLevel, n = 1) %>%
  ungroup()

# --------------------------------------------- METADATA

metadata <- read_tsv(paste0(base_dir, '/path/to/metadata.tsv')) %>%
  # select relevant columns
  select(sample_id, cohort, cancer_type, cancer_type_code) %>%
  # rename 'Hartwig' to 'Hartwig'
  mutate(cohort = if_else(cohort == 'HMF', 'Hartwig', cohort))

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

# --------------------------------------------- ACTIONABLE VARIANTS
actionable_variants_clean <- read_tsv(file = paste0(base_dir, '/path/to/',
                                                    actionable_variants_results_date ,'/actionable_variants.tsv')) %>%
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
  mutate(on_off_label = if_else(cancer_type == cancerType, 'On_label', 'Off_label')) %>%
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

################ complete the AV table to include all sample - label combinations
complete_cohort <- metadata %>% select(-cancer_type_code)

complete_tier_label_table <- actionable_variants_cleantt %>%
  group_by(sample_id, cohort, cancer_type, ampCapAscoEvidenceLevel, on_off_label) %>%
  count() %>%
  ungroup() %>%
  right_join(complete_cohort, by = c('sample_id', 'cohort', 'cancer_type')) %>%
  tidyr::replace_na(., replace = list(ampCapAscoEvidenceLevel = 'A', on_off_label = 'Off_label', n = 0)) %>%
  complete(sample_id, ampCapAscoEvidenceLevel, on_off_label, fill = list(n = 0)) %>%
  select(-cohort, -cancer_type) %>%
  inner_join(metadata, by = 'sample_id') %>%
  relocate(., n, .after = last_col()) %>%
  unite(., col = 'label', ampCapAscoEvidenceLevel, on_off_label)

################ calculate fisher p-value per cancer type and label between cohorts

# get all cancer type names and labels
all_labels <- unique(complete_tier_label_table$label)
all_cancers <- unique(sample_size$cancer_type)

### global
# initiate list
global_contingency_list <- list()

# calculate Fishers exact test p-value, FDR adjusted p-value (per cancer type, per driver mutation type),
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

# filter out cases where there are less than 4 mutated samples in both primary or metastatic group
min_freq_idx <- pmax(global_contingency_df$n_PCAWG, global_contingency_df$n_Hartwig) > 3
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
  mutate(global_significance = if_else(fisher_pval_adj < q_val_cutoff, '*', NA_character_),
         cohort = 'Hartwig') %>%
  mutate(prop_Hartwig = n_Hartwig / (n_Hartwig + n_non_Hartwig),
         prop_PCAWG = n_PCAWG / (n_PCAWG + n_non_PCAWG)) %>%
  mutate(fold_change = round((prop_Hartwig - prop_PCAWG) / prop_PCAWG, 1),
         fold_change = if_else(is.na(global_significance), NA_real_, fold_change),
         fold_change_label = if_else(fold_change > 0, str_c('+', fold_change), as.character(fold_change))) %>%
  mutate(abs_diff = round((prop_Hartwig - prop_PCAWG), 1),
         abs_diff = if_else(is.na(global_significance), NA_real_, abs_diff),
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

# filter out cases where there are less than 4 mutated samples in both primary or metastatic group
min_freq_idx <- pmax(contingency_df$n_PCAWG, contingency_df$n_Hartwig) > 3
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

## ggplot custom theme
resize_factor <- 1.2
theme_bigfont_barplot <- theme(plot.title = element_text(size=16 * resize_factor),
                               plot.subtitle = element_text(size = 14 * resize_factor),
                               axis.text.x= element_text(size=12 * resize_factor),
                               axis.text.y= element_text(size=11 * resize_factor), 
                               axis.title=element_text(size=18 * resize_factor),
                               legend.title = element_text(size = 16 * resize_factor),
                               legend.text = element_text(size = 14 * resize_factor),
                               strip.text.x = element_text(size = 11 * resize_factor),
                               strip.text.y = element_text(size = 11 * resize_factor))

# define cancer type order
source(paste0(base_dir, '/path/to/r_objects/cancer_type_order.R'))
cancer_type_order <- c('Pancancer', cancer_type_order)

# define ASCO evidence tier level colors
label_colors <- c('#1e6220', '#2F913F', '#35E358', '#C4FFD0')
names(label_colors) <- c('A_On_label', 'A_Off_label', 'B_On_label', 'B_Off_label')

# get the correct cancer type label order
cancer_type_label_order <- sample_size %>%
  arrange(factor(cancer_type, levels = cancer_type_order)) %>%
  pull(cancer_type_label) %>%
  unique()

# impose order on cancer type label & tier level, filter NA level
actionable_variants_clean_plottt <- actionable_variants_clean_plot %>%
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
  mutate(fold_change = round((prop_Hartwig - prop_PCAWG) / prop_PCAWG, 1)) %>%
  mutate(abs_diff = round((prop_Hartwig - prop_PCAWG), 1)) %>%
  filter(fold_change != Inf) %>%
  group_by(label) %>%
  summarise(mean_fold_change = mean(fold_change),
            sd_fold_change = sd(fold_change),
            mean_abs_diff = mean(abs_diff),
            sd_abs_diff = sd(abs_diff))

# # get the alternating index number for the grey background rectangles
rect_idx <- seq(2, 23, by = 2)
rect_idx_vec <- unique(actionable_variants_clean_plottt$cancer_type)[order(match(unique(actionable_variants_clean_plottt$cancer_type), cancer_type_order))]
rect_idx <- rect_idx_vec[rect_idx]

# plot stacked barplot: Proportion of A & B On/Off-label AVs, facetted by cancer type
# export: 900 x 900
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
    data = actionable_variants_clean_plottt %>%
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
    subtitle = paste0('global p < ', p_val_cutoff),
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

ggsave(
  path = output_dir,
  filename = paste0(current_date, '_', response_filter, '_av_barplot_noleuk.png'),
  plot = global_av_plot,
  device = 'png',
  width = 10,
  height = 10,
  dpi = 600
)

pdf(
  file = paste0(output_dir, '/', current_date, '_', response_filter, '_av_barplot_noleuk.pdf'),
  width = 10,
  height = 10,
  useDingbats = FALSE,
  compress = FALSE
)
print(global_av_plot)
dev.off()

# save the plot as a SVG file
svglite(
  filename = paste0(output_dir, '/', current_date, '_', response_filter, '_av_barplot_noleuk.svg'),
  width = 10,
  height = 10
)
global_av_plot
dev.off()

################################# Tier level events per sample (Total, A & B On/Off-label)

################ avg number of actionable variants per sample by Tier level and On/Off-label
# calculate total number of AVs per sample + add this total to the original table
av_per_sample_tier_label <- complete_tier_label_table %>%
  group_by(sample_id, cohort, cancer_type, cancer_type_code) %>%
  summarise(n = sum(n), .groups = 'drop') %>%
  mutate(label = 'Total')

av_per_sample_tier_label <- complete_tier_label_table %>%
  union(., av_per_sample_tier_label)

############# calculate mean of mean AVs per sample to get an idea on pancancer level
av_per_sample_tier_label_mean_of_mean <- av_per_sample_tier_label %>%
  # filter(label == 'Total') %>%
  group_by(cancer_type, cohort, label) %>%
  summarise(mean_av = mean(n), .groups = 'drop') %>%
  group_by(cohort, label) %>%
  summarise(mean_av = round(mean(mean_av), 1), .groups = 'drop') %>%
  pivot_wider(., names_from = 'cohort', values_from = 'mean_av') %>%
  mutate(mean_diff = round(Hartwig - PCAWG, 2)) %>%
  mutate(mean_diff_label = if_else(mean_diff > 0, paste0('+ ', mean_diff), paste0('- ', abs(mean_diff)))) %>%
  pivot_longer(., cols = c('Hartwig', 'PCAWG'), names_to = 'cohort', values_to = 'mean_av')
############

# duplicate the data, but add 'Pancancer' as the cancer type, then add this duplicated table to the original table
pancan_av_per_sample_tier_label <- av_per_sample_tier_label %>%
  mutate(cancer_type = 'Pancancer', cancer_type_code = 'Pancancer')

av_per_sample_tier_label <- av_per_sample_tier_label %>%
  union(., pancan_av_per_sample_tier_label)

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

sig_total <- c(
  'Pancancer\n(1916 vs. 3835)',
  'Breast cancer\n(201 vs. 787)',
  'Hepatocellular carcinoma\n(288 vs. 53)',
  'Pancreas neuroendocrine\n(81 vs. 38)',
  'Ovarian cancer\n(110 vs. 168)',
  'Kidney clear cell carcinoma\n(109 vs. 129)', 
  'Non small cell lung cancer\n(83 vs. 511)',
  'Prostate carcinoma\n(153 vs. 404)'
)

sig_a_on_label <- c(
  'Pancancer\n(1916 vs. 3835)',
  'Breast cancer\n(201 vs. 787)',
  'Non small cell lung cancer\n(83 vs. 511)'
)

sig_a_off_label <- c(
  'Pancancer\n(1916 vs. 3835)',
  'Pancreas neuroendocrine\n(81 vs. 38)',
  'Kidney clear cell carcinoma\n(109 vs. 129)',
  'Prostate carcinoma\n(153 vs. 404)'
)

sig_b_on_label <- c(
  'Pancancer\n(1916 vs. 3835)',
  'Breast cancer\n(201 vs. 787)',
  'Non small cell lung cancer\n(83 vs. 511)',
  'Prostate carcinoma\n(153 vs. 404)',
  'Liposarcoma\n(17 vs. 25)'
)

sig_b_off_label <- c(
  'Pancancer\n(1916 vs. 3835)',
  'Breast cancer\n(201 vs. 787)',
  'Hepatocellular carcinoma\n(288 vs. 53)',
  'Pancreas neuroendocrine\n(81 vs. 38)',
  'Ovarian cancer\n(110 vs. 168)',
  'Kidney clear cell carcinoma\n(109 vs. 129)',
  'Prostate carcinoma\n(153 vs. 404)'
)

if (response_filter == 'sensitive') {
  total_y_limit <- 9
  a_on_label_limit <- 6
  a_off_label_limit <- 8
  b_on_label_limit <- 6
  b_off_label_limit <- 8
} else {
  total_y_limit <- 5
  a_on_label_limit <- 5
  a_off_label_limit <- 5
  b_on_label_limit <- 5
  b_off_label_limit <- 4
}

# get the alternating index number for the grey background rectangles
rect_idx <- seq(2, 23, by = 2)
rect_idx <- cancer_type_order[rect_idx]

# import cohort color
source(paste0(base_dir, '/analysis/linx_driver/r_objects/color_codes.R'))

# label_filter <- 'B_Off_label'

# av per sample violin plots
# total AVs plot: export 600 x 800
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
  geom_text(data = mean_av_sample_tier_label %>% filter(cohort == 'Hartwig', 
                                                        cancer_type_label %in% sig_total,
                                                        label == 'Total'),
            aes(label = mean_diff_label, y = (total_y_limit - 2) * 1.3), nudge_x = 0.5) +
  geom_point(data = mean_av_sample_tier_label %>% filter(label == 'Total'), aes(y = log2(mean_av + 1))) +
  stat_compare_means(aes(group = cohort), method = 'wilcox.test', label.y = (total_y_limit - 2) * 1.1,
                     vjust = 0, size = 4, label = 'p.signif', hide.ns = TRUE,
                     symnum.args = list(cutpoints = c(0, p_val_cutoff, 1), symbols = c('*', 'ns'))) +
  facet_wrap(cancer_type_label ~ ., ncol = 1, strip.position = 'left') +
  coord_flip(
    ylim = c(0, total_y_limit - 1.9),
    clip = 'off'
  ) +
  scale_y_continuous(
    limits = c(0, Inf),
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
                      ' actionable variants | p < ', p_val_cutoff),
    x = '',
    y = 'Actionable variants per sample',
    color = ''
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = 'right',
    panel.spacing = unit(0, 'cm'),
    plot.margin = unit(c(1,0,1,0), 'cm'))

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
  geom_text(data = mean_av_sample_tier_label %>% filter(cohort == 'Hartwig', 
                                                        cancer_type_label %in% sig_a_on_label,
                                                        label == 'A_On_label'),
            aes(label = mean_diff_label, y = (a_on_label_limit - 2) * 1.3), nudge_x = 0.5) +
  stat_summary(aes(group = cohort), fun=mean, geom="point", shape=16, size=2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  stat_compare_means(aes(group = cohort), method = 'wilcox.test', label.y = (a_on_label_limit - 2) * 1.1,
                     vjust = 0, size = 4, label = 'p.signif', hide.ns = TRUE,
                     symnum.args = list(cutpoints = c(0, p_val_cutoff, 1), symbols = c('*', 'ns'))) +
  facet_wrap(cancer_type_label ~ ., ncol = 1, strip.position = 'left') +
  coord_flip(
    ylim = c(0, a_on_label_limit - 1.9),
    clip = 'off'
  ) +
  scale_y_continuous(
    limits = c(0, Inf), 
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
                      ' actionable variants | p < ', p_val_cutoff),
    x = '',
    y = 'Actionable variants per sample',
    color = ''
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = 'right',
    panel.spacing = unit(0, 'cm'),
    plot.margin = unit(c(1,0,1,0), 'cm'))

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
  geom_text(data = mean_av_sample_tier_label %>% filter(cohort == 'Hartwig', 
                                                        cancer_type_label %in% sig_a_off_label,
                                                        label == 'A_Off_label'),
            aes(label = mean_diff_label, y = (a_off_label_limit - 2) * 1.3), nudge_x = 0.5) +
  stat_summary(aes(group = cohort), fun=mean, geom="point", shape=16, size=2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  stat_compare_means(aes(group = cohort), method = 'wilcox.test', label.y = (a_off_label_limit - 2) * 1.1,
                     vjust = 0, size = 4, label = 'p.signif', hide.ns = TRUE,
                     symnum.args = list(cutpoints = c(0, p_val_cutoff, 1), symbols = c('*', 'ns'))) +
  facet_wrap(cancer_type_label ~ ., ncol = 1, strip.position = 'left') +
  coord_flip(
    ylim = c(0, a_off_label_limit - 1.9),
    clip = 'off'
  ) +
  scale_y_continuous(
    limits = c(0, Inf), 
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
                      ' actionable variants | p < ', p_val_cutoff),
    x = '',
    y = 'Actionable variants per sample',
    color = ''
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = 'right',
    panel.spacing = unit(0, 'cm'),
    plot.margin = unit(c(1,0,1,0), 'cm'))

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
  geom_text(data = mean_av_sample_tier_label %>% filter(cohort == 'Hartwig', 
                                                        cancer_type_label %in% sig_b_on_label,
                                                        label == 'B_On_label'),
            aes(label = mean_diff_label, y = (b_on_label_limit - 2) * 1.3), nudge_x = 0.5) +
  stat_summary(aes(group = cohort), fun=mean, geom="point", shape=16, size=2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  stat_compare_means(aes(group = cohort), method = 'wilcox.test', label.y = (b_on_label_limit - 2) * 1.1,
                     vjust = 0, size = 4, label = 'p.signif', hide.ns = TRUE,
                     symnum.args = list(cutpoints = c(0, p_val_cutoff, 1), symbols = c('*', 'ns'))) +
  facet_wrap(cancer_type_label ~ ., ncol = 1, strip.position = 'left') +
  coord_flip(
    ylim = c(0, b_on_label_limit - 1.9),
    clip = 'off'
  ) +
  scale_y_continuous(
    limits = c(0, Inf),
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
                      ' actionable variants | p < ', p_val_cutoff),
    x = '',
    y = 'Actionable variants per sample',
    color = ''
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = 'right',
    panel.spacing = unit(0, 'cm'),
    plot.margin = unit(c(1,0,1,0), 'cm'))

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
  geom_text(data = mean_av_sample_tier_label %>% filter(cohort == 'Hartwig',
                                                        cancer_type_label %in% sig_b_off_label,
                                                        label == 'B_Off_label'),
            aes(label = mean_diff_label, y = (b_off_label_limit - 2) * 1.3), nudge_x = 0.5) +
  stat_summary(aes(group = cohort), fun=mean, geom="point", shape=16, size=2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  stat_compare_means(aes(group = cohort), method = 'wilcox.test', label.y = (b_off_label_limit - 2) * 1.1,
                     vjust = 0, size = 4, label = 'p.signif', hide.ns = TRUE,
                     symnum.args = list(cutpoints = c(0, p_val_cutoff, 1), symbols = c('*', 'ns'))) +
  facet_wrap(cancer_type_label ~ ., ncol = 1, strip.position = 'left') +
  coord_flip(
    ylim = c(0, b_off_label_limit - 1.9),
    clip = 'off'
  ) + # coord_flip() needs to be before scale_y_continuous()!
  scale_y_continuous(
    limits = c(0, b_off_label_limit), # needed so that the jitter points are not plotted below 0
    breaks = c(0, log2(2 ^ seq(0, b_off_label_limit - 2, by = 1) + 1)), 
    labels = c(0, 2 ^ seq(0, b_off_label_limit - 2, by = 1))) +
  scale_color_manual(
    values = rev(cohort_colors),
    guide = guide_legend(title = '', reverse = TRUE, override.aes = list(fill = cohort_colors))
    ) +
  labs(
    title = 'B Off-label',
    subtitle = paste0(str_to_title(response_filter), 
                      ' actionable variants | p < ', p_val_cutoff),
    x = '',
    y = 'Actionable variants per sample',
    color = ''
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = 'right',
    panel.spacing = unit(0, 'cm'),
    plot.margin = unit(c(1,0,1,0), 'cm'))

# save a PNG and PDF of the plots
ggsave(
  path = output_dir,
  filename = paste0(current_date, '_', response_filter, '_total_avps_noleuk.png'),
  plot = total_avps_plot,
  device = 'png',
  width = 6,
  height = 7,
  dpi = 600
)

pdf(
  file = paste0(output_dir, '/', current_date, '_', response_filter, '_total_avps_noleuk.pdf'),
  width = 6,
  height = 7,
  useDingbats = FALSE,
  compress = FALSE
)
print(total_avps_plot)
dev.off()

ggsave(
  path = output_dir,
  filename = paste0(current_date, '_', response_filter, '_a_on_label_avps_noleuk.png'),
  plot = a_on_label_avps_plot,
  device = 'png',
  width = 6,
  height = 7,
  dpi = 600
)

pdf(
  file = paste0(output_dir, '/', current_date, '_', response_filter, '_a_on_label_avps_noleuk.pdf'),
  width = 6,
  height = 7,
  useDingbats = FALSE,
  compress = FALSE
)
print(a_on_label_avps_plot)
dev.off()

ggsave(
  path = output_dir,
  filename = paste0(current_date, '_', response_filter, '_a_off_label_avps_noleuk.png'),
  plot = a_off_label_avps_plot,
  device = 'png',
  width = 6,
  height = 7,
  dpi = 600
)

pdf(
  file = paste0(output_dir, '/', current_date, '_', response_filter, '_a_off_label_avps_noleuk.pdf'),
  width = 6,
  height = 7,
  useDingbats = FALSE,
  compress = FALSE
)
print(a_off_label_avps_plot)
dev.off()

ggsave(
  path = output_dir,
  filename = paste0(current_date, '_', response_filter, '_b_on_label_avps_noleuk.png'),
  plot = b_on_label_avps_plot,
  device = 'png',
  width = 6,
  height = 7,
  dpi = 600
)

pdf(
  file = paste0(output_dir, '/', current_date, '_', response_filter, '_b_on_label_avps_noleuk.pdf'),
  width = 6,
  height = 7,
  useDingbats = FALSE,
  compress = FALSE
)
print(b_on_label_avps_plot)
dev.off()

ggsave(
  path = output_dir,
  filename = paste0(current_date, '_', response_filter, '_b_off_label_avps_noleuk.png'),
  plot = b_off_label_avps_plot,
  device = 'png',
  width = 6,
  height = 7,
  dpi = 600
)

pdf(
  file = paste0(output_dir, '/', current_date, '_', response_filter, '_b_off_label_avps_noleuk.pdf'),
  width = 6,
  height = 7,
  useDingbats = FALSE,
  compress = FALSE
)
print(b_off_label_avps_plot)
dev.off()

# save the plots as an SVG file
svglite(
  filename = paste0(output_dir, '/', current_date, '_', response_filter, '_total_avps_noleuk.svg'),
  width = 6,
  height = 7
)
total_avps_plot
dev.off()

svglite(
  filename = paste0(output_dir, '/', current_date, '_', response_filter, '_a_on_label_avps_noleuk.svg'),
  width = 6,
  height = 7
)
a_on_label_avps_plot
dev.off()

svglite(
  filename = paste0(output_dir, '/', current_date, '_', response_filter, '_a_off_label_avps_noleuk.svg'),
  width = 6,
  height = 7
)
a_off_label_avps_plot
dev.off()

svglite(
  filename = paste0(output_dir, '/', current_date, '_', response_filter, '_b_on_label_avps_noleuk.svg'),
  width = 6,
  height = 7
)
b_on_label_avps_plot
dev.off()

svglite(
  filename = paste0(output_dir, '/', current_date, '_', response_filter, '_b_off_label_avps_noleuk.svg'),
  width = 6,
  height = 7
)
b_off_label_avps_plot
dev.off()

##################### Which actionable variants contribute to the significant increase in the different tier levels?

# read in treatment information to include as a label in the figure
treatment <- read_tsv(file = paste0(base_dir, '/path/to/hartwig_actionability.tsv')) %>%
  ### test del: remove liquid tumors because they are weird
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

# define the cancer types which show a significatly different AV number change in PCAWG vs Hartwig
a_on_label_cancer_types <- c('Breast cancer', 'Non small cell lung cancer')
a_off_label_cancer_types <- c('Pancreas neuroendocrine', 'Kidney clear cell carcinoma', 'Prostate carcinoma')
b_on_label_cancer_types <- c('Breast cancer', 'Non small cell lung cancer', 'Prostate carcinoma', 'Liposarcoma')
b_off_label_cancer_types <- c('Breast cancer', 'Hepatocellular carcinoma', 'Pancreas neuroendocrine', 'Ovarian cancer', 'Kidney clear cell carcinoma', 'Prostate carcinoma')

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
min_freq_idx <- pmax(variant_freq_df2$n_PCAWG, variant_freq_df2$n_Hartwig) > 3
variant_freq_df2 <- variant_freq_df2 %>%
  filter(min_freq_idx)

# source the functions to calculate the fishers exact test and wilcoxon test p-values and effect sizes
# from test matrices
source(paste0(base_dir, '/path/to/0_func/fisher_wilcoxon_matrix_functions.R'))

### calculate fisher's exact test p-value
# create fisher matrix
variant_freq_mat <- variant_freq_df2 %>%
  mutate(cancer_type_gene_variant = paste0(cancer_type_code, '_', tier_level, '_', gene_variant)) %>%
  select(cancer_type_gene_variant, n_PCAWG, n_non_PCAWG, n_Hartwig, n_non_Hartwig) %>%
  column_to_rownames(., var = 'cancer_type_gene_variant') %>%
  as.matrix()

# calculate exact test p-value
variant_freq_df2$fisher_pval <- fisherTest.matrix(variant_freq_mat)

# adjust p-value
variant_freq_df2 <- variant_freq_df2 %>%
  group_by(cancer_type) %>%
  mutate(fisher_pval_adj = p.adjust(fisher_pval, method = 'BH')) %>%
  ungroup() %>%
  mutate(significance_level = if_else(fisher_pval_adj < q_val_cutoff, '*', NA_character_))

# prepare ylim for plot
pcawg_freq_limit <- -100
hartwig_freq_limit <- 100

# define ASCO evidence tier level colors
label_colors <- c('#1e6220', '#2F913F', '#35E358', '#C4FFD0')
tier_label_levels <- c('A_On_label', 'A_Off_label', 'B_On_label', 'B_Off_label')
names(label_colors) <- tier_label_levels

# define cancer type label order
source(paste0(base_dir, '/path/to/r_objects/cancer_type_order.R'))

# arrange data for plotting
variant_freq_df_final <- variant_freq_df2  %>%
  filter(pmax(freq_Hartwig, freq_PCAWG) > 5) %>%
  filter(abs(freq_diff) > 5) %>%
  select(-n_Hartwig, -n_PCAWG, -sample_size_Hartwig, -sample_size_PCAWG) %>%
  pivot_longer(., cols = c(freq_Hartwig, freq_PCAWG), names_to = 'cohort', names_prefix = 'freq_', values_to = 'freq') %>%
  mutate(freq = if_else(cohort == 'PCAWG', freq * -1, freq)) %>%
  mutate(cancer_type_code = factor(cancer_type_code, levels = cancer_type_code_order)) %>%
  mutate(cancer_type_gene_variant = paste0(cancer_type_code, ': ', str_replace(gene_variant, pattern = '_', replacement = ' '))) %>%
  mutate(tier_level = factor(tier_level, levels = tier_label_levels)) %>%
  arrange(factor(cancer_type, levels = cancer_type_order), cancer_type_gene_variant)

# get the order of the camcer type gene variants to match general cancer type order we agreed upon
cancer_type_gene_variant_order <- unique(variant_freq_df_final$cancer_type_gene_variant)

variant_freq_df_final <- variant_freq_df_final %>%
  # impose order on the x axis variable
  mutate(cancer_type_gene_variant = factor(cancer_type_gene_variant, levels = rev(cancer_type_gene_variant_order)))

# prepare treatment labels for plotting
treatment <- treatment %>%
  inner_join(variant_freq_df_final, by = c('cancer_type_label', 'tier_level', 'gene_variant')) %>%
  filter(cohort == 'Hartwig') %>%
  mutate(tier_level = factor(tier_level, levels = tier_label_levels))

# export: 1000 x 1000
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
  geom_label_repel(data = treatment %>% filter(!is.na(significance_level)), aes(y = freq + 20, label = treatment_label),
                   ylim = c(110, NA), label.padding = 0.2, point.padding = 1, size = 4, min.segment.length = 0) +
  facet_grid(tier_level ~ ., 
             scales = 'free_y', space = 'free_y') +
  coord_flip(clip = 'off') +
  scale_y_continuous(
    limits = c(pcawg_freq_limit - 5, hartwig_freq_limit + 5),
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
    strip.text.y = element_blank())

# save the plot as PNG, PDF and SVG
ggsave(
  path = output_dir,
  filename = paste0(current_date, '_', response_filter, '_av_frequency_q', q_val_cutoff, '_noleuk',  '.png'),
  plot = av_frequency_plot,
  device = 'png',
  width = 11,
  height = 15,
  dpi = 600
)

pdf(
  file = paste0(output_dir, '/', current_date, '_', response_filter, '_av_frequency_q', q_val_cutoff, '_noleuk', '.pdf'),
  width = 11,
  height = 15,
  useDingbats = FALSE,
  compress = FALSE
)
print(av_frequency_plot)
dev.off()

# save the plots as an SVG file
svglite(
  filename = paste0(output_dir, '/', current_date, '_', response_filter, '_av_frequency_q', q_val_cutoff,  '.svg'),
  width = 11,
  height = 15
)
av_frequency_plot
dev.off()

######################################### SUPPLEMENTARY TABLES
# refresh token for gsheets
gs4_auth(email = 'your_email@gmail.com')

# get the correct IDs for supp tables (HMF IDs for HMF samples)
hmf_ids <- read_tsv(paste0(base_dir, '/path/to/metadata.tsv')) %>%
  # select relevant columns
  select(sample_id, supp_table_id, cohort, cancer_type)

# 1. number of actionable variants per tier label per sample
av_per_sample_tier_label_supp <- av_per_sample_tier_label %>%
  filter(cancer_type != 'Pancancer') %>%
  select(sample_id, label, n) %>%
  pivot_wider(., names_from = 'label', values_from = 'n', names_prefix = 'n_') %>%
  inner_join(hmf_ids, by = 'sample_id') %>%
  select(supp_table_id, cancer_type, n_A_On_label, n_A_Off_label, n_B_On_label, n_B_Off_label, n_Total) %>%
  rename(sample_id = supp_table_id)

# push to gsheets
googlesheets4::write_sheet(av_per_sample_tier_label_supp,
                           ss = 'https://docs.google.com/spreadsheets/...',
                           sheet = 'av_per_sample_and_tier')

# 2. global fisher test results
global_fisher_results_supp <- global_contingency_df_supp %>%
  mutate(prop_Hartwig = n_Hartwig / (n_Hartwig + n_non_Hartwig),
         prop_PCAWG = n_PCAWG / (n_PCAWG + n_non_PCAWG)) %>%
  mutate(fold_change = round((prop_Hartwig - prop_PCAWG) / prop_PCAWG, 1),
         fold_change_label = if_else(fold_change > 0, str_c('+', fold_change), as.character(fold_change))) %>%
  rename_all(~ gsub('_non', '_without_av', .x)) %>%
  rename(n_with_av_PCAWG = n_PCAWG,
         n_with_av_Hartwig = n_Hartwig) %>%
  select(cancer_type:fisher_pval_adj, fold_change)

# push to gsheets
googlesheets4::write_sheet(global_fisher_results_supp,
                           ss = 'https://docs.google.com/spreadsheets/...',
                           sheet = 'statistical_comparison_results')

# 3. av frequencies
av_frequency_supp <- variant_freq_df2 %>%
  select(cancer_type, cancer_type_code, gene_variant, tier_level, n_PCAWG, n_non_PCAWG, n_Hartwig, n_non_Hartwig, fisher_pval, fisher_pval_adj) %>%
  rename_all(~ gsub('_non', '_without_av', .x)) %>%
  rename(n_with_av_PCAWG = n_PCAWG,
         n_with_av_Hartwig = n_Hartwig) %>%
  group_by(tier_level) %>%
  group_split()

googlesheets4::write_sheet(av_frequency_supp[[1]],
                           ss = 'https://docs.google.com/spreadsheets/...',
                           sheet = paste0('av_frequency_', unique(av_frequency_supp[[1]]$tier_level)))

googlesheets4::write_sheet(av_frequency_supp[[2]],
                           ss = 'https://docs.google.com/spreadsheets/...',
                           sheet = paste0('av_frequency_', unique(av_frequency_supp[[2]]$tier_level)))

googlesheets4::write_sheet(av_frequency_supp[[3]],
                           ss = 'https://docs.google.com/spreadsheets/...',
                           sheet = paste0('av_frequency_', unique(av_frequency_supp[[3]]$tier_level)))

googlesheets4::write_sheet(av_frequency_supp[[4]],
                           ss = 'https://docs.google.com/spreadsheets/...',
                           sheet = paste0('av_frequency_', unique(av_frequency_supp[[4]]$tier_level)))
