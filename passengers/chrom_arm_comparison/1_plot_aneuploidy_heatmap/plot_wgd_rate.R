############### Plot whole genome duplication rate histograms ###############
# author: remy (sascha)
# date: 30/08/2021
# last updated: 28/03/2022

### Description
# This script creates the barplots that compare the WGD rate stratified by a certain strata (e.g. cohort) that is specified
# by the user under 'Splitting the dataset'. These barplots are combined with the mean arm ploidy heatmap.

# libraries
library(tidyverse)
library(naniar)

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

# define variables
exclude_hypermutators <- FALSE

## ggplot custom theme
theme_bigfont_barplot <- theme(plot.title = element_text(size=16),
                               plot.subtitle = element_text(size = 14),
                               axis.text.x= element_text(size=15),
                               axis.text.y= element_text(size=13), 
                               axis.title=element_text(size=18),
                               legend.title = element_text(size = 16),
                               legend.text = element_text(size = 14),
                               strip.text.x = element_text(size = 15),
                               strip.text.y = element_text(size = 15))

# source the functions to calculate the fishers exact test and wilcoxon test p-values and effect sizes
# from test matrices
source(paste0(base_dir, '/path/to/0_func/fisher_wilcoxon_matrix_functions.R'))
# define cancer type order
source(paste0(base_dir, '/path/to/r_objects/cancer_type_order.R'))
# load cohort colors
source(paste0(base_dir, '/path/to/r_objects/color_codes.R'))

# read metadata in
metadata <- read_tsv(paste0(base_dir, '/path/to/metadata.tsv')) %>%
  # select relevant columns
  select(sample_id, cohort, is_hypermutated, is_metastatic, cancer_type, whole_genome_duplication)

# get the 22 cancer types with enough samples in both PCAWG and HMF cohort
cancer_type_abundance <- metadata %>%
  select(sample_id, cancer_type, cohort) %>%
  group_by(cancer_type, cohort) %>%
  count(name = 'n_per_stadium') %>% 
  ungroup() %>%
  pivot_wider(names_from = cohort, 
              values_from = n_per_stadium) %>%
  filter(!is.na(PCAWG),
         !is.na(HMF)) %>%
  # filter out any cancer with less than 10 samples in either primary or metastatic group
  filter(pmin(HMF, PCAWG) > 14) %>%
  mutate(n_total = PCAWG + HMF,
         primary_ratio = PCAWG / HMF) %>%
  arrange(desc(n_total))
cancers_with_enough_samples <- cancer_type_abundance$cancer_type
rm(cancer_type_abundance)

# filter the metadata for those 22 cancer types
metadata <- metadata %>%
  filter(cancer_type %in% cancers_with_enough_samples)

# exclude hypermutators?
if (exclude_hypermutators) {
  metadata <- metadata %>%
    filter(is_hypermutated == FALSE)
}

# --------------- Splitting the dataset
# in order to stratify the dataset and perform the tests by a binary strata, this step transforms the strata column into a 'cohort'
# column which is then subsequently used by functions downstream
# DEFAULT strata is the PCAWG vs HMF 'cohort' column
metadata <- metadata %>%
  mutate(cohort = if_else(
    cohort == 'HMF', ### <---------------------------------- CHANGE HERE 
    'HMF', 'PCAWG')) # create a new cohort column

# specify the strata you created (used in plot subtitle)
# strata <- 'cohort'

# ---------------

# count the number of WGD samples per cancer type and stadium
metadata <- metadata %>%
  group_by(cancer_type, cohort, whole_genome_duplication) %>%
  mutate(n_wgd = n()) %>%
  ungroup()

metadata_wide <- metadata %>%
  select(cancer_type, cohort, whole_genome_duplication, n_wgd) %>%
  distinct() %>%
  pivot_wider(names_from = whole_genome_duplication, names_prefix = 'wgd_', values_from = n_wgd)

# pivot wider per wgd_TRUE and wgd_FALSE, then join the two table together to form the final table
# that is used to calculate the fishers exact test pvalue
metadata_wider_wgd_true <- metadata_wide %>%
  select(-wgd_FALSE) %>%
  pivot_wider(names_from = cohort, values_from = wgd_TRUE, names_prefix = 'wgd_true_')

metadata_wider_wgd_false <- metadata_wide %>%
  select(-wgd_TRUE) %>%
  pivot_wider(names_from = cohort, values_from = wgd_FALSE, names_prefix = 'wgd_false_')

metadata_wider_wgd_combined <- metadata_wider_wgd_true %>%
  inner_join(metadata_wider_wgd_false, by = 'cancer_type') %>%
  rename(metastatic_wgd = 'wgd_true_HMF',
         metastatic_non_wgd = 'wgd_false_HMF',
         primary_wgd = 'wgd_true_PCAWG',
         primary_non_wgd = 'wgd_false_PCAWG') %>%
  select(cancer_type, primary_wgd, primary_non_wgd, metastatic_wgd, metastatic_non_wgd)

# prepare matrix for fisher pval calculation
metadata_wider_wgd_combined_mat <- metadata_wider_wgd_combined %>% 
  column_to_rownames(var = 'cancer_type') %>% 
  as.matrix()

# calculate fisher, odds ratio and cramers v + add these as separate columns
metadata_wider_wgd_combined$fisher_pval <- fisherTest.matrix(metadata_wider_wgd_combined_mat)

# get the alternating index number for the grey background rectangles
rect_idx <- seq(2, 22, by = 2)
rect_idx <- cancer_type_order[rect_idx]


metadata_wider_wgd_combined_plot <- metadata_wider_wgd_combined %>%
  mutate(primary_wgd_fraction = primary_wgd / (primary_wgd + primary_non_wgd),
         metastatic_wgd_fraction = metastatic_wgd / (metastatic_wgd + metastatic_non_wgd)) %>%
  mutate(cancer_type = factor(cancer_type, levels = cancer_type_order)) %>%
  pivot_longer(cols = primary_wgd_fraction:metastatic_wgd_fraction, names_to = 'stadium_wgd', values_to = 'fraction') %>%
  mutate(significance_level = if_else(fisher_pval < 0.01 & stadium_wgd == 'primary_wgd_fraction', '*', NA_character_))

# plot the histogram facets
wgd_histogram <- metadata_wider_wgd_combined_plot %>%
  ggplot(., aes(x = stadium_wgd, y = fraction)) +
  geom_rect(data = . %>% 
              filter(cancer_type %in% rect_idx) %>%
              distinct(cancer_type), aes(fill = cancer_type),
            inherit.aes = FALSE,
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#a0a0a0') +
  geom_col(aes(fill = stadium_wgd), color = 'black') +
  geom_text(aes(label = significance_level, y = 0.9), size = 8, nudge_x = -1) +
  facet_wrap(cancer_type ~ ., ncol = 1, strip.position = 'left') +
  scale_x_discrete(
    labels = c('Hartwig', 'PCAWG'),
    position = 'top',
    expand = c(0.5,0.5)) + # 'expand' is used to center the histogram in the panel
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1), expand = c(0,0)) +
  scale_fill_manual(values = rev(cohort_colors)) +
  labs(title = 'WGD\nfraction',
       x = '',
       y = '') +
  coord_flip() +
  theme_classic() +
  theme_bigfont_barplot +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    strip.background = element_blank(),
    strip.text.y.left = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    legend.position = 'none',
    panel.spacing = unit(0, 'cm'))

############## SUPPLEMENTARY TABLE
supp_wgd <- metadata %>%
  select(sample_id, whole_genome_duplication)
