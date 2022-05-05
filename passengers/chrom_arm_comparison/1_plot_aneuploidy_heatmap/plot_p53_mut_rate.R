############### Plot p53 mutated rate histograms ###############
# author: remy (sascha)
# date: 05/11/2021
# last updated: 28/03/2022

### Description
# This script creates the barplots that compare the p53 mutation rate stratified by a certain strata (e.g. cohort) that is specified
# by the user under the 'Splitting the dataset' section. These barplots are combined with the mean arm ploidy heatmap.

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

mutated_driver <- 'TP53'

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
  select(sample_id, cohort, is_hypermutated, whole_genome_duplication, cancer_type)

# read in pretreatment info
pretreatment <- read_tsv(paste0(base_dir, '/path/to/pretreatment.tsv')) %>%
  select(sample_id, is_pretreated)

# join the pretreatment table to the metadata table
metadata <- metadata %>%
  left_join(pretreatment, by = 'sample_id')

# driver results date
driver_results_date <- '2021_12_03' # old: 2021_07_22, new: 2021_10_25, new pannel == 70: 2021_12_03

# read in linx drivers
linx_p53_drivers <- read_tsv(paste0(base_dir, '/path/to/', driver_results_date ,'/linx_drivers.tsv'),
                         col_names = TRUE) %>%
  filter(gene == mutated_driver) %>%
  select(sample_id, gene) %>%
  distinct()

# join the p53 table to the metadata table, then add a new column 'is_mutated'
metadata <- metadata %>%
  left_join(linx_p53_drivers, by = 'sample_id') %>%
  mutate(is_mutated = if_else(is.na(gene), FALSE, TRUE))

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

# count the number of samples with mutated p53 per cancer type and stadium
metadata <- metadata %>%
  group_by(cancer_type, cohort, is_mutated) %>%
  mutate(n_mutated = n()) %>%
  ungroup()

metadata_wide <- metadata %>%
  select(cancer_type, cohort, is_mutated, n_mutated) %>%
  distinct() %>%
  pivot_wider(names_from = is_mutated, names_prefix = 'is_mutated_', values_from = n_mutated) %>%
  # replace NAs with 0 in the numeric columns
  replace_na(., replace = list(is_mutated_TRUE = 0, is_mutated_FALSE = 0))

# pivot wider per is_mutated_TRUE and is_mutated_FALSE, then join the two table together to form the final table
# that is used to calculate the fishers exact test pvalue
metadata_wider_is_mutated_true <- metadata_wide %>%
  select(-is_mutated_FALSE) %>%
  pivot_wider(names_from = cohort, values_from = is_mutated_TRUE, names_prefix = 'is_mutated_true_')

metadata_wider_is_mutated_false <- metadata_wide %>%
  select(-is_mutated_TRUE) %>%
  pivot_wider(names_from = cohort, values_from = is_mutated_FALSE, names_prefix = 'is_mutated_false_')

metadata_wider_is_mutated_combined <- metadata_wider_is_mutated_true %>%
  inner_join(metadata_wider_is_mutated_false, by = 'cancer_type') %>%
  rename(metastatic_is_mutated = 'is_mutated_true_HMF',
         metastatic_non_mutated = 'is_mutated_false_HMF',
         primary_is_mutated = 'is_mutated_true_PCAWG',
         primary_non_mutated = 'is_mutated_false_PCAWG') %>%
  select(cancer_type, primary_is_mutated, primary_non_mutated, metastatic_is_mutated, metastatic_non_mutated)

# prepare matrix for fisher pval calculation
metadata_wider_is_mutated_combined_mat <- metadata_wider_is_mutated_combined %>% 
  column_to_rownames(var = 'cancer_type') %>% 
  as.matrix()

# calculate fisher, odds ratio and cramers v + add these as separate columns
metadata_wider_is_mutated_combined$fisher_pval <- fisherTest.matrix(metadata_wider_is_mutated_combined_mat)

# get the alternating index number for the grey background rectangles
rect_idx <- seq(2, 22, by = 2)
rect_idx <- cancer_type_order[rect_idx]

# finalize data to plot it
metadata_wider_is_mutated_combined_plot <- metadata_wider_is_mutated_combined %>%
  mutate(primary_is_mutated_fraction = primary_is_mutated / (primary_is_mutated + primary_non_mutated),
         metastatic_is_mutated_fraction = metastatic_is_mutated / (metastatic_is_mutated + metastatic_non_mutated)) %>%
  mutate(cancer_type = factor(cancer_type, levels = cancer_type_order)) %>%
  pivot_longer(cols = primary_is_mutated_fraction:metastatic_is_mutated_fraction, names_to = 'stadium_is_mutated', values_to = 'fraction') %>%
  mutate(significance_level = if_else(fisher_pval < 0.05 & stadium_is_mutated == 'primary_is_mutated_fraction', '*', NA_character_))

# plot the histogram facets
mutated_histogram <- metadata_wider_is_mutated_combined_plot %>%
  ggplot(., aes(x = stadium_is_mutated, y = fraction)) +
  geom_rect(data = . %>% 
              filter(cancer_type %in% rect_idx) %>%
              distinct(cancer_type), aes(fill = cancer_type),
            inherit.aes = FALSE,
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#a0a0a0') +
  geom_col(aes(fill = stadium_is_mutated), color = 'black') +
  geom_text(aes(label = significance_level, y = 1.1), size = 8, nudge_x = -1) +
  facet_wrap(cancer_type ~ ., ncol = 1, strip.position = 'left') +
  scale_x_discrete(
    labels = c('Hartwig', 'PCAWG'),
    position = 'top',
    expand = c(0.5,0.5)) + # 'expand' is used to center the histogram in the panel
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1.2), expand = c(0,0)) +
  scale_fill_manual(values = rev(cohort_colors)) +
  labs(title = paste0('Mutated ', mutated_driver, '\nfraction'),
       subtitle = '*Fisher\'s exact test\np-value < 0.05',
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
    legend.position = 'none',
    panel.spacing = unit(0, 'cm')
    )

################## SUPPLEMENTARY TABLE
supp_p53_mut <- metadata %>%
  select(sample_id, is_mutated) %>%
  rename(p53_mut = 'is_mutated')
