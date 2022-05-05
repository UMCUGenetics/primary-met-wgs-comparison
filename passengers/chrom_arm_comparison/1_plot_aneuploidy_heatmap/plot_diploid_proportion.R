############### Plot diploid proportion histograms ###############
# author: remy (sascha)
# date: 09/11/2021
# last updated: 28/03/2022

### Description
# This script creates the barplots that compare the diploid proportion from diploid_proportion files stratified by a certain strata (e.g. cohort) that is specified
# by the user under the 'Splitting the dataset' section. A Mann Whitney test is performed to compare the diploid proportion between strata.
# These barplots are combined with the mean arm ploidy heatmap.

# libraries
library(tidyverse) # data maniplutation ans plotting
library(broom) # tidy model output

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
                               axis.text.y= element_text(size=12), 
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

# read in the diploid_proportion table (from purity files)
gz_file <- gzfile(paste0(base_dir, '/path/to/purple-purity-fraction-changed.txt.gz'))
diploid_proportion <- read_tsv(file = gz_file)

# join diploid_proportion table to metadata
metadata <- metadata %>%
  left_join(diploid_proportion, by = 'sample_id') %>%
# calculate the LOH proportion of diploid genomes
mutate(loh_proportion = 1 - diploidProportion)

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

# custom function to calculate the Mann Whitney test p-value between diploid proportion of cohorts per cancer type
get_wilcox_table <- function(df_input, cancer_type_input, is_wgd = TRUE) {
  
  # create a list of two dfs that represent the samples in the cohorts for one cancer type
  wilcox_list <- df_input %>%
    filter(cancer_type == cancer_type_input) %>%
    filter(whole_genome_duplication == is_wgd) %>%
    select(sample_id, cohort, loh_proportion) %>%
    group_by(cohort) %>%
    group_split()
  
  # calculate the MW p-value between diploid proportions of the cohorts
  wilcox_table <- tidy(wilcox.test(wilcox_list[[1]]$loh_proportion, 
                                   wilcox_list[[2]]$loh_proportion, 
                                   alternative = 'two.sided',
                                   paired = FALSE,
                                   mu = 0)) %>%
    mutate(cancer_type = cancer_type_input,
           significance_level = if_else(p.value < 0.05, '*', NA_character_)) %>%
    select(cancer_type, method, p.value, significance_level)
  
  return(wilcox_table)
  
}

# get the cancer type names into a vector
cancer_type_vec <- unique(metadata$cancer_type)

# calculate the Mann whitney p-value for all cancer types
wgd_wilcox_table <- map_df(cancer_type_vec, ~get_wilcox_table(metadata, .x, is_wgd = TRUE)) %>%
  mutate(whole_genome_duplication = TRUE)
non_wgd_wilcox_table <- map_df(cancer_type_vec, ~get_wilcox_table(metadata, .x, is_wgd = FALSE)) %>%
  mutate(whole_genome_duplication = FALSE)
wilcox_table <- union(wgd_wilcox_table, non_wgd_wilcox_table)

# calculate median diploid proportion per cancer type and cohort then add this info to the wilcox_table (for point layer in plot)
median_diploid <- metadata %>%
  group_by(cancer_type, cohort, whole_genome_duplication) %>%
  summarise(median_loh_proportion = median(loh_proportion), .groups = 'drop')
wilcox_table <- wilcox_table %>%
  left_join(median_diploid, by = c('cancer_type', 'whole_genome_duplication'))

# get the alternating index number for the grey background rectangles
rect_idx <- seq(2, 22, by = 2)
rect_idx <- cancer_type_order[rect_idx]

# get the cancer types in the right order
# NOTE: HAVE TO DO THIS FOR ALL TABLES THAT ARE USED IN THE PLOT
metadata <- metadata %>%
  mutate(cancer_type = factor(cancer_type, levels = cancer_type_order),
         stadium = if_else(cohort == 'HMF', 'HMF', 'PCAWG'),
         stadium = factor(stadium, levels = c('PCAWG', 'HMF')))
wilcox_table <- wilcox_table %>%
  mutate(cancer_type = factor(cancer_type, levels = cancer_type_order),
         stadium = if_else(cohort == 'HMF', 'HMF', 'PCAWG'),
         stadium = factor(stadium, levels = c('PCAWG', 'HMF')))

# plot the violin facets
wgd_diploid_violin <- metadata %>%
  filter(whole_genome_duplication) %>%
  ggplot(., aes(x = stadium, y = loh_proportion)) +
  geom_rect(data = . %>% 
              filter(cancer_type %in% rect_idx) %>%
              distinct(cancer_type), aes(fill = cancer_type),
            inherit.aes = FALSE,
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#a0a0a0') +
  geom_violin(aes(fill = cohort), color = 'black', scale = 'width') +
  geom_point(data = wilcox_table %>% filter(whole_genome_duplication), 
             aes(y = median_loh_proportion), size = 2) +
  geom_text(data = wilcox_table %>% filter(cohort == 'HMF') %>% filter(whole_genome_duplication), 
            aes(label = significance_level, y = 1.1), size = 8, nudge_x = 0) +
  facet_wrap(cancer_type ~ ., ncol = 1, strip.position = 'left') +
  scale_x_discrete(
    labels = c('Hartwig', 'PCAWG'),
    position = 'top',
    expand = c(0.5,0.5)) + # 'expand' is used to center the histogram in the panel
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1.2), expand = c(0,0)) +
  scale_fill_manual(values = cohort_colors) +
  labs(title = 'Genome diploid\nproportion of\nWGD samples',
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

non_wgd_diploid_violin <- metadata %>%
  filter(!whole_genome_duplication) %>%
  ggplot(., aes(x = stadium, y = loh_proportion)) +
  geom_rect(data = . %>% 
              filter(cancer_type %in% rect_idx) %>%
              distinct(cancer_type), aes(fill = cancer_type),
            inherit.aes = FALSE,
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#a0a0a0') +
  geom_violin(aes(fill = cohort), color = 'black', scale = 'width') +
  geom_point(data = wilcox_table %>% filter(!whole_genome_duplication), 
             aes(y = median_loh_proportion), size = 2) +
  geom_text(data = wilcox_table %>% filter(cohort == 'HMF') %>% filter(!whole_genome_duplication), 
            aes(label = significance_level, y = 1.1), size = 8, nudge_x = 0) +
  facet_wrap(cancer_type ~ ., ncol = 1, strip.position = 'left') +
  scale_x_discrete(
    labels = c('Hartwig', 'PCAWG'),
    position = 'top',
    expand = c(0.5,0.5)) + # 'expand' is used to center the histogram in the panel
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1.2), expand = c(0,0)) +
  scale_fill_manual(values = rev(cohort_colors)) +
  labs(title = 'LOH proportion in\ndiploid tumors',
       subtitle = '*Mann Whitney\ntest p-value < 0.05',
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

################### SUPPLEMENTARY TABLE
supp_loh_proportion <- diploid_proportion %>%
  mutate(loh_proportion = 1 - diploidProportion) %>%
  select(-diploidProportion)
  