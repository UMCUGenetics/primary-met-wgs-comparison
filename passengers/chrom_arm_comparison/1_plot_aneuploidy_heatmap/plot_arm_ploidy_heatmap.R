############### Plot chromosome arm ploidy comparison ###############
# author: remy (sascha)
# date: 17/08/2021
# last updated: 21/04/2022

### Description
# This script prepares the arm ploidy matrix to plot the mean arm gains and losses, and compares these mean ploidy levels
# by a certain binary strata which can be specified (e.g. cohort or treatment). Meanwhile a Mann Whitney U test is performed
# per cancer type and arm between the two groups in the specified strata. P-values are FDR adjusted per cancer type.
# Significant q-values (< 0.05) are plotted on top of the heatmap tiles to highlight sig. differences. Additionally,
# a barplot which shows the number of cancer types with a sig. gain / loss in that particular arm is plotted on top of the heatmap.
# This heatmap is then arranged together with 90 degrees rotated barplots (e.g. WGD rate).

# libraries
library(tidyverse) # data manipulation and plotting
library(broom) # tidy model output
library(naniar) # handling NAs
library(patchwork) # plot grids
library(svglite) # handle the issue with SVGs in Inkscape
library(googlesheets4) # push supplementary tables to gsheets

#========= Path prefixes =========#
heatmap_base_dir <- list(
  hpc='/base/path',
  mnt='/base/path',
  umc='/base/path'
)

for(i in heatmap_base_dir){
  if(dir.exists(i)){
    heatmap_base_dir <- i
    break
  }
}

output_dir <- list(
  local='/output/path/',
  umc='/output/path'
)

for(i in output_dir){
  if(dir.exists(i)){
    output_dir <- i
    break
  }
}

current_date <- format(Sys.Date(), '%Y_%m_%d')

# define variables
exclude_hypermutators <- FALSE

## ggplot custom theme
theme_bigfont <- theme(plot.title = element_text(size=18),
                       plot.subtitle = element_text(size = 16),
                       axis.text.x= element_text(size=15),
                       axis.text.y= element_text(size=15), 
                       axis.title=element_text(size=16),
                       legend.title = element_text(size = 16),
                       legend.text = element_text(size = 14),
                       strip.text.x = element_text(size = 15),
                       strip.text.y = element_text(size = 15))

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
source(paste0(heatmap_base_dir, '/path/to/0_func/fisher_wilcoxon_matrix_functions.R'))
# load cohort colors
source(paste0(heatmap_base_dir, '/path/to/r_objects/color_codes.R'))
# define cancer type order
source(paste0(heatmap_base_dir, '/path/to/r_objects/cancer_type_order.R'))

# read metadata in
metadata <- read_tsv(paste0(heatmap_base_dir, '/path/to/metadata.tsv')) %>%
  # select relevant columns
  select(sample_id, supp_table_id, cohort, is_metastatic, is_hypermutated, cancer_type, cancer_type_code)

# join the pretreatment table to the metadata table
metadata <- metadata %>%
  left_join(pretreatment, by = 'sample_id')

# read arm matrix into R
arm_matrix_input <- gzfile(paste0(heatmap_base_dir, '/path/to/arm_ploidy_matrix.txt.gz'))
arm_matrix <- read.table(arm_matrix_input, header = TRUE)

# change weird colnames to something more meaningful
colnames(arm_matrix) <- sub(pattern = 'X(\\d)', replacement = 'chr\\1', colnames(arm_matrix))
colnames(arm_matrix) <- sub(pattern = '(X)', replacement = 'chr\\1', colnames(arm_matrix)) # X chrom

# subset matrix, arm and genome ploidy
arm_ploidy_mat <- arm_matrix[, 1:ncol(arm_matrix)-1]
genome_ploidy_vec <- arm_matrix[, ncol(arm_matrix)]

# subtract genome ploidy level from arm level, row by row
arm_ploidy_mat <- sweep(arm_ploidy_mat, MARGIN = 1, FUN = '-', genome_ploidy_vec)

# transform into dataframe, add rownames as a separate column
arm_ploidy_df1 <- as.data.frame(arm_ploidy_mat)
arm_ploidy_df1 <- arm_ploidy_df1 %>%
  mutate(sample_id = rownames(arm_ploidy_df1), .before = chr1p)

# join metadata to that long table
arm_ploidy_df <- arm_ploidy_df1 %>%
  inner_join(metadata, by = 'sample_id')

# filter hypermutators?
if (exclude_hypermutators) {
  arm_ploidy_df <- arm_ploidy_df %>%
    filter(is_hypermutated == FALSE)
}

# --------------- Splitting the dataset
# in order to stratify the dataset and perform the tests by strata, this step transforms the strata column into a 'cohort'
# column which is then subsequently used by functions downstream
# default strata is the PCAWG vs HMF cohort column
arm_ploidy_df <- arm_ploidy_df %>%
  mutate(cohort = if_else(
    cohort == 'HMF', ### <---------------------------------- CHANGE HERE 
    'HMF', 'PCAWG')) # create a new cohort column

# specify the strata you created (used in plot subtitle)
strata <- 'cohort' ### <---------------------------------- CHANGE HERE

# ---------------

# number of total samples per cancer type
sample_size <- arm_ploidy_df %>%
  group_by(cancer_type, cancer_type_code, cohort) %>%
  summarise(sample_size_tot = n()) %>%
  pivot_wider(names_from = cohort, values_from = sample_size_tot, names_prefix = 'n_') %>%
  mutate(cancer_type_label = str_c(cancer_type, ' (', n_PCAWG, ' vs. ', n_HMF, ')')) %>%
  mutate(cancer_type_code_label = str_c(cancer_type_code, ' (', n_PCAWG, ' vs. ', n_HMF, ')')) %>%
  ungroup()
  
# get the axis in the right orders to compare the plot to Nguyen et al. 2021 Fig. 2a
arm_order <- c('chr1p', 'chr1q', 'chr2p', 'chr2q',
               'chr3p', 'chr3q', 'chr4p', 'chr4q',
               'chr5p', 'chr5q', 'chr6p', 'chr6q',
               'chr7p', 'chr7q', 'chr8p', 'chr8q',
               'chr9p', 'chr9q', 'chr10p', 'chr10q',
               'chr11p', 'chr11q', 'chr12p', 'chr12q',
               'chr13q', 'chr14q',
               'chr15q', 'chr16p', 'chr16q',
               'chr17p', 'chr17q', 'chr18p', 'chr18q',
               'chr19p', 'chr19q', 'chr20p', 'chr20q',
               'chr21q', 'chr22q')

# ------------ Prepare matrix for statistical testing: Wilcoxon
# rename var
arm_ploidy_df_wilcox <- arm_ploidy_df

# filter hypermutators?
if (exclude_hypermutators) {
  arm_ploidy_df_wilcox <- arm_ploidy_df_wilcox %>%
    filter(is_hypermutated == FALSE)
}

# custom function to calculate the wilcoxon test pvalue and cliff delta per cancer_type for each 
# chromosome arm (primary vs. metastatic)
get_wilcoxon_table <- function(df_input, cancer_type_input) {
 
  # filter the table for the desired cancer_type
  arm_ploidy_df_wilcox <- df_input %>%
    filter(cancer_type == cancer_type_input)
  
  # prime the output table
  wilcox_table <- arm_ploidy_df_wilcox %>%
    select(cancer_type, chr1p:chr22q) %>%
    pivot_longer(cols = chr1p:chr22q, names_to = 'arm', values_to = 'value') %>%
    select(cancer_type, arm) %>%
    distinct()
  
  wilcox_list <- arm_ploidy_df_wilcox %>%
    group_by(cohort) %>%
    group_split()
  
  # set the rownames, select only the arms, and transform both dfs into matrices
  # [[1]] contains the cohort == 'HMF' (metastatic)
  # [[2]] contains the cohort == 'PCAWG' (primary)
  wilcox_list <- map(wilcox_list, ~column_to_rownames(.x, var = 'sample_id'))
  wilcox_list <- map(wilcox_list, ~select(.x, chr1p:chr22q))
  wilcox_list <- map(wilcox_list, ~as.matrix(.x))
  
  # run columnwise wilcoxon test and cliffs delta calculation, then adjust the p-value 
  wilcox_table$wilcoxon_pval <- wilcoxTest.matrix(wilcox_list[[1]], wilcox_list[[2]], alternative = 'two.sided')
  wilcox_table$wilcoxon_pval_adj <- p.adjust(wilcox_table$wilcoxon_pval, method = 'fdr')
  wilcox_table$cliff_delta <- cliffDelta.matrix(wilcox_list[[1]], wilcox_list[[2]])
  
  return(wilcox_table)
   
}

# get the cancer type names into a vector
cancer_type_vec <- unique(metadata$cancer_type)

# perform mann whitney U test for comparison of mean ploidy between cohorts
cancer_wilcoxon_list <- map(cancer_type_vec, ~get_wilcoxon_table(arm_ploidy_df_wilcox, cancer_type_input = .x))
# combine the individual list elements into on table
cancer_wilcoxon_df <- do.call(rbind, cancer_wilcoxon_list)

# only keep significant p-values
cancer_wilcoxon_df_sig <- cancer_wilcoxon_df %>%
  filter(wilcoxon_pval_adj < 0.05)

# add the significance * to each data point (for the plot annotation later)
cancer_wilcoxon_df_sig <- cancer_wilcoxon_df_sig %>%
  mutate(significance_level = case_when(
    wilcoxon_pval_adj < 0.05 ~ '*',
    TRUE ~ 'none')) %>%
  replace_with_na(list(significance_level = 'none'))

# prepare the data for plotting: calculate the mean arm ploidy
arm_ploidy_df_wilcox_mean <- arm_ploidy_df_wilcox %>%
  group_by(cancer_type, cohort) %>%
  summarise(across(chr1p:chr22q, ~mean(.x)), .groups = 'drop') %>%
  filter(cancer_type %in% cancer_type_vec)

arm_ploidy_df_wilcox_mean <- arm_ploidy_df_wilcox_mean %>%
  pivot_longer(cols = chr1p:chr22q, names_to = 'arm', values_to = 'mean_arm_ploidy')

# join the sample size to the arm df
mean_ploidy_plot_wilcox <- sample_size %>%
  select(cancer_type, cancer_type_label, cancer_type_code_label) %>%
  inner_join(arm_ploidy_df_wilcox_mean, by = c('cancer_type'))

# get the right cancer type label order
cancer_type_order <- mean_ploidy_plot_wilcox %>%
  select(cancer_type, cancer_type_label) %>%
  distinct() %>%
  mutate(cancer_type = factor(cancer_type, levels = cancer_type_order)) %>%
  arrange(cancer_type)
cancer_type_order <- cancer_type_order$cancer_type_label

# get cancer type label to finalize the sig table for manual plot annotation
cancer_wilcoxon_df_sig <- mean_ploidy_plot_wilcox %>%
  select(cancer_type, arm, cancer_type_label) %>%
  distinct()  %>%
  right_join(cancer_wilcoxon_df_sig, by = c('cancer_type', 'arm')) %>%
  mutate(arm = factor(arm, levels = arm_order),
         cancer_type_label = factor(cancer_type_label, levels = cancer_type_order))

# plot subtitle conditional
if (exclude_hypermutators) {
  subtitle_text <- 'Excluding hypermutators'
} else {
  subtitle_text <- 'Including hypermutators'
}

# colors
tile_colors <- c('#000080' ,'#009ACD','#E5E5E5','#FF758C', '#FF007F', '#A8152E')
breaks <- seq(-1, 2.5, by = 0.25)

# plot the nguyen et al. 2021 replicated plot: wilcoxon method
# export: 2000 x 1100
mean_ploidy_heatmap <- mean_ploidy_plot_wilcox %>%
  mutate(arm = factor(arm, levels = arm_order),
         cancer_type_label = factor(cancer_type_label, levels = cancer_type_order),
         stadium = if_else(cohort == 'HMF', 'Metastatic', 'Primary'),
         stadium = factor(stadium, levels = c('Metastatic', 'Primary'))) %>%
  ggplot(., aes(x = arm, y = stadium, fill = mean_arm_ploidy)) +
  geom_tile(color = '#FFFFFF', size = 0.4) +
  # wilcoxon
  geom_text(data = cancer_wilcoxon_df_sig, aes(y = 1.5, fill = NULL, label = significance_level), size = 8, nudge_y = -0.5) +
  facet_grid(cancer_type_label ~ ., switch = 'y') +
  scale_fill_stepsn(colours = tile_colors, na.value = '#000000',
                    breaks = breaks,
                    limits = c(-1.5, 2.5),
                    guide = guide_colorsteps(title = 'Arm ploidy gains / losses relative to\nmean genome ploidy', title.position = 'top',
                                             barwidth = 25, show.limits = TRUE)) +
  scale_x_discrete(labels = str_replace(arm_order, pattern = 'chr(.+)', replacement = '\\1')) +
  scale_y_discrete(position = 'right') +
  labs(
    x = 'Chromosome arm',
    y = '') +
  theme_bigfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y.right = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    legend.position = 'bottom',
    legend.text = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.margin = margin(r = 0, l = 0, unit = 'mm')
  )

# create arm level sig difference barplot
total_arm_df <- data.frame(arm = arm_order)

# barplot of number of sig arm per cancer type
sig_diff_barplot <- cancer_wilcoxon_df_sig %>%
  inner_join(arm_ploidy_df_wilcox_mean, by = c('cancer_type', 'arm')) %>%
  pivot_wider(names_from = 'cohort', values_from = 'mean_arm_ploidy') %>%
  mutate(sig_gain_loss = case_when(HMF > PCAWG ~ 'gain',
                                   HMF < PCAWG ~ 'loss',
                                   TRUE ~ 'nix')) %>%
  group_by(arm, sig_gain_loss) %>%
  summarise(n_sig = n(), .groups = 'drop') %>%
  right_join(total_arm_df, by = 'arm') %>%
  complete(., arm, sig_gain_loss, fill = list(n_sig = 0)) %>%
  drop_na(., sig_gain_loss) %>%
  mutate(arm = factor(arm, levels = arm_order)) %>%
  ggplot(., aes(x = arm, y = n_sig)) +
  geom_col(aes(fill = sig_gain_loss), color = 'black', position = 'stack') +
  scale_y_continuous(expand = c(0.01,0), breaks = seq(1, 10, by = 2)) +
  scale_fill_manual(name = '', values = c('#FF007F', '#009ACD')) +
  labs(
    x = '',
    y = ''
  ) +
  theme_classic() +
  theme_bigfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    axis.text.x = element_blank(),
    axis.title.y = element_text(angle = 0),
    axis.text.y = element_text(size = 13),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    panel.grid.major.y = element_line(size = 1, color = '#CCCCCC'),
    plot.margin = unit(c(0,0,0,0), units = 'cm'),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )

# calc number of sig gains and losses
sig_gl_count <- cancer_wilcoxon_df_sig %>%
  inner_join(arm_ploidy_df_wilcox_mean, by = c('cancer_type', 'arm')) %>%
  pivot_wider(names_from = 'cohort', values_from = 'mean_arm_ploidy') %>%
  mutate(sig_gain_loss = case_when(HMF > PCAWG ~ 'gain',
                                   HMF < PCAWG ~ 'loss',
                                   TRUE ~ 'nix')) %>%
  group_by(arm, sig_gain_loss) %>%
  summarise(n_sig = n(), .groups = 'drop') %>%
  right_join(total_arm_df, by = 'arm') %>%
  complete(., arm, sig_gain_loss, fill = list(n_sig = 0)) %>%
  drop_na(., sig_gain_loss) %>%
  mutate(arm = factor(arm, levels = arm_order)) %>%
  group_by(sig_gain_loss) %>%
  summarise(sum_gl = sum(n_sig))

# plot number of significantly different arm ploidies per cancer type
total_cancer_types <- data.frame(cancer_type_label = cancer_type_order)

sig_diff_barplot_cancer_type <- cancer_wilcoxon_df_sig %>%
  group_by(cancer_type_label) %>%
  summarise(n_sig = n(), .groups = 'drop') %>%
  right_join(total_cancer_types, by = c('cancer_type_label')) %>%
  complete(., cancer_type_label, fill = list(n_sig = 0)) %>%
  mutate(cancer_type_label = factor(cancer_type_label, levels = rev(cancer_type_order))) %>%
  ggplot(., aes(x = cancer_type_label, y = n_sig)) +
  geom_col(color = 'black', fill = 'black', position = 'stack', width = 0.75) +
  scale_y_continuous(expand = c(0.01,0), breaks = seq(0, 20, by = 10)) +
  scale_fill_manual(name = '', values = c('#FF007F', '#009ACD')) +
  labs(
    x = '',
    y = ''
  ) +
  coord_flip() +
  theme_classic() +
  theme_bigfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_text(angle = 0),
    axis.text.y.right = element_text(size = 15, angle = 90, hjust = 0.5),
    panel.grid.major.x = element_line(size = 1, color = '#CCCCCC'),
    plot.margin = unit(c(0.2,0,0,0), units = 'cm'),
    plot.background = element_rect(fill = "transparent",colour = NA))

# ------------- Aneuploidy score plot

# custom function to calculate the Mann Whitney test p-value between diploid proportion of cohorts per cancer type
get_wilcox_table_aneuploidy_score <- function(df_input, cancer_type_input) {
  
  # create a list of two dfs that represent the samples in the cohorts for one cancer type
  wilcox_list <- df_input %>%
    filter(cancer_type == cancer_type_input) %>%
    select(sample_id, cohort, aneuploidy_score) %>%
    group_by(cohort) %>%
    group_split()
  
  # calculate the MW p-value between aneuploidy scores of the cohorts
  wilcox_table <- tidy(wilcox.test(wilcox_list[[1]]$aneuploidy_score, 
                                   wilcox_list[[2]]$aneuploidy_score, 
                                   alternative = 'two.sided',
                                   paired = FALSE,
                                   mu = 0)) %>%
    mutate(cancer_type = cancer_type_input,
           significance_level = if_else(p.value < 0.05, '*', NA_character_)) %>%
    select(cancer_type, method, p.value, significance_level)
  
  return(wilcox_table)
  
}

# prepare data: make ploidy df long format
aneuploidy_score <- arm_ploidy_df %>%
  pivot_longer(., cols = c(chr1p:chrXq), names_to = 'arm', values_to = 'gain_loss') %>%
  mutate(gain_loss = case_when(gain_loss > 0 ~ 1,
                               gain_loss < 0 ~ 1,
                               TRUE ~ 0)) %>%
  group_by(sample_id) %>%
  summarise(aneuploidy_score = sum(gain_loss), .groups = 'drop') %>%
  inner_join(metadata, by = 'sample_id') %>%
  inner_join(sample_size, by = c('cancer_type', 'cancer_type_code')) %>%
  # calc median aneuploidy score per cancer type and cohort
  group_by(cohort, cancer_type) %>%
  mutate(median_aneuploidy_score = median(aneuploidy_score)) %>%
  ungroup()

# calculate the Mann whitney p-value for all cancer types + join this data to the aneuploidy_score table
aneuploidy_wilcox_table <- map_df(cancer_type_vec, ~get_wilcox_table_aneuploidy_score(aneuploidy_score, .x)) %>%
  select(-method)

aneuploidy_score <- aneuploidy_score %>%
  inner_join(aneuploidy_wilcox_table, by = 'cancer_type')

# get the alternating index number for the grey background rectangles
rect_idx_aneuploidy <- seq(2, 22, by = 2)
rect_idx_aneuploidy <- cancer_type_order[rect_idx_aneuploidy]

# plot aneuploidy score violins
aneuploidy_score_violin <- aneuploidy_score %>%
  mutate(cancer_type_label = factor(cancer_type_label, levels = cancer_type_order)) %>%
  ggplot(., aes(x = cohort, y = aneuploidy_score)) +
  geom_rect(data = . %>% 
              filter(cohort == 'HMF') %>%
              filter(cancer_type_label %in% rect_idx_aneuploidy) %>%
              distinct(cancer_type_label), aes(fill = cancer_type_label), 
            inherit.aes = FALSE,
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#CCCCCC') +
  geom_violin(aes(fill = cohort), color = 'black', scale = 'width') +
  geom_point(aes(y = median_aneuploidy_score), size = 2) +
  geom_text(data = . %>% filter(cohort == 'HMF'), aes(label = significance_level, y = 36), size = 8, nudge_x = 0) +
  facet_wrap(cancer_type_label ~ ., ncol = 1, strip.position = 'left') +
  scale_x_discrete(
    labels = c('Hartwig', 'PCAWG'),
    position = 'top',
    expand = c(0.5, 0.5)) +
  scale_y_continuous(limits = c(0, 39), breaks = c(0, 10, 20, 30, 39)) +
  scale_fill_manual(values = rev(cohort_colors)) +
  labs(
    title = 'Aneuploidy\nscore',
    x = '',
    y = ''
  ) +
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

# get the WGD, diploid, median polyclonal and mutated p53 fraction histograms / violins
source(paste0(heatmap_base_dir, '/path/to/plot_wgd_rate.R'))
source(paste0(heatmap_base_dir, '/path/to/plot_p53_mut_rate.R'))
source(paste0(heatmap_base_dir, '/path/to/plot_diploid_proportion.R'))

# plot all plots in one window
assembled_plot <- mean_ploidy_heatmap + plot_spacer() + aneuploidy_score_violin +
  non_wgd_diploid_violin + wgd_histogram +
  mutated_histogram + plot_layout(widths = c(4,0.25,0.5,0.5,0.5,0.5), nrow = 1)

# save SVG files of the plots
svglite(
  filename = paste0(output_dir, '/', current_date, '_sig_gl_barplot.svg'),
  width = 12,
  height = 1
)
sig_diff_barplot
dev.off()

svglite(
  filename = paste0(output_dir, '/', current_date, '_sig_cancer_type_barplot.svg'),
  width = 2,
  height = 12
)
sig_diff_barplot_cancer_type
dev.off()

svglite(
  filename = paste0(output_dir, '/', current_date, '_chromosome_arm_heatmap_wilcoxon_mean_gl_qval0.05_cohort.svg'),
  width = 26,
  height = 12
)
assembled_plot
dev.off()

################################### SUPPLEMENTARY TABLES
# refresh token for gsheets
gs4_auth(email = 'your_email@gmail.com')

# get the correct IDs for supp tables (HMF IDs for HMF samples)
hmf_ids <- read_tsv(paste0(heatmap_base_dir, '/path/to/metadata.tsv')) %>%
  # select relevant columns
  select(sample_id, supp_table_id, cohort, cancer_type)

# 1: primary karyotype
primary_karyotype <- arm_matrix %>%
  rownames_to_column(., var = 'sample_id') %>%
  inner_join(hmf_ids, by = 'sample_id') %>%
  filter(cohort == 'PCAWG') %>%
  select(supp_table_id, genome, everything(), -sample_id, -cohort, -cancer_type) %>%
  rename(sample_id = 'supp_table_id')

# push to gsheets
googlesheets4::write_sheet(primary_karyotype,
                           ss = 'https://docs.google.com/spreadsheets/...',
                           sheet = 'primary_karyotype')

# 2. metastatic karytype
metastatic_karytype <- arm_matrix %>%
  rownames_to_column(., var = 'sample_id') %>%
  inner_join(hmf_ids, by = 'sample_id') %>%
  filter(cohort == 'HMF') %>%
  select(supp_table_id, genome, everything(), -sample_id, -cohort, -cancer_type) %>%
  rename(sample_id = 'supp_table_id')

# push to gsheets
googlesheets4::write_sheet(metastatic_karytype,
                           ss = 'https://docs.google.com/spreadsheets/...',
                           sheet = 'metastatic_karyotype')

# 3. computed diff per cancer type and arm, wilcoxon p-val and cliff's delta
computed_differences <- cancer_wilcoxon_df

# push to gsheets
googlesheets4::write_sheet(computed_differences,
                           ss = 'https://docs.google.com/spreadsheets/...',
                           sheet = 'mean_arm_gl')

# 4. genome instability measures
supp_aneuploidy_score <- aneuploidy_score %>%
  select(sample_id, aneuploidy_score)

genomic_instability_measures <- hmf_ids %>%
  inner_join(supp_aneuploidy_score, by = 'sample_id') %>%
  inner_join(supp_loh_proportion, by = 'sample_id') %>%
  inner_join(supp_wgd, by = 'sample_id') %>%
  inner_join(supp_p53_mut, by = 'sample_id') %>%
  select(-sample_id, -cohort) %>%
  rename(sample_id = 'supp_table_id')

# push to gsheets
googlesheets4::write_sheet(genomic_instability_measures,
                           ss = 'https://docs.google.com/spreadsheets/...',
                           sheet = 'genomic_instability_measures')
