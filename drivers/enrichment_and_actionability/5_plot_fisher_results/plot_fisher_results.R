############### Investigate LINX driver output ###############
# author: remy (sascha)
# date: 21/07/2021
# last updated: 06/04/2022

### Description
# This script creates the dotplot heatmap and volcano plot for the driver prevalence analysis.

# libraries
library(tidyverse) # data manipulation and plotting
library(ggthemes) # plot theme
library(ggrepel) # handling labels
library(ggnewscale) # mutliple plot color scales
library(naniar) # replace with NA function
library(svglite) # handle inkscape SVG issue
library(googlesheets4) # push supp table to gsheets

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
  local='output/path',
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
theme_bigfont_heatmap <- theme(plot.title = element_text(size=22),
                               plot.subtitle = element_text(size = 16),
                               axis.text.x= element_text(size=15),
                               axis.text.y= element_text(size=10), 
                               axis.title=element_text(size=18),
                               legend.text = element_text(size = 12),
                               strip.text.x = element_text(size = 15))

theme_bigfont_volcano <- theme(plot.title = element_text(size=22),
                               plot.subtitle = element_text(size = 16),
                               axis.text.x= element_text(size=15),
                               axis.text.y= element_text(size=15), 
                               axis.title=element_text(size=18),
                               legend.text = element_text(size = 12),
                               strip.text.x = element_text(size = 15))

theme_smallfont <- theme(plot.title = element_text(size=17),
                         plot.subtitle = element_text(size = 11),
                         axis.text.x= element_text(size=10),
                         axis.text.y= element_text(size=8), 
                         axis.title=element_text(size=13),
                         legend.text = element_text(size = 12),
                         strip.text.x = element_text(size = 10))


#  ----------------------- Custom functions

# source the functions to calculate the fishers exact test and effect sizes (cramer's V and odds Ratio)
# from test matrices
source(paste0(base_dir, '/0_func/fisher_wilcoxon_matrix_functions.R'))

#  ----------------------------------------

# define variables
contingency_matrix_results_date <- '2021_12_06' # old: 2021_07_30, new (with intogen filter): 2021_11_18, PANEL == 70: 2021_12_06, PANEL == 90: 2021_12_07
fusion_contingency_matrix_results_date <- '2021_12_06' # old: 2021_08_04, new: 2021_10_26, PANEL == 70: 2021_12_06, PANEL == 90: 2021_12_07
q_val_cutoff <- 0.01 # any value between 0 - 1
exclude_hypermutators <- FALSE # T/F
plot_clonality <- FALSE # T/F
clonality <- 'clonal' # clonal/subclonal
strata <- 'cohort' # basically any binary strata column you have in the dataset

# ---------------------- METADATA

metadata <- read_tsv(paste0(base_dir, '/path/to/metadata.tsv')) %>%
  # select relevant columns
  select(sample_id, cohort, is_metastatic, is_hypermutated, cancer_type, cancer_type_code)

# ---------------------- Sample size stratified, calculated from metadata

# calculate sample size by strata
sample_size <- metadata %>%
  group_by(cancer_type, cancer_type_code, across(all_of(strata))) %>%
  summarise(sample_size = n(), .groups = 'drop') %>%
  pivot_wider(names_from = strata, values_from = 'sample_size') %>%
  filter(!is.na(PCAWG),
         !is.na(HMF)) %>%
  # filter out any cancer with less than 15 samples in either group
  filter(pmin(HMF, PCAWG) > 14)

################################################################################################################
################################### CONTINGENCY MATRIX ONLY ####################################################

if (exclude_hypermutators) {
  # read in the generated contingency matrices
  drivers_contingency_matrix_df <- read_tsv(paste0(base_dir, '/path/to/', contingency_matrix_results_date ,'/', 
                                                   contingency_matrix_results_date,  '_', clonality, '_driver_contingency_matrix_no_hyp.tsv'),
                                            col_names = TRUE)
  
  fusion_contingency_matrix_df <- read_tsv(paste0(base_dir, '/path/to/', fusion_contingency_matrix_results_date ,'/', 
                                                  fusion_contingency_matrix_results_date,  '_fusion_contingency_matrix_no_hyp.tsv'),
                                           col_names = TRUE)
} else {
  if (plot_clonality) {
    # read in the generated contingency matrices
    drivers_contingency_matrix_df <- read_tsv(paste0(base_dir, '/path/to/', contingency_matrix_results_date ,'/', 
                                                     contingency_matrix_results_date,  '_', clonality, '_driver_contingency_matrix_with_hyp.tsv'),
                                              col_names = TRUE) 
    
    fusion_contingency_matrix_df <- read_tsv(paste0(base_dir, '/path/to/', fusion_contingency_matrix_results_date ,'/', 
                                                    fusion_contingency_matrix_results_date,  '_fusion_contingency_matrix_with_hyp.tsv'),
                                             col_names = TRUE) 
  } else {
    # read in the generated contingency matrices
    drivers_contingency_matrix_df <- read_tsv(paste0(base_dir, '/path/to/', contingency_matrix_results_date ,'/', 
                                                     contingency_matrix_results_date, '_driver_contingency_matrix_with_hyp.tsv'),
                                              col_names = TRUE) 
    
    fusion_contingency_matrix_df <- read_tsv(paste0(base_dir, '/path/to/', fusion_contingency_matrix_results_date ,'/', 
                                                    fusion_contingency_matrix_results_date,  '_fusion_contingency_matrix_with_hyp.tsv'),
                                             col_names = TRUE)
  }
}

# combine both matrices
contingency_matrix_df <- rbind(drivers_contingency_matrix_df, fusion_contingency_matrix_df)

# filter out cases where there are less than 4 mutated samples in both primary or metastatic group
min_freq_idx <- pmax(contingency_matrix_df$mut_false_group, contingency_matrix_df$mut_true_group) > 3
contingency_matrix_df <- contingency_matrix_df %>%
  filter(min_freq_idx)

# create the contingency matrix with only the four columns needed in the correct order
contingency_matrix <- contingency_matrix_df %>%
  select(mut_false_group, wt_false_group, mut_true_group, wt_true_group) %>%
  as.matrix(.)

# calculate the fisher p-value, the odds ratio and cramers V from that matrix, row by row
contingency_matrix_df$cm_fisher_pval <- fisherTest.matrix(contingency_matrix)
contingency_matrix_df$cm_odds_ratio <- oddsRatio.matrix(contingency_matrix)
contingency_matrix_df$cm_cramers_v <- cramerV.matrix(contingency_matrix)

# adjust p-value for multiple testing across all the rows using the FDR method
contingency_matrix_df <- contingency_matrix_df %>%
  group_by(cancer_type, driver) %>%
  mutate(cm_fisher_pval_adj = p.adjust(cm_fisher_pval, method = 'fdr'), .after = cm_fisher_pval) %>%
  ungroup()

# define cancer type order
source(paste0(base_dir, '/path/to/r_objects/cancer_type_order.R'))

# add the sample size as a separate column (for the dotplot heatmap)
contingency_matrix_df <- contingency_matrix_df %>%
  inner_join(sample_size, by = 'cancer_type') %>%
  mutate(cancer_type_label = paste0(cancer_type, ' (', PCAWG, ' vs. ', HMF, ')'),
         sample_size = PCAWG + HMF) %>%
  mutate(cancer_type = factor(cancer_type, levels = cancer_type_order)) %>%
  # multiply cramers V by -1 so that the metastatic cancers get the positive numbers
  mutate(cm_cramers_v = -1 * cm_cramers_v)

# filter for significantly enriched drivers (either in primary or metastatic group)
# q_val_cutoff <- 0.01
contingency_matrix_df_sig <- contingency_matrix_df %>%
  filter(cm_fisher_pval_adj < q_val_cutoff)

# prepare the df for plotting
contingency_matrix_df_sig <- contingency_matrix_df_sig %>%
  # add exclusivity column to show which drivers are exclusive to either prim or met
  mutate(exclusivity = case_when(mut_false_group == 0 ~ 'Hartwig only',
                                 mut_true_group == 0 ~ 'PCAWG only',
                                 TRUE ~ 'Both'))

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

if (exclude_hypermutators) {
  plot_subtitle2 <- 'Excluding hypermutators'
} else {
  plot_subtitle2 <- 'Including hypermutators'
}

# import cohort color
source(paste0(base_dir, '/path/to/r_objects/color_codes.R'))

########### calculate basic statistics about the results
freq_driver_muts <- contingency_matrix_df_sig %>%
  mutate(freq_pcawg = mut_false_group / wt_false_group * 100,
         freq_hartwig = mut_true_group / wt_true_group * 100, .after = wt_true_group)

enrichement_stats <- contingency_matrix_df_sig %>%
  mutate(enriched_in = if_else(cm_cramers_v > 0, 'enriched in hartwig', 'enriched in PCAWG'))

total_enrichment <- nrow(enrichement_stats)
unique_enriched_driver_genes <- n_distinct(enrichement_stats$gene)
unique_cancer_types <- n_distinct(enrichement_stats$cancer_type)

recurrent_drivers <- enrichement_stats %>%
  group_by(gene) %>%
  summarise(n = n_distinct(cancer_type), .groups = 'drop') %>%
  mutate(percent_cancer_type = n / n_distinct(enrichement_stats$cancer_type) * 100)

percent_enriched_in <- enrichement_stats %>%
  group_by(enriched_in) %>%
  count() %>%
  mutate(percent = n / total_enrichment * 100)

percent_exclusive <- enrichement_stats %>%
  group_by(exclusivity) %>%
  count() %>%
  mutate(percent = n / total_enrichment * 100)

percent_mut_type_cohort <- enrichement_stats %>%
  group_by(driver, enriched_in) %>%
  count() %>%
  group_by(driver) %>%
  mutate(n_per_driver = sum(n)) %>%
  mutate(percent = n / n_per_driver * 100)
########### 

### dotplot heatmap
# export q < 0.01 as: 1200 x 800
# export q < 0.05 as: 1600 x 800
# export q < 0.1 as: 1200 x 800
dotplot_heatmap <- contingency_matrix_df_sig %>%
  ggplot(., aes(x = cancer_type, y = fct_reorder(gene, desc(gene)))) +
  geom_point(aes(size = sample_size,
                 shape = exclusivity,
                 fill = cm_cramers_v), color = 'lightgrey') +
  facet_grid(~driver, scales = 'free_x', space = 'free_x', labeller = as_labeller(facet_labels)) +
  scale_fill_gradient2(low = cohort_colors[1], mid = '#FFFFFF', high = cohort_colors[2], breaks = c(0.5, 0.4, 0.3, 0.2, 0.1, 0, -0.1, -0.2, -0.3, -0.4),
                       guide = guide_colorbar(title = "Eff size (Cramer's V)", order = 1,
                                              show.limits = TRUE, label.position = 'left', barheight = 10)) +
  scale_size_continuous(range = c(2, 6), breaks = c(50, 250, 500, 1000), limits = c(30, 1000),
                        guide = guide_legend(title = 'Sample size', order = 2, override.aes = list(color = 'black'))) + # for q < 0.1
  scale_shape_manual(values = c(21, 24, 25), guide = guide_legend(title = '', order = 3, override.aes = list(size = 4, fill = 'black'))) +
  labs(subtitle = paste0(plot_subtitle, ' | ', plot_subtitle2, ' | Stratified by ', strata, ' | q < ', q_val_cutoff),
       x = '',
       y = '') +
  theme_bw() +
  theme_bigfont_heatmap +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = unit(c(1,1,1,1), 'cm'))

### volcano plot
  
# prepare df for volcano plot
contingency_matrix_df_volcano <- contingency_matrix_df %>%
 { if (q_val_cutoff == 0.01)  filter(., driver != 'FUSION') else . } %>% # optinal filter for prettier display of low q-val plot
  mutate(significance = if_else(cm_fisher_pval_adj < q_val_cutoff, 'yes', 'no'))

# color palette that is assigning a color to a cancer type
source(paste0(base_dir, '/analysis/linx_driver/r_objects/color_codes.R'))
# names(cancer_type_palette) <- cancer_type_order

contingency_matrix_df_sig <- contingency_matrix_df_sig %>%
  mutate(gene_label = paste0(gene, ' (', cancer_type_code, ')'))

# give correct names to colors
names(cancer_type_palette) <- cancer_type_order

# volcano plot: odds ratio
# export q < 0.01 as: 1600 x 800
# export q < 0.05 as: 1600 x 1600
# export q < 0.1 as: 1600 x 1000
# the log2(OR) is multiplied by -1 so the drivers which are more prevalent in metastatic have positive log2(odds)
volcano_plot_odds_ratio <- ggplot(contingency_matrix_df_volcano, aes(x = log2(cm_odds_ratio) * -1, y = -log10(cm_fisher_pval_adj))) +
  geom_point(aes(color = significance, size = abs(cm_cramers_v)),
             ) +
  geom_vline(xintercept = 0, color = '#B3B3B3', linetype = 'dashed') +
  geom_hline(yintercept = -log10(q_val_cutoff), color = '#B3B3B3', linetype = 'dashed') +
  scale_color_manual(values = c('#E5E5E5', '#000000'), guide = 'none') +
  labs(color = 'Significance') +
  new_scale_color() +
  geom_point(
    data = contingency_matrix_df_sig,
    aes(fill = cancer_type, size = abs(cm_cramers_v)),
    color = 'black',
    shape = 21,
    stroke = 1.1
  ) +
  geom_text_repel(data = contingency_matrix_df_sig, 
                  aes(label = gene_label, color = cancer_type), key_glyph = 'rect', size = 5,
                  xlim = c(-8, 8), 
                  max.overlaps = Inf,
                  force = 9, force_pull = 0, 
                  max.time = 10, direction = 'both', 
                  min.segment.length = 0, point.padding = 0.5, box.padding = 0.5) +
  facet_wrap(~driver, labeller = as_labeller(facet_labels)) +
  scale_x_continuous(breaks = seq(-10, 10, by = 2), limits = c(-8, 8)) +
  scale_y_continuous(breaks = c(2, 5, 10, 15, 20), limits = c(0, 20)) +
  scale_color_manual(values = cancer_type_palette) +
  scale_fill_manual(values = cancer_type_palette) +
  labs(subtitle = paste0(plot_subtitle, ' | ', plot_subtitle2, ' | Stratified by ', strata, ' | q < ', q_val_cutoff),
       x = 'Log2(Odds ratio Hartwig / PCAWG)',
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

# volcano plot; cramer's v
volcano_plot_cramersv <- ggplot(contingency_matrix_df_volcano, aes(x = cm_cramers_v, y = -log10(cm_fisher_pval_adj))) +
  geom_point(aes(x = cm_cramers_v, color = significance, size = abs(cm_cramers_v)),
  ) +
  geom_vline(xintercept = 0, color = '#B3B3B3', linetype = 'dashed') +
  geom_hline(yintercept = -log10(q_val_cutoff), color = '#B3B3B3', linetype = 'dashed') +
  scale_color_manual(values = c('#E5E5E5', '#000000'), guide = 'none') +
  labs(color = 'Significance') +
  new_scale_color() +
  geom_point(
    data = contingency_matrix_df_sig,
    aes(fill = cancer_type, size = abs(cm_cramers_v)),
    color = 'black',
    shape = 21,
    stroke = 1.1
  ) +
  geom_text_repel(data = contingency_matrix_df_sig, 
                  aes(label = gene_label, color = cancer_type), key_glyph = 'rect', size = 5,
                  xlim = c(-1.1, 1.1), 
                  max.overlaps = Inf,
                  force = 9, force_pull = 0, 
                  max.time = 10, direction = 'both', 
                  min.segment.length = 0, point.padding = 0.5, box.padding = 0.5) +
  facet_wrap(~driver, labeller = as_labeller(facet_labels)) +
  scale_x_continuous(breaks = seq(-1, 1, by = 0.2), limits = c(-1, 1)) +
  scale_y_continuous(breaks = c(2, 5, 10, 15, 20), limits = c(0, 20)) +
  scale_color_manual(values = cancer_type_palette) +
  scale_fill_manual(values = cancer_type_palette) +
  labs(subtitle = paste0(plot_subtitle, ' | ', plot_subtitle2, ' | Stratified by ', strata, ' | q < ', q_val_cutoff),
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

# save SVG files of the plots
svglite(
  filename = paste0(output_dir, '/', current_date, '_driver_fusion_dotplot_heatmap_minfreq4_q', q_val_cutoff, '_with_hyp.svg'),
  width = 16,
  height = 8
)
dotplot_heatmap
dev.off()

svglite(
  filename = paste0(output_dir, '/', current_date, '_driver_volcanoplot_minfreq4_q', q_val_cutoff, '_with_hyp.svg'),
  width = 16, 
  height = 8
)
volcano_plot_cramersv
dev.off()


################################ SUPPLEMENTARY TABLES
# refresh token for gsheets
gs4_auth(email = 'your_email@gmail.com')

# results from fishers exact test with effect size
fisher_results_supp <- contingency_matrix_df_volcano %>%
  # mutate(cm_cramers_v = cm_cramers_v * -1) %>%
  rename_all(~ gsub('_false_group$', '_pcawg', .x)) %>%
  rename_all(~ gsub('_true_group', '_hartwig', .x)) %>%
  rename_all(~ gsub('^cm_', '', .x)) %>%
  select(cancer_type, cancer_type_code, driver, gene, mut_pcawg, wt_pcawg, mut_hartwig, wt_hartwig, fisher_pval, fisher_pval_adj, cramers_v) %>%
  arrange(cancer_type, driver, gene)


# push to gsheets
googlesheets4::write_sheet(fisher_results_supp,
                           ss = 'https://docs.google.com/spreadsheets/...',
                           sheet = 'statistical_comparison_results')
