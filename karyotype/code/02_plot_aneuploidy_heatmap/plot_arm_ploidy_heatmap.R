############### Plot chromosome arm ploidy comparison ############### 
# author: remy (sascha)
# date: 17/08/2021
# last updated: 23/11/2022

### Description
# This script prepares the arm ploidy matrix to plot the normalized mean arm gains and losses, and compares these mean ploidy levels
# by a certain binary strata which can be specified (e.g. cohort). A pairwise Mann Whitney U test is performed
# per cancer type and arm between the two groups in the specified strata. P-values are FDR adjusted per cancer type.
# Significant q-values are plotted on top of the heatmap tiles to highlight sig. differences. Additionally,
# a barplot which shows the number of cancer types with a sig. gain / loss in that particular arm is plotted on top of the heatmap.
# This heatmap is then stitched together barplots and violins that compare genomic instability measures between the specified groups 
# (WGD rate, mut TP53 rate, aneuploidy score).

### Input
# metadata.tsv
# arm ploidy matrix incl estimated genome ploidy

### Output
# Figure 1d, Figure 1e

# libraries
here::set_here(path='../../')
source(paste0(here::here(), '/code/r_objects/libs.R'))

#========= Path prefixes =========#
heatmap_base_dir <- list(
  path=paste0(here::here(), '/data')
)

for(i in heatmap_base_dir){
  if(dir.exists(i)){
    heatmap_base_dir <- i
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

# get current date
current_date <- format(Sys.Date(), '%Y_%m_%d')

### define variables
## general
filter_wgd <- F # filter data based on whole genome duplication?
exclude_wgd <- F # T: filter out WGD tumors, F: filter out non-WGD tumors, filter_wgd must be TRUE
by_subtype <- F # compare primary vs met by tumor subtype
by_subtype_met_location <- F # compare prim vs. met by subtype and met location, by_subtype must be TRUE
by_progression <- F # compare primary progression subtypes (stable, progressive/relapse) to metastatic lesions
by_met_location <- F # compare primary vs met location
q_val_cutoff <- 0.01 # significance threshold
## for p53 mutation rate plot
mutated_driver <- 'TP53' # mutated gene (can be any gene from the LINX driver catalog)
driver_results_date <- '2021_12_03' # the date you generated the LINX driver catalog in step 2 of the driver analysis.

# checks
if (sum(c(by_subtype, by_progression, by_met_location)) > 1) {
  stop('Multiple conditionals are TRUE, choose only either by_subtype, by_progression or by_met_location to be TRUE while leaving the others FALSE.')
}

## custom font
font_add_google(
  name = 'Inter',
  family = 'Helvetica'
)

## ggplot custom theme
theme_bigfont <- theme(plot.title = element_text(size = 18),
                       plot.subtitle = element_text(size = 16),
                       plot.caption = element_text(size = 14),
                       axis.text.x = element_text(size = 15),
                       axis.text.y = element_text(size = 15),
                       axis.title = element_text(size = 16),
                       legend.title = element_text(size = 16),
                       legend.text = element_text(size = 14),
                       strip.text.x = element_text(size = 15),
                       strip.text.y = element_text(size = 15))

theme_bigfont_barplot <- theme(plot.title = element_text(size = 16),
                               plot.subtitle = element_text(size = 14),
                               axis.text.x = element_text(size = 15),
                               axis.text.y = element_text(size = 12), 
                               axis.title = element_text(size = 18),
                               legend.title = element_text(size = 16),
                               legend.text = element_text(size = 14),
                               strip.text.x = element_text(size = 15),
                               strip.text.y = element_text(size = 15))

# ----------------- Custom objects and functions

# source the functions to calculate the fishers exact test and wilcoxon test p-values and effect sizes
# from test matrices
source(paste0(here::here(), '/code/00_func/fisher_wilcoxon_matrix_functions.R'))

# load cohort colors
source(paste0(here::here(), '/code/r_objects/color_codes.R'))

# define cancer type order
source(paste0(here::here(), '/code/r_objects/cancer_type_order.R'))

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

# ----------------- METADATA

metadata <- read_tsv(paste0(heatmap_base_dir, '/processed/metadata/metadata.tsv'),
                     col_types = 'icccccccccccccccllllcccdiiildclccllllllliiiiiiiiiii') %>%
  # select relevant columns
  select(supp_table_id, cohort, cancer_type, cancer_type_code, 
         cancer_subtype, whole_genome_duplication, progression_status_code, is_blacklisted_subtype, metastatic_location) %>%
  # add has_subtype column for easier filtering
  mutate(has_subtype = if_else(cancer_type_code %in% c('BRCA', 'COREAD', 'UCEC'), TRUE, FALSE)) %>%
  rename(sample_id = 'supp_table_id')

# ----------------- SUBTYPES

# do analysis by cancer subtype
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
}

# ---------------------- PRIMARY PROGRESSION

# do analysis by primary progression
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
  
}

# ----------------- ARM PLOIDY MATRIX

arm_matrix <- read_tsv(
  file = paste0(heatmap_base_dir, '/processed/arm_ploidy_matrix.tsv')) %>%
  column_to_rownames(., var = 'sample_id') %>%
  select(chr1p:chrXq, genome) %>%
  as.matrix()

# subset matrix, arm and genome ploidy
arm_ploidy_mat <- arm_matrix[, 1:(ncol(arm_matrix)-1)]
genome_ploidy_vec <- arm_matrix[, ncol(arm_matrix)]

# subtract genome ploidy level from arm level, row by row
# old matrix & new matrix:
arm_ploidy_mat <- sweep(arm_ploidy_mat, MARGIN = 1, FUN = '-', genome_ploidy_vec)
# binarize gains and losses, column wise
arm_ploidy_mat <- apply(arm_ploidy_mat, MARGIN = 2, FUN = function(x) ifelse(x >= 0, ifelse(x == 0, 0, 1), -1))

# transform into dataframe, add rownames as a separate column
arm_ploidy_df1 <- as.data.frame(arm_ploidy_mat)
arm_ploidy_df1 <- arm_ploidy_df1 %>%
  mutate(sample_id = rownames(arm_ploidy_df1), .before = chr1p)

# join metadata to that long table
arm_ploidy_df <- arm_ploidy_df1 %>%
  inner_join(metadata, by = 'sample_id')

# --------------- Splitting the dataset
# in order to stratify the dataset and perform the tests by strata, this step transforms the strata column into a 'cohort'
# column which is then subsequently used by functions downstream
# default strata is the PCAWG vs HMF cohort column
metadata <- metadata %>%
  mutate(cohort = if_else(
    cohort == 'Hartwig',
    'Hartwig', 'PCAWG'))

# specify the strata you created (used in plot subtitle)
strata <- 'cohort'

# calculate sample size
if (filter_wgd) {
  if (exclude_wgd) {
    sample_size <- arm_ploidy_df %>%
      filter(!whole_genome_duplication) %>%
      group_by(cancer_type, cancer_type_code, cohort) %>%
      summarise(sample_size_tot = n(), .groups = 'drop') %>%
      pivot_wider(names_from = c('cohort'), values_from = 'sample_size_tot', names_prefix = 'n_') %>%
      replace_na(., replace = list(n_Hartwig = 0, n_PCAWG = 0)) %>%
      mutate(cancer_type_label = str_c(cancer_type, ' (', n_PCAWG, ' vs. ', n_Hartwig, ')')) %>%
      mutate(cancer_type_code_label = str_c(cancer_type_code, ' (', n_PCAWG, ' vs. ', n_Hartwig, ')')) %>%
      ungroup()
  } else {
    sample_size <- arm_ploidy_df %>%
      filter(whole_genome_duplication) %>%
      group_by(cancer_type, cancer_type_code, cohort) %>%
      summarise(sample_size_tot = n(), .groups = 'drop') %>%
      pivot_wider(names_from = c('cohort'), values_from = 'sample_size_tot', names_prefix = 'n_') %>%
      replace_na(., replace = list(n_Hartwig = 0, n_PCAWG = 0)) %>%
      mutate(cancer_type_label = str_c(cancer_type, ' (', n_PCAWG, ' vs. ', n_Hartwig, ')')) %>%
      mutate(cancer_type_code_label = str_c(cancer_type_code, ' (', n_PCAWG, ' vs. ', n_Hartwig, ')')) %>%
      ungroup()
  }
} else {
  sample_size <- arm_ploidy_df %>%
    group_by(cancer_type, cancer_type_code, cohort) %>%
    summarise(sample_size_tot = n(), .groups = 'drop') %>%
    pivot_wider(names_from = c('cohort'), values_from = 'sample_size_tot', names_prefix = 'n_') %>%
    replace_na(., replace = list(n_Hartwig = 0, n_PCAWG = 0)) %>%
    mutate(cancer_type_label = str_c(cancer_type, ' (', n_PCAWG, ' vs. ', n_Hartwig, ')')) %>%
    mutate(cancer_type_code_label = str_c(cancer_type_code, ' (', n_PCAWG, ' vs. ', n_Hartwig, ')')) %>%
    ungroup()
}

# sample size filter
sample_size_filter <- sample_size %>%
  filter(pmin(n_Hartwig, n_PCAWG) > 4) %>%
  pull(cancer_type)

if (!filter_wgd) {
  arm_ploidy_df <- arm_ploidy_df %>%
    filter(cancer_type %in% sample_size_filter) 
}

# ----------------- Prepare matrix for statistical testing: Wilcoxon
# rename
arm_ploidy_df_wilcox <- arm_ploidy_df

# filter (non-)WGD tumors
if(filter_wgd) {
  if (exclude_wgd) {
    arm_ploidy_df_wilcox <- arm_ploidy_df_wilcox %>%
      filter(!whole_genome_duplication)
  } else {
    arm_ploidy_df_wilcox <- arm_ploidy_df_wilcox %>%
      filter(whole_genome_duplication)
  }
}

# custom function to calculate the wilcoxon test pvalue and cliff delta per cancer_type for each 
# chromosome arm (primary vs. metastatic)
get_wilcoxon_table <- function(df_input, cancer_type_input) {
  
  print(paste0('processing...', cancer_type_input))
 
  # filter the table for the desired cancer_type
  arm_ploidy_df_wilcox <- arm_ploidy_df_wilcox %>%
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
  # [[1]] contains the cohort == 'Hartwig' (metastatic)
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

# get the cancer type names into a vector, but exclude one of the prostate primary 
cancer_type_vec <- unique(arm_ploidy_df$cancer_type)

# perform mann-whitney U test for comparison of mean ploidy between cohorts
cancer_wilcoxon_list <- map(cancer_type_vec, ~get_wilcoxon_table(arm_ploidy_df_wilcox, 
                                                                 cancer_type_input = .x))
# combine the individual list elements into on table
cancer_wilcoxon_df <- do.call(rbind, cancer_wilcoxon_list)

# only keep significant p-values
cancer_wilcoxon_df_sig <- cancer_wilcoxon_df %>%
  filter(wilcoxon_pval_adj < q_val_cutoff)

# prepare the data for plotting: calculate the mean arm ploidy
arm_ploidy_df_wilcox_mean <- arm_ploidy_df_wilcox %>%
  group_by(cancer_type, cohort) %>%
  summarise(across(chr1p:chr22q, ~mean(.x)), .groups = 'drop') %>%
  filter(cancer_type %in% cancer_type_vec)

arm_ploidy_df_wilcox_mean <- arm_ploidy_df_wilcox_mean %>%
  pivot_longer(cols = chr1p:chr22q, names_to = 'arm', values_to = 'mean_arm_ploidy') %>%
  pivot_wider(., names_from = 'cohort', values_from = 'mean_arm_ploidy') %>%
  mutate(mean_arm_ploidy_diff = abs(Hartwig - PCAWG)) %>%
  pivot_longer(., cols = c('Hartwig', 'PCAWG'), names_to = 'cohort', values_to = 'mean_arm_ploidy')

# join the sample size to the arm df
mean_ploidy_plot_wilcox <- sample_size %>%
  select(cancer_type, cancer_type_label, cancer_type_code_label) %>%
  inner_join(arm_ploidy_df_wilcox_mean, by = c('cancer_type'))

# get the right cancer type label order
cancer_type_label_order <- mean_ploidy_plot_wilcox %>%
  select(cancer_type, cancer_type_label) %>%
  distinct() %>%
  { if (sum(by_subtype, by_met_location, by_progression) == 0) { arrange(., factor(cancer_type, levels = cancer_type_order)) } else { arrange(., cancer_type) } }
cancer_type_label_order <- cancer_type_label_order$cancer_type_label

# get cancer type label to finalize the sig table for manual plot annotation
cancer_wilcoxon_df_sig <- mean_ploidy_plot_wilcox %>%
  select(cancer_type, arm, cancer_type_label, mean_arm_ploidy_diff) %>%
  distinct()  %>%
  right_join(cancer_wilcoxon_df_sig, by = c('cancer_type', 'arm')) %>%
  mutate(arm = factor(arm, levels = arm_order),
         cancer_type_label = factor(cancer_type_label, levels = cancer_type_label_order))

# add the significance * to each data point (for the plot annotation later)
cancer_wilcoxon_df_sig <- cancer_wilcoxon_df_sig %>%
  mutate(significance_level = case_when(
    wilcoxon_pval_adj < q_val_cutoff & mean_arm_ploidy_diff > 0.25 & cliff_delta > 0 ~ '*',
    wilcoxon_pval_adj < q_val_cutoff & mean_arm_ploidy_diff > 0.25 & cliff_delta < 0 ~ '*',
    TRUE ~ 'none')) %>%
  # make sure to get rid of non-significant rows
  filter(significance_level != 'none')

# colors
tile_colors <- c('#000080','#009ACD', 
                 '#B2DFEE',
                 '#FFFFFF',
                 '#FF758C','#FF007F','#A8152E')

# breaks
breaks <- c(seq(-1, -0.1, by = 0.1), seq(0.1, 1, by = 0.1))

# plot karyotype heatmap
mean_ploidy_heatmap <- mean_ploidy_plot_wilcox %>%
  mutate(arm = factor(arm, levels = arm_order),
         cancer_type_label = factor(cancer_type_label, levels = cancer_type_label_order),
         stadium = if_else(cohort == 'Hartwig', 'Metastatic', 'Primary'),
         stadium = factor(stadium, levels = c('Metastatic', 'Primary'))) %>%
  # define neutral bin as NA values for mean_arm_ploidy between -0.1 and 0.1 (important because color stepsn function keeps adjusting gradient color)
  mutate(mean_arm_ploidy_na = if_else(mean_arm_ploidy >= -0.1 & mean_arm_ploidy <= 0.1, NA_real_, mean_arm_ploidy)) %>%
  ggplot(., aes(x = arm, y = stadium, fill = mean_arm_ploidy_na)) +
  geom_tile(color = '#FFFFFF', size = 0.4) +
  # wilcoxon
  { if (nrow(cancer_wilcoxon_df_sig) != 0)
    geom_text(data = cancer_wilcoxon_df_sig,
              aes(y = 1.5, fill = NULL, label = significance_level),
              fontface = 'bold', size = 8, nudge_y = -0.25)
  } +
  facet_grid(cancer_type_label ~ ., switch = 'y') +
  scale_fill_stepsn(colours = tile_colors, na.value = '#e5e5e5',
                    breaks = breaks,
                    limits = c(-1, 1),
                    guide = guide_colorsteps(
                      title = 'Normalized arm ploidy gains / losses relative to\nestimated genome ploidy',
                      title.position = 'top',
                      barwidth = 35,
                      show.limits = TRUE)
                    ) +
  scale_x_discrete(labels = str_replace(arm_order, pattern = 'chr(.+)', replacement = '\\1')) +
  scale_y_discrete(position = 'right') +
  labs(
    caption = '',
    x = 'Chromosome arm',
    y = '') +
  theme_bigfont +
  theme(
    text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y.right = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0,margin = margin(l = 15, unit = 'mm'), hjust = 1),
    legend.position = 'bottom',
    legend.text = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.margin = margin(r = 0, l = 0, unit = 'mm')
  )

if (nrow(cancer_wilcoxon_df_sig) != 0) {
  # create arm level sig difference barplot
  total_arm_df <- data.frame(arm = arm_order)
  
  sig_diff_barplot <- cancer_wilcoxon_df_sig %>%
    inner_join(arm_ploidy_df_wilcox_mean, by = c('cancer_type', 'arm')) %>%
    pivot_wider(names_from = 'cohort', values_from = 'mean_arm_ploidy') %>%
    mutate(sig_gain_loss = case_when(Hartwig > PCAWG & Hartwig > 0 ~ 'gain',
                                     Hartwig < PCAWG & Hartwig < 0 ~ 'loss',
                                     TRUE ~ 'nix')) %>%
    filter(sig_gain_loss != 'nix') %>%
    group_by(arm, sig_gain_loss) %>%
    summarise(n_sig = n(), .groups = 'drop') %>%
    right_join(total_arm_df, by = 'arm') %>%
    complete(., arm, sig_gain_loss, fill = list(n_sig = 0)) %>%
    drop_na(., sig_gain_loss) %>%
    mutate(arm = factor(arm, levels = arm_order))
  
  height <- max(sig_diff_barplot$n_sig) + 1.2
  
  # plot
  sig_diff_barplot <- sig_diff_barplot %>%
    ggplot(., aes(x = arm, y = n_sig)) +
    geom_col(aes(fill = sig_gain_loss), color = 'black', position = 'stack') +
    scale_y_continuous(expand = c(0.01,0), limits = c(0,height), breaks = seq(1, 10, by = 2)) +
    scale_fill_manual(name = '', values = c(`gain` = '#FF007F', `loss` = '#009ACD')) +
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
  
  # calc number of sig gains and losses (for the text)
  sig_gl_count <- cancer_wilcoxon_df_sig %>%
    inner_join(arm_ploidy_df_wilcox_mean, by = c('cancer_type', 'arm')) %>%
    pivot_wider(names_from = 'cohort', values_from = 'mean_arm_ploidy') %>%
    mutate(sig_gain_loss = case_when(Hartwig > PCAWG & Hartwig > 0 ~ 'gain',
                                     Hartwig < PCAWG & Hartwig < 0 ~ 'loss',
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
  total_cancer_types <- data.frame(cancer_type = unique(mean_ploidy_plot_wilcox$cancer_type))
  
  sig_diff_barplot_cancer_type <- cancer_wilcoxon_df_sig %>%
    group_by(cancer_type) %>%
    summarise(n_sig = n(), .groups = 'drop') %>%
    right_join(total_cancer_types, by = c('cancer_type')) %>%
    complete(., cancer_type, fill = list(n_sig = 0)) %>%
    { if (sum(by_subtype, by_met_location, by_progression) == 0) mutate(., cancer_type = factor(cancer_type, levels = rev(cancer_type_order))) else arrange(., cancer_type) } %>%
    ggplot(., aes(x = cancer_type, y = n_sig)) +
    geom_col(color = 'black', fill = 'black', position = 'stack', width = 0.75) +
    scale_y_continuous(expand = c(0.01,0), breaks = seq(0, 20, by = 3)) +
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
      plot.background = element_rect(fill = "transparent", colour = NA))
}

# ----------------- Aneuploidy score plot

# custom function to calculate the Mann Whitney test p-value between diploid proportion of cohorts per cancer type
get_wilcox_table_aneuploidy_score <- function(df_input, cancer_type_input) {
  
  print(paste0('processing...', cancer_type_input))
  
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
           significance_level = if_else(p.value < 0.01, '*', NA_character_)) %>%
    select(cancer_type, method, p.value, significance_level)
  
  return(wilcox_table)
  
}

# prepare data: make ploidy df long format
aneuploidy_score <- arm_ploidy_df_wilcox %>%
  pivot_longer(., cols = c(chr1p:chrXq), names_to = 'arm', values_to = 'gain_loss') %>%
  # filter for relevant arms only
  filter(arm %in% arm_order) %>%
  mutate(gain_loss = case_when(gain_loss > 0 ~ 1,
                               gain_loss < 0 ~ 1,
                               TRUE ~ 0)) %>%
  group_by(sample_id, cancer_type) %>%
  summarise(aneuploidy_score = sum(gain_loss), .groups = 'drop') %>%
  inner_join(metadata, by = c('sample_id', 'cancer_type')) %>%
  inner_join(sample_size, by = c('cancer_type', 'cancer_type_code')) %>%
  # calc median aneuploidy score per cancer type and cohort
  group_by(cohort, cancer_type) %>%
  mutate(median_aneuploidy_score = median(aneuploidy_score)) %>%
  ungroup()

# calculate the Mann whitney (FDR adjusted) p-values for all cancer types + join this data to the aneuploidy_score table
aneuploidy_wilcox_table <- map_df(cancer_type_vec, ~get_wilcox_table_aneuploidy_score(aneuploidy_score, .x)) %>%
  select(-method) %>%
  mutate(pval_adj = p.adjust(p.value, method = 'fdr'),
         significance_level = if_else(pval_adj < q_val_cutoff, '*', NA_character_))

aneuploidy_score <- aneuploidy_score %>%
  inner_join(aneuploidy_wilcox_table, by = 'cancer_type')

# get the alternating index number for the grey background rectangles
if (by_subtype | by_met_location) {
  rect_idx_aneuploidy <- seq(2, n_distinct(aneuploidy_score$cancer_type_label), by = 2)
  rect_idx_aneuploidy <- sort(unique(aneuploidy_score$cancer_type_label))[rect_idx_aneuploidy]
} else {
  rect_idx_aneuploidy <- seq(2, 22, by = 2)
  rect_idx_aneuploidy <- cancer_type_label_order[rect_idx_aneuploidy]
}

# plot aneuploidy violins
aneuploidy_score_violin <- aneuploidy_score %>%
  mutate(cancer_type_label = factor(cancer_type_label, levels = cancer_type_label_order)) %>%
  ggplot(., aes(x = cohort, y = aneuploidy_score)) +
  geom_rect(data = . %>% 
              filter(cohort == 'Hartwig') %>%
              filter(cancer_type_label %in% rect_idx_aneuploidy) %>%
              distinct(cancer_type_label), aes(fill = cancer_type_label), 
            inherit.aes = FALSE,
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#a0a0a0') +
  geom_violin(aes(fill = cohort), color = 'black', scale = 'width') +
  geom_point(
    data = . %>%
      distinct(cancer_type_label, cohort, median_aneuploidy_score),
    aes(y = median_aneuploidy_score), size = 2) +
  geom_text(data = . %>% filter(cohort == 'Hartwig') %>% distinct(cancer_type_label, cohort, significance_level), 
            aes(label = significance_level, y = 36), size = 8, nudge_x = 0) +
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
source(paste0(here::here(), '/code/02_plot_aneuploidy_heatmap/plot_wgd_rate.R'))
source(paste0(here::here(), '/code/02_plot_aneuploidy_heatmap/plot_p53_mut_rate.R'))
source(paste0(here::here(), '/code/02_plot_aneuploidy_heatmap/plot_loh_proportion.R'))

# plot all plots in one window
assembled_plot <- mean_ploidy_heatmap + plot_spacer() + aneuploidy_score_violin +
  loh_violin + wgd_histogram +
  mutated_histogram + plot_layout(widths = c(4,0.25,0.5,0.5,0.5,0.5), nrow = 1)

if (filter_wgd) {
  assembled_plot <- mean_ploidy_heatmap + plot_spacer() + aneuploidy_score_violin +
    loh_violin + mutated_histogram + plot_layout(widths = c(4,0.25,0.5,0.5,0.5), nrow = 1)
}

# create output dir
output_dir <- paste0(output_dir, '/02_plot_aneuploidy_heatmap/', current_date, '/results/plots')
dir.create(path = paste0(output_dir), recursive = TRUE)

# save plots
# as PDF
if (nrow(cancer_wilcoxon_df_sig) != 0) {
  pdf(
    file = paste0(output_dir, '/', current_date, '_sig_gl_barplot',
                  if (by_subtype) { 
                    if (by_subtype_met_location) { '_subtype_and_met_location.pdf' } 
                    else { '_subtype.pdf' } } 
                  else if (by_met_location) { '_met_location.pdf' } 
                  else if (by_progression) { '_progression.pdf' }
                  else { if (filter_wgd & exclude_wgd) { '_nonWGD.pdf' }
                    else if (filter_wgd) { '_WGD.pdf' }
                    else { '.pdf' } }),
    width = 12,
    height = 1,
    useDingbats = FALSE,
    compress = FALSE
  )
  print(sig_diff_barplot)
  dev.off()
  
  pdf(
    file = paste0(output_dir, '/', current_date, '_sig_cancer_type_barplot',
                  if (by_subtype) {
                    if (by_subtype_met_location) { '_subtype_and_met_location.pdf' }
                    else { '_subtype.pdf' } }
                  else if (by_met_location) { '_met_location.pdf' } 
                  else if (by_progression) { '_progression.pdf' }
                  else { if (filter_wgd & exclude_wgd) { '_nonWGD.pdf' }
                    else if (filter_wgd) { '_WGD.pdf' }
                    else { '.pdf' } }),
    width = 2,
    height = 12,
    useDingbats = FALSE,
    compress = FALSE
  )
  print(sig_diff_barplot_cancer_type)
  dev.off()
}

pdf(
  file = paste0(output_dir, '/', current_date, '_chromosome_arm_heatmap_wilcoxon_mean_gl_qval', q_val_cutoff, '_cohort',
                    if (by_subtype) { 
                      if (by_subtype_met_location) { '_subtype_and_met_location.pdf'} 
                      else { '_subtype.pdf' } } 
                    else if (by_met_location) { '_met_location.pdf' }
                    else if (by_progression) { '_progression.pdf' }
                    else { if (filter_wgd & exclude_wgd) { '_nonWGD.pdf' }
                      else if (filter_wgd) { '_WGD.pdf' }
                      else { '.pdf' } }),
  width = 26,
  height = if (by_subtype) { if (by_subtype_met_location) { 14 } else { 10 } }
  else if (by_met_location) { 30 }
  else if (by_progression) { 6 }
  else { 12 },
  useDingbats = FALSE,
  compress = FALSE
)
print(assembled_plot)
dev.off()