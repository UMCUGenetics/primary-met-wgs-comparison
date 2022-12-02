############### Plot diploid proportion histograms ###############
# author: remy (sascha)
# date: 09/11/2021
# last updated: 02/11/2022

### Description
# This script creates the barplots that compare the diploid proportion from loh_proportion files stratified by a certain strata (e.g. cohort) that is specified
# by the user under the 'Splitting the dataset' section. A Mann Whitney test is performed to compare the diploid proportion between strata.
# These barplots are combined with the mean arm ploidy heatmap.

### Input
# see plot_arm_ploidy_heatmap.R

### Output
# see plot_arm_ploidy_heatmap.R

#========= Path prefixes =========#
base_dir <- list(
  path=paste0(here::here(), '/data')
)

for(i in base_dir){
  if(dir.exists(i)){
    base_dir <- i
    break
  }
}

## ggplot custom theme
theme_bigfont_barplot <- theme(plot.title = element_text(size = 16),
                               plot.subtitle = element_text(size = 14),
                               axis.text.x = element_text(size = 15),
                               axis.text.y = element_text(size = 12), 
                               axis.title = element_text(size = 18),
                               legend.title = element_text(size = 16),
                               legend.text = element_text(size = 14),
                               strip.text.x = element_text(size = 15),
                               strip.text.y = element_text(size = 15))

# ----------------- METADATA

metadata <- read_tsv(paste0(base_dir, '/processed/metadata/metadata.tsv'),
                     col_types = 'icccccccccccccccllllcccdiiildclccllllllliiiiiiiiiii') %>%
  # select relevant columns
  select(supp_table_id, cohort, cancer_type, cancer_type_code, 
         cancer_subtype, whole_genome_duplication, progression_status_code, is_blacklisted_subtype, metastatic_location) %>%
  # add has_subtype column for easier filtering
  mutate(has_subtype = if_else(cancer_type_code %in% c('BRCA', 'COREAD', 'UCEC'), TRUE, FALSE)) %>%
  rename(sample_id = 'supp_table_id')

# WGD filter
if (filter_wgd) {
  if (exclude_wgd) {
    metadata <- metadata %>%
      filter(!whole_genome_duplication)
  } else {
    metadata <- metadata %>%
      filter(whole_genome_duplication)
  }
}

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

# sample size filtwer
if (!filter_wgd) {
  metadata <- metadata %>%
    filter(cancer_type %in% sample_size_filter) 
}


# ----------------- LOH PROPORTION (parsed by Arne)

loh_proportion <- read_tsv(
  file = paste0(base_dir, '/processed/LOH.tsv'))

# join loh_proportion table to metadata
metadata <- metadata %>%
  inner_join(loh_proportion, by = 'sample_id')

# --------------- Splitting the dataset
# in order to stratify the dataset and perform the tests by a binary strata, this step transforms the strata column into a 'cohort'
# column which is then subsequently used by functions downstream
# DEFAULT strata is the PCAWG vs Hartwig 'cohort' column
metadata <- metadata %>%
  mutate(cohort = if_else(
    cohort == 'Hartwig',
    'Hartwig', 'PCAWG')) # create a new cohort column

# custom function to calculate the Mann Whitney test p-value between diploid proportion of cohorts per cancer type
get_wilcox_table <- function(df_input, cancer_type_input) {
  
  print(paste0('processing...', cancer_type_input))
  
  # create a list of two dfs that represent the samples in the cohorts for one cancer type
  wilcox_list <- df_input %>%
    filter(cancer_type == cancer_type_input) %>%
    select(sample_id, cohort, loh_proportion) %>%
    group_by(cohort) %>%
    group_split()
  
  if (n_distinct(wilcox_list) == 1) {
    wilcox_table <- data.frame(
      cancer_type = cancer_type_input,
      method = 'Wilcoxon rank sum test with continuity correction',
      p.value = NA_real_,
      significance_level = NA_character_
    )
    
    return(wilcox_table)
  }
  
  # calculate the MW p-value between diploid proportions of the cohorts
  wilcox_table <- tidy(wilcox.test(wilcox_list[[1]]$loh_proportion, 
                                   wilcox_list[[2]]$loh_proportion, 
                                   alternative = 'two.sided',
                                   paired = FALSE,
                                   mu = 0)) %>%
    mutate(cancer_type = cancer_type_input,
           significance_level = if_else(p.value < 0.01, '*', NA_character_)) %>%
    select(cancer_type, method, p.value, significance_level)
  
  return(wilcox_table)
  
}

# get the cancer type names into a vector
cancer_type_vec <- unique(metadata$cancer_type)

# calculate the Mann whitney p-value for all cancer types
wilcox_table <- map_df(cancer_type_vec, ~get_wilcox_table(metadata, .x))

# adjust p-values
wilcox_table <- wilcox_table %>%
  mutate(pval_adj = p.adjust(p.value, method = 'fdr')) %>%
  mutate(significance_level = if_else(pval_adj < q_val_cutoff, '*', NA_character_))
  
# calculate median diploid proportion per cancer type and cohort then add this info to the wilcox_table (for point layer in plot)
median_diploid <- metadata %>%
  group_by(cancer_type, cohort) %>%
  summarise(median_loh_proportion = median(loh_proportion, na.rm = T), .groups = 'drop')
wilcox_table <- wilcox_table %>%
  left_join(median_diploid, by = c('cancer_type'))

# get the alternating index number for the grey background rectangles
if (by_subtype | by_met_location | by_progression) {
  rect_idx <- seq(2, n_distinct(metadata$cancer_type), by = 2)
  rect_idx <- sort(unique(metadata$cancer_type))[rect_idx]
} else {
  rect_idx <- seq(2, 22, by = 2)
  rect_idx <- cancer_type_order[rect_idx]
}

# get the cancer types in the right order
# NOTE: HAVE TO DO THIS FOR ALL TABLES THAT ARE USED IN THE PLOT
metadata <- metadata %>%
  { if (sum(by_subtype, by_met_location, by_progression) == 0) mutate(., cancer_type = factor(cancer_type, levels = cancer_type_order)) else arrange(., cancer_type) } %>%
  mutate(stadium = if_else(cohort == 'Hartwig', 'Hartwig', 'PCAWG'),
         stadium = factor(stadium, levels = c('Hartwig', 'PCAWG')))
wilcox_table <- wilcox_table %>%
  { if (sum(by_subtype, by_met_location, by_progression) == 0) mutate(., cancer_type = factor(cancer_type, levels = cancer_type_order)) else arrange(., cancer_type) } %>%
  mutate(stadium = if_else(cohort == 'Hartwig', 'Hartwig', 'PCAWG'),
         stadium = factor(stadium, levels = c('Hartwig', 'PCAWG')))

# plot the LOH violins
loh_violin <- metadata %>%
  ggplot(., aes(x = stadium, y = loh_proportion)) +
  geom_rect(data = . %>% 
              filter(cancer_type %in% rect_idx) %>%
              distinct(cancer_type), aes(fill = cancer_type),
            inherit.aes = FALSE,
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#a0a0a0') +
  geom_violin(aes(fill = cohort), color = 'black', scale = 'width') +
  geom_point(data = wilcox_table %>% distinct(cancer_type, stadium, median_loh_proportion), 
             aes(y = median_loh_proportion), size = 2) +
  geom_text(data = wilcox_table %>% filter(cohort == 'Hartwig'), 
            aes(label = significance_level, y = 1.1), size = 8, nudge_x = 0) +
  facet_wrap(cancer_type ~ ., ncol = 1, strip.position = 'left') +
  scale_x_discrete(
    labels = c('Hartwig', 'PCAWG'),
    position = 'top',
    expand = c(0.5,0.5)) + # 'expand' is used to center the histogram in the panel
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1.2), expand = c(0,0)) +
  scale_fill_manual(values = rev(cohort_colors)) +
  labs(title = 'LOH\nproportion',
       subtitle = paste0('*Mann Whitney\ntest q-value < ', q_val_cutoff),
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