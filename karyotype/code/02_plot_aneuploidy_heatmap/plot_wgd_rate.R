############### Plot whole genome duplication rate histograms ###############
# author: remy (sascha)
# date: 30/08/2021
# last updated: 02/11/2022

### Description
# This script creates the barplots that compare the WGD rate stratified by a certain strata (e.g. cohort) that is specified
# by the user under 'Splitting the dataset'. These barplots are combined with the mean arm ploidy heatmap.

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
                               axis.text.y = element_text(size = 13), 
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

# sample size filter
metadata <- metadata %>%
  filter(cancer_type %in% sample_size_filter)

# --------------- Splitting the dataset
# in order to stratify the dataset and perform the tests by a binary strata, this step transforms the strata column into a 'cohort'
# column which is then subsequently used by functions downstream
# DEFAULT strata is the PCAWG vs Hartwig 'cohort' column
metadata <- metadata %>%
  mutate(cohort = if_else(
    cohort == 'Hartwig',
    'Hartwig', 'PCAWG'))

# count the number of WGD samples per cancer type and stadium
metadata <- metadata %>%
  group_by(cancer_type, cohort, whole_genome_duplication) %>%
  mutate(n_wgd = n()) %>%
  ungroup()

metadata_wide <- metadata %>%
  select(cancer_type, cohort, whole_genome_duplication, n_wgd) %>%
  distinct() %>%
  pivot_wider(names_from = whole_genome_duplication, names_prefix = 'wgd_', values_from = n_wgd) %>%
  replace_na(., replace = list(wgd_TRUE = 0,
                               wgd_FALSE = 0))

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
  rename(metastatic_wgd = 'wgd_true_Hartwig',
         metastatic_non_wgd = 'wgd_false_Hartwig',
         primary_wgd = 'wgd_true_PCAWG',
         primary_non_wgd = 'wgd_false_PCAWG') %>%
  select(cancer_type, primary_wgd, primary_non_wgd, metastatic_wgd, metastatic_non_wgd)

# prepare matrix for fisher pval calculation
metadata_wider_wgd_combined_mat <- metadata_wider_wgd_combined %>% 
  column_to_rownames(var = 'cancer_type') %>% 
  as.matrix()

# calculate fisher, odds ratio and cramers v + add these as separate columns
metadata_wider_wgd_combined$fisher_pval <- fisherTest.matrix(metadata_wider_wgd_combined_mat)
metadata_wider_wgd_combined$fisher_pval_adj <- p.adjust(metadata_wider_wgd_combined$fisher_pval, method = 'fdr')

# get the alternating index number for the grey background rectangles
if (by_subtype | by_met_location | by_progression) {
  rect_idx <- seq(2, n_distinct(metadata_wider_wgd_combined$cancer_type), by = 2)
  rect_idx <- sort(unique(metadata_wider_wgd_combined$cancer_type))[rect_idx]
} else {
  rect_idx <- seq(2, 22, by = 2)
  rect_idx <- cancer_type_order[rect_idx]
}

metadata_wider_wgd_combined_plot <- metadata_wider_wgd_combined %>%
  mutate(primary_wgd_fraction = primary_wgd / (primary_wgd + primary_non_wgd),
         metastatic_wgd_fraction = metastatic_wgd / (metastatic_wgd + metastatic_non_wgd)) %>%
  { if (sum(by_subtype, by_met_location, by_progression) == 0) mutate(., cancer_type = factor(cancer_type, levels = cancer_type_order)) else arrange(., cancer_type) } %>%
  pivot_longer(cols = primary_wgd_fraction:metastatic_wgd_fraction, names_to = 'stadium_wgd', values_to = 'fraction') %>%
  mutate(significance_level = if_else(fisher_pval_adj < q_val_cutoff & stadium_wgd == 'primary_wgd_fraction', '*', NA_character_))

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