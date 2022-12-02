############### Plot drivers per sample plots ###############
# author: remy (sascha)
# date: 16/11/2021
# last updated: 24/10/2022

### Description
# This script creates all the drivers per samples plots. 
# The user can choose whether to perform the analysis:
# 1. At cancer type level (by_subtype == FALSE; by_progression == FALSE, by_met_location == FALSE)
# 2. by cancer subtypes (by_subtype == TRUE; by_progression == FALSE, by_met_location == FALSE)
# 3. by primary progression (by_subtype == FALSE; by_progression == TRUE, by_met_location == FALSE)
# 4. by metastatic location (by_subtype == FALSE; by_progression == FALSE, by_met_location == TRUE).

### Input
# metadata.tsv
# complete LINX drivers table from step 6. (optional: filtered for resistance drivers)

### Output
# Figure 5a

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

# create output dir
dir.create(path = paste0(output_dir, '/06_plot_drivers_per_sample'), recursive = TRUE)
output_dir <- paste0(output_dir, '/06_plot_drivers_per_sample')

# get current date
current_date <- format(Sys.Date(), '%Y_%m_%d')

## ggplot custom themes
font_add_google(
  name = 'Inter',
  family = 'Helvetica'
)

## ggplot custom themes
resize_factor <- 1

theme_bigfont <- theme(plot.title = element_text(size = 22),
                       plot.subtitle = element_text(size = 16),
                       axis.text.x = element_text(size = 15),
                       axis.text.y = element_text(size = 15), 
                       axis.title = element_text(size = 18),
                       legend.text = element_text(size = 14),
                       strip.text.x = element_text(size = 15))

theme_smallfont <- theme(plot.title = element_text(size = 17 * resize_factor),
                         plot.subtitle = element_text(size = 11 * resize_factor),
                         axis.text.x = element_text(size = 10 * resize_factor),
                         axis.text.y = element_text(size = 10 * resize_factor), 
                         axis.title = element_text(size = 13 * resize_factor),
                         legend.text = element_text(size = 12 * resize_factor),
                         strip.text = element_text(size = 10 * resize_factor))

#  ----------------------- Custom functions & objects

# define cancer type label order
source(paste0(base_dir, '/code/r_objects/cancer_type_order.R'))

# import cohort color
source(paste0(base_dir, '/code/r_objects/color_codes.R'))

#  ---------------------- define dates of the input files & other variables

driver_results_date <- '2021_12_03' # old: 2021_07_22, new: 2021_10_25, PANEL == 70: 2021_12_03, PANEL == 90: 2021_12_07
fusion_results_date <- '2021_12_03' # old: 2021_08_03, new: 2021_10_25, PANEL == 70: 2021_12_03, PANEL == 90: 2021_12_07
q_val_cutoff <- 0.01 # significance threshold for FDR corrected p-values
plot_clonality <- FALSE # TRUE: include either clonal or subclonal drivers (specify below), FALSE: include clonal & subclonal driver in analysis
clonality_filter <- 'clonal' # can be either 'clonal' or 'subclonal', plot_clonality must be TRUE
filter_resistance <- FALSE # TRUE: filter resistance drivers, FALSE: include resistance drivers
filter_sig <- FALSE # TRUE: only show effect size labels for significantly different cases, FALSE: show all labels
plot_cancer_type_only <- FALSE # TRUE: exclude the violin plot that shows pancancer trends, FLALSE: include pancancer violin plot
by_subtype <- FALSE # TRUE: perform analysis by cancer subtype, FALSE: dont
by_subtype_met_location <- FALSE # TRUE: perform analysis by cancer subtype and metastatic location, FALSE: dont, by_subtype must be TRUE!
by_progression <- FALSE # TRUE: perform analysis by primary progression type, FALSE: dont
by_met_location <- FALSE # TRUE: perform analysis by metastatic location, FALSE: dont

# ---------------------- METADATA

metadata <- read_tsv(paste0(base_dir, '/data/processed/metadata/metadata.tsv'),
                     col_types = 'icccccccccccccccllllcccdiiildclccllllllliiiiiiiiiii') %>%
  # select relevant columns
  select(patient_id, sample_id, cohort, cancer_type, cancer_type_code, 
         cancer_subtype, progression_status_code, is_blacklisted_subtype, metastatic_location) %>%
  # add has_subtype column for easier filtering
  mutate(has_subtype = if_else(cancer_type_code %in% c('BRCA', 'COREAD', 'UCEC'), TRUE, FALSE)) %>%
  mutate(cancer_maintype = cancer_type, .before = cancer_type) %>%
  mutate(cancer_maintype_code = cancer_type_code)

# ---------------------- COMPLETE LINX DRIVER TABLE

if (filter_resistance) {
  linx_drivers_nores <- read_tsv(file = paste0(base_dir, '/results/02_combine_driver_catalogs/', driver_results_date ,'/driver/results/linx_drivers_nores.tsv'))
} else {
  linx_drivers <- read_tsv(file = paste0(base_dir, '/results/02_combine_driver_catalogs/', driver_results_date ,'/driver/results/linx_drivers.tsv'))
}

linx_drivers <- linx_drivers %>%
  rename(clonality = 'clonality2')

# make sure only those samples that are in metadata sheet are included
linx_drivers <- linx_drivers %>%
  semi_join(metadata, by = 'sample_id')

# ---------------------- SUBTYPES

# do analysis by subtype
if (by_subtype) {
  
  # duplicate the rows of samples with subtype annotation, fuse together cancer_type and cancer_subtype column
  # do analysis by subtype and metastatic location
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
    
    # combine artificial pcawg duplicate rows
    subtypes <- subtypes_met_loc %>%
      union(., pcawg_local) %>%
      union(., pcawg_distant) %>%
      union(., pcawg_lymph)
    
    # filter out those cancer subtypes with low hartwig samples in a certain metastatic location category (<5 samples)
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

# do analysis by comparing different primary progression phenotypes to metastatic samples
if (by_progression) {
  
  progression_cancer_types <- metadata %>%
    filter(cancer_type_code %in% c('PRAD', 'PANET')) %>%
    mutate(cancer_type = cancer_type_code)
  
  # create new contrast groups (primary progression stable vs. relapse, relapse vs. met, stable vs. met) for all cancer types with annotation
  progression <- metadata %>%
    filter(progression_status_code %in% c('PROG/RL', 'ST/RM')) %>%
    mutate(cancer_type = paste0(cancer_type_code, '_', progression_status_code))
  
  # duplicate metastatic rows and change the cancer type accordingly to create artificial cancer types contrasts
  
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
    mutate(cancer_type = paste0(cancer_type, '_', metastatic_location)) %>%
    mutate(cancer_type_code = paste0(cancer_type_code, '_', metastatic_location))
  
  # duplicate PCAWG samples rows to create artificial metastatic location contrast groups
  pcawg_local <- metadata %>%
    filter(cohort == 'PCAWG') %>%
    mutate(cancer_type = paste0(cancer_type, '_Local')) %>%
    mutate(cancer_type_code = paste0(cancer_type_code, '_Local'))
  
  pcawg_distant <- metadata %>%
    filter(cohort == 'PCAWG') %>%
    mutate(cancer_type = paste0(cancer_type, '_Distant')) %>%
    mutate(cancer_type_code = paste0(cancer_type_code, '_Distant'))
  
  pcawg_lymph <- metadata %>%
    filter(cohort == 'PCAWG') %>%
    mutate(cancer_type = paste0(cancer_type, '_Lymph')) %>%
    mutate(cancer_type_code = paste0(cancer_type_code, '_Lymph'))
  
  # filter the dataset for cancer types that actually have metastatic location annotation
  metadata <- metadata %>%
    filter(cancer_type %in% met_loc_cancer_types)
  
  # add artificial pcawg duplicate rows to metadata
  metadata <- metadata %>%
    union(., met_loc) %>%
    union(., pcawg_local) %>%
    union(., pcawg_distant) %>%
    union(., pcawg_lymph)
  
  # filter out those cancer_types with low hartwig samples in a certain metastatic location category
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

# calculate sample size by strata
sample_size <- metadata %>%
  group_by(cancer_type, cohort) %>%
  summarise(sample_size = n(), .groups = 'drop') %>%
  pivot_wider(names_from = 'cohort', values_from = 'sample_size') %>%
  mutate(cancer_type_label = paste0(cancer_type, ' (', PCAWG, ' vs. ', Hartwig, ')'))

# calculate cohort size
cohort_size <- metadata %>%
  group_by(cohort) %>%
  summarise(cohort_size = n(), .groups = 'drop') %>%
  pivot_wider(names_from = 'cohort', values_from = 'cohort_size') %>%
  mutate(cancer_type_label = paste0('Pancancer (', PCAWG, ' vs. ', Hartwig, ')')) %>%
  mutate(cancer_type = 'Pancancer', cancer_type_code = 'Pancancer')

# left join metadata and sample size to the linx driver table (some samples have 0 drivers)
linx_drivers <- linx_drivers %>%
  left_join(metadata, by = 'sample_id') %>%
  left_join(sample_size, by = 'cancer_type')

# get the correct cancer types order
label_order <- sample_size %>% 
  { if (by_subtype | by_progression | by_met_location) arrange(., factor(cancer_type)) else arrange(., factor(cancer_type, levels = cancer_type_order)) } %>%
  pull(cancer_type_label)

# get the alternating index number for the grey background rectangles
rect_idx <- seq(2, length(label_order), by = 2)
rect_idx <- label_order[rect_idx]

# impose the correct order on cancer types, and 
linx_drivers <- linx_drivers %>%
  mutate(cancer_type = factor(cancer_type_label, levels = label_order)) %>%
  mutate(strata = if_else(cohort == 'Hartwig', 'Hartwig', 'PCAWG'),
         strata = factor(strata, levels = c('Hartwig', 'PCAWG')))

# calculate total drivers per sample (for total driver plot)
linx_drivers_final <- linx_drivers %>%
  { if (plot_clonality) filter(., clonality == clonality_filter | clonality == 'none') else . } %>%
  group_by(sample_id, cancer_type) %>%
  mutate(total_drivers = sum(driver_per_clonality)) %>%
  ungroup()

if(plot_clonality) {
  # rename the driver_per_clonality column in clonality mode
  linx_drivers_final <- linx_drivers_final %>%
    rename(driver_per_driver_type = 'driver_per_clonality')
} else {
  # sum up clonal and subclonal MUTATION driver for total dps plots (necessary for facetted plot and mean calculation)
  linx_drivers_final <- linx_drivers_final %>%
    group_by(sample_id, cancer_type, cohort, driver) %>%
    mutate(driver_per_driver_type = sum(driver_per_clonality)) %>%
    ungroup()
}

# ---------------------- calculate the mean of mean dps (difference on pancancer level); for text

# total
linx_drivers_mean_of_mean1 <- linx_drivers_final %>%
  group_by(cancer_type, cohort) %>%
  summarise(mean_driver = mean(total_drivers), .groups = 'drop') %>%
  group_by(cohort) %>%
  summarise(mean_driver = mean(mean_driver), .groups = 'drop') %>%
  pivot_wider(., names_from = 'cohort', values_from = 'mean_driver') %>%
  mutate(mean_diff = round(Hartwig - PCAWG, 2)) %>%
  mutate(mean_diff_label = if_else(mean_diff > 0, paste0('+ ', mean_diff), paste0('- ', abs(mean_diff)))) %>%
  pivot_longer(., cols = c('Hartwig', 'PCAWG'), names_to = 'cohort', values_to = 'mean_driver')

# per driver type
linx_drivers_mean_of_mean2 <- linx_drivers_final %>%
  group_by(., cancer_type, cohort, driver) %>%
  summarise(., mean_driver = mean(driver_per_driver_type), .groups = 'drop') %>%
  group_by(cohort, driver) %>%
  summarise(mean_driver = round(mean(mean_driver),1), .groups = 'drop') %>%
  pivot_wider(., names_from = 'cohort', values_from = 'mean_driver') %>%
  mutate(mean_diff = round(Hartwig - PCAWG, 2)) %>%
  mutate(mean_diff_label = if_else(mean_diff > 0, paste0('+ ', mean_diff), paste0('- ', abs(mean_diff)))) %>%
  pivot_longer(., cols = c('Hartwig', 'PCAWG'), names_to = 'cohort', values_to = 'mean_driver')

# calculate driver per sample difference on avg per cancer type between PCAWG / Hartwig and on pancancer level
# pancancer
mean_driver_difference_total_pancan <- linx_drivers_final %>%
  group_by(cohort) %>%
  summarise(mean_driver = mean(total_drivers), .groups = 'drop') %>%
  pivot_wider(., names_from = 'cohort', values_from = 'mean_driver') %>%
  mutate(mean_diff = round(Hartwig - PCAWG, 2)) %>%
  mutate(mean_diff_label = if_else(mean_diff > 0, paste0('+ ', mean_diff), paste0('- ', abs(mean_diff)))) %>%
  pivot_longer(., cols = c('Hartwig', 'PCAWG'), names_to = 'cohort', values_to = 'mean_driver') %>%
  select(-mean_driver)

# per cancer type
mean_driver_difference_total_percan <- linx_drivers_final %>%
  group_by(cancer_type, cohort) %>%
  summarise(mean_driver = mean(total_drivers), .groups = 'drop') %>%
  pivot_wider(., names_from = 'cohort', values_from = 'mean_driver') %>%
  mutate(mean_diff = round(Hartwig - PCAWG, 2)) %>%
  mutate(mean_diff_label = if_else(mean_diff > 0, paste0('+ ', mean_diff), paste0('- ', abs(mean_diff)))) %>%
  pivot_longer(., cols = c('Hartwig', 'PCAWG'), names_to = 'cohort', values_to = 'mean_driver') %>%
  select(-mean_driver)

# ---------------------- custom labels for plots

# custom cancer type labels because the names get way too long for some subgroups..
if (by_subtype) {
  if (by_subtype_met_location) {
    cancer_type_label_df <- linx_drivers_final %>%
      mutate(cancer_type_label = paste0(cancer_maintype_code, if_else(!is.na(cancer_subtype), paste0(' ', cancer_subtype, ': Primary vs. ', metastatic_location, ' (', PCAWG, ' vs. ', Hartwig, ')'), paste0(' (', PCAWG, ' vs. ', Hartwig, ')')))) %>%
      distinct(cancer_type, cancer_type_label)
  } else {
    cancer_type_label_df <- linx_drivers_final %>%
      mutate(cancer_type_label = paste0(str_replace(cancer_type, '_', ': '))) %>%
      distinct(cancer_type, cancer_type_label)
  }
} else if (by_progression) {
  cancer_type_label_df <- linx_drivers_final %>%
    mutate(cancer_type_label = paste0(str_replace(cancer_type, '_([A-Z]+\\/[A-Z]+)(.+$)', ': \\1 vs. Met \\2'))) %>%
    distinct(cancer_type, cancer_type_label)
} else if (by_met_location) {
  cancer_type_label_df <- linx_drivers_final %>%
    mutate(cancer_type_label = paste0(str_replace(cancer_type_code, '_', ': Primary vs. '), ' (', PCAWG, ' vs. ', Hartwig, ')')) %>%
    distinct(cancer_type, cancer_type_label)
} else {
  cancer_type_label_df <- linx_drivers_final
}

cancer_type_label <- as.character(cancer_type_label_df$cancer_type_label)
names(cancer_type_label) <- as.character(cancer_type_label_df$cancer_type)

# new custom facet labels
if (plot_clonality) {
  facet_labels <- c(
    `AMP` = 'AMP',
    `DEL` = 'DEL',
    `FUSION` = 'FUS',
    `MUTATION` = paste0('MUT (', clonality_filter, ')') 
  )
} else {
  facet_labels <- c(
    `AMP` = 'AMP',
    `DEL` = 'DEL',
    `FUSION` = 'FUS',
    `MUTATION` = 'MUT'
  )
}

# ----------------------  plot total drivers per sample

# define plot limits and label positions iu plot
mean_diff_label_pos <- 1.15
asterisk_pos <- 0.95
legend_x_pos <- 1.8
total_y_limit <- 25

# pancancer
dps_total_pancan <- linx_drivers_final %>%
  inner_join(mean_driver_difference_total_pancan, by = c('cohort')) %>%
  mutate(cancer_type = 'Pancancer') %>%
  select(-cancer_type_label, -Hartwig, -PCAWG) %>%
  inner_join(cohort_size, by = c('cancer_type')) %>%
  distinct(sample_id, cancer_type, cancer_type_label, strata, total_drivers) %>%
  ggplot(., aes(x = strata, y = total_drivers)) +
  geom_rect(
    data = . %>%
      distinct(cancer_type),
    inherit.aes = FALSE, # fixed the inheritance error
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#a0a0a0') +
  geom_jitter(color = 'gray40', size = 0.1, alpha = 0.8) +
  geom_violin(aes(color = strata), fill = NA, position = position_dodge(1), size = 1.01, 
              scale = 'width', width = 0.5, adjust = 2.5) +
  geom_text(data = . %>% 
              compare_means(
                data = ., formula = total_drivers ~ strata,
                method = 'wilcox.test',
                symnum.args = list(cutpoints = c(0, q_val_cutoff, 1), symbols = c('*', 'ns'))) %>% 
              inner_join(mean_driver_difference_total_pancan, by = c('group1' = 'cohort')) %>%
              rename(strata = 'group1') %>%
              { if (filter_sig) filter(., p.signif != 'ns') else . }, 
            aes(label = p.signif, y = total_y_limit * asterisk_pos), size = 3, nudge_x = 0.5) +
  geom_text(data = . %>% 
              compare_means(
                data = ., formula = total_drivers ~ strata,
                method = 'wilcox.test',
                symnum.args = list(cutpoints = c(0, q_val_cutoff, 1), symbols = c('*', 'ns'))) %>% 
              inner_join(mean_driver_difference_total_pancan, by = c('group1' = 'cohort')) %>%
              rename(strata = 'group1') %>%
              { if (filter_sig) filter(., p.signif != 'ns') else . },
            aes(label = mean_diff_label, y = total_y_limit * mean_diff_label_pos), size = 3, nudge_x = 0.5) +
  stat_summary(aes(group = strata), fun = mean, geom = 'point', shape = 16, size = 2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  facet_wrap(cancer_type_label ~ ., ncol = 1, strip.position = 'left') +
  scale_color_manual(
    values = rev(cohort_colors),
    guide = guide_legend(title = '', reverse = TRUE, override.aes = list(fill = cohort_colors))
  ) +
  coord_flip(
    ylim = c(0, total_y_limit),
    clip = 'off'
  ) +
  labs(
    title = 'All drivers',
    x = '',
    y = 'Drivers per sample') +
  theme_classic() +
  theme_smallfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = 'none',
    panel.spacing = unit(0, 'cm')
    )

# per cancer type
dps_total_percan <- linx_drivers_final %>%
  inner_join(mean_driver_difference_total_percan, by = c('cancer_type', 'cohort')) %>%
  distinct(sample_id, cancer_type, cancer_type_label, strata, total_drivers) %>%
  ggplot(., aes(x = strata, y = total_drivers)) +
  geom_rect(data = linx_drivers_final %>% 
              filter(cohort == 'Hartwig') %>%
              filter(cancer_type %in% rect_idx) %>%
              distinct(cancer_type),
            inherit.aes = FALSE, # fixed the inheritance error
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#a0a0a0') +
  geom_jitter(color = 'gray40', size = 0.1, alpha = 0.15) +
  geom_violin(aes(color = strata), fill = NA, position = position_dodge(1), size = 1.01, 
              scale = 'width', width = 0.5, adjust = 2.5) +
  geom_text(data = . %>% 
              compare_means(
                data = ., formula = total_drivers ~ strata,
                method = 'wilcox.test', group.by = 'cancer_type', p.adjust.method = 'fdr',
                symnum.args = list(cutpoints = c(0, q_val_cutoff, 1), symbols = c('*', 'ns'))) %>% 
              inner_join(mean_driver_difference_total_percan, by = 'cancer_type') %>%
              rename(strata = 'group1') %>%
              { if (filter_sig) filter(., p.signif != 'ns') else . }, 
            aes(label = p.signif, y = total_y_limit * asterisk_pos), size = 3, nudge_x = 0.5) +
  geom_text(data = . %>% 
              compare_means(
                data = ., formula = total_drivers ~ strata,
                method = 'wilcox.test', group.by = 'cancer_type', p.adjust.method = 'fdr',
                symnum.args = list(cutpoints = c(0, q_val_cutoff, 1), symbols = c('*', 'ns'))) %>% 
              inner_join(mean_driver_difference_total_percan, by = 'cancer_type') %>%
              rename(strata = 'group1') %>%
              { if (filter_sig) filter(., p.signif != 'ns') else . },
            aes(label = mean_diff_label, y = total_y_limit * mean_diff_label_pos), size = 3, nudge_x = 0.5) +
  stat_summary(aes(group = strata), fun = mean, geom = 'point', shape = 16, size = 2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  facet_wrap(cancer_type ~ ., ncol = 1, strip.position = 'left', labeller =  as_labeller(cancer_type_label)) +
  scale_color_manual(
    values = rev(cohort_colors),
    guide = guide_legend(title = '', reverse = TRUE, override.aes = list(fill = cohort_colors))
    ) +
  coord_flip(
    ylim = c(0, total_y_limit),
    clip = 'off'
    ) +
  labs(
    title = 'All drivers',
    x = '',
    y = 'Drivers per sample') +
  theme_classic() +
  theme_smallfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.title = if (!by_subtype & !by_met_location & !by_progression & !plot_cancer_type_only) { element_blank() } else { },
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = c(legend_x_pos, 0.5),
    panel.spacing = unit(0, 'cm'),
    plot.margin = margin(1,5,1,0, unit = 'cm')
    )

# ---------------------- calculate the mean driver difference per driver type

# pancancer
mean_driver_difference_per_driver_pancan <- linx_drivers_final %>%
  group_by(cohort, driver) %>%
  summarise(mean_driver = mean(driver_per_driver_type), .groups = 'drop') %>%
  pivot_wider(names_from = 'cohort', values_from = 'mean_driver') %>%
  mutate(mean_diff = round(Hartwig - PCAWG, 2)) %>%
  mutate(mean_diff_label = if_else(mean_diff > 0, paste0('+ ', mean_diff), paste0('- ', abs(mean_diff)))) %>%
  pivot_longer(., cols = c('Hartwig', 'PCAWG'), names_to = 'cohort', values_to = 'mean_driver') %>%
  select(-mean_driver)

# per cancer type
mean_driver_difference_per_driver <- linx_drivers_final %>%
  group_by(cancer_type, cohort, driver) %>%
  summarise(mean_driver = mean(driver_per_driver_type), .groups = 'drop') %>%
  pivot_wider(names_from = 'cohort', values_from = 'mean_driver') %>%
  mutate(mean_diff = round(Hartwig - PCAWG, 2)) %>%
  mutate(mean_diff_label = if_else(mean_diff > 0, paste0('+ ', mean_diff), paste0('- ', abs(mean_diff)))) %>%
  pivot_longer(., cols = c('Hartwig', 'PCAWG'), names_to = 'cohort', values_to = 'mean_driver') %>%
  select(-mean_driver)

# join mean driver difference to the complete table
# pancancer
linx_drivers_final_pancan <- linx_drivers_final %>%
  inner_join(mean_driver_difference_per_driver_pancan, by = c('cohort', 'driver')) %>%
  mutate(cancer_type = 'Pancancer') %>%
  # get the cohort size in there
  select(-cancer_type_label, -Hartwig, -PCAWG) %>%
  inner_join(cohort_size, by = c('cancer_type'))

# per cancer type
linx_drivers_final_percan <- linx_drivers_final %>%
  inner_join(mean_driver_difference_per_driver, by = c('cancer_type', 'cohort', 'driver'))

# define ylim for each subplot (to align labels in plot at the same position relative to ylim)
amp_y_limit <- 16
del_y_limit <- 20
fus_y_limit <- 4
mut_y_limit <- 20

# ----------------------  plot amplification dps violin plot

# pancancer
amp_violin_pancan <- linx_drivers_final_pancan %>%
  filter(driver == 'AMP') %>%
  distinct(sample_id, strata, cancer_type, cancer_type_label, driver, driver_per_driver_type) %>%
  ggplot(., aes(x = strata, y = driver_per_driver_type)) +
  geom_rect(
    data = . %>%
      distinct(cancer_type),
    inherit.aes = FALSE, # fixed the inheritance error
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#a0a0a0') +
  geom_jitter(color = 'gray40', size = 0.1, alpha = 0.25) +
  geom_violin(aes(color = strata), fill = NA, position = position_dodge(1), size = 1.01, 
              scale = 'width', width = 0.5, adjust = 2.5) +
  geom_text(data = . %>% 
              compare_means(
                data = ., formula = driver_per_driver_type ~ strata,
                method = 'wilcox.test',
                symnum.args = list(cutpoints = c(0, q_val_cutoff, 1), symbols = c('*', 'ns'))) %>% 
              inner_join(mean_driver_difference_per_driver_pancan, by = c('group1' = 'cohort')) %>%
              filter(driver == 'AMP') %>%
              rename(strata = 'group1') %>%
              { if (filter_sig) filter(., p.signif != 'ns') else . }, 
            aes(label = p.signif, y = amp_y_limit * asterisk_pos), size = 3, nudge_x = 0.5) +
  geom_text(data = . %>% 
              compare_means(
                data = ., formula = driver_per_driver_type ~ strata,
                method = 'wilcox.test',
                symnum.args = list(cutpoints = c(0, q_val_cutoff, 1), symbols = c('*', 'ns'))) %>% 
              inner_join(mean_driver_difference_per_driver_pancan, by = c('group1' = 'cohort')) %>%
              filter(driver == 'AMP') %>%
              rename(strata = 'group1') %>%
              { if (filter_sig) filter(., p.signif != 'ns') else . },
            aes(label = mean_diff_label, y = amp_y_limit * mean_diff_label_pos), size = 3, nudge_x = 0.5) +
  stat_summary(aes(group = strata), fun = mean, geom = 'point', shape = 16, size = 2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  facet_wrap(cancer_type_label ~ ., ncol = 1, strip.position = 'left') +
  scale_color_manual(
    values = rev(cohort_colors),
    guide = guide_legend(title = '', reverse = TRUE, override.aes = list(fill = cohort_colors))
  ) +
  coord_flip(
    ylim = c(0, amp_y_limit),
    clip = 'off'
  ) +
  labs(
    title = 'Amplifications',
    x = '',
    y = 'Drivers per sample') +
  theme_classic() +
  theme_smallfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = 'none',
    panel.spacing = unit(0, 'cm')
    )

# per cancer type
amp_violin <- linx_drivers_final_percan %>%
  filter(driver == 'AMP') %>%
  distinct(sample_id, strata, cancer_type, cancer_type_label, driver, driver_per_driver_type) %>%
  ggplot(., aes(x = strata, y = driver_per_driver_type)) +
  geom_rect(data = linx_drivers_final %>% 
              filter(cohort == 'Hartwig') %>%
              filter(cancer_type %in% rect_idx) %>%
              distinct(cancer_type),
            inherit.aes = FALSE, # fixed the inheritance error
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#a0a0a0') +
  geom_jitter(color = 'gray40', size = 0.1, alpha = 0.25) +
  geom_violin(aes(color = strata), fill = NA, position = position_dodge(1), size = 1.01, 
              scale = 'width', width = 0.5, adjust = 2.5) +
  geom_text(data = . %>% 
              compare_means(
                data = ., formula = driver_per_driver_type ~ strata,
                method = 'wilcox.test', group.by = 'cancer_type', p.adjust.method = 'fdr',
                symnum.args = list(cutpoints = c(0, q_val_cutoff, 1), symbols = c('*', 'ns'))) %>% 
              inner_join(mean_driver_difference_per_driver, by = 'cancer_type') %>%
              filter(driver == 'AMP', cohort == 'Hartwig') %>%
              rename(strata = 'group1') %>%
              { if (filter_sig) filter(., p.signif != 'ns') else . }, 
            aes(label = p.signif, y = amp_y_limit * asterisk_pos), size = 3, nudge_x = 0.5) +
  geom_text(data = . %>%
              compare_means(
                data = ., formula = driver_per_driver_type ~ strata,
                method = 'wilcox.test', group.by = 'cancer_type', p.adjust.method = 'fdr',
                symnum.args = list(cutpoints = c(0, q_val_cutoff, 1), symbols = c('*', 'ns'))) %>%
              inner_join(mean_driver_difference_per_driver, by = c('cancer_type')) %>%
              filter(driver == 'AMP', cohort == 'Hartwig') %>%
              rename(strata = 'group1') %>%
              { if (filter_sig) filter(., p.signif != 'ns') else . },
            aes(label = mean_diff_label, y = amp_y_limit * mean_diff_label_pos), size = 3, nudge_x = 0.5) +
  stat_summary(aes(group = strata), fun = mean, geom = 'point', shape = 16, size = 2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  facet_wrap(cancer_type ~ ., ncol = 1, strip.position = 'left', labeller =  as_labeller(cancer_type_label)) +
  scale_color_manual(
    values = rev(cohort_colors),
    guide = guide_legend(title = '', reverse = TRUE, override.aes = list(fill = cohort_colors))
  ) +
  coord_flip(
    ylim = c(0, amp_y_limit),
    clip = 'off'
  ) +
  labs(
    title = 'Amplifications',
    x = '',
    y = 'Drivers per sample') +
  theme_classic() +
  theme_smallfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.title = if (!by_subtype & !by_met_location & !by_progression & !plot_cancer_type_only) { element_blank() } else { },
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = c(legend_x_pos, 0.5),
    panel.spacing = unit(0, 'cm'),
    plot.margin = margin(1,5,1,0, unit = 'cm')
    )

# ----------------------  plot deletion dps violin plot

# pancancer
del_violin_pancan <- linx_drivers_final_pancan %>%
  filter(driver == 'DEL') %>%
  distinct(sample_id, strata, cancer_type, cancer_type_label, driver, driver_per_driver_type) %>%
  ggplot(., aes(x = strata, y = driver_per_driver_type)) +
  geom_rect(
    data = . %>%
      distinct(cancer_type),
    inherit.aes = FALSE, # fixed the inheritance error
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#a0a0a0') +
  geom_jitter(color = 'gray40', size = 0.1, alpha = 0.25) +
  geom_violin(aes(color = strata), fill = NA, position = position_dodge(1), size = 1.01, 
              scale = 'width', width = 0.5, adjust = 2.5) +
  geom_text(data = . %>% 
              compare_means(
                data = ., formula = driver_per_driver_type ~ strata,
                method = 'wilcox.test',
                symnum.args = list(cutpoints = c(0, q_val_cutoff, 1), symbols = c('*', 'ns'))) %>% 
              inner_join(mean_driver_difference_per_driver_pancan, by = c('group1' = 'cohort')) %>%
              filter(driver == 'DEL') %>%
              rename(strata = 'group1') %>%
              { if (filter_sig) filter(., p.signif != 'ns') else . }, 
            aes(label = p.signif, y = del_y_limit * asterisk_pos), size = 3, nudge_x = 0.5) +
  geom_text(data = . %>% 
              compare_means(
                data = ., formula = driver_per_driver_type ~ strata,
                method = 'wilcox.test',
                symnum.args = list(cutpoints = c(0, q_val_cutoff, 1), symbols = c('*', 'ns'))) %>% 
              inner_join(mean_driver_difference_per_driver_pancan, by = c('group1' = 'cohort')) %>%
              filter(driver == 'DEL') %>%
              rename(strata = 'group1') %>%
              { if (filter_sig) filter(., p.signif != 'ns') else . },
            aes(label = mean_diff_label, y = del_y_limit * mean_diff_label_pos), size = 3, nudge_x = 0.5) +
  stat_summary(aes(group = strata), fun = mean, geom = 'point', shape = 16, size = 2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  facet_wrap(cancer_type_label ~ ., ncol = 1, strip.position = 'left') +
  scale_color_manual(
    values = rev(cohort_colors),
    guide = guide_legend(title = '', reverse = TRUE, override.aes = list(fill = cohort_colors))
  ) +
  coord_flip(
    ylim = c(0, del_y_limit),
    clip = 'off'
  ) +
  labs(
    title = 'Deletions',
    x = '',
    y = 'Drivers per sample') +
  theme_classic() +
  theme_smallfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = 'none',
    panel.spacing = unit(0, 'cm')
  )

# per cancer type
del_violin <- linx_drivers_final_percan %>%
  filter(driver == 'DEL') %>%
  distinct(sample_id, strata, cancer_type, cancer_type_label, driver, driver_per_driver_type) %>%
  ggplot(., aes(x = strata, y = driver_per_driver_type)) +
  geom_rect(data = linx_drivers_final %>% 
              filter(cohort == 'Hartwig') %>%
              filter(cancer_type %in% rect_idx) %>%
              distinct(cancer_type),
            inherit.aes = FALSE, # fixed the inheritance error
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#a0a0a0') +
  geom_jitter(color = 'gray40', size = 0.1, alpha = 0.25) +
  geom_violin(aes(color = strata), fill = NA, position = position_dodge(1), size = 1.01, 
              scale = 'width', width = 0.5, adjust = 2.5) +
  geom_text(data = . %>% 
              compare_means(
                data = ., formula = driver_per_driver_type ~ strata,
                method = 'wilcox.test', group.by = 'cancer_type', p.adjust.method = 'fdr',
                symnum.args = list(cutpoints = c(0, q_val_cutoff, 1), symbols = c('*', 'ns'))) %>% 
              inner_join(mean_driver_difference_per_driver, by = 'cancer_type') %>%
              filter(driver == 'DEL', cohort == 'Hartwig') %>%
              rename(strata = 'group1') %>%
              { if (filter_sig) filter(., p.signif != 'ns') else . }, 
            aes(label = p.signif, y = del_y_limit * asterisk_pos), size = 3, nudge_x = 0.5) +
  geom_text(data = . %>% 
              compare_means(
                data = ., formula = driver_per_driver_type ~ strata,
                method = 'wilcox.test', group.by = 'cancer_type', p.adjust.method = 'fdr',
                symnum.args = list(cutpoints = c(0, q_val_cutoff, 1), symbols = c('*', 'ns'))) %>% 
              inner_join(mean_driver_difference_per_driver, by = 'cancer_type') %>%
              filter(driver == 'DEL', cohort == 'Hartwig') %>%
              rename(strata = 'group1') %>%
              { if (filter_sig) filter(., p.signif != 'ns') else . },
            aes(label = mean_diff_label, y = del_y_limit * mean_diff_label_pos), size = 3, nudge_x = 0.5) +
  stat_summary(aes(group = strata), fun = mean, geom = 'point', shape = 16, size = 2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  facet_wrap(cancer_type ~ ., ncol = 1, strip.position = 'left', labeller =  as_labeller(cancer_type_label)) +
  scale_color_manual(
    values = rev(cohort_colors),
    guide = guide_legend(title = '', reverse = TRUE, override.aes = list(fill = cohort_colors))
  ) +
  coord_flip(
    ylim = c(0, del_y_limit),
    clip = 'off'
  ) +
  labs(
    title = 'Deletions',
    x = '',
    y = 'Drivers per sample') +
  theme_classic() +
  theme_smallfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.title = if (!by_subtype & !by_met_location & !by_progression & !plot_cancer_type_only) { element_blank() } else { },
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = c(legend_x_pos, 0.5),
    panel.spacing = unit(0, 'cm'),
    plot.margin = margin(1,5,1,0, unit = 'cm')
    )

# ----------------------  plot fusion dps violin plot

# pancancer
fus_violin_pancan <- linx_drivers_final_pancan %>%
  filter(driver == 'FUSION') %>%
  distinct(sample_id, strata, cancer_type, cancer_type_label, driver, driver_per_driver_type) %>%
  ggplot(., aes(x = strata, y = driver_per_driver_type)) +
  geom_rect(
    data = . %>%
      distinct(cancer_type),
    inherit.aes = FALSE, # fixed the inheritance error
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#a0a0a0') +
  geom_jitter(color = 'gray40', size = 0.1, alpha = 0.25) +
  geom_violin(aes(color = strata), fill = NA, position = position_dodge(1), size = 1.01, 
              scale = 'width', width = 0.5, adjust = 2.5) +
  geom_text(data = . %>% 
              compare_means(
                data = ., formula = driver_per_driver_type ~ strata,
                method = 'wilcox.test',
                symnum.args = list(cutpoints = c(0, q_val_cutoff, 1), symbols = c('*', 'ns'))) %>% 
              inner_join(mean_driver_difference_per_driver_pancan, by = c('group1' = 'cohort')) %>%
              filter(driver == 'FUSION') %>%
              rename(strata = 'group1') %>%
              { if (filter_sig) mutate(., p.signif = if_else(p.signif == 'ns', NA_character_, p.signif)) else . }, 
            aes(label = p.signif, y = fus_y_limit * asterisk_pos), size = 3, nudge_x = 0.5) +
  geom_text(data = . %>% 
              compare_means(
                data = ., formula = driver_per_driver_type ~ strata,
                method = 'wilcox.test',
                symnum.args = list(cutpoints = c(0, q_val_cutoff, 1), symbols = c('*', 'ns'))) %>% 
              inner_join(mean_driver_difference_per_driver_pancan, by = c('group1' = 'cohort')) %>%
              filter(driver == 'FUSION') %>%
              rename(strata = 'group1') %>%
              { if (filter_sig) mutate(., mean_diff_label = if_else(p.signif == 'ns', NA_character_, p.signif)) else . },
            aes(label = mean_diff_label, y = fus_y_limit * mean_diff_label_pos), size = 3, nudge_x = 0.5) +
  stat_summary(aes(group = strata), fun = mean, geom = 'point', shape = 16, size = 2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  facet_wrap(cancer_type_label ~ ., ncol = 1, strip.position = 'left') +
  scale_color_manual(
    values = rev(cohort_colors),
    guide = guide_legend(title = '', reverse = TRUE, override.aes = list(fill = cohort_colors))
  ) +
  coord_flip(
    ylim = c(0, fus_y_limit),
    clip = 'off'
  ) +
  scale_y_continuous(limits = c(0, fus_y_limit * legend_x_pos), oob = scales::censor) + # for some reason ggplot plots points below 0 for fusions, so need to add this after coord_flip()
  labs(
    title = 'Fusions',
    x = '',
    y = 'Drivers per sample') +
  theme_classic() +
  theme_smallfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = 'none',
    panel.spacing = unit(0, 'cm')
  )

# per cancer type
fus_violin <- linx_drivers_final_percan %>%
  filter(driver == 'FUSION') %>%
  distinct(sample_id, strata, cancer_type, cancer_type_label, driver, driver_per_driver_type) %>%
  ggplot(., aes(x = strata, y = driver_per_driver_type)) +
  geom_rect(data = linx_drivers_final %>% 
              filter(cohort == 'Hartwig') %>%
              filter(cancer_type %in% rect_idx) %>%
              distinct(cancer_type),
            inherit.aes = FALSE, # fixed the inheritance error
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#a0a0a0') +
  geom_jitter(color = 'gray40', size = 0.1, alpha = 0.25) +
  geom_violin(aes(color = strata), fill = NA, position = position_dodge(1), size = 1.01, 
              scale = 'width', width = 0.5, adjust = 2.5) +
  geom_text(data = . %>% 
              compare_means(
                data = ., formula = driver_per_driver_type ~ strata,
                method = 'wilcox.test', group.by = 'cancer_type', p.adjust.method = 'fdr',
                symnum.args = list(cutpoints = c(0, q_val_cutoff, 1), symbols = c('*', 'ns'))) %>% 
              inner_join(mean_driver_difference_per_driver, by = 'cancer_type') %>%
              filter(driver == 'FUSION', cohort == 'Hartwig') %>%
              rename(strata = 'group1') %>%
              { if (filter_sig) mutate(., p.signif = if_else(p.signif == 'ns', NA_character_, p.signif)) else . }, 
            aes(label = p.signif, y = fus_y_limit * asterisk_pos), size = 3, nudge_x = 0.5) +
  geom_text(data = . %>%
              compare_means(
                data = ., formula = driver_per_driver_type ~ strata,
                method = 'wilcox.test', group.by = 'cancer_type', p.adjust.method = 'fdr',
                symnum.args = list(cutpoints = c(0, q_val_cutoff, 1), symbols = c('*', 'ns'))) %>%
              inner_join(mean_driver_difference_per_driver, by = 'cancer_type') %>%
              filter(driver == 'FUSION', cohort == 'Hartwig') %>%
              rename(strata = 'group1') %>%
              { if (filter_sig) mutate(., mean_diff_label = if_else(p.signif == 'ns', NA_character_, p.signif)) else . },
            aes(label = mean_diff_label, y = fus_y_limit * mean_diff_label_pos), size = 3, nudge_x = 0.5) +
  stat_summary(aes(group = strata), fun = mean, geom = 'point', shape = 16, size = 2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  facet_wrap(cancer_type ~ ., ncol = 1, strip.position = 'left', labeller =  as_labeller(cancer_type_label)) +
  scale_color_manual(
    values = rev(cohort_colors),
    guide = guide_legend(title = '', reverse = TRUE, override.aes = list(fill = cohort_colors))
  ) +
  coord_flip(
    ylim = c(0, fus_y_limit),
    clip = 'off'
  ) +
  scale_y_continuous(limits = c(0, fus_y_limit * legend_x_pos), oob = scales::censor) + # for some reason ggplot plots points below 0 for fusions, so need to add this after coord_flip()
  labs(
    title = 'Fusions',
    x = '',
    y = 'Drivers per sample') +
  theme_classic() +
  theme_smallfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.title = if (!by_subtype & !by_met_location & !by_progression & !plot_cancer_type_only) { element_blank() } else { },
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = c(legend_x_pos, 0.5),
    panel.spacing = unit(0, 'cm'),
    plot.margin = margin(1,5,1,0, unit = 'cm')
    )

# ----------------------  plot (clonal/subclonal) point mutation dps violin plot

# pancancer
mut_violin_pancan <- linx_drivers_final_pancan %>%
  filter(driver == 'MUTATION') %>%
  distinct(sample_id, strata, cancer_type, cancer_type_label, driver, driver_per_driver_type) %>%
  ggplot(., aes(x = strata, y = driver_per_driver_type)) +
  geom_rect(
    data = . %>%
      distinct(cancer_type),
    inherit.aes = FALSE, # fixed the inheritance error
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#a0a0a0') +
  geom_jitter(color = 'gray40', size = 0.1, alpha = 0.25) +
  geom_violin(aes(color = strata), fill = NA, position = position_dodge(1), size = 1.01, 
              scale = 'width', width = 0.5, adjust = 2.5) +
  geom_text(data = . %>% 
              compare_means(
                data = ., formula = driver_per_driver_type ~ strata,
                method = 'wilcox.test',
                symnum.args = list(cutpoints = c(0, q_val_cutoff, 1), symbols = c('*', 'ns'))) %>% 
              inner_join(mean_driver_difference_per_driver_pancan, by = c('group1' = 'cohort')) %>%
              filter(driver == 'MUTATION') %>%
              rename(strata = 'group1') %>%
              { if (filter_sig) filter(., p.signif != 'ns') else . }, 
            aes(label = p.signif, y = mut_y_limit * asterisk_pos), size = 3, nudge_x = 0.5) +
  geom_text(data = . %>% 
              compare_means(
                data = ., formula = driver_per_driver_type ~ strata,
                method = 'wilcox.test',
                symnum.args = list(cutpoints = c(0, q_val_cutoff, 1), symbols = c('*', 'ns'))) %>% 
              inner_join(mean_driver_difference_per_driver_pancan, by = c('group1' = 'cohort')) %>%
              filter(driver == 'MUTATION') %>%
              rename(strata = 'group1') %>%
              { if (filter_sig) filter(., p.signif != 'ns') else . },
            aes(label = mean_diff_label, y = mut_y_limit * mean_diff_label_pos), size = 3, nudge_x = 0.5) +
  stat_summary(aes(group = strata), fun = mean, geom = 'point', shape = 16, size = 2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  facet_wrap(cancer_type_label ~ ., ncol = 1, strip.position = 'left') +
  scale_color_manual(
    values = rev(cohort_colors),
    guide = guide_legend(title = '', reverse = TRUE, override.aes = list(fill = cohort_colors))
  ) +
  coord_flip(
    ylim = c(0, mut_y_limit),
    clip = 'off'
  ) +
  labs(
    title = 'Mutations',
    x = '',
    y = 'Drivers per sample') +
  theme_classic() +
  theme_smallfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = 'none',
    panel.spacing = unit(0, 'cm')
  )

# per cancer type
mut_violin <- linx_drivers_final_percan %>%
  filter(driver == 'MUTATION') %>%
  distinct(sample_id, strata, cancer_type, cancer_type_label, driver, driver_per_driver_type) %>%
  ggplot(., aes(x = strata, y = driver_per_driver_type)) +
  geom_rect(data = linx_drivers_final %>% 
              filter(cohort == 'Hartwig') %>%
              filter(cancer_type %in% rect_idx) %>%
              distinct(cancer_type),
            inherit.aes = FALSE, # fixed the inheritance error
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#a0a0a0') +
  geom_jitter(color = 'gray40', size = 0.1, alpha = 0.25) +
  geom_violin(aes(color = strata), fill = NA, position = position_dodge(1), size = 1.01, 
              scale = 'width', width = 0.5, adjust = 2.5) +
  geom_text(data = . %>% 
              compare_means(
                data = ., formula = driver_per_driver_type ~ strata,
                method = 'wilcox.test', group.by = 'cancer_type', p.adjust.method = 'fdr',
                symnum.args = list(cutpoints = c(0, q_val_cutoff, 1), symbols = c('*', 'ns'))) %>% 
              inner_join(mean_driver_difference_per_driver, by = 'cancer_type') %>%
              filter(driver == 'MUTATION', cohort == 'Hartwig') %>%
              rename(strata = 'group1') %>%
              { if (filter_sig) filter(., p.signif != 'ns') else . }, 
            aes(label = p.signif, y = mut_y_limit * asterisk_pos), size = 3, nudge_x = 0.5) +
  geom_text(data = . %>% 
              compare_means(
                data = ., formula = driver_per_driver_type ~ strata,
                method = 'wilcox.test', group.by = 'cancer_type', p.adjust.method = 'fdr',
                symnum.args = list(cutpoints = c(0, q_val_cutoff, 1), symbols = c('*', 'ns'))) %>% 
              inner_join(mean_driver_difference_per_driver, by = 'cancer_type') %>%
              filter(driver == 'MUTATION', cohort == 'Hartwig') %>%
              rename(strata = 'group1') %>%
              { if (filter_sig) filter(., p.signif != 'ns') else . },
            aes(label = mean_diff_label, y = mut_y_limit * mean_diff_label_pos), size = 3, nudge_x = 0.5) +
  stat_summary(aes(group = strata), fun = mean, geom = 'point', shape = 16, size = 2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  facet_wrap(cancer_type ~ ., ncol = 1, strip.position = 'left', labeller =  as_labeller(cancer_type_label)) +
  scale_color_manual(
    values = rev(cohort_colors),
    guide = guide_legend(title = '', reverse = TRUE, override.aes = list(fill = cohort_colors))
  ) +
  coord_flip(
    ylim = c(0, mut_y_limit),
    clip = 'off'
  ) +
  labs(
    title = 'Mutations',
    x = '',
    y = 'Drivers per sample') +
  theme_classic() +
  theme_smallfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.title = if (!by_subtype & !by_met_location & !by_progression & !plot_cancer_type_only) { element_blank() } else { },
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = c(legend_x_pos, 0.5),
    panel.spacing = unit(0, 'cm'),
    plot.margin = margin(1,5,1,0, unit = 'cm')
    )

# assemble individual plots
total_assembled_plot <- dps_total_pancan + dps_total_percan + plot_layout(nrow = 2, ncol = 1, heights = c(1,length(label_order)))
amp_assembled_plot <- amp_violin_pancan + amp_violin + plot_layout(nrow = 2, ncol = 1, heights = c(1,length(label_order)))
del_assembled_plot <- del_violin_pancan + del_violin + plot_layout(nrow = 2, ncol = 1, heights = c(1,length(label_order)))
fus_assembled_plot <- fus_violin_pancan + fus_violin + plot_layout(nrow = 2, ncol = 1, heights = c(1,length(label_order)))
mut_assembled_plot <- mut_violin_pancan + mut_violin + plot_layout(nrow = 2, ncol = 1, heights = c(1,length(label_order)))

# save plots
# as PDF
# total
pdf(
  file = paste0(output_dir, '/', current_date, '_total_dps_q', q_val_cutoff,
                if (by_subtype) { if (by_subtype_met_location) { '_by_subtype_and_met_location.pdf' } else { '_by_subtype.pdf' } } 
                else if (by_met_location) { '_by_met_location.pdf' } 
                else if (by_progression) { '_by_progression.pdf' } 
                else { '.pdf' }),
  width = 7, 
  height = if (by_subtype & !by_subtype_met_location) { 5 } 
  else if (by_subtype_met_location) { 11 } 
  else if (by_met_location) { 16 }
  else if (by_progression) { 3 }
  else { 7 },
  useDingbats = FALSE,
  compress = FALSE
)
if (by_subtype | by_met_location | by_progression | plot_cancer_type_only) { print(dps_total_percan) } else { print(total_assembled_plot) }
dev.off()

# amplification
pdf(
  file = paste0(output_dir, '/', current_date, '_amp_dps_q', q_val_cutoff,
                if (by_subtype) { if (by_subtype_met_location) { '_by_subtype_and_met_location.pdf' } else { '_by_subtype.pdf' } } 
                else if (by_met_location) { '_by_met_location.pdf' } 
                else if (by_progression) { '_by_progression.pdf' } 
                else { '.pdf' }),
  width = 7, 
  height = if (by_subtype & !by_subtype_met_location) { 5 } 
  else if (by_subtype_met_location) { 11 } 
  else if (by_met_location) { 16 }
  else if (by_progression) { 3 }
  else { 7 },
  useDingbats = FALSE,
  compress = FALSE
)
if (by_subtype | by_met_location | by_progression | plot_cancer_type_only) { print(amp_violin) } else { print(amp_assembled_plot) }
dev.off()

# deletion
pdf(
  file = paste0(output_dir, '/', current_date, '_del_dps_q', q_val_cutoff,
                if (by_subtype) { if (by_subtype_met_location) { '_by_subtype_and_met_location.pdf' } else { '_by_subtype.pdf' } } 
                else if (by_met_location) { '_by_met_location.pdf' } 
                else if (by_progression) { '_by_progression.pdf' } 
                else { '.pdf' }),
  width = 7, 
  height = if (by_subtype & !by_subtype_met_location) { 5 } 
  else if (by_subtype_met_location) { 11 } 
  else if (by_met_location) { 16 }
  else if (by_progression) { 3 }
  else { 7 },
  useDingbats = FALSE,
  compress = FALSE
)
if (by_subtype | by_met_location | by_progression | plot_cancer_type_only) { print(del_violin) } else { print(del_assembled_plot) }
dev.off()

# fusion
pdf(
  file = paste0(output_dir, '/', current_date, '_fus_dps_q', q_val_cutoff,
                if (by_subtype) { if (by_subtype_met_location) { '_by_subtype_and_met_location.pdf' } else { '_by_subtype.pdf' } } 
                else if (by_met_location) { '_by_met_location.pdf' } 
                else if (by_progression) { '_by_progression.pdf' } 
                else { '.pdf' }),
  width = 7, 
  height = if (by_subtype & !by_subtype_met_location) { 5 } 
  else if (by_subtype_met_location) { 11 } 
  else if (by_met_location) { 16 }
  else if (by_progression) { 3 }
  else { 7 },
  useDingbats = FALSE,
  compress = FALSE
)
if (by_subtype | by_met_location | by_progression | plot_cancer_type_only) { print(fus_violin) } else { print(fus_assembled_plot) }
dev.off()

# point mutation
pdf(
  file = paste0(output_dir, '/', current_date, '_mut_dps_q', q_val_cutoff,
                if (by_subtype) { if (by_subtype_met_location) { '_by_subtype_and_met_location.pdf' } else { '_by_subtype.pdf' } } 
                else if (by_met_location) { '_by_met_location.pdf' } 
                else if (by_progression) { '_by_progression.pdf' } 
                else { '.pdf' }),
  width = 7, 
  height = if (by_subtype & !by_subtype_met_location) { 5 } 
  else if (by_subtype_met_location) { 11 } 
  else if (by_met_location) { 16 }
  else if (by_progression) { 3 }
  else { 7 },
  useDingbats = FALSE,
  compress = FALSE
)
if (by_subtype | by_met_location | by_progression | plot_cancer_type_only) { print(mut_violin) } else { print(mut_assembled_plot) }
dev.off()