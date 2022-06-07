############### Plot drivers per sample plots ###############
# author: remy (sascha)
# date: 16/11/2021
# last updated: 08/03/2022

### Description
# This script creates all the drivers per samples plots

# libraries
library(tidyverse) # data manipulation and plotting
library(ggpubr) # wilcoxon test for plots
library(naniar) # replace with NA function
library(patchwork) # stitch plots together
library(svglite) # handle the issue with SVGs in Inkscape

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
resize_factor <- 1

theme_bigfont <- theme(plot.title = element_text(size=22),
                       plot.subtitle = element_text(size = 16),
                       axis.text.x = element_text(size=15),
                       axis.text.y = element_text(size=15), 
                       axis.title = element_text(size=18),
                       legend.text = element_text(size = 14),
                       strip.text.x = element_text(size = 15))

theme_smallfont <- theme(plot.title = element_text(size=17 * resize_factor),
                         plot.subtitle = element_text(size = 11 * resize_factor),
                         axis.text.x = element_text(size=10 * resize_factor),
                         axis.text.y = element_text(size=10 * resize_factor), 
                         axis.title = element_text(size=13 * resize_factor),
                         legend.text = element_text(size = 12 * resize_factor),
                         strip.text = element_text(size = 10 * resize_factor))

# define dates of the input files
driver_results_date <- '2021_12_03' 
fusion_results_date <- '2021_12_03' 
p_val_cutoff <- 0.01
exclude_hypermutators <- FALSE

# read in metadata
metadata <- read_tsv(paste0(base_dir, '/path/to/metadata.tsv')) %>%
  # select relevant columns
  select(sample_id, cohort, is_hypermutated, cancer_type)

# exclude hypermutators?
if (exclude_hypermutators) {
  metadata <-  metadata %>%
    filter(is_hypermutated == FALSE)
}

# calculate sample size by strata
sample_size <- metadata %>%
  group_by(cancer_type, cohort) %>%
  summarise(sample_size = n(), .groups = 'drop') %>%
  pivot_wider(names_from = 'cohort', values_from = 'sample_size') %>%
  mutate(cancer_type_label = paste0(cancer_type, ' (', PCAWG, ' vs. ', HMF, ')'))

# calc cohort size
cohort_size <- metadata %>%
  group_by(cohort) %>%
  summarise(cohort_size = n(), .groups = 'drop') %>%
  pivot_wider(names_from = 'cohort', values_from = 'cohort_size') %>%
  mutate(cancer_type_label = paste0('Pancancer (', PCAWG, ' vs. ', HMF, ')')) %>%
  mutate(cancer_type = 'Pancancer', cancer_type_code = 'Pancancer')

# read in the completed table
linx_drivers_high_plot <- read_tsv(file = paste0(base_dir, '/path/to/', 
                                                 driver_results_date, '/full_drivers_per_sample_table.tsv'))


# --------------------------------------- Priestley et al. 2019 comparison

### Fig. 4a drivers per sample plot

if (exclude_hypermutators) {
  plot_subtitle <- 'Excluding hypermutators'
} else {
  plot_subtitle <- 'Including hypermutators'
}

# left_join metadata and sample size
linx_drivers_high_plot <- linx_drivers_high_plot %>%
  left_join(metadata, by = 'sample_id') %>%
  left_join(sample_size, by = 'cancer_type')

# define cancer type label order
source(paste0(base_dir, '/path/to/r_objects/cancer_type_order.R'))

label_order <- sample_size %>% arrange(factor(cancer_type, levels = cancer_type_order)) %>% pull(cancer_type_label)

# get the alternating index number for the grey background rectangles
rect_idx <- seq(2, 22, by = 2)
rect_idx <- label_order[rect_idx]

# finalize data to plot it, define the strata (stadium here)
linx_drivers_high_plot <- linx_drivers_high_plot %>%
  mutate(cancer_type = factor(cancer_type_label, levels = label_order)) %>%
  mutate(stadium = if_else(cohort == 'HMF', 'Hartwig', 'PCAWG'),
         stadium = factor(stadium, levels = c('Hartwig', 'PCAWG')))


#### plot clonal, subclonal or total drivers per sample
plot_clonality <- FALSE ### <-------------- CHANGE: TRUE = plot for clonality defined below will be plotted, FALSE = total drivers per sample will be plotted
clonality <- 'subclonal' ### <--------------- CHANGE: can be either 'clonal' or 'subclonal', plot_clonality must be TRUE

# calculate total drivers per sample (for total driver plot)
linx_drivers_high_plot_final <- linx_drivers_high_plot %>%
  { if (plot_clonality) filter(., clonality2 == clonality | clonality2 == 'none') else . } %>%
  group_by(sample_id) %>%
  mutate(total_drivers = sum(driver_per_clonality)) %>%
  ungroup()

# sum up clonal and subclonal MUTATION driver for total dps plots (necessary for facetted plot and mean calculation)
if(plot_clonality == FALSE) {
  linx_drivers_high_plot_final <- linx_drivers_high_plot_final %>%
    group_by(sample_id, cancer_type, cohort, driver) %>%
    mutate(driver_per_clonality = sum(driver_per_clonality)) %>%
    ungroup()
}

### NOTE: ONLY USED FOR FACETTED PLOT
if (plot_clonality) {
  
  # summarize the mean drivers per sample by cancer_type and per clonality
  linx_drivers_high_plot_mean <- linx_drivers_high_plot_final %>%
    mutate(cancer_type = factor(cancer_type_label)) %>%
    group_by(cancer_type, driver, cohort, clonality2) %>%
    summarise(mean_drivers = mean(driver_per_clonality),
              n = n(), .groups = 'drop') %>%
    filter(n > 1) %>%
    filter(clonality2 == clonality | clonality2 == 'none')
  
} else {
  
  # summarize the mean drivers per sample by cancer_type
  linx_drivers_high_plot_mean <- linx_drivers_high_plot_final %>%
    mutate(cancer_type = factor(cancer_type_label)) %>%
    group_by(cancer_type, cohort, driver) %>% # for facetted total plot
    summarise(mean_drivers = mean(driver_per_clonality),
              n = n(), .groups = 'drop') %>%
    filter(n > 1)
  
}

############ calculate the mean of mean dps across cancer types to get and idea on pancancer level
# total
linx_drivers_high_plot_mean_of_mean <- linx_drivers_high_plot_final %>%
  group_by(cancer_type, cohort) %>%
  summarise(mean_driver = mean(total_drivers), .groups = 'drop') %>%
  group_by(cohort) %>%
  summarise(mean_driver = mean(mean_driver), .groups = 'drop') %>%
  pivot_wider(., names_from = 'cohort', values_from = 'mean_driver') %>%
  mutate(mean_diff = round(HMF - PCAWG, 2)) %>%
  mutate(mean_diff_label = if_else(mean_diff > 0, paste0('+ ', mean_diff), paste0('- ', abs(mean_diff)))) %>%
  pivot_longer(., cols = c('HMF', 'PCAWG'), names_to = 'cohort', values_to = 'mean_driver')

# per driver type
linx_drivers_high_plot_mean_of_mean2 <- linx_drivers_high_plot_final %>%
  group_by(cancer_type, cohort, driver) %>%
  summarise(mean_driver = mean(driver_per_clonality), .groups = 'drop') %>%
  group_by(cohort, driver) %>%
  summarise(mean_driver = round(mean(mean_driver),1 ), .groups = 'drop') %>%
  pivot_wider(., names_from = 'cohort', values_from = 'mean_driver') %>%
  mutate(mean_diff = round(HMF - PCAWG, 2)) %>%
  mutate(mean_diff_label = if_else(mean_diff > 0, paste0('+ ', mean_diff), paste0('- ', abs(mean_diff)))) %>%
  pivot_longer(., cols = c('HMF', 'PCAWG'), names_to = 'cohort', values_to = 'mean_driver')
############

# new custom facet labels
if (plot_clonality) {
  facet_labels <- c(
    `AMP` = 'AMP',
    `DEL` = 'DEL',
    `FUSION` = 'FUS',
    `MUTATION` = paste0('MUT (', clonality, ')') 
  )
} else {
  facet_labels <- c(
    `AMP` = 'AMP',
    `DEL` = 'DEL',
    `FUSION` = 'FUS',
    `MUTATION` = 'MUT'
  )
}

# calculate driver per sample difference on avg per cancer type between PCAWG / Hartwig and on pancancer level
mean_driver_difference_overall <- linx_drivers_high_plot_final %>%
  group_by(cancer_type, cohort) %>%
  summarise(mean_driver = mean(total_drivers), .groups = 'drop') %>%
  pivot_wider(., names_from = 'cohort', values_from = 'mean_driver') %>%
  mutate(mean_diff = round(HMF - PCAWG, 2)) %>%
  mutate(mean_diff_label = if_else(mean_diff > 0, paste0('+ ', mean_diff), paste0('- ', abs(mean_diff)))) %>%
  pivot_longer(., cols = c('HMF', 'PCAWG'), names_to = 'cohort', values_to = 'mean_driver') %>%
  select(-mean_driver)

mean_driver_difference_overall_pancan <- linx_drivers_high_plot_final %>%
  group_by(cohort) %>%
  summarise(mean_driver = mean(total_drivers), .groups = 'drop') %>%
  pivot_wider(., names_from = 'cohort', values_from = 'mean_driver') %>%
  mutate(mean_diff = round(HMF - PCAWG, 2)) %>%
  mutate(mean_diff_label = if_else(mean_diff > 0, paste0('+ ', mean_diff), paste0('- ', abs(mean_diff)))) %>%
  pivot_longer(., cols = c('HMF', 'PCAWG'), names_to = 'cohort', values_to = 'mean_driver') %>%
  select(-mean_driver)

# test del: only keep significant mean increase / decrease
sig_total <- c(
  'Breast cancer (201 vs. 787)',
  'Esophagus cancer (85 vs. 140)',
  'Stomach cancer (55 vs. 43)',
  'Cholangiocarcinoma (26 vs. 66)',
  'Hepatocellular carcinoma (288 vs. 53)',
  'Pancreas neuroendocrine (81 vs. 38)',
  'Uterus carcinoma (42 vs. 33)',
  'Kidney clear cell carcinoma (109 vs. 129)',
  'Non small cell lung cancer (83 vs. 511)',
  'Prostate carcinoma (153 vs. 404)',
  'Liposarcoma (17 vs. 25)',
  'Thyroid cancer (44 vs. 22)'
  )

sig_amp <- c(
  'Breast cancer (201 vs. 787)',
  'Esophagus cancer (85 vs. 140)',
  'Hepatocellular carcinoma (288 vs. 53)',
  'Prostate carcinoma (153 vs. 404)',
  'Skin melanoma (104 vs. 303)'
)

sig_del <- c(
  'Pancreas neuroendocrine (81 vs. 38)',
  'Kidney clear cell carcinoma (109 vs. 129)',
  'Prostate carcinoma (153 vs. 404)',
  'Thyroid cancer (44 vs. 22)'
)

sig_fus <- c(
  'Thyroid cancer (44 vs. 22)'
)

sig_mut <- c(
  'Breast cancer (201 vs. 787)',
  'Colorectum carcinoma (51 vs. 628)',
  'Hepatocellular carcinoma (288 vs. 53)',
  'Pancreas neuroendocrine (81 vs. 38)',
  'Kidney clear cell carcinoma (109 vs. 129)',
  'Non small cell lung cancer (83 vs. 511)',
  'Prostate carcinoma (153 vs. 404)',
  'Liposarcoma (17 vs. 25)',
  'Thyroid cancer (44 vs. 22)'
)

# import cohort color
source(paste0(base_dir, '/path/to/r_objects/color_codes.R'))

# ------------------------------------ plot drivers per sample overall

total_y_limit <- 25

# pancancer
dps_overall_pancan <- linx_drivers_high_plot_final %>%
  inner_join(mean_driver_difference_overall_pancan, by = c('cohort')) %>%
  mutate(cancer_type = 'Pancancer') %>%
  select(-cancer_type_label, -HMF, -PCAWG) %>%
  inner_join(cohort_size, by = c('cancer_type')) %>%
  { if (plot_clonality) filter(., clonality2 == clonality | clonality2 == 'none') else . } %>% # optional clonality filter
  ggplot(., aes(x = stadium, y = total_drivers)) +
  geom_rect(
    data = . %>%
      distinct(cancer_type),
    inherit.aes = FALSE, # fixed the inheritance error
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#a0a0a0') +
  geom_jitter(color = 'gray40', size = 0.1, alpha = 0.15) +
  geom_violin(aes(color = stadium), fill = NA, position = position_dodge(1), size = 1.01, scale = 'width', width = 0.5, adjust = 2.5) +
  geom_text(data = . %>% filter(cohort == 'HMF') %>% distinct(stadium, cancer_type_label, mean_diff_label), aes(label = mean_diff_label, y = total_y_limit * 1.3), size = 3, nudge_x = 0.5) +
  stat_summary(aes(group = stadium), fun=mean, geom='point', shape=16, size=2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  stat_compare_means(aes(group = stadium), method = 'wilcox.test', label.y = total_y_limit * 1.08,
                     vjust = 0, size = 4, label = 'p.signif', hide.ns = TRUE,
                     symnum.args = list(cutpoints = c(0, p_val_cutoff, 1), symbols = c('*', 'ns'))) +
  facet_wrap(cancer_type_label ~ ., ncol = 1, strip.position = 'left') +
  scale_color_manual(
    values = rev(cohort_colors),
    guide = guide_legend(title = '', reverse = TRUE, override.aes = list(fill = cohort_colors))
  ) +
  coord_flip(
    ylim = c(0,total_y_limit),
    clip = 'off'
  ) +
  labs(
    title = 'Total',
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

# export: 800 x 800
dps_overall_percan <- linx_drivers_high_plot_final %>%
  inner_join(mean_driver_difference_overall, by = c('cancer_type', 'cohort')) %>%
  { if (plot_clonality) filter(., clonality2 == clonality | clonality2 == 'none') else . } %>% # optional clonality filter
  ggplot(., aes(x = stadium, y = total_drivers)) +
  geom_rect(data = linx_drivers_high_plot_final %>% 
              filter(cohort == 'HMF') %>%
              filter(cancer_type %in% rect_idx) %>%
              distinct(cancer_type),
            inherit.aes = FALSE, # fixed the inheritance error
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#a0a0a0') +
  geom_jitter(color = 'gray40', size = 0.1, alpha = 0.15) +
  geom_violin(aes(color = stadium), fill = NA, position = position_dodge(1), size = 1.01, scale = 'width', width = 0.5, adjust = 2.5) +
  geom_text(data = . %>% filter(cohort == 'HMF') %>% filter(cancer_type %in% sig_total) %>% distinct(stadium, cancer_type, mean_diff_label), 
            aes(label = mean_diff_label, y = total_y_limit * 1.3), size = 3, nudge_x = 0.5) +
  stat_summary(aes(group = stadium), fun=mean, geom='point', shape=16, size=2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  stat_compare_means(aes(group = stadium), method = 'wilcox.test', label.y = total_y_limit * 1.08,
                     vjust = 0, size = 4, label = 'p.signif', hide.ns = TRUE,
                     symnum.args = list(cutpoints = c(0, p_val_cutoff, 1), symbols = c('*', 'ns'))) +
  facet_wrap(cancer_type ~ ., ncol = 1, strip.position = 'left') +
  scale_color_manual(
    values = rev(cohort_colors),
    guide = guide_legend(title = '', reverse = TRUE, override.aes = list(fill = cohort_colors))
    ) +
  coord_flip(
    ylim = c(0,total_y_limit),
    clip = 'off'
    ) +
  labs(
    x = '',
    y = 'Drivers per sample') +
  theme_classic() +
  theme_smallfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = 'right',
    panel.spacing = unit(0, 'cm'),
    plot.margin = unit(c(1,0,1,0), 'cm')
    )

# ------------------------------------

# calculate the mean drivers per driver type & join to the final table
mean_driver_difference_per_driver <- linx_drivers_high_plot_final %>%
  group_by(cancer_type, cohort, driver) %>%
  summarise(mean_driver = mean(driver_per_clonality), .groups = 'drop') %>%
  pivot_wider(., names_from = 'cohort', values_from = 'mean_driver') %>%
  mutate(mean_diff = round(HMF - PCAWG, 2)) %>%
  mutate(mean_diff_label = if_else(mean_diff > 0, paste0('+ ', mean_diff), paste0('- ', abs(mean_diff)))) %>%
  pivot_longer(., cols = c('HMF', 'PCAWG'), names_to = 'cohort', values_to = 'mean_driver') %>%
  select(-mean_driver)

mean_driver_difference_per_driver_pancan <- linx_drivers_high_plot_final %>%
  group_by(cohort, driver) %>%
  summarise(mean_driver = mean(driver_per_clonality), .groups = 'drop') %>%
  pivot_wider(., names_from = 'cohort', values_from = 'mean_driver') %>%
  mutate(mean_diff = round(HMF - PCAWG, 2)) %>%
  mutate(mean_diff_label = if_else(mean_diff > 0, paste0('+ ', mean_diff), paste0('- ', abs(mean_diff)))) %>%
  pivot_longer(., cols = c('HMF', 'PCAWG'), names_to = 'cohort', values_to = 'mean_driver') %>%
  select(-mean_driver)

linx_drivers_high_plot_final_percan <- linx_drivers_high_plot_final %>%
  inner_join(mean_driver_difference_per_driver, by = c('cancer_type', 'cohort', 'driver'))

linx_drivers_high_plot_final_pancan <- linx_drivers_high_plot_final %>%
  inner_join(mean_driver_difference_per_driver_pancan, by = c('cohort', 'driver')) %>%
  mutate(cancer_type = 'Pancancer') %>%
  # get the cohort size in there
  select(-cancer_type_label, -HMF, -PCAWG) %>%
  inner_join(cohort_size, by = c('cancer_type'))

# define ylim for each subplot (to align labels in plot at the same position relative to ylim)
amp_y_limit <- 16
del_y_limit <- 20
fus_y_limit <- 4
mut_y_limit <- 20

# ------------------------------------ plot amplification dps violin plot

# pancancer
amp_violin_pancan <- linx_drivers_high_plot_final_pancan %>%
  filter(driver == 'AMP') %>%
  { if (plot_clonality) filter(., clonality2 == clonality | clonality2 == 'none') else . } %>% # optional clonality filter
  ggplot(., aes(x = stadium, y = driver_per_clonality)) +
  geom_rect(
    data = . %>%
      distinct(cancer_type),
    inherit.aes = FALSE, # fixed the inheritance error
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#a0a0a0') +
  geom_jitter(color = 'gray40', size = 0.1, alpha = 0.25) +
  geom_violin(aes(color = stadium), fill = NA, position = position_dodge(1), size = 1.01, scale = 'width', width = 0.5, adjust = 2.5) +
  geom_text(data = . %>% filter(cohort == 'HMF') %>% distinct(stadium, cancer_type_label, mean_diff_label),aes(label = mean_diff_label, y = amp_y_limit * 1.3), nudge_x = 0.5) +
  stat_summary(aes(group = stadium), fun=mean, geom='point', shape=16, size=2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  stat_compare_means(aes(group = stadium), method = 'wilcox.test', label.y = amp_y_limit * 1.08,
                     vjust = 0, size = 4, label = 'p.signif', hide.ns = TRUE,
                     symnum.args = list(cutpoints = c(0, p_val_cutoff, 1), symbols = c('*', 'ns'))) +
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
    title = 'AMP',
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

amp_violin <- linx_drivers_high_plot_final_percan %>%
  filter(driver == 'AMP') %>%
  { if (plot_clonality) filter(., clonality2 == clonality | clonality2 == 'none') else . } %>% # optional clonality filter
  ggplot(., aes(x = stadium, y = driver_per_clonality)) +
  geom_rect(data = linx_drivers_high_plot_final %>% 
              filter(cohort == 'HMF') %>%
              filter(cancer_type %in% rect_idx) %>%
              distinct(cancer_type),
            inherit.aes = FALSE, # fixed the inheritance error
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#a0a0a0') +
  geom_jitter(color = 'gray40', size = 0.1, alpha = 0.25) +
  geom_violin(aes(color = stadium), fill = NA, position = position_dodge(1), size = 1.01, scale = 'width', width = 0.5, adjust = 2.5) +
  geom_text(data = . %>% filter(cohort == 'HMF') %>% filter(cancer_type %in% sig_amp) %>% distinct(stadium, cancer_type, mean_diff_label),
            aes(label = mean_diff_label, y = amp_y_limit * 1.3), nudge_x = 0.5) +
  stat_summary(aes(group = stadium), fun=mean, geom='point', shape=16, size=2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  stat_compare_means(aes(group = stadium), method = 'wilcox.test', label.y = amp_y_limit * 1.08,
                     vjust = 0, size = 4, label = 'p.signif', hide.ns = TRUE,
                     symnum.args = list(cutpoints = c(0, p_val_cutoff, 1), symbols = c('*', 'ns'))) +
  facet_wrap(cancer_type ~ ., ncol = 1, strip.position = 'left') +
  scale_color_manual(
    values = rev(cohort_colors),
    guide = guide_legend(title = '', reverse = TRUE, override.aes = list(fill = cohort_colors))
  ) +
  coord_flip(
    ylim = c(0, amp_y_limit),
    clip = 'off'
  ) +
  labs(
    title = 'AMP',
    x = '',
    y = 'Drivers per sample') +
  theme_classic() +
  theme_smallfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.title = element_blank(),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = 'right',
    panel.spacing = unit(0, 'cm'),
    plot.margin = unit(c(1,0,1,0), 'cm'))

# ------------------------------------ plot deletion dps violin plot

# pancancer
del_violin_pancan <- linx_drivers_high_plot_final_pancan %>%
  filter(driver == 'DEL') %>%
  { if (plot_clonality) filter(., clonality2 == clonality | clonality2 == 'none') else . } %>% # optional clonality filter
  ggplot(., aes(x = stadium, y = driver_per_clonality)) +
  geom_rect(
    data = . %>%
      distinct(cancer_type),
    inherit.aes = FALSE, # fixed the inheritance error
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#a0a0a0') +
  geom_jitter(color = 'gray40', size = 0.1, alpha = 0.25) +
  geom_violin(aes(color = stadium), fill = NA, position = position_dodge(1), size = 1.01, scale = 'width', width = 0.5, adjust = 2.5) +
  geom_text(data = . %>% filter(cohort == 'HMF') %>% distinct(stadium, cancer_type_label, mean_diff_label),aes(label = mean_diff_label, y = del_y_limit * 1.3), nudge_x = 0.5) +
  stat_summary(aes(group = stadium), fun=mean, geom='point', shape=16, size=2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  stat_compare_means(aes(group = stadium), method = 'wilcox.test', label.y = del_y_limit * 1.08,
                     vjust = 0, size = 4, label = 'p.signif', hide.ns = TRUE,
                     symnum.args = list(cutpoints = c(0, p_val_cutoff, 1), symbols = c('*', 'ns'))) +
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
    title = 'DEL',
    x = '',
    y = 'Drivers per sample') +
  theme_classic() +
  theme_smallfont +
  theme(
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

del_violin <- linx_drivers_high_plot_final_percan %>%
  filter(driver == 'DEL') %>%
  { if (plot_clonality) filter(., clonality2 == clonality | clonality2 == 'none') else . } %>% # optional clonality filter
  ggplot(., aes(x = stadium, y = driver_per_clonality)) +
  geom_rect(data = linx_drivers_high_plot_final %>% 
              filter(cohort == 'HMF') %>%
              filter(cancer_type %in% rect_idx) %>%
              distinct(cancer_type),
            inherit.aes = FALSE, # fixed the inheritance error
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#a0a0a0') +
  geom_jitter(color = 'gray40', size = 0.1, alpha = 0.25) +
  geom_violin(aes(color = stadium), fill = NA, position = position_dodge(1), size = 1.01, scale = 'width', width = 0.5, adjust = 2.5) +
  geom_text(data = . %>% filter(cohort == 'HMF') %>% filter(cancer_type %in% sig_del) %>% distinct(stadium, cancer_type, mean_diff_label),aes(label = mean_diff_label, y = del_y_limit * 1.3), nudge_x = 0.5) +
  stat_summary(aes(group = stadium), fun=mean, geom='point', shape=16, size=2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  stat_compare_means(aes(group = stadium), method = 'wilcox.test', label.y = del_y_limit * 1.08,
                     vjust = 0, size = 4, label = 'p.signif', hide.ns = TRUE,
                     symnum.args = list(cutpoints = c(0, p_val_cutoff, 1), symbols = c('*', 'ns'))) +
  facet_wrap(cancer_type ~ ., ncol = 1, strip.position = 'left') +
  scale_color_manual(
    values = rev(cohort_colors),
    guide = guide_legend(title = '', reverse = TRUE, override.aes = list(fill = cohort_colors))
  ) +
  coord_flip(
    ylim = c(0, del_y_limit),
    clip = 'off'
  ) +
  labs(
    title = 'DEL',
    x = '',
    y = 'Drivers per sample') +
  theme_classic() +
  theme_smallfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.title = element_blank(),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = 'right',
    panel.spacing = unit(0, 'cm'),
    plot.margin = unit(c(1,0,1,0), 'cm'))

# ------------------------------------ plot fusion dps violin plot

# pancancer
fus_violin_pancan <- linx_drivers_high_plot_final_pancan %>%
  filter(driver == 'FUSION') %>%
  { if (plot_clonality) filter(., clonality2 == clonality | clonality2 == 'none') else . } %>% # optional clonality filter
  ggplot(., aes(x = stadium, y = driver_per_clonality)) +
  geom_rect(
    data = . %>%
      distinct(cancer_type),
    inherit.aes = FALSE, # fixed the inheritance error
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#a0a0a0') +
  geom_jitter(color = 'gray40', size = 0.1, alpha = 0.25) +
  geom_violin(aes(color = stadium), fill = NA, position = position_dodge(1), size = 1.01, scale = 'width', width = 0.5, adjust = 2.5) +
  stat_summary(aes(group = stadium), fun=mean, geom='point', shape=16, size=2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  stat_compare_means(aes(group = stadium), method = 'wilcox.test', label.y = fus_y_limit * 1.08,
                     vjust = 0, size = 4, label = 'p.signif', hide.ns = TRUE,
                     symnum.args = list(cutpoints = c(0, p_val_cutoff, 1), symbols = c('*', 'ns'))) +
  facet_wrap(cancer_type_label ~ ., ncol = 1, strip.position = 'left') +
  scale_color_manual(
    values = rev(cohort_colors),
    guide = guide_legend(title = '', reverse = TRUE, override.aes = list(fill = cohort_colors))
  ) +
  coord_flip(
    ylim = c(0, fus_y_limit),
    clip = 'off'
  ) +
  scale_y_continuous(limits = c(0, fus_y_limit), oob = scales::censor) + # for some reason ggplot plots points below 0 for fusions, so need to add this after coord_flip()
  labs(
    title = 'FUS',
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

fus_violin <- linx_drivers_high_plot_final_percan %>%
  filter(driver == 'FUSION') %>%
  { if (plot_clonality) filter(., clonality2 == clonality | clonality2 == 'none') else . } %>% # optional clonality filter
  ggplot(., aes(x = stadium, y = driver_per_clonality)) +
  geom_rect(data = linx_drivers_high_plot_final %>% 
              filter(cohort == 'HMF') %>%
              filter(cancer_type %in% rect_idx) %>%
              distinct(cancer_type),
            inherit.aes = FALSE, # fixed the inheritance error
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#a0a0a0') +
  geom_jitter(color = 'gray40', size = 0.1, alpha = 0.25) +
  geom_violin(aes(color = stadium), fill = NA, position = position_dodge(1), size = 1.01, scale = 'width', width = 0.5, adjust = 2.5) +
  stat_summary(aes(group = stadium), fun=mean, geom='point', shape=16, size=2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  stat_compare_means(aes(group = stadium), method = 'wilcox.test', label.y = fus_y_limit * 1.08,
                     vjust = 0, size = 4, label = 'p.signif', hide.ns = TRUE,
                     symnum.args = list(cutpoints = c(0, p_val_cutoff, 1), symbols = c('*', 'ns'))) +
  facet_wrap(cancer_type ~ ., ncol = 1, strip.position = 'left') +
  scale_color_manual(
    values = rev(cohort_colors),
    guide = guide_legend(title = '', reverse = TRUE, override.aes = list(fill = cohort_colors))
  ) +
  coord_flip(
    ylim = c(0, fus_y_limit),
    clip = 'off'
  ) +
  scale_y_continuous(limits = c(0, fus_y_limit), oob = scales::censor) + # for some reason ggplot plots points below 0 for fusions, so need to add this after coord_flip()
  labs(
    title = 'FUS',
    x = '',
    y = 'Drivers per sample') +
  theme_classic() +
  theme_smallfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.title = element_blank(),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = 'right',
    panel.spacing = unit(0, 'cm'),
    plot.margin = unit(c(1,0,1,0), 'cm'))

# ------------------------------------ plot (clonal/subclonal) point mutation dps violin plot

# pancancer
mut_violin_pancan <- linx_drivers_high_plot_final_pancan %>%
  filter(driver == 'MUTATION') %>%
  { if (plot_clonality) filter(., clonality2 == clonality | clonality2 == 'none') else . } %>% # optional clonality filter
  ggplot(., aes(x = stadium, y = driver_per_clonality)) +
  geom_rect(
    data = . %>%
      distinct(cancer_type),
    inherit.aes = FALSE, # fixed the inheritance error
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#a0a0a0') +
  geom_jitter(color = 'gray40', size = 0.1, alpha = 0.25) +
  geom_violin(aes(color = stadium), fill = NA, position = position_dodge(1), size = 1.01, scale = 'width', width = 0.5, adjust = 2.5) +
  geom_text(data = . %>% filter(cohort == 'HMF') %>% distinct(stadium, cancer_type_label, mean_diff_label),aes(label = mean_diff_label, y = mut_y_limit * 1.3), nudge_x = 0.5) +
  stat_summary(aes(group = stadium), fun=mean, geom='point', shape=16, size=2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  stat_compare_means(aes(group = stadium), method = 'wilcox.test', label.y = mut_y_limit * 1.08,
                     vjust = 0, size = 4, label = 'p.signif', hide.ns = TRUE,
                     symnum.args = list(cutpoints = c(0, p_val_cutoff, 1), symbols = c('*', 'ns'))) +
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
    title = 'MUT',
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

mut_violin <- linx_drivers_high_plot_final_percan %>%
  filter(driver == 'MUTATION') %>%
  { if (plot_clonality) filter(., clonality2 == clonality | clonality2 == 'none') else . } %>% # optional clonality filter
  ggplot(., aes(x = stadium, y = driver_per_clonality)) +
  geom_rect(data = linx_drivers_high_plot_final %>% 
              filter(cohort == 'HMF') %>%
              filter(cancer_type %in% rect_idx) %>%
              distinct(cancer_type),
            inherit.aes = FALSE, # fixed the inheritance error
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#a0a0a0') +
  geom_jitter(color = 'gray40', size = 0.1, alpha = 0.25) +
  geom_violin(aes(color = stadium), fill = NA, position = position_dodge(1), size = 1.01, scale = 'width', width = 0.5, adjust = 2.5) +
  geom_text(data = . %>% filter(cohort == 'HMF') %>% filter(cancer_type %in% sig_mut) %>% distinct(stadium, cancer_type, mean_diff_label),aes(label = mean_diff_label, y = mut_y_limit * 1.3), nudge_x = 0.5) +
  stat_summary(aes(group = stadium), fun=mean, geom='point', shape=16, size=2, color = 'black',
               position = position_dodge(1), show.legend = FALSE) + # for unfacetted plot
  stat_compare_means(aes(group = stadium), method = 'wilcox.test', label.y = mut_y_limit * 1.08,
                     vjust = 0, size = 4, label = 'p.signif', hide.ns = TRUE,
                     symnum.args = list(cutpoints = c(0, p_val_cutoff, 1), symbols = c('*', 'ns'))) +
  facet_wrap(cancer_type ~ ., ncol = 1, strip.position = 'left') +
  scale_color_manual(
    values = rev(cohort_colors),
    guide = guide_legend(title = '', reverse = TRUE, override.aes = list(fill = cohort_colors))
  ) +
  coord_flip(
    ylim = c(0, mut_y_limit),
    clip = 'off'
  ) +
  labs(
    title = 'MUT',
    x = '',
    y = 'Drivers per sample') +
  theme_classic() +
  theme_smallfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.title = element_blank(),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = 'right',
    panel.spacing = unit(0, 'cm'),
    plot.margin = unit(c(1,0,1,0), 'cm'))

# plot individually
total_assembled_plot <- dps_overall_pancan + dps_overall_percan + plot_layout(nrow = 2, ncol = 1, heights = c(1,22))
amp_assembled_plot <- amp_violin_pancan + amp_violin + plot_layout(nrow = 2, ncol = 1, heights = c(1,22))
del_assembled_plot <- del_violin_pancan + del_violin + plot_layout(nrow = 2, ncol = 1, heights = c(1,22))
fus_assembled_plot <- fus_violin_pancan + fus_violin + plot_layout(nrow = 2, ncol = 1, heights = c(1,22))
mut_assembled_plot <- mut_violin_pancan + mut_violin + plot_layout(nrow = 2, ncol = 1, heights = c(1,22))

### save PNG file of the plots
# total
ggsave(
  path = output_dir,
  filename = paste0(current_date, '_total_dps_q', p_val_cutoff, '.png'),
  plot = total_assembled_plot,
  device = 'png',
  width = 6,
  height = 7,
  dpi = 600
)

pdf(
  file = paste0(output_dir, '/', current_date, '_total_dps_q', p_val_cutoff, '.pdf'),
  width = 6, 
  height = 7,
  useDingbats = FALSE,
  compress = FALSE
)
print(total_assembled_plot)
dev.off()

# amplification
ggsave(
  path = output_dir,
  filename = paste0(current_date, '_amp_dps_q', p_val_cutoff, '.png'),
  plot = amp_assembled_plot,
  device = 'png',
  width = 6,
  height = 7,
  dpi = 600
)

pdf(
  file = paste0(output_dir, '/', current_date, '_amp_dps_q', p_val_cutoff, '.pdf'),
  width = 6, 
  height = 7,
  useDingbats = FALSE,
  compress = FALSE
)
amp_assembled_plot
dev.off()

# deletion
ggsave(
  path = output_dir,
  filename = paste0(current_date, '_del_dps_q', p_val_cutoff, '.png'),
  plot = del_assembled_plot,
  device = 'png',
  width = 6,
  height = 7,
  dpi = 600
)

pdf(
  file = paste0(output_dir, '/', current_date, '_del_dps_q', p_val_cutoff, '.pdf'),
  width = 6, 
  height = 7,
  useDingbats = FALSE,
  compress = FALSE
)
del_assembled_plot
dev.off()

# fusion
ggsave(
  path = output_dir,
  filename = paste0(current_date, '_fus_dps_q', p_val_cutoff, '.png'),
  plot = fus_assembled_plot,
  device = 'png',
  width = 6,
  height = 7,
  dpi = 600
)

pdf(
  file = paste0(output_dir, '/', current_date, '_fus_dps_q', p_val_cutoff, '.pdf'),
  width = 6, 
  height = 7,
  useDingbats = FALSE,
  compress = FALSE
)
fus_assembled_plot
dev.off()

# point mutation
ggsave(
  path = output_dir,
  filename = paste0(current_date, '_mut_dps_q', p_val_cutoff, '.png'),
  plot = mut_assembled_plot,
  device = 'png',
  width = 6,
  height = 7,
  dpi = 600
)

pdf(
  file = paste0(output_dir, '/', current_date, '_mut_dps_q', p_val_cutoff, '.pdf'),
  width = 6, 
  height = 7,
  useDingbats = FALSE,
  compress = FALSE
)
mut_assembled_plot
dev.off()

# save SVG files for all the individual plots
svglite(
  filename = paste0(output_dir, '/', current_date, '_total_dps_q', p_val_cutoff, '.svg'),
  width = 6, 
  height = 7
)
total_assembled_plot
dev.off()

svglite(
  filename = paste0(output_dir, '/', current_date, '_amp_dps_q', p_val_cutoff, '.svg'),
  width = 8, 
  height = 8
)
amp_assembled_plot
dev.off()

svglite(
  filename = paste0(output_dir, '/', current_date, '_del_dps_q', p_val_cutoff, '.svg'),
  width = 8, 
  height = 8
)
del_assembled_plot
dev.off()

svglite(
  filename = paste0(output_dir, '/', current_date, '_fus_dps_q', p_val_cutoff, '.svg'),
  width = 8, 
  height = 8
)
fus_assembled_plot
dev.off()

svglite(
  filename = paste0(output_dir, '/', current_date, '_mut_dps_q', p_val_cutoff, '.svg'),
  width = 8, 
  height = 8
)
mut_assembled_plot
dev.off()
