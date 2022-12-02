############### Compare purity values per cancer types between primary and metastatic ###############
# author: remy (sascha)
# date: 20/09/2022
# last updated: 27/10/2022

### Description
# This script compares the purity values estimated by PURPLE between primary and metastatic samples
# for every cancer type.

### Input
# - metadata.tsv

### Output
# purity violin plots

### Usage
# Just run it

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

# get todays date
current_date <- format(Sys.Date(), '%Y_%m_%d')

# define variables
q_val_cutoff <- 0.01 # significance threshold for the FDR corrected p-value

## custom font
font_add_google(
  name = 'Inter',
  family = 'Helvetica'
)

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

# -------------------------- Custom objects and functions

# define cancer type order
source(paste0(base_dir, '/code/r_objects/cancer_type_order.R'))

# load cohort colors
source(paste0(base_dir, '/code/r_objects/color_codes.R'))

# custom function to calculate the Mann Whitney test p-value between diploid proportion of cohorts per cancer type
get_wilcox_table <- function(df_input, cancer_type_input) {
  
  print(paste0('processing...', cancer_type_input))
  
  # create a list of two dfs that represent the samples in the cohorts for one cancer type
  wilcox_list <- df_input %>%
    filter(cancer_type == cancer_type_input) %>%
    select(sample_id, cohort, tumor_purity) %>%
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
  wilcox_table <- tidy(wilcox.test(wilcox_list[[1]]$tumor_purity, 
                                   wilcox_list[[2]]$tumor_purity, 
                                   alternative = 'two.sided',
                                   paired = FALSE,
                                   mu = 0)) %>%
    mutate(cancer_type = cancer_type_input,
           significance_level = if_else(p.value < 0.01, '*', NA_character_)) %>%
    select(cancer_type, method, p.value, significance_level)
  
  return(wilcox_table)
  
}

# -------------------------- METADATA

metadata <- read_tsv(file = paste0(base_dir, '/data/processed/metadata/metadata.tsv')) %>%
  select(sample_id, cohort, cancer_type, tumor_purity)

# calc sample size
sample_size <- metadata %>%
  group_by(cancer_type, cohort) %>%
  summarise(sample_size_tot = n(), .groups = 'drop') %>%
  pivot_wider(names_from = c('cohort'), values_from = sample_size_tot, names_prefix = 'n_') %>%
  replace_na(., replace = list(n_Hartwig = 0, n_PCAWG = 0)) %>%
  mutate(cancer_type_label = str_c(cancer_type, ' (', n_PCAWG, ' vs. ', n_Hartwig, ')')) %>%
  ungroup()

# join sample size to metadata
metadata <- metadata %>%
  inner_join(sample_size, by = 'cancer_type')

# get the cancer type names into a vector
cancer_type_vec <- unique(metadata$cancer_type)

# calculate the Mann whitney p-value for all cancer types
wilcox_table <- map_df(cancer_type_vec, ~get_wilcox_table(metadata, .x))

# adjust p-values
wilcox_table <- wilcox_table %>%
  mutate(pval_adj = p.adjust(p.value, method = 'fdr'), .before = significance_level) %>%
  mutate(significance_level = if_else(pval_adj < q_val_cutoff, '*', NA_character_))

# calculate median diploid proportion per cancer type and cohort then add this info to the wilcox_table (for point layer in plot)
median_purity <- metadata %>%
  group_by(cancer_type, cancer_type_label, cohort) %>%
  summarise(median_purity = median(tumor_purity), .groups = 'drop')
wilcox_table <- wilcox_table %>%
  left_join(median_purity, by = c('cancer_type'))

# get the alternating index number for the grey background rectangles
rect_idx <- seq(2, 22, by = 2)
rect_idx <- cancer_type_order[rect_idx]

# get the cancer type label order
cancer_type_label_order <-metadata %>%
  arrange(factor(cancer_type, levels = cancer_type_order)) %>%
  distinct(cancer_type_label) %>%
  pull(cancer_type_label)

# get the cancer types in the right order
# NOTE: HAVE TO DO THIS FOR ALL TABLES THAT ARE USED IN THE PLOT
metadata <- metadata %>%
  mutate(cancer_type_label = factor(cancer_type_label, levels = cancer_type_label_order)) %>%
  mutate(stadium = if_else(cohort == 'Hartwig', 'Hartwig', 'PCAWG'),
         stadium = factor(stadium, levels = c('Hartwig', 'PCAWG')))
wilcox_table <- wilcox_table %>%
  mutate(cancer_type_label = factor(cancer_type_label, levels = cancer_type_label_order)) %>%
  mutate(stadium = if_else(cohort == 'Hartwig', 'Hartwig', 'PCAWG'),
         stadium = factor(stadium, levels = c('Hartwig', 'PCAWG')))

# plot the violin facets
purity_violin <- metadata %>%
  ggplot(., aes(x = stadium, y = tumor_purity)) +
  geom_rect(data = . %>% 
              filter(cancer_type %in% rect_idx) %>%
              distinct(cancer_type_label), aes(fill = cancer_type_label),
            inherit.aes = FALSE,
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.3, fill = '#a0a0a0') +
  geom_violin(aes(fill = cohort), color = 'black', scale = 'width') +
  geom_point(data = wilcox_table, 
             aes(y = median_purity), size = 2) +
  geom_text(data = wilcox_table %>% filter(cohort == 'Hartwig'), 
            aes(label = significance_level, y = 1.1), size = 8, nudge_x = 0) +
  facet_wrap(cancer_type_label ~ ., ncol = 1, strip.position = 'left') +
  scale_x_discrete(
    labels = c('Hartwig', 'PCAWG'),
    position = 'top',
    expand = c(0.5,0.5)) + # 'expand' is used to center the histogram in the panel
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1.2), expand = c(0,0)) +
  scale_fill_manual(values = rev(cohort_colors)) +
  labs(title = 'Purity',
       x = '',
       y = '') +
  coord_flip() +
  theme_classic() +
  theme_bigfont_barplot +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    legend.position = 'none',
    panel.spacing = unit(0, 'cm'))

# create output dir
output_dir <- paste0(output_dir, '/12_purity/03_plot_pure_samples/', current_date, '/results/plots')
dir.create(path = paste0(output_dir), recursive = TRUE)

# save as PDF
pdf(
  file = paste0(output_dir, '/', current_date, '_purity_comparison.pdf'),
  width = 8,
  height = 10,
  useDingbats = FALSE,
  compress = FALSE
)
print(purity_violin)
dev.off()
