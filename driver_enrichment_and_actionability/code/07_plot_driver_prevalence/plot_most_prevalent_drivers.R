############### Plot the most prevalent somatic mutated driver genes ###############
# author: remy (sascha)
# date: 12/11/2021
# last updated: 24/10/2022

### Description
# This script creates all the plots according to Fig.3 in Priestley et al. 2019.
# Uses the driver prevalence table from create_driver_prevalence_table.R as input.

### Input
# driver prevalence table created in step 7.

### Output
# Supplementary Figure 5b

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
dir.create(path = paste0(output_dir, '/07_plot_driver_prevalence/plots'), recursive = TRUE)

# get todays date
current_date <- format(Sys.Date(), '%Y_%m_%d')

## custom fonts
font_add_google(
  name = 'Inter',
  family = 'Helvetica'
)

## ggplot custom theme
theme_smallfont <- theme(plot.title = element_text(size = 17),
                         plot.subtitle = element_text(size = 11),
                         axis.text.x = element_text(size = 12),
                         axis.text.y = element_text(size = 12), 
                         axis.title = element_text(size = 13),
                         legend.text = element_text(size = 12),
                         strip.text.y = element_text(size = 12))

#  ----------------------- Custom functions & objects

# load cohort color codes
source(paste0(base_dir, '/code/r_objects/color_codes.R'))

# ----------------------- DRIVER PREVALENCE

driver_prevalence <- read_tsv(file = paste0(base_dir, '/results/07_plot_driver_prevalence/driver_prevalence.tsv'))

# calculate percentage of most prevalent driver genes across cohorts
most_prevalent_drivers_across_cohorts <- driver_prevalence %>%
  distinct(sample_id, gene, .keep_all = TRUE) %>%
  group_by(cohort, gene) %>%
  mutate(n_cohort = n()) %>%
  ungroup() %>%
  mutate(percent_cohort = n_cohort / cohort_size * 100) %>%
  distinct(cohort, gene, category, percent_cohort, .keep_all = TRUE) %>%
  arrange(category, cohort, desc(percent_cohort))

# most prevalent drivers per cancer type
most_prevalent_drivers_across_cancer_types <- driver_prevalence %>%
  distinct(sample_id, gene, .keep_all = TRUE) %>%
  group_by(cancer_type, cancer_type_label, cohort, gene) %>%
  mutate(n_cancer_type = n()) %>%
  summarise(percent_cancer_type = n_cancer_type / sample_size * 100, .groups = 'drop') %>%
  distinct()

# label cosmetics
most_prevalent_drivers_across_cancer_types <- most_prevalent_drivers_across_cancer_types %>%
  mutate(cancer_type_label = str_replace(cancer_type_label, pattern = regex('\\s\\('), replacement = '\n('))

##############################################################################

# get the names of the top 20 onco / tsg in each cohort, then use those to filter 
# the most_prevalent_drivers_across_cohorts_split_up table
top_20_genes_onco <- most_prevalent_drivers_across_cohorts %>%
  filter(cohort == 'Hartwig', category == 'ONCO') %>%
  slice(1:20) %>% pull(gene)
top_20_genes_tsg <- most_prevalent_drivers_across_cohorts %>%
  filter(cohort == 'Hartwig', category == 'TSG') %>%
  slice(1:20) %>% pull(gene)

##############################################################################

### Rows split by cohort like karyotype heatmap

onco_split_heatmap <- most_prevalent_drivers_across_cancer_types %>%
  mutate(label_col = if_else(percent_cancer_type > 70, 'black', 'black')) %>%
  filter(gene %in% top_20_genes_onco) %>%
  mutate(gene = factor(gene, levels = top_20_genes_onco)) %>%
  ggplot(., aes(x = cohort, y = gene)) +
  geom_tile(aes(fill = percent_cancer_type)) +
  geom_text(aes(label = round(percent_cancer_type, 1), color = label_col)) +
  facet_grid(cancer_type_label ~ ., switch = 'y') +
  scale_color_manual(values = c('black', 'white')) +
  scale_fill_gradient(low = 'gray99', high = 'darkorange1') +
  scale_x_discrete(position = 'top') +
  coord_flip() +
  labs(
    title = paste0('Top 20 Hartwig ONCO drivers'),
    x = '',
    y = '',
    fill = '% cohort') +
  theme_classic() +
  theme_smallfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    plot.margin = unit(c(5,5,5,0), 'mm'),
    legend.position = 'none'
  )


tsg_split_heatmap <- most_prevalent_drivers_across_cancer_types %>%
  mutate(label_col = if_else(percent_cancer_type > 70, 'white', 'black')) %>%
  filter(gene %in% top_20_genes_tsg) %>%
  mutate(gene = factor(gene, levels = top_20_genes_tsg)) %>%
  ggplot(., aes(x = cohort, y = gene)) +
  geom_tile(aes(fill = percent_cancer_type)) +
  geom_text(aes(label = round(percent_cancer_type, 1), color = label_col)) +
  facet_grid(cancer_type_label ~ ., switch = 'y') +
  scale_color_manual(values = c('black', 'white')) +
  scale_fill_gradient(low = 'gray99', high = 'darkorchid3') +
  scale_x_discrete(position = 'top') +
  coord_flip() +
  labs(
    title = paste0('Top 20 Hartwig TSG drivers'),
    x = '',
    y = '',
    fill = '% cohort') +
  theme_classic() +
  theme_smallfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    strip.background = element_blank(),
    strip.text.y.left = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.5),
    plot.margin = unit(c(5,5,5,0), 'mm'),
    legend.position = 'none'
  )

# stitch plots together
assembled_plot <- onco_split_heatmap + tsg_split_heatmap + plot_layout(ncol = 2, nrow = 1, byrow = TRUE)

# save plot
# as PDF
pdf(
  file = paste0(output_dir, '/07_plot_driver_prevalence/plots/', current_date, '_top20drivers_combined_heatmap.pdf'),
  width = 24, 
  height = 13,
  useDingbats = FALSE,
  compress = FALSE
)
print(assembled_plot)
dev.off()
