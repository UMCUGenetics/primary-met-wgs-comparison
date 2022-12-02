############### PCAWG vs. Hartwig ###########################################################
############### Analyse the variants that are private and shared between pipelines ##########
# author: remy (sascha)
# date: 02/11/2021
# last updated: 25/04/2022

### Description
# This script investigates and plots the amount of private and shared variants 
# called for PCAWG samples by the Hartwig and the PCAWG pipeline

### Input
# The variant privacy table you generated in step 12.

### Output
# Supplementary Figure 3A in supplementary note 1

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

# create output dir
output_dir <- paste0(output_dir, '/10_private_vs_shared_mutations/', current_date, '/results/plots')
dir.create(path = paste0(output_dir), recursive = TRUE)

## ggplot custom theme
theme_bigfont <- theme(plot.title = element_text(size = 22),
                       plot.subtitle = element_text(size = 16),
                       axis.text.x = element_text(size = 15),
                       axis.text.y = element_text(size = 15), 
                       axis.title = element_text(size = 18),
                       legend.title = element_text(size = 16),
                       legend.text = element_text(size = 14),
                       strip.text.x = element_text(size = 15))

# def vars
coverage_info_results_date <- '2022_02_07' # date name of your PCAWG coverage table you generated in step 7
variant_privacy_results_date <- '2022_02_17' # date name of your extract variants from step 10

# ----------------------------- VARIANT PRIVACY

variant_privacy <- read_tsv(paste0(base_dir, '/results/10_private_vs_shared_mutations/', 
                                   variant_privacy_results_date,'/results/variant_privacy.tsv')) %>%
  # recode 'HMF_only' to 'Hartwig_only'
  mutate(privacy = if_else(privacy == 'Hartwig_only', 'Hartwig_only', privacy))

# ----------------------------- METADATA

metadata <- read_tsv(paste0(base_dir, '/data/processed/metadata/metadata.tsv'),
                     col_types = 'icccccccccccccccllllcccdiiildclccllllllliiiiiiiiiii')

# ----------------------------- PCAWG COVERAGE

coverage <- read_tsv(file = paste0(base_dir, '/results/7_create_coverage_tables/', coverage_info_results_date, '/results/pcawg_coverage_info.tsv')) %>%
  select(sample_id, cohort, mean_coverage) %>%
  mutate(coverage_label = case_when(
    mean_coverage < 49 ~ '38X',
    mean_coverage >= 49 ~ '60X'
  ))

# join metadata to the coverage table
coverage <- coverage %>%
  inner_join(metadata, by = c('sample_id', 'cohort')) %>%
  select(1:4, cancer_type_code)

### filled barplot of coverage proportion across cancer types
coverage_proportion_plot <- coverage %>%
  group_by(cancer_type_code, coverage_label) %>%
  summarise(n = n()) %>%
  drop_na() %>%
  ggplot(., aes(x = cancer_type_code, y = n, fill = coverage_label)) +
  geom_bar(
    color = 'black',
    position = 'fill', 
    stat = 'identity'
  ) +
  scale_y_continuous(expand = c(0,0.02)) +
  labs(
    title = '',
    x = '',
    y = 'Proportion of PCAWG samples',
    fill = ''
  ) +
  theme_classic() +
  theme_bigfont +
  theme(plot.title = element_text(size = 18),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# save as pdf
pdf(
  file = paste0(output_dir, '/', current_date, '_coverage_proportion.pdf'),
  width = 11,
  height = 5,
  useDingbats = FALSE,
  compress = FALSE
)
print(coverage_proportion_plot)
dev.off()

# ------------------ FOR COVERAGE GROUPS
# join coverage info to the variant_privacy table
coverage <- coverage %>%
  select(sample_id, coverage_label)
variant_privacy <- variant_privacy %>%
  left_join(coverage, by = 'sample_id')
# ------------------

# choose relevant columns
metadata <- metadata %>%
  select(sample_id, cancer_type, cancer_type_code)

# join the metadata to the variant privacy table
variant_privacy <- variant_privacy %>%
  left_join(metadata, by = 'sample_id') %>%
  # drop samples that have no cancer_type label
  drop_na(cancer_type) %>%
  # exclude MNVs bc these are weird
  filter(mut_type != 'mnv')
  
# summarise amount of private and shared variants
variant_privacy_summed <- variant_privacy %>%
  group_by(mut_type, privacy) %>%
  summarise(n = sum(n), .groups = 'drop')

### stacked histogram across samples
variant_privacy_summed_plot <- variant_privacy_summed %>%
  ### optional coverage filter
  filter(!privacy %in% c('Hartwig_only_NA', 'shared_NA', 'Hartwig_only_90X', 'shared_90X')) %>%
  ggplot(., aes(x = mut_type, y = n, fill = privacy)) +
  geom_bar(
    color = 'black',
    position = 'fill', 
    stat = 'identity'
  ) +
  scale_y_continuous(expand = c(0,0.02)) +
  labs(
    title = 'Amount of private and shared variants called by each pipeline overall',
    x = '',
    y = 'Rel. contrib. to total amount',
    fill = ''
  ) +
  theme_classic() +
  theme_bigfont +
  theme(plot.title = element_text(size = 18))

# save as PDF
pdf(
  file = paste0(output_dir, '/', current_date, '_manuscript_violin_plot.pdf'),
  width = 10,
  height = 8,
  useDingbats = FALSE,
  compress = FALSE
)
print(variant_privacy_summed_plot)
dev.off()

# manipulate privacy column to segment it more finely by coverage_label
variant_privacy <- variant_privacy %>%
  # split shared and Hartwig only mutations into clonal subclonal
  mutate(privacy = if_else(privacy == 'shared', paste0(privacy, '_', coverage_label), privacy)) %>%
  mutate(privacy = if_else(privacy == 'Hartwig_only', paste0(privacy, '_', coverage_label), privacy)) %>%
  # drop the 'shared_none' privacy groups, because these are a technical artefact from the clonality method
  filter(privacy != 'shared_none')
  
# summarise amount of private and shared variants per cancer type
variant_privacy_summed_percan <- variant_privacy %>%
  group_by(cancer_type_code, mut_type, privacy) %>%
  summarise(n = sum(n), .groups = 'drop')

# order of colors
privacy_order <- c('Hartwig_only_60X', 'Hartwig_only_30X', 'PCAWG_only', 'shared_60X', 'shared_30X')

### stacked histogram across samples
variant_privacy_summed_percan %>%
  ### optional filter for clonality split
  # filter(privacy != 'shared_none') %>%
  ### optional coverage filter
  filter(!privacy %in% c('Hartwig_only_NA', 'shared_NA', 'Hartwig_only_90X', 'shared_90X')) %>%
  ### optional coverage order for comparability
  # mutate(privacy = factor(privacy, levels = privacy_order)) %>%
  ggplot(., aes(x = mut_type, y = n, fill = privacy)) +
  geom_bar(
    color = 'black',
    position = 'fill', 
    stat = 'identity'
  ) +
  scale_y_continuous(expand = c(0,0.02)) +
  scale_fill_manual(values = c('#F8766D', '#D674FD', '#7CAE00', '#00C0BE', '#00BCF4', '#22A3FF')) +
  facet_wrap(~ cancer_type_code, nrow = 2) +
  labs(
    title = 'Amount of private and shared variants called by each pipeline per cancer type',
    subtitle = 'Excluding MSI samples',
    x = '',
    y = 'Rel. contrib to total amount',
    fill = ''
  ) +
  theme_classic() +
  theme_bigfont +
  theme(
    strip.background = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank()
  )