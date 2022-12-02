############### PCAWG vs. Hartwig #######################################################
############### Analyse the cosine similarity and total mutational burden ###############
############### for the SNV, DBS, Indel & SV profiles ###################################
############### Plots arranged for manuscript ###########################################
# author: remy (sascha)
# date: 19/10/2021
# last updated: 26/04/2022

### Description
# This script plots the total mutational burden ratio vs the cosine similarity for the
# 96 SNV mutational profiles, Indel profiles, DBS and SV (Both with and without small length cutoff) profiles.
# This script depends on the four scrits that each plot and save the output plots to a specified folder:

# slope_plot.R, ratio_plot.R, waterfall_plot.R

# The above scripts should be located in the same folder as this script for it to work.

### Input
# - The cosine similarity and mutationa burden tables you generated in step 4 for all the mutation types.
# - metadata.tsv file.

### Output
# Supplementary Figure 7 in supplementary note 1

### Usage
# Open this script in a new R session and just run it

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

# get today's date
current_date <- format(Sys.Date(), '%Y_%m_%d')

# create output dir
output_dir <- paste0(output_dir, '/5_plot_cosine_similarity_vs_mutational_burden/', current_date, '/results/plots')
dir.create(path = paste0(output_dir), recursive = TRUE)

## ggplot custom theme
theme_bigfont <- theme(plot.title = element_text(size = 21),
                       plot.subtitle = element_text(size = 16),
                       axis.text.x = element_text(size = 15),
                       axis.text.y = element_text(size = 15), 
                       axis.title = element_text(size = 18),
                       legend.title = element_text(size = 16),
                       legend.text = element_text(size = 14),
                       strip.text.x = element_text(size = 15))

# define variables
tmb_input_results_date <- '2022_02_16' # date name of your tumor mutational burden results from step 4
cos_input_results_date <- '2022_02_16' # date name of your cosine similarity results from step 4

# ------------------ TOTAL MUTATIONAL BURDEN

snv_total_mut_burden <- read_tsv(paste0(base_dir, '/results/4_calculate_total_mutational_burden/', tmb_input_results_date, '/results/snv_total_mutational_burden.tsv'))
dbs_total_mut_burden <- read_tsv(paste0(base_dir, '/results/4_calculate_total_mutational_burden/', tmb_input_results_date, '/results/dbs_total_mutational_burden.tsv'))
indel_total_mut_burden <- read_tsv(paste0(base_dir, '/results/4_calculate_total_mutational_burden/', tmb_input_results_date, '/results/indel_total_mutational_burden.tsv')) %>%
  mutate(indel_pcawg_total_mutational_burden = as.numeric(indel_pcawg_total_mutational_burden))
sv_total_mut_burden <- read_tsv(paste0(base_dir, '/results/4_calculate_total_mutational_burden/', tmb_input_results_date, '/results/sv_total_mutational_burden.tsv'))
sv_total_mut_burden_len_cutoff <- read_tsv(paste0(base_dir, '/results/4_calculate_total_mutational_burden/', tmb_input_results_date, '/results/sv_total_mutational_burden_len_cutoff.tsv')) %>%
  mutate(sv_pcawg_total_mutational_burden_len_cutoff = as.numeric(sv_pcawg_total_mutational_burden_len_cutoff))

# ------------------ COSINE SIMILARITY

snv_cosine_similarity <- read_tsv(paste0(base_dir, '/results/4_calculate_cosine_similarity/', cos_input_results_date, '/results/snv_cosine_similarity.tsv'),
                              col_names = TRUE) %>%
  mutate(snv_cosine_similarity = as.numeric(snv_cosine_similarity))
dbs_cosine_similarity <- read_tsv(paste0(base_dir, '/results/4_calculate_cosine_similarity/', cos_input_results_date, '/results/dbs_cosine_similarity.tsv'),
                                  col_names = TRUE) %>%
  mutate(dbs_cosine_similarity = as.numeric(dbs_cosine_similarity))
indel_cosine_similarity <- read_tsv(paste0(base_dir, '/results/4_calculate_cosine_similarity/', cos_input_results_date, '/results/indel_cosine_similarity.tsv'),
                                  col_names = TRUE) %>%
  mutate(indel_cosine_similarity = as.numeric(indel_cosine_similarity))
sv_cosine_similarity <- read_tsv(paste0(base_dir, '/results/4_calculate_cosine_similarity/', cos_input_results_date, '/results/sv_cosine_similarity.tsv'),
                                  col_names = TRUE) %>%
  mutate(sv_cosine_similarity = as.numeric(sv_cosine_similarity))
sv_cosine_similarity_len_cufoff <- read_tsv(paste0(base_dir, '/results/4_calculate_cosine_similarity/', cos_input_results_date, '/results/sv_cosine_similarity_len_cutoff.tsv'),
                                 col_names = TRUE) %>%
  mutate(sv_cosine_similarity_len_cutoff = as.numeric(sv_cosine_similarity_len_cutoff))

# ------------------ METADATA

metadata <- read_tsv(paste0(base_dir, '/data/processed/metadata/metadata.tsv'),
                     col_types = 'icccccccccccccccllllcccdiiildclccllllllliiiiiiiiiii')

# filter for microsatellite stable (MSS) PCAWG samples only
metadata <- metadata %>%
  filter(msi_status == 'MSS') %>%
  filter(cohort == 'PCAWG')

# ------------------ prepare data

# join the two tables together, calculate mutational burden ratio and add column
# that indicates whether this observation has a cosine similarity over 0.90
# filter the one tumor with a ratio over 2

# SNV
snv_cos_tot_mut_combined <- snv_cosine_similarity %>%
  inner_join(snv_total_mut_burden, by = 'sample_id') %>%
  mutate(snv_pcawg_total_mutational_burden = as.numeric(snv_pcawg_total_mutational_burden)) %>%
  mutate(mut_burden_ratio = snv_hartwig_total_mutational_burden / snv_pcawg_total_mutational_burden,
         over_90 = ifelse(snv_cosine_similarity > 0.9, TRUE, FALSE))

# join the analysis table to the metadata table to get the cancer type and hypermutation information
snv_cos_tot_mut_combined <- snv_cos_tot_mut_combined %>%
  inner_join(metadata, by = c('sample_id')) %>%
  select(1:6, cancer_type)

# DBS
dbs_cos_tot_mut_combined <- dbs_cosine_similarity %>%
  inner_join(dbs_total_mut_burden, by = 'sample_id') %>%
  mutate(dbs_pcawg_total_mutational_burden = as.numeric(dbs_pcawg_total_mutational_burden),
         mut_burden_ratio = dbs_hartwig_total_mutational_burden / dbs_pcawg_total_mutational_burden,
         over_90 = ifelse(dbs_cosine_similarity > 0.9, TRUE, FALSE))

dbs_cos_tot_mut_combined <- dbs_cos_tot_mut_combined %>%
  inner_join(metadata, by = c('sample_id')) %>%
  select(1:6, cancer_type)

# Indel
indel_cos_tot_mut_combined <- indel_cosine_similarity %>%
  inner_join(indel_total_mut_burden, by = 'sample_id') %>%
  mutate(indel_pcawg_total_mutational_burden = as.numeric(indel_pcawg_total_mutational_burden),
         mut_burden_ratio = indel_hartwig_total_mutational_burden / indel_pcawg_total_mutational_burden,
         over_90 = ifelse(indel_cosine_similarity > 0.9, TRUE, FALSE))

indel_cos_tot_mut_combined <- indel_cos_tot_mut_combined %>%
  inner_join(metadata, by = c('sample_id')) %>%
  select(1:6, cancer_type)

# SV
sv_cos_tot_mut_combined <- sv_cosine_similarity %>%
  inner_join(sv_cosine_similarity_len_cufoff, by = 'sample_id') %>%
  inner_join(sv_total_mut_burden, by = 'sample_id') %>%
  inner_join(sv_total_mut_burden_len_cutoff, by = 'sample_id') %>%
  mutate(sv_pcawg_total_mutational_burden = as.numeric(sv_pcawg_total_mutational_burden)) %>%
  mutate(mut_burden_ratio = sv_hartwig_total_mutational_burden / as.numeric(sv_pcawg_total_mutational_burden),
         mut_burden_ratio_len_cutoff = sv_hartwig_total_mutational_burden_len_cutoff / as.numeric(sv_pcawg_total_mutational_burden_len_cutoff),
         over_90 = ifelse(sv_cosine_similarity > 0.9, TRUE, FALSE),
         over_90_len_cutoff = ifelse(sv_cosine_similarity_len_cutoff > 0.90, TRUE, FALSE))

# join the analysis table to the metadata table to get the cancer type and hypermutation information
sv_cos_tot_mut_combined <- sv_cos_tot_mut_combined %>%
  inner_join(metadata, by = c('sample_id')) %>%
  select(1:11, cancer_type)

# define sv length cutoff you used (will show up in plot title)
sv_length_cutoff <- 500

############################## BOXPLOT ###############################################

source('boxplot.R')

############################## SLOPE PLOT ############################################

# source('slope_plot.R')

############################## RATIO PLOT ############################################

source('ratio_plot.R')

############################## WATERFALL PLOT ########################################

source('waterfall_plot.R')

# combine all plots into one final plot
waterfall_ratio_plots <- ggarrange(waterfall_plots,
                                   ratio_plots,
                                   ncol = 2, nrow = 1)

technical_comparison_plot <- ggarrange(boxplot,
                                       waterfall_ratio_plots,
                                       heights = c(1,4),
                                       ncol = 1, nrow = 2)

# save a PDF file of the plots
pdf(
  file = paste0(output_dir, '/', current_date, '_boxplot.pdf'),
  width = 20,
  height = 5,
  useDingbats = FALSE,
  compress = FALSE
)
boxplot
dev.off()

pdf(
  file = paste0(output_dir, '/', current_date, '_technical_plot.pdf'),
  width = 20,
  height = 25,
  useDingbats = FALSE,
  compress = FALSE
)
technical_comparison_plot
dev.off()