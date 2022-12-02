############### Compare purity values per cancer types between primary and metastatic ###############
# author: remy (sascha)
# date: 02/11/2022
# last updated: 02/11/2022

### Description
# This script plots the 

### Input
# metadata.tsv
# the parsed VAF valu from variants of samples with 100% purity (step 12.2.)

### Output
# VAF distribution of 100% pure tumor samples

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

# define variables
vaf_results_date <- '2022_11_11' # date name of your parsed vaf values from step 12.2
subcl_thresh <- 0.1 # the threshold at which to consider a variant subclonal (the higher the more likely for it to be subclonal)
filter_cohort <- F # TRUE: Filter for samples in either PCAWG or Hartwig cohort (specified below)
cohort_filter <- 'PCAWG' # 'PCAWG' or 'Hartwig'

# -------------------------- METADATA

metadata <- read_tsv(
  file = paste0(base_dir, '/data/processed/metadata/metadata.tsv')
) %>%
  select(sample_id, cohort, cancer_type, tumor_purity, whole_genome_duplication, msi_status) %>%
  filter(near(tumor_purity, 1))

# -------------------------- VAF

vaf_input <- list.files(paste0(base_dir, '/results/12_purity/02_parse_vaf/', vaf_results_date, '/results/'),
                        full.names = T)
if (!exists('vaf')) {
  vaf <- lapply(vaf_input, read_tsv)
  vaf <- do.call(rbind, vaf)
  vaf <- vaf %>%
    mutate(
      vaf = round(vaf, 3)
    ) %>%
    replace_na(., replace = list(subclonal_likelihood = 0))
}

vaf_filt <- vaf %>%
  filter(
    !is.na(gene),
    subclonal_likelihood < subcl_thresh
    ) %>%
  group_by(sample_id, gene) %>%
  slice_max(order_by = vaf) %>%
  ungroup() %>%
  distinct(sample_id, gene, vaf)

purity100_vaf_histogram <- vaf_filt %>%
  inner_join(metadata, by = c('sample_id')) %>%
  { if (filter_cohort) filter(., cohort == cohort_filter) else . } %>%
  ggplot(., aes(x = vaf)) +
  geom_histogram(binwidth = 0.01) +
  labs(
    title = paste0('VAF dist. of point mutations in 100% pure tumor samples'),
    subtitle = if (filter_cohort) { paste0('Subclonal likelihood < ', subcl_thresh, ' | ', cohort_filter) }
    else { paste0('Subclonal likelihood < ', subcl_thresh) },
    x = 'VAF',
    y = 'Count'
  ) +
  theme_bw() +
  theme_bigfont_barplot +
  theme(
    plot.margin = margin(0.5,0.5,0.5,0.5, unit = 'cm')
  )

# create output dir
output_dir <- paste0(output_dir, '/12_purity/03_plot_pure_samples/', current_date, '/results/plots')
dir.create(path = paste0(output_dir), recursive = TRUE)

# save as PDF
pdf(
  file = paste0(output_dir, '/', current_date, '_purity100_vaf',
                if (filter_cohort) { paste0('_', cohort_filter) }
                  else { '.pdf' }),
  width = 7.5,
  height = 5,
  useDingbats = FALSE,
  compress = FALSE
)
print(purity100_vaf_histogram)
dev.off()
