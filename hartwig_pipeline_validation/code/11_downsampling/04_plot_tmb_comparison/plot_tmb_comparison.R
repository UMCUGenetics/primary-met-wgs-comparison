### BOXPLOTS: Plot the TMB of downsampled vs normal coverage Hartwig samples
# author: remy (sascha)
# date: 04/04/2022
# last updated: 03/05/2022

### Description
# This script plots TMB comparison of postprocessed (non-)SAGE filtered VCFs that were downsampled in silico to 38X and 60X coverage against
# original 100X coverage VCFs as a scatterplot and draws a regression line through the points.

### Input
# The combined context matrix you generated in step 11.3. (downsampling suppl. table)

### Output
# Supplementary Figure 2 from supplementary note 1.

### Usage
# Just run it

# global options
options(stringsAsFactors = F, scipen = 999)

# libs
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
context_results_date <- '2022_11_11' # date name of your merged SMNV/SV context matrix results you generate in step 11.3
plot_60X <- TRUE # TRUE: plot 60X downsampled vs 109X, FALSE: plot 38X downsampled vs 109X

# -------------------- TMB per mut type

# import SMNV TMB of downsampled Hartwig tumors (parsed by Luan)
gz_file <- gzfile(paste0(base_dir, '/results/11_downsampling/02_extract_mutational_signatures_from_vcfs/', context_results_date, 
                         '/results/matrices/smnv_contexts.txt.gz'))
smnv_tmb <- read.table(file = gz_file,
                       row.names = 1, # first col are the rownames
                       check.names = FALSE # avoid changes to colnames with special characters
                       )
smnv_tmb <- smnv_tmb %>%
  rownames_to_column(., var = 'sample_id') %>%
  # add the coverage label
  mutate(coverage_label = case_when(str_detect(sample_id, pattern = regex('AT?$')) ~ 'cov_60X',
                                    str_detect(sample_id, pattern = regex('BT?$')) ~ 'cov_30X',
                                    str_detect(sample_id, pattern = regex('CT?$')) ~ 'cov_38X',
                                    TRUE ~ 'cov_100X'), .after = sample_id) %>%
  mutate(sample_id_2 = if_else(str_detect(sample_id, pattern = regex('[ABC]T?$')), NA_character_, sample_id), .after = sample_id) %>%
  fill(., sample_id_2, .direction = 'up')

# impoort SV TMB of downsampled Hartwig tumors
sv_tmb_init <- read_tsv(file = paste0(base_dir, '/results/11_downsampling/02_extract_mutational_signatures_from_vcfs/', context_results_date, 
                                      '/results/matrices/sv_contexts.txt'))

#########################################################################################
### SNV #################################################################################

# prepare SNV data
snv_tmb <- smnv_tmb %>%
  select(sample_id_2, coverage_label, contains('['))

# calc snv tmb
snv_tmb <- snv_tmb %>%
  rowwise() %>%
  summarise(
    sample_id = sample_id_2,
    coverage_label = coverage_label,
    tmb = sum(c_across(where(is.numeric))))

# define ylim for SNV plot
snv_y_limit <- 6

# plot as boxplot with median fold change as label
if (plot_60X) {
  
  snv_tmb_plot <- snv_tmb %>%
    filter(coverage_label %in% c('cov_60X', 'cov_100X'))
  
  snv_p1 <- ggplot(snv_tmb_plot, aes(x = coverage_label, y = log10(tmb + 1))) +
    geom_jitter(color = '#C0C0C0') +
    geom_boxplot(alpha = 0.8) +
    geom_label(
      data = . %>% group_by(coverage_label) %>% summarise(median_tmb = median(tmb), .groups = 'drop'),
      aes(y = log10(median_tmb), label = round(median_tmb, 0))
    ) +
    geom_label(
      data = . %>% group_by(coverage_label) %>%
        summarise(median_tmb = median(tmb), .groups = 'drop') %>%
        arrange(coverage_label) %>%
        mutate(fold_change = paste0(round((median_tmb - lead(median_tmb)) / lead(median_tmb), 2), 'x')) %>%
        filter(fold_change != 'NAx'),
      aes(y = log10(median_tmb), label = fold_change),nudge_x = 0.5
    ) +
    scale_x_discrete(labels = c(`cov_100X` = '109X', `cov_60X` = '60X')) +
    scale_y_continuous(
      breaks = seq(0,snv_y_limit, by = 1),
      labels = 10^seq(0,snv_y_limit, by = 1), 
      limits = c(0, snv_y_limit)) +
    labs(title = 'SNV',
         x = '',
         y = 'Log10(Mut. load + 1)') +
    theme_bw() +
    theme_bigfont +
    theme(
      axis.text = element_text(color = '#000000', family = 'Helvetica'),
      axis.line = element_line(color = '#000000'),
      plot.margin = margin(t= 5, l = 5, b = 5, r = 15))
  
} else {
  
  snv_tmb_plot <- snv_tmb %>%
    filter(coverage_label %in% c('cov_38X', 'cov_100X'))
  
  snv_p1 <- ggplot(snv_tmb_plot, aes(x = coverage_label, y = log10(tmb + 1))) +
    geom_jitter(color = '#C0C0C0') +
    geom_boxplot(alpha = 0.8) +
    geom_label(
      data = . %>% group_by(coverage_label) %>% summarise(median_tmb = median(tmb), .groups = 'drop'),
      aes(y = log10(median_tmb), label = round(median_tmb, 0))
    ) +
    geom_label(
      data = . %>% group_by(coverage_label) %>%
        summarise(median_tmb = median(tmb), .groups = 'drop') %>%
        arrange(coverage_label) %>%
        mutate(fold_change = paste0(round((median_tmb - lead(median_tmb)) / lead(median_tmb), 2), 'x')) %>%
        filter(fold_change != 'NAx'),
      aes(y = log10(median_tmb), label = fold_change),nudge_x = 0.5
    ) +
    scale_x_discrete(labels = c(`cov_100X` = '109X', `cov_38X` = '38X')) +
    scale_y_continuous(
      breaks = seq(0,snv_y_limit, by = 1),
      labels = 10^seq(0,snv_y_limit, by = 1), 
      limits = c(0, snv_y_limit)) +
    labs(title = 'SNV',
         x = '',
         y = 'Log10(Mut. load + 1)') +
    theme_bw() +
    theme_bigfont +
    theme(
      axis.text = element_text(color = '#000000', family = 'Helvetica'),
      axis.line = element_line(color = '#000000'),
      plot.margin = margin(t= 5, l = 5, b = 5, r = 15))
}

#########################################################################################
### DBS #################################################################################

# prepare DBS data
dbs_tmb <- smnv_tmb %>%
  select(sample_id_2, coverage_label, matches('^[A-Z]{2}\\>'))

# calc dbs tmb
dbs_tmb <- dbs_tmb %>%
  rowwise() %>%
  summarise(
    sample_id = sample_id_2,
    coverage_label = coverage_label,
    tmb = sum(c_across(where(is.numeric))))

# define DBS ylim for plot
dbs_y_limit <- 4

# plot as boxplot with median fold change as label
if (plot_60X) {
  
  dbs_tmb_plot <- dbs_tmb %>%
    filter(coverage_label %in% c('cov_60X', 'cov_100X'))
  
  dbs_p1 <- ggplot(dbs_tmb_plot, aes(x = coverage_label, y = log10(tmb + 1))) +
    geom_jitter(color = '#C0C0C0') +
    geom_boxplot(alpha = 0.8) +
    geom_label(
      data = . %>% group_by(coverage_label) %>% summarise(median_tmb = median(tmb), .groups = 'drop'),
      aes(y = log10(median_tmb), label = round(median_tmb, 0))
    ) +
    geom_label(
      data = . %>% group_by(coverage_label) %>%
        summarise(median_tmb = median(tmb), .groups = 'drop') %>%
        arrange(coverage_label) %>%
        mutate(fold_change = paste0(round((median_tmb - lead(median_tmb)) / lead(median_tmb), 2), 'x')) %>%
        filter(fold_change != 'NAx'),
      aes(y = log10(median_tmb), label = fold_change),nudge_x = 0.5
    ) +
    scale_x_discrete(labels = c(`cov_100X` = '109X', `cov_60X` = '60X')) +
    scale_y_continuous(
      breaks = seq(0,dbs_y_limit, by = 1),
      labels = 10^seq(0,dbs_y_limit, by = 1), 
      limits = c(0, dbs_y_limit)) +
    labs(title = 'DBS',
         x = '',
         y = 'Log10(Mut. load + 1)') +
    theme_bw() +
    theme_bigfont +
    theme(
      axis.text = element_text(color = '#000000', family = 'Helvetica'),
      axis.line = element_line(color = '#000000'),
      plot.margin = margin(t= 5, l = 5, b = 5, r = 15))
  
} else {
  
  dbs_tmb_plot <- dbs_tmb %>%
    filter(coverage_label %in% c('cov_38X', 'cov_100X'))
  
  dbs_p1 <- ggplot(dbs_tmb_plot, aes(x = coverage_label, y = log10(tmb + 1))) +
    geom_jitter(color = '#C0C0C0') +
    geom_boxplot(alpha = 0.8) +
    geom_label(
      data = . %>% group_by(coverage_label) %>% summarise(median_tmb = median(tmb), .groups = 'drop'),
      aes(y = log10(median_tmb), label = round(median_tmb, 0))
    ) +
    geom_label(
      data = . %>% group_by(coverage_label) %>%
        summarise(median_tmb = median(tmb), .groups = 'drop') %>%
        arrange(coverage_label) %>%
        mutate(fold_change = paste0(round((median_tmb - lead(median_tmb)) / lead(median_tmb), 2), 'x')) %>%
        filter(fold_change != 'NAx'),
      aes(y = log10(median_tmb), label = fold_change),nudge_x = 0.5
    ) +
    scale_x_discrete(labels = c(`cov_100X` = '109X', `cov_38X` = '38X')) +
    scale_y_continuous(
      breaks = seq(0,dbs_y_limit, by = 1),
      labels = 10^seq(0,dbs_y_limit, by = 1), 
      limits = c(0, dbs_y_limit)) +
    labs(title = 'DBS',
         x = '',
         y = 'Log10(Mut. load + 1)') +
    theme_bw() +
    theme_bigfont +
    theme(
      axis.text = element_text(color = '#000000', family = 'Helvetica'),
      axis.line = element_line(color = '#000000'),
      plot.margin = margin(t= 5, l = 5, b = 5, r = 15))
}

#########################################################################################
### Indel ###############################################################################

# prepare INDEL data
indel_tmb <- smnv_tmb %>%
  select(sample_id_2, coverage_label, matches('^del|^ins'))

# calc dbs tmb
indel_tmb <- indel_tmb %>%
  rowwise() %>%
  summarise(
    sample_id = sample_id_2,
    coverage_label = coverage_label,
    tmb = sum(c_across(where(is.numeric))))

# # filter out the 3 extreme outliers
# indel_tmb <- indel_tmb %>%
#   filter(tmb < 100000)

# define INDEL ylim for plot
indel_y_limit <- 7

# plot as boxplot with median fold change as label
if (plot_60X) {
  
  indel_tmb_plot <- indel_tmb %>%
    filter(coverage_label %in% c('cov_60X', 'cov_100X'))
  
  indel_p1 <- ggplot(indel_tmb_plot, aes(x = coverage_label, y = log10(tmb + 1))) +
    geom_jitter(color = '#C0C0C0') +
    geom_boxplot(alpha = 0.8) +
    geom_label(
      data = . %>% group_by(coverage_label) %>% summarise(median_tmb = median(tmb), .groups = 'drop'),
      aes(y = log10(median_tmb), label = round(median_tmb, 0))
    ) +
    geom_label(
      data = . %>% group_by(coverage_label) %>%
        summarise(median_tmb = median(tmb), .groups = 'drop') %>%
        arrange(coverage_label) %>%
        mutate(fold_change = paste0(round((median_tmb - lead(median_tmb)) / lead(median_tmb), 2), 'x')) %>%
        filter(fold_change != 'NAx'),
      aes(y = log10(median_tmb), label = fold_change),nudge_x = 0.5
    ) +
    scale_x_discrete(labels = c(`cov_100X` = '109X', `cov_60X` = '60X')) +
    scale_y_continuous(
      breaks = seq(0,indel_y_limit, by = 1),
      labels = 10^seq(0,indel_y_limit, by = 1), 
      limits = c(0, indel_y_limit)) +
    labs(title = 'INDEL',
         x = '',
         y = 'Log10(Mut. load + 1)') +
    theme_bw() +
    theme_bigfont +
    theme(
      axis.text = element_text(color = '#000000', family = 'Helvetica'),
      axis.line = element_line(color = '#000000'),
      plot.margin = margin(t= 5, l = 5, b = 5, r = 15))
  
} else {
  
  indel_tmb_plot <- indel_tmb %>%
    filter(coverage_label %in% c('cov_38X', 'cov_100X'))
  
  indel_p1 <- ggplot(indel_tmb_plot, aes(x = coverage_label, y = log10(tmb + 1))) +
    geom_jitter(color = '#C0C0C0') +
    geom_boxplot(alpha = 0.8) +
    geom_label(
      data = . %>% group_by(coverage_label) %>% summarise(median_tmb = median(tmb), .groups = 'drop'),
      aes(y = log10(median_tmb), label = round(median_tmb, 0))
    ) +
    geom_label(
      data = . %>% group_by(coverage_label) %>%
        summarise(median_tmb = median(tmb), .groups = 'drop') %>%
        arrange(coverage_label) %>%
        mutate(fold_change = paste0(round((median_tmb - lead(median_tmb)) / lead(median_tmb), 2), 'x')) %>%
        filter(fold_change != 'NAx'),
      aes(y = log10(median_tmb), label = fold_change),nudge_x = 0.5
    ) +
    scale_x_discrete(labels = c(`cov_100X` = '109X', `cov_38X` = '38X')) +
    scale_y_continuous(
      breaks = seq(0,indel_y_limit, by = 1),
      labels = 10^seq(0,indel_y_limit, by = 1), 
      limits = c(0, indel_y_limit)) +
    labs(title = 'INDEL',
         x = '',
         y = 'Log10(Mut. load + 1)') +
    theme_bw() +
    theme_bigfont +
    theme(
      axis.text = element_text(color = '#000000', family = 'Helvetica'),
      axis.line = element_line(color = '#000000'),
      plot.margin = margin(t= 5, l = 5, b = 5, r = 15))
}

#########################################################################################
### SV ##################################################################################

# prepare SV data
sv_tmb <- sv_tmb_init %>%
  mutate(sample = str_replace(sample, pattern = regex('([ABC])T$'), replacement = '\\1')) %>%
  inner_join(smnv_tmb, by = c('sample' = 'sample_id')) %>%
  select(sample_id_2, coverage_label, SV_LOAD) %>%
  rename(tmb = 'SV_LOAD',
         sample_id = 'sample_id_2')

# define SV ylim for plot
sv_y_limit <- 4

# plot as boxplot with median fold change as label
if (plot_60X) {
  
  sv_tmb_plot <- sv_tmb %>%
    filter(coverage_label %in% c('cov_60X', 'cov_100X'))
  
  sv_p1 <- ggplot(sv_tmb_plot, aes(x = coverage_label, y = log10(tmb + 1))) +
    geom_jitter(color = '#C0C0C0') +
    geom_boxplot(alpha = 0.8) +
    geom_label(
      data = . %>% group_by(coverage_label) %>% summarise(median_tmb = median(tmb), .groups = 'drop'),
      aes(y = log10(median_tmb), label = round(median_tmb, 0))
    ) +
    geom_label(
      data = . %>% group_by(coverage_label) %>%
        summarise(median_tmb = median(tmb), .groups = 'drop') %>%
        arrange(coverage_label) %>%
        mutate(fold_change = paste0(round((median_tmb - lead(median_tmb)) / lead(median_tmb), 2), 'x')) %>%
        filter(fold_change != 'NAx'),
      aes(y = log10(median_tmb), label = fold_change),nudge_x = 0.5
    ) +
    scale_x_discrete(labels = c(`cov_100X` = '109X', `cov_60X` = '60X')) +
    scale_y_continuous(
      breaks = seq(0,sv_y_limit, by = 1),
      labels = 10^seq(0,sv_y_limit, by = 1), 
      limits = c(0, sv_y_limit)) +
    labs(title = 'SV',
         x = '',
         y = 'Log10(Mut. load + 1)') +
    theme_bw() +
    theme_bigfont +
    theme(
      axis.text = element_text(color = '#000000', family = 'Helvetica'),
      axis.line = element_line(color = '#000000'),
      plot.margin = margin(t= 5, l = 5, b = 5, r = 15))
  
} else {
  
  sv_tmb_plot <- sv_tmb %>%
    filter(coverage_label %in% c('cov_38X', 'cov_100X'))
  
  sv_p1 <- ggplot(sv_tmb_plot, aes(x = coverage_label, y = log10(tmb + 1))) +
    geom_jitter(color = '#C0C0C0') +
    geom_boxplot(alpha = 0.8) +
    geom_label(
      data = . %>% group_by(coverage_label) %>% summarise(median_tmb = median(tmb), .groups = 'drop'),
      aes(y = log10(median_tmb), label = round(median_tmb, 0))
    ) +
    geom_label(
      data = . %>% group_by(coverage_label) %>%
        summarise(median_tmb = median(tmb), .groups = 'drop') %>%
        arrange(coverage_label) %>%
        mutate(fold_change = paste0(round((median_tmb - lead(median_tmb)) / lead(median_tmb), 2), 'x')) %>%
        filter(fold_change != 'NAx'),
      aes(y = log10(median_tmb), label = fold_change),nudge_x = 0.5
    ) +
    scale_x_discrete(labels = c(`cov_100X` = '109X', `cov_38X` = '38X')) +
    scale_y_continuous(
      breaks = seq(0,sv_y_limit, by = 1),
      labels = 10^seq(0,sv_y_limit, by = 1), 
      limits = c(0, sv_y_limit)) +
    labs(title = 'SV',
         x = '',
         y = 'Log10(Mut. load + 1)') +
    theme_bw() +
    theme_bigfont +
    theme(
      axis.text = element_text(color = '#000000', family = 'Helvetica'),
      axis.line = element_line(color = '#000000'),
      plot.margin = margin(t= 5, l = 5, b = 5, r = 15))
}

# stitch plots together
boxplots <- ggarrange(snv_p1,
                      dbs_p1,
                      indel_p1,
                      sv_p1,
                      labels = NULL,
                      ncol = 4, nrow = 1,
                      align = 'hv',
                      font.label = list(size = 8, color = 'black', face = 'bold', family = NULL, position = 'top'))

# create output dir
output_dir <- paste0(output_dir, '/11_downsampling/04_plot_tmb_comparison/', current_date, '/results/plots')
dir.create(path = paste0(output_dir), recursive = TRUE)

# save as pdf
pdf(
  file = paste0(output_dir, '/', current_date, 
                if (plot_60X) { '_60X' } else { '_38X' }, '_smnv_slope_plot.pdf'),
  width = 20, 
  height = 5,
  useDingbats = FALSE,
  compress = FALSE
)
print(boxplots)
dev.off()