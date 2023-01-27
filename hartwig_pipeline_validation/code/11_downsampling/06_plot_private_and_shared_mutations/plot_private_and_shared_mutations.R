############### Plot variant call overlaps in downsampled hartwig samples ###############
# author: remy (sascha)
# date: 24/01/2023
# last updated: 24/01/2023

### Description
# This script calculates the percentage of variants that are recapitulated in downsampled hartwig samples (38X, 60X)
# versus the original coverage samples (109X). 
# This rate is stratified by purity adjusted VAF (bins of 0.1 VAF) of the variant in the original sample (109X).
# Further, the rate are shown separately for every mutation type (SNV, DBS, INDEL).
# A sample is only shown in a specific bin, if the total number of mutations (recap and non-recap) is at least 10.

### Input
# list of extracted variants tables of downsampled and normal hartwig samples.

### Output
# Supplementary Figure 2C in supplementary note 1.

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

## ggpelot custom theme
theme_bigfont <- theme(plot.title = element_text(size = 21),
                       plot.subtitle = element_text(size = 16),
                       axis.text.x = element_text(size = 15),
                       axis.text.y = element_text(size = 15), 
                       axis.title = element_text(size = 18),
                       legend.title = element_text(size = 16),
                       legend.text = element_text(size = 14),
                       strip.text = element_text(size = 15))

# define variables
filtered_variants_results_date <- '2023_01_25' # date name of the folder where you stored your extracted VCF variants of the downsampled samples.
use_sage_filtered_vcfs <- TRUE # TRUE: use SAGE VCF data, FALSE: non SAGE filtered VCF data

# -------------------- METADATA DOWNSAMPLED

downsampled_samples <- read_tsv(paste0(base_dir, '/data/processed/downsampling/cov_downsampling_samples.tsv')) %>%
  rename_with(tolower) %>%
  mutate(cov_30x_id = paste0(sample_id, 'B'),
         cov_38x_id = paste0(sample_id, 'C'),
         cov_60x_id = paste0(sample_id, 'A'),
         cov_109x_id = sample_id) %>%
  select(-sample_id) %>%
  pivot_longer(., cols = c(cov_30x_id:cov_109x_id), names_to = 'coverage', values_to = 'sample_id')

# -------------------- FILTERED VARIANTS

filtered_variants <- list.files(
  path = paste0(base_dir, '/results/11_downsampling/06_plot_private_and_shared_mutations/', 
                filtered_variants_results_date, if (use_sage_filtered_vcfs) { '/sage_filtered' } else { '/original' }, '/results'),
  full.names = TRUE
)

# remove 30X
filtered_variants <- filtered_variants[!str_detect(filtered_variants, pattern = regex('CPCT[0-9]+B'))]

# read all the table into R and transform to df
filtered_variants <- map_df(filtered_variants, function(x) {
  read_tsv(
    file = x
  )
})

# drop NA VAF, or unreasonably high VAF
filtered_variants <- filtered_variants %>%
  filter(!is.na(Purple_af), Purple_af < 2)

# annotate mutation type
filtered_variants$MutType <- mutSigExtractor::detSmnvType(filtered_variants$Ref, filtered_variants$Alt)

# create purple VAF bins, add coverage annotation
filtered_variants <- filtered_variants %>%
  mutate(
    Coverage = case_when(
      str_detect(SampleId, pattern = regex('CPCT[0-9]+C')) ~ 'cov_38X',
      str_detect(SampleId, pattern = regex('CPCT[0-9]+A')) ~ 'cov_60X',
      TRUE ~ 'cov_109X'
    )
  )

# add coverage annotation, nest data based on strata
filtered_variantstt <- filtered_variants %>%
  group_by(Coverage, SampleId) %>%
  nest()

# add sample_id to variants table, make df wide
filtered_variantstt <- downsampled_samples %>%
  inner_join(filtered_variantstt, by = c('sample_id' = 'SampleId')) %>%
  select(-sample_id) %>%
  pivot_wider(., names_from = c('Coverage'), values_from = 'data')

# determine the number of variants in 109X sample, that are recapitulated in the 
filtered_variantstt2 <- filtered_variantstt %>%
  mutate(
    Overlap_38_109 = map2(cov_109X, cov_38X, function(df1, df2) semi_join(df1, df2, by = c('Chrom', 'Pos', 'Ref', 'Alt'))),
    Overlap_60_109 = map2(cov_109X, cov_60X, function(df1, df2) semi_join(df1, df2, by = c('Chrom', 'Pos', 'Ref', 'Alt'))),
  ) %>%
  mutate(
    Overlap_38_109 = map2(Overlap_38_109, cov_109X, function(df1, df2) {
      df1 %>% 
        mutate(Recap = TRUE) %>% 
        right_join(., df2, by = colnames(df2)) %>%
        replace_na(., replace = list('Recap' = FALSE)) %>%
        mutate(
          Purple_af_bin = cut(Purple_af, breaks = c(seq(0,0.9, by = 0.1), 2), labels = str_c('vaf_', seq(0,0.9, by = 0.1), "_", seq(0.1,1, by = 0.1)))
        ) %>%
        group_by(MutType, Purple_af_bin, Recap) %>%
        summarise(Count = n(), .groups = 'drop') %>%
        complete(., MutType, Purple_af_bin, Recap, fill = list(Count = 0)) %>%
        group_by(MutType, Purple_af_bin) %>%
        mutate(Sum = sum(Count)) %>%
        ungroup() %>%
        filter(Recap) %>%
        mutate(Percent = (Count / Sum) * 100)
        }),
    Overlap_60_109 = map2(Overlap_60_109, cov_109X, function(df1, df2) {
      df1 %>% 
        mutate(Recap = TRUE) %>% 
        right_join(., df2, by = colnames(df2)) %>%
        replace_na(., replace = list('Recap' = FALSE)) %>%
        mutate(
          Purple_af_bin = cut(Purple_af, breaks = c(seq(0,0.9, by = 0.1), 2), labels = str_c('vaf_', seq(0,0.9, by = 0.1), "_", seq(0.1,1, by = 0.1)))
        ) %>%
        group_by(MutType, Purple_af_bin, Recap) %>%
        summarise(Count = n(), .groups = 'drop') %>%
        complete(., MutType, Purple_af_bin, Recap, fill = list(Count = 0)) %>%
        group_by(MutType, Purple_af_bin) %>%
        mutate(Sum = sum(Count)) %>%
        ungroup() %>%
        filter(Recap) %>%
        mutate(Percent = (Count / Sum) * 100)
    }) 
  )

# unnest for plotting
filtered_variants_plot <- filtered_variantstt2 %>% 
  select(sample_id, starts_with('Overlap')) %>%
  pivot_longer(., c('Overlap_38_109', 'Overlap_60_109'), names_to = 'Group', values_to = 'unnest_col') %>%
  unnest(., cols = unnest_col) %>%
  filter(Sum > 10) %>%
  mutate(MutType = str_to_upper(MutType)) %>%
  mutate(Group = recode(Group, Overlap_38_109 = '38X vs. 109X', Overlap_60_109 = '60X vs. 109X'))

### plots
# color palette
group_colors <- c('#1874CD', '#FFD700')

# histogram
hist <- filtered_variants_plot %>%
  group_by(MutType, Group, Purple_af_bin) %>%
  summarise(Sum = sum(Count), .groups = 'drop')

# snv
snv_p1 <- hist %>%
  filter(MutType == 'SNV') %>%
  ggplot(., aes(
    x = Group,
    y = Sum,
    fill = Group
  )) +
  geom_col(color = '#000000') +
  facet_grid(MutType ~ Purple_af_bin) +
  scale_y_continuous(
    expand = c(0,0),
    limits = c(0, 8e5),
    breaks = seq(0, 7e5, by = 1e5)
  ) +
  scale_fill_manual(values = group_colors) +
  labs(
    x = '',
    y = 'Number of mutations'
  ) +
  theme_bigfont +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = 'none',
    plot.margin = margin(0,2,0,0, unit = 'cm')
  )

# dbs
dbs_p1 <- hist %>%
  filter(MutType == 'DBS') %>%
  ggplot(., aes(
    x = Group,
    y = Sum,
    fill = Group
  )) +
  geom_col(color = '#000000') +
  facet_grid(MutType ~ Purple_af_bin) +
  scale_y_continuous(
    expand = c(0,0),
    limits = c(0, 8500),
    breaks = seq(0, 8000, by = 1000)
  ) +
  scale_fill_manual(values = group_colors) +
  labs(
    x = '',
    y = 'Number of mutations'
  ) +
  theme_bigfont +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_blank(),
    legend.position = 'none',
    plot.margin = margin(0,2,0,0, unit = 'cm')
  )

# indel
indel_p1 <- hist %>%
  filter(MutType == 'INDEL') %>%
  ggplot(., aes(
    x = Group,
    y = Sum,
    fill = Group
  )) +
  geom_col(color = '#000000') +
  facet_grid(MutType ~ Purple_af_bin) +
  scale_y_continuous(
    expand = c(0,0),
    limits = c(0, 4.1 * 1e5),
    breaks = seq(0, 4e5, by = 1e5)
  ) +
  scale_fill_manual(values = group_colors) +
  labs(
    x = '',
    y = 'Number of mutations'
  ) +
  theme_bigfont +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_blank(),
    legend.position = 'none',
    plot.margin = margin(0,2,0,0, unit = 'cm')
  )

# boxplot
snv_p2 <- filtered_variants_plot %>%
  filter(MutType == 'SNV') %>%
  ggplot(., aes(
    x = Group,
    y = Percent,
    fill = Group
  )) +
  geom_boxplot() +
  geom_jitter() +
  facet_grid(MutType ~ Purple_af_bin) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20)
  ) +
  scale_fill_manual(values = group_colors) +
  labs(
    x = '',
    y = 'Percent overlap'
  ) +
  theme_bigfont +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_blank(),
    legend.position = 'none',
    plot.margin = margin(0,2,0.1,0, unit = 'cm')
  )

dbs_p2 <- filtered_variants_plot %>%
  filter(MutType == 'DBS') %>%
  ggplot(., aes(
    x = Group,
    y = Percent,
    fill = Group
  )) +
  geom_boxplot() +
  geom_jitter() +
  facet_grid(MutType ~ Purple_af_bin) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20)
  ) +
  scale_fill_manual(values = group_colors) +
  labs(
    x = '',
    y = 'Percent overlap'
  ) +
  theme_bigfont +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_blank(),
    legend.position = 'none',
    plot.margin = margin(0,2,0.1,0, unit = 'cm')
  )

indel_p2 <- filtered_variants_plot %>%
  filter(MutType == 'INDEL') %>%
  ggplot(., aes(
    x = Group,
    y = Percent,
    fill = Group
  )) +
  geom_boxplot() +
  geom_jitter() +
  facet_grid(MutType ~ Purple_af_bin) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20)
  ) +
  scale_fill_manual(values = group_colors) +
  labs(
    x = '',
    y = 'Percent overlap'
  ) +
  theme_bigfont +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    strip.text.x = element_blank(),
    legend.position = 'none',
    plot.margin = margin(-0.1,2,0.1,0, unit = 'cm')
  )

# combine plots
combined_plots <- snv_p1 + snv_p2 + dbs_p1 + dbs_p2 + indel_p1 + indel_p2 + plot_layout(
  nrow = 6,
  ncol = 1,
)

# create output dir
output_dir <- paste0(output_dir, '/11_downsampling/06_plot_private_and_shared_mutations/', current_date, if (use_sage_filtered_vcfs) { '/sage_filtered' } else { '/non_sage_filtered' }, '/results/plots')
dir.create(path = paste0(output_dir), recursive = TRUE)

# save figure as pdf
pdf(
  file = paste0(output_dir, '/', current_date, '_downsampled_overlap_', if (use_sage_filtered_vcfs) { 'sage_filtered' } else { 'non_sage_filtered' }, '.pdf'),
  width = 18, 
  height = 22,
  useDingbats = FALSE,
  compress = FALSE
)
print(combined_plots)
dev.off()

### table: calculate median values
median_overlap <- filtered_variants_plot %>%
  filter(MutType != 'MNV') %>%
  group_by(Group, MutType, Purple_af_bin) %>%
  summarise(MedianOverlap = median(Percent), .groups = 'drop') %>%
  arrange(factor(MutType, levels = c('SNV', 'DBS', 'INDEL')), Purple_af_bin)

write_tsv(
  median_overlap,
  file = paste0(output_dir, '/', current_date, '_median_overlap_', if (use_sage_filtered_vcfs) { 'sage_filtered' } else { 'non_sage_filtered' }, '.tsv')
)