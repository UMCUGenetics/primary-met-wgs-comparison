### Plot sensitivity loss
# author: remy (sascha)
# date: 28/03/2022
# last updated: 21/04/2022

### Description
# This script performs stratified sampling of 25 semi-randomly selected Hartwig samples that were downsampled to 38X and 60X coverage. 
# The stratified sampling is done at different probabilities (0, 0.2, 0.4, 0.6, 0.8, 1) for 38X samples to achieve a realistic proportion of 38X samples similar to proportions there are in the PCAWG dataset
# which are calculated thereafter. After that a mean sensitivity loss is calculated by subtracting the downsampled TMB from the 100X original TMB and then averaging this loss across the 25 samples for each probability. 
# This process is replicated 100 times to get a realistic estimate of the mean sensitivity loss (with 95% confidence interval) for each 38X proportion.

### Input
# metadata.tsv

### Output
# Supplementary Figure 3b in supplementary note 1.

### Usage
# Just run it

# global options
options(stringsAsFactors = F)

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
                       strip.text.x = element_text(size = 15))

# source color codes and order for cancer types
source(paste0(base_dir, '/code/r_objects/color_codes.R'))
source(paste0(base_dir, '/code/r_objects/cancer_type_order.R'))

# define variables
coverage_results_date <- '2022_11_11'# date name of your PCAWG coverage table you generate in step 7
context_results_date <- '2022_11_11' # date name of your merged SMNV/SV context matrix results you generate in step 11.3

# -------------------- METADATA

metadata <- read_tsv(file = paste0(base_dir, '/data/processed/metadata/metadata.tsv')) %>%
  filter(cohort == 'PCAWG') %>%
  select(sample_id, cancer_type, cancer_type_code)

# -------------------- PCAWG COVERAGE INFO

pcawg_coverage <- read_tsv(file = paste0(base_dir, '/results/07_create_coverage_tables/', coverage_results_date, '/results/pcawg_coverage_info.tsv')) %>%
  # change coverage label to split at 49X: '38X' = <= 49X, '60X' + > 49X
  mutate(coverage_label = if_else(mean_coverage <= 49, 'cov_38X', 'cov_60X')) %>%
  inner_join(metadata, by = 'sample_id') %>%
  # calculate the proportion of 38X samples per cancer type
  group_by(cancer_type, cancer_type_code, coverage_label) %>%
  mutate(n_per_cov_label = n()) %>%
  group_by(cancer_type) %>%
  mutate(n_per_cancer_type = n()) %>%
  ungroup() %>%
  distinct(., cancer_type, cancer_type_code, coverage_label, n_per_cov_label, n_per_cancer_type) %>%
  complete(., cancer_type, coverage_label, fill = list(n_per_cov_label = 0)) %>%
  group_by(cancer_type) %>%
  fill(., c(cancer_type_code, n_per_cancer_type), .direction = 'updown') %>%
  ungroup() %>%
  # filter for 38X proportion only
  filter(coverage_label == 'cov_38X') %>%
  mutate(cov_38X_prop = n_per_cov_label / n_per_cancer_type)
  
  
# -------------------- TMB per mut type for SMNVs and SVs

gz_file <- gzfile(paste0(base_dir, '/results/11_downsampling/02_extract_mutational_signatures_from_vcfs/', context_results_date, 
                         '/results/matrices/smnv_contexts.txt.gz'))
smnv_tmb <- read.table(file = gz_file,
                       row.names = 1, # first col are the rownames
                       check.names = F # avoid changes to colnames with special characters
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

sv_tmb_init <- read_tsv(file = paste0(base_dir, '/results/11_downsampling/02_extract_mutational_signatures_from_vcfs/', context_results_date, 
                                      '/results/matrices/sv_contexts.txt'))

#########################################################################################
### STRATIFIED SAMPLING #################################################################

# def vars
n_replicates <- 100

# custom function to perform stratified sampling
do_stratified_sampling <- function(df_input, cov_38X_prop) {
  
  stratified_samples <- df_input %>%
    filter(coverage_label %in% c('cov_38X', 'cov_60X')) %>%
    arrange(sample_id, coverage_label) %>%
    group_by(sample_id) %>%
    slice_sample(., n = 1, weight_by = c(cov_38X_prop, 1 - cov_38X_prop), replace = FALSE) %>%
    ungroup()
  
  return(stratified_samples)
  
}

#########################################################################################
### 95% CONFIDENCE INTERVAL CALCULATION VARIABLES #######################################

# variables for 95% confint calculation
sample_size <- nrow(smnv_tmb) / 4
alpha <- 0.05
degrees_of_freedom <- sample_size - 1
t_score <- qt(p = alpha/2, df = degrees_of_freedom, lower.tail = FALSE)

#########################################################################################
### SNV #################################################################################

snv_tmb <- smnv_tmb %>%
  select(sample_id_2, coverage_label, contains('['))

# calc snv tmb
snv_tmb <- snv_tmb %>%
  rowwise() %>%
  summarise(
    sample_id = sample_id_2,
    coverage_label = coverage_label,
    tmb = sum(c_across(where(is.numeric))))

snv_original <- snv_tmb %>%
  filter(coverage_label == 'cov_100X') %>%
  select(-coverage_label) %>%
  rename(snv_tmb_original = tmb)

calc_sensitivity_loss <- function(df_input, original_df_input, cov_38X_prop) {
  
  df_stratified <- df_input %>%
    do_stratified_sampling(., cov_38X_prop = cov_38X_prop) %>%
    select(-coverage_label) %>%
    rename(snv_tmb_downsampled = tmb)
  
  df_out <- df_stratified %>%
    inner_join(original_df_input, by = 'sample_id') %>%
    summarise(median_tmb_downsampled = median(snv_tmb_downsampled),
              median_tmb_original = median(snv_tmb_original)) %>%
    mutate(sensitivity_loss = (median_tmb_original - median_tmb_downsampled) / median_tmb_downsampled) %>%
    mutate(cov_38X_prop = cov_38X_prop, 
           cov_60X_prop = 1 - cov_38X_prop,
           cov_38X_prop_label = paste0(cov_38X_prop * 100, '%'),
           cov_60X_prop_label = paste0((1 - cov_38X_prop) * 100, '%'))
  
  return(df_out)
  
}

# initiate vector of cov_38X sampling weights
cov_38X_sampling_weights <- seq(0, 1, by = 0.2)

# replicate stratified sampling N times
replication_list <- replicate(n_replicates, map_dfr(cov_38X_sampling_weights, ~calc_sensitivity_loss(df_input = snv_tmb, 
                                                                        original_df_input = snv_original, 
                                                                        cov_38X_prop = .x)), simplify = F)

# add a replicate ID column per replicate
snv_sensitivity_loss <- map_dfr(replication_list, ~mutate(., replicate_id = uuid::UUIDgenerate()))

# calculate mean, sd, sem, and lower/upper CI of mean sensitivity loss estimate
snv_sensitivity_loss <- snv_sensitivity_loss %>%
  group_by(cov_38X_prop, cov_60X_prop, cov_38X_prop_label, cov_60X_prop_label) %>%
  summarise(mean_sensitivity_loss = mean(sensitivity_loss),
            sd_sensitivity_loss = sd(sensitivity_loss),
            sem_sensitivity_loss = sd_sensitivity_loss / sqrt(sample_size),
            lower_ci = mean_sensitivity_loss - t_score * sem_sensitivity_loss,
            upper_ci = mean_sensitivity_loss + t_score * sem_sensitivity_loss, .groups = 'drop')

# create linear model for SNVs
snv_model <- lm(mean_sensitivity_loss ~ cov_38X_prop, data = snv_sensitivity_loss)

# use model to predict mean_sensitivity_loss at calculated 38X sample proportions for each cancer type
pcawg_coverage$snv_sensitivity_loss <-  predict(snv_model, newdata = pcawg_coverage)

# plot sensitivity loss estimate plus standard deviation
snv_p1 <- snv_sensitivity_loss %>%
  mutate(sd_sensitivity_loss = if_else(sd_sensitivity_loss == 0, NA_real_, sd_sensitivity_loss)) %>%
  ggplot(., aes(x = cov_38X_prop, y = mean_sensitivity_loss)) +
  geom_errorbar(aes(ymin = mean_sensitivity_loss + sd_sensitivity_loss, ymax = mean_sensitivity_loss - sd_sensitivity_loss), width = 0.1) +
  geom_point() +
  geom_line(size = 1) +
  geom_point(data = pcawg_coverage %>% mutate(cancer_type_code = factor(cancer_type_code, levels = cancer_type_code_order)),
             aes(x = cov_38X_prop, y = snv_sensitivity_loss, color = cancer_type_code), key_glyph = draw_key_rect, size = 2) +
  geom_text_repel(
    data = pcawg_coverage %>% mutate(cancer_type_code = factor(cancer_type_code, levels = cancer_type_code_order)),
    aes(x = cov_38X_prop, y = snv_sensitivity_loss, label = cancer_type_code),
    min.segment.length = 0, ylim = c(-Inf, Inf), box.padding = 0.5
  ) +
  scale_x_continuous(
    breaks = cov_38X_sampling_weights,
    labels = paste0((cov_38X_sampling_weights * 100), '%')
  ) +
  scale_y_continuous(
    limits = c(-0.1,1),
    breaks = seq(-0.1, 1, by = 0.1),
    labels = paste0(seq(-0.1, 1, by = 0.1), 'x')
  ) +
  scale_color_manual(values = cancer_type_palette) +
  labs(
    title = 'SNV',
    x = 'Proportion of 38X samples',
    y = 'Sensitivity loss [fold change]'
  ) +
  theme_bw() +
  theme_bigfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    legend.position = 'none'
  )

#########################################################################################
### DBS #################################################################################

dbs_tmb <- smnv_tmb %>%
  select(sample_id_2, coverage_label, matches('^[A-Z]{2}\\>'))

# calc dbs tmb
dbs_tmb <- dbs_tmb %>%
  rowwise() %>%
  summarise(
    sample_id = sample_id_2,
    coverage_label = coverage_label,
    tmb = sum(c_across(where(is.numeric))))

dbs_original <- dbs_tmb %>%
  filter(coverage_label == 'cov_100X') %>%
  select(-coverage_label) %>%
  rename(dbs_tmb_original = tmb)

calc_sensitivity_loss <- function(df_input, original_df_input, cov_38X_prop) {
  
  df_stratified <- df_input %>%
    do_stratified_sampling(., cov_38X_prop = cov_38X_prop) %>%
    select(-coverage_label) %>%
    rename(dbs_tmb_downsampled = tmb)
  
  df_out <- df_stratified %>%
    inner_join(original_df_input, by = 'sample_id') %>%
    summarise(median_tmb_downsampled = median(dbs_tmb_downsampled),
              median_tmb_original = median(dbs_tmb_original)) %>%
    mutate(sensitivity_loss = (median_tmb_original - median_tmb_downsampled) / median_tmb_downsampled) %>%
    mutate(cov_38X_prop = cov_38X_prop, 
           cov_60X_prop = 1 - cov_38X_prop,
           cov_38X_prop_label = paste0(cov_38X_prop * 100, '%'),
           cov_60X_prop_label = paste0((1 - cov_38X_prop) * 100, '%'))
  
  return(df_out)
  
}

# initiate vector of cov_38X sampling weights
cov_38X_sampling_weights <- seq(0, 1, by = 0.2)

# replicate stratified sampling N times
replication_list <- replicate(n_replicates, map_dfr(cov_38X_sampling_weights, ~calc_sensitivity_loss(df_input = dbs_tmb, 
                                                                                            original_df_input = dbs_original, 
                                                                                            cov_38X_prop = .x)),simplify = F)

# add a replicate ID column per replicate
dbs_sensitivity_loss <- map_dfr(replication_list, ~mutate(., replicate_id = uuid::UUIDgenerate()))

# calculate mean, sd, sem, and lower/upper CI of mean sensitivity loss estimate
dbs_sensitivity_loss <- dbs_sensitivity_loss %>%
  group_by(cov_38X_prop, cov_60X_prop, cov_38X_prop_label, cov_60X_prop_label) %>%
  summarise(mean_sensitivity_loss = mean(sensitivity_loss),
            sd_sensitivity_loss = sd(sensitivity_loss),
            sem_sensitivity_loss = sd_sensitivity_loss / sqrt(sample_size),
            lower_ci = mean_sensitivity_loss - t_score * sem_sensitivity_loss,
            upper_ci = mean_sensitivity_loss + t_score * sem_sensitivity_loss, .groups = 'drop')

# create linear model for DBS
dbs_model <- lm(mean_sensitivity_loss ~ cov_38X_prop, data = dbs_sensitivity_loss)

# use model to predict mean_sensitivity_loss at calculated 38X sample proportions for each cancer type
pcawg_coverage$dbs_sensitivity_loss <-  predict(dbs_model, newdata = pcawg_coverage)

# plot sensitivity loss estimate plus 95% CI
dbs_p1 <- dbs_sensitivity_loss %>%
  mutate(sd_sensitivity_loss = if_else(sd_sensitivity_loss == 0, NA_real_, sd_sensitivity_loss)) %>%
  ggplot(., aes(x = cov_38X_prop, y = mean_sensitivity_loss)) +
  # geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.1) +
  geom_errorbar(aes(ymin = mean_sensitivity_loss + sd_sensitivity_loss, ymax = mean_sensitivity_loss - sd_sensitivity_loss), width = 0.1) +
  geom_point() +
  geom_line(size = 1) +
  geom_point(data = pcawg_coverage %>% mutate(cancer_type_code = factor(cancer_type_code, levels = cancer_type_code_order)),
             aes(x = cov_38X_prop, y = dbs_sensitivity_loss, color = cancer_type_code), key_glyph = draw_key_rect, size = 2) +
  geom_text_repel(
    data = pcawg_coverage %>% mutate(cancer_type_code = factor(cancer_type_code, levels = cancer_type_code_order)),
    aes(x = cov_38X_prop, y = dbs_sensitivity_loss, label = cancer_type_code),
    min.segment.length = 0, ylim = c(-Inf, Inf), box.padding = 0.5
  ) +
  scale_x_continuous(
    breaks = cov_38X_sampling_weights,
    labels = paste0((cov_38X_sampling_weights * 100), '%')
  ) +
  scale_y_continuous(
    limits = c(-0.1,1),
    breaks = seq(-0.1, 1, by = 0.1),
    labels = paste0(seq(-0.1, 1, by = 0.1), 'x')
  ) +
  scale_color_manual(values = cancer_type_palette) +
  labs(
    title = 'DBS',
    x = 'Proportion of 38X samples',
    y = 'Sensitivity loss [fold change]'
  ) +
  theme_bw() +
  theme_bigfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    legend.position = 'none'
  )

#########################################################################################
### Indel ###############################################################################

indel_tmb <- smnv_tmb %>%
  select(sample_id_2, coverage_label, matches('^del|^ins'))

# calc dbs tmb
indel_tmb <- indel_tmb %>%
  rowwise() %>%
  summarise(
    sample_id = sample_id_2,
    coverage_label = coverage_label,
    tmb = sum(c_across(where(is.numeric))))

indel_original <- indel_tmb %>%
  filter(coverage_label == 'cov_100X') %>%
  select(-coverage_label) %>%
  rename(indel_tmb_original = tmb)

calc_sensitivity_loss <- function(df_input, original_df_input, cov_38X_prop) {
  
  df_stratified <- df_input %>%
    do_stratified_sampling(., cov_38X_prop = cov_38X_prop) %>%
    select(-coverage_label) %>%
    rename(indel_tmb_downsampled = tmb)
  
  df_out <- df_stratified %>%
    inner_join(original_df_input, by = 'sample_id') %>%
    summarise(median_tmb_downsampled = median(indel_tmb_downsampled),
              median_tmb_original = median(indel_tmb_original)) %>%
    mutate(sensitivity_loss = (median_tmb_original - median_tmb_downsampled) / median_tmb_downsampled) %>%
    mutate(cov_38X_prop = cov_38X_prop, 
           cov_60X_prop = 1 - cov_38X_prop,
           cov_38X_prop_label = paste0(cov_38X_prop * 100, '%'),
           cov_60X_prop_label = paste0((1 - cov_38X_prop) * 100, '%'))
  
  return(df_out)
  
}

# initiate vector of cov_38X sampling weights
cov_38X_sampling_weights <- seq(0, 1, by = 0.2)

# replicate stratified sampling N times
replication_list <- replicate(n_replicates, map_dfr(cov_38X_sampling_weights, ~calc_sensitivity_loss(df_input = indel_tmb, 
                                                                                            original_df_input = indel_original, 
                                                                                            cov_38X_prop = .x)),simplify = F)

# add a replicate ID column per replicate
indel_sensitivity_loss <- map_dfr(replication_list, ~mutate(., replicate_id = uuid::UUIDgenerate()))

# calculate mean, sd, sem, and lower/upper CI of mean sensitivity loss estimate
indel_sensitivity_loss <- indel_sensitivity_loss %>%
  group_by(cov_38X_prop, cov_60X_prop, cov_38X_prop_label, cov_60X_prop_label) %>%
  summarise(mean_sensitivity_loss = mean(sensitivity_loss),
            sd_sensitivity_loss = sd(sensitivity_loss),
            sem_sensitivity_loss = sd_sensitivity_loss / sqrt(sample_size),
            lower_ci = mean_sensitivity_loss - t_score * sem_sensitivity_loss,
            upper_ci = mean_sensitivity_loss + t_score * sem_sensitivity_loss, .groups = 'drop')

# create linear model for INDELs
indel_model <- lm(mean_sensitivity_loss ~ cov_38X_prop, data = indel_sensitivity_loss)

# use model to predict mean_sensitivity_loss at calculated 38X sample proportions for each cancer type
pcawg_coverage$indel_sensitivity_loss <-  predict(indel_model, newdata = pcawg_coverage)

# plot sensitivity loss estimate plus 95% CI
indel_p1 <- indel_sensitivity_loss %>%
  mutate(sd_sensitivity_loss = if_else(sd_sensitivity_loss == 0, NA_real_, sd_sensitivity_loss)) %>%
  ggplot(., aes(x = cov_38X_prop, y = mean_sensitivity_loss)) +
  # geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.1) +
  geom_errorbar(aes(ymin = mean_sensitivity_loss + sd_sensitivity_loss, ymax = mean_sensitivity_loss - sd_sensitivity_loss), width = 0.1) +
  geom_point() +
  geom_line(size = 1) +
  geom_point(data = pcawg_coverage %>% mutate(cancer_type_code = factor(cancer_type_code, levels = cancer_type_code_order)),
             aes(x = cov_38X_prop, y = indel_sensitivity_loss, color = cancer_type_code), key_glyph = draw_key_rect, size = 2) +
  geom_text_repel(
    data = pcawg_coverage %>% mutate(cancer_type_code = factor(cancer_type_code, levels = cancer_type_code_order)),
    aes(x = cov_38X_prop, y = indel_sensitivity_loss, label = cancer_type_code),
    min.segment.length = 0, ylim = c(-Inf, Inf), box.padding = 0.5
  ) +
  scale_x_continuous(
    breaks = cov_38X_sampling_weights,
    labels = paste0((cov_38X_sampling_weights * 100), '%')
  ) +
  scale_y_continuous(
    limits = c(-0.1,1),
    breaks = seq(-0.1, 1, by = 0.1),
    labels = paste0(seq(-0.1, 1, by = 0.1), 'x')
  ) +
  scale_color_manual(values = cancer_type_palette) +
  labs(
    title = 'INDEL',
    x = 'Proportion of 38X samples',
    y = 'Sensitivity loss [fold change]'
  ) +
  theme_bw() +
  theme_bigfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    legend.position = 'none'
  )

#########################################################################################
### SV ##################################################################################

sv_tmb <- sv_tmb_init %>%
  mutate(sample = str_replace(sample, pattern = regex('([ABC])T$'), replacement = '\\1')) %>%
  inner_join(smnv_tmb, by = c('sample' = 'sample_id')) %>%
  select(sample_id_2, coverage_label, SV_LOAD) %>%
  rename(
    sample_id = sample_id_2,
    tmb = SV_LOAD)

sv_original <- sv_tmb %>%
  filter(coverage_label == 'cov_100X') %>%
  select(-coverage_label) %>%
  rename(sv_tmb_original = tmb)

calc_sensitivity_loss <- function(df_input, original_df_input, cov_38X_prop) {
  
  df_stratified <- df_input %>%
    do_stratified_sampling(., cov_38X_prop = cov_38X_prop) %>%
    select(-coverage_label) %>%
    rename(sv_tmb_downsampled = tmb)
  
  df_out <- df_stratified %>%
    inner_join(original_df_input, by = 'sample_id') %>%
    summarise(median_tmb_downsampled = median(sv_tmb_downsampled),
              median_tmb_original = median(sv_tmb_original)) %>%
    mutate(sensitivity_loss = (median_tmb_original - median_tmb_downsampled) / median_tmb_downsampled) %>%
    mutate(cov_38X_prop = cov_38X_prop, 
           cov_60X_prop = 1 - cov_38X_prop,
           cov_38X_prop_label = paste0(cov_38X_prop * 100, '%'),
           cov_60X_prop_label = paste0((1 - cov_38X_prop) * 100, '%'))
  
  return(df_out)
  
}

# initiate vector of cov_38X sampling weights
cov_38X_sampling_weights <- seq(0, 1, by = 0.2)

# replicate stratified sampling N times
replication_list <- replicate(n_replicates, map_dfr(cov_38X_sampling_weights, ~calc_sensitivity_loss(df_input = sv_tmb, 
                                                                                            original_df_input = sv_original, 
                                                                                            cov_38X_prop = .x)),simplify = F)

# add a replicate ID column per replicate
sv_sensitivity_loss <- map_dfr(replication_list, ~mutate(., replicate_id = uuid::UUIDgenerate()))

# calculate mean, sd, sem, and lower/upper CI of mean sensitivity loss estimate
sv_sensitivity_loss <- sv_sensitivity_loss %>%
  group_by(cov_38X_prop, cov_60X_prop, cov_38X_prop_label, cov_60X_prop_label) %>%
  summarise(mean_sensitivity_loss = mean(sensitivity_loss),
            sd_sensitivity_loss = sd(sensitivity_loss),
            sem_sensitivity_loss = sd_sensitivity_loss / sqrt(sample_size),
            lower_ci = mean_sensitivity_loss - t_score * sem_sensitivity_loss,
            upper_ci = mean_sensitivity_loss + t_score * sem_sensitivity_loss, .groups = 'drop')

# create linear model for SVs
sv_model <- lm(mean_sensitivity_loss ~ cov_38X_prop, data = sv_sensitivity_loss)

# use model to predict mean_sensitivity_loss at calculated 38X sample proportions for each cancer type
pcawg_coverage$sv_sensitivity_loss <-  predict(sv_model, newdata = pcawg_coverage)

# plot sensitivity loss estimate plus 95% CI
sv_p1 <- sv_sensitivity_loss %>%
  mutate(sd_sensitivity_loss = if_else(sd_sensitivity_loss == 0, NA_real_, sd_sensitivity_loss)) %>%
  ggplot(., aes(x = cov_38X_prop, y = mean_sensitivity_loss)) +
  geom_errorbar(aes(ymin = mean_sensitivity_loss + sd_sensitivity_loss, ymax = mean_sensitivity_loss - sd_sensitivity_loss), width = 0.1) +
  geom_point() +
  geom_line(size = 1) +
  geom_point(data = pcawg_coverage %>% mutate(cancer_type_code = factor(cancer_type_code, levels = cancer_type_code_order)),
             aes(x = cov_38X_prop, y = sv_sensitivity_loss, color = cancer_type_code), key_glyph = draw_key_rect, size = 2) +
  geom_text_repel(
    data = pcawg_coverage %>% mutate(cancer_type_code = factor(cancer_type_code, levels = cancer_type_code_order)),
    aes(x = cov_38X_prop, y = sv_sensitivity_loss, label = cancer_type_code),
    min.segment.length = 0, ylim = c(-Inf, Inf), box.padding = 0.5
  ) +
  scale_x_continuous(
    breaks = cov_38X_sampling_weights,
    labels = paste0((cov_38X_sampling_weights * 100), '%')
  ) +
  scale_y_continuous(
    limits = c(-0.1,1),
    breaks = seq(-0.1, 1, by = 0.1),
    labels = paste0(seq(-0.1, 1, by = 0.1), 'x')
  ) +
  scale_color_manual(values = cancer_type_palette) +
  labs(
    title = 'SV',
    x = 'Proportion of 38X samples',
    y = 'Sensitivity loss [fold change]'
  ) +
  theme_bw() +
  theme_bigfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    legend.position = 'none'
  )

# stitch plots together
slope_plots <- ggarrange(snv_p1,
                         dbs_p1,
                         indel_p1,
                         sv_p1,
                         labels = NULL,
                         ncol = 2, nrow = 2,
                         align = "hv",
                         font.label = list(size = 8, color = "black", face = "bold", family = NULL, position = "top"))

# create output dir
output_dir <- paste0(output_dir, '/11_downsampling/05_plot_sensitivity_loss/', current_date, '/results/plots')
dir.create(path = paste0(output_dir), recursive = TRUE)

# save figure as pdf
pdf(
  file = paste0(output_dir, '/', current_date, '_smnv_sensitivity_loss.pdf'),
  width = 18, 
  height = 12,
  useDingbats = FALSE,
  compress = FALSE
)
print(slope_plots)
dev.off()
