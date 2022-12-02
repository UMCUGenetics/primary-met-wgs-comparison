### BOXPLOTS: Plot the TMB of PCAWG vs Hartwig pipeline processed PCAWG samples

options(scipen = 999) # disable scientific notation of numbers in plots

#########################################################################################
### SNV #################################################################################

# data prep
snv_cos_tot_mut_combined_box <- snv_cos_tot_mut_combined %>%
  rename(PCAWG = 'snv_pcawg_total_mutational_burden',
         Hartwig = 'snv_hartwig_total_mutational_burden') %>%
  pivot_longer(., cols = c(Hartwig, PCAWG), names_to = 'pipeline', values_to = 'mut_burden')

# get the number of rows
snv_n_samples <- nrow(snv_cos_tot_mut_combined)

snv_y_limit <- 6

# plot a box plot with median tmb as a label and fold change from PCAWG to Hartwig
snv_p1 <- ggplot(snv_cos_tot_mut_combined_box, aes(x = pipeline, y = log10(mut_burden + 1))) +
  geom_jitter(color = '#C0C0C0') +
  geom_boxplot(alpha = 0.8) +
  geom_label(
    data = . %>% group_by(pipeline) %>% summarise(median_tmb = median(mut_burden), .groups = 'drop'),
    aes(y = log10(median_tmb), label = round(median_tmb, 0))
  ) +
  geom_label(
    data = . %>% group_by(pipeline) %>%
      summarise(median_tmb = median(mut_burden), .groups = 'drop') %>%
      arrange(pipeline) %>%
      mutate(fold_change = paste0(round((median_tmb - lead(median_tmb)) / lead(median_tmb), 2), 'x')) %>%
      filter(fold_change != 'NAx'),
    aes(y = log10(median_tmb), label = fold_change),nudge_x = 0.5
  ) +
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

#########################################################################################
### DBS #################################################################################

# data prep
dbs_cos_tot_mut_combined_box <- dbs_cos_tot_mut_combined %>%
  rename(PCAWG = 'dbs_pcawg_total_mutational_burden',
         Hartwig = 'dbs_hartwig_total_mutational_burden') %>%
  pivot_longer(., cols = c(Hartwig, PCAWG), names_to = 'pipeline', values_to = 'mut_burden')

# get the number of rows
dbs_n_samples <- nrow(dbs_cos_tot_mut_combined)

dbs_y_limit <- 5

# plot a box plot with median tmb as a label and fold change from PCAWG to Hartwig
dbs_p1 <- ggplot(dbs_cos_tot_mut_combined_box, aes(x = pipeline, y = log10(mut_burden + 1))) +
  geom_jitter(color = '#C0C0C0') +
  geom_boxplot(alpha = 0.8) +
  geom_label(
    data = . %>% group_by(pipeline) %>% summarise(median_tmb = median(mut_burden), .groups = 'drop'),
    aes(y = log10(median_tmb), label = round(median_tmb, 0))
  ) +
  geom_label(
    data = . %>% group_by(pipeline) %>%
      summarise(median_tmb = median(mut_burden), .groups = 'drop') %>%
      arrange(pipeline) %>%
      mutate(fold_change = paste0(round((median_tmb - lead(median_tmb)) / lead(median_tmb), 2), 'x')) %>%
      filter(fold_change != 'NAx'),
    aes(y = log10(median_tmb), label = fold_change),nudge_x = 0.5
  ) +
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

#########################################################################################
### Indel ###############################################################################

# data prep
indel_cos_tot_mut_combined_box <- indel_cos_tot_mut_combined %>%
  rename(PCAWG = 'indel_pcawg_total_mutational_burden',
         Hartwig = 'indel_hartwig_total_mutational_burden') %>%
  pivot_longer(., cols = c(Hartwig, PCAWG), names_to = 'pipeline', values_to = 'mut_burden')

# get the number of rows
indel_n_samples <- nrow(indel_cos_tot_mut_combined)

indel_y_limit <- 5

# plot a box plot with median tmb as a label and fold change from PCAWG to Hartwig
indel_p1 <- ggplot(indel_cos_tot_mut_combined_box, aes(x = pipeline, y = log10(mut_burden + 1))) +
  geom_jitter(color = '#C0C0C0') +
  geom_boxplot(alpha = 0.8) +
  geom_label(
    data = . %>% group_by(pipeline) %>% summarise(median_tmb = median(mut_burden), .groups = 'drop'),
    aes(y = log10(median_tmb), label = round(median_tmb, 0))
  ) +
  geom_label(
    data = . %>% group_by(pipeline) %>%
      summarise(median_tmb = median(mut_burden), .groups = 'drop') %>%
      arrange(pipeline) %>%
      mutate(fold_change = paste0(round((median_tmb - lead(median_tmb)) / lead(median_tmb), 2), 'x')) %>%
      filter(fold_change != 'NAx'),
    aes(y = log10(median_tmb), label = fold_change),nudge_x = 0.5
  ) +
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

#########################################################################################
### SV ##################################################################################

# ------------------------------------ SV with length cutoff

# data prep
sv_total_mut_burden_len_cutoff_box <- sv_total_mut_burden_len_cutoff %>%
  rename(PCAWG = 'sv_pcawg_total_mutational_burden_len_cutoff',
         Hartwig = 'sv_hartwig_total_mutational_burden_len_cutoff') %>%
  pivot_longer(., cols = c(Hartwig, PCAWG), names_to = 'pipeline', values_to = 'mut_burden')

# get the number of rows
sv_n_samples <- nrow(sv_cos_tot_mut_combined)

sv_y_limit <- 4

# plot a box plot with median tmb as a label and fold change from PCAWG to Hartwig
sv_p1_len_cutoff <- ggplot(sv_total_mut_burden_len_cutoff_box, aes(x = pipeline, y = log10(mut_burden + 1))) +
  geom_jitter(color = '#C0C0C0') +
  geom_boxplot(alpha = 0.8) +
  geom_label(
    data = . %>% group_by(pipeline) %>% summarise(median_tmb = median(mut_burden), .groups = 'drop'),
    aes(y = log10(median_tmb), label = round(median_tmb, 0))
  ) +
  geom_label(
    data = . %>% group_by(pipeline) %>%
      summarise(median_tmb = median(mut_burden), .groups = 'drop') %>%
      arrange(pipeline) %>%
      mutate(fold_change = paste0(round((median_tmb - lead(median_tmb)) / lead(median_tmb), 2), 'x')) %>%
      filter(fold_change != 'NAx'),
    aes(y = log10(median_tmb), label = fold_change),nudge_x = 0.5
  ) +
  scale_y_continuous(
    breaks = seq(0,sv_y_limit, by = 1),
    labels = 10^seq(0,sv_y_limit, by = 1), 
    limits = c(0, sv_y_limit)) +
  labs(title = paste0('SV length >', sv_length_cutoff, 'bp'),
       x = '',
       y = 'Log10(Mut. load + 1)') +
  theme_bw() +
  theme_bigfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.margin = margin(t= 5, l = 5, b = 5, r = 15)) 

# stitch plots together
boxplot <- ggarrange(snv_p1,
                     dbs_p1,
                     indel_p1,
                     sv_p1_len_cutoff,
                     labels = NULL,
                     ncol = 4, nrow = 1,
                     align = "hv",
                     font.label = list(size = 8, color = "black", face = "bold", family = NULL, position = "top"))

boxplot <- annotate_figure(boxplot,
                           fig.lab = 'A', fig.lab.size = 24, fig.lab.face = 'bold')