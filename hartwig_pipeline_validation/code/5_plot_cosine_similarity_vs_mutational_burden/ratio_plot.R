### RATIO PLOTS

#########################################################################################
### SNV #################################################################################

# get the percentage of tumors that have a consine similarity over 0.90
snv_percent_over_90 <- snv_cos_tot_mut_combined %>%
  summarise(percentage_over_90 = mean(over_90, na.rm = TRUE)) %>% 
  pull()

# get the reverse value (percentage under 0.90)
snv_n_below_90 <- floor((1 - snv_percent_over_90) * snv_n_samples)

# plot the mutational burden ratio against the cosine similarity of 96 SNV context profiles
snv_p2 <- ggplot(snv_cos_tot_mut_combined, aes(x = snv_cosine_similarity, y = log2(mut_burden_ratio), color = over_90)) +
  geom_point(alpha = 0.4, shape = 16, size = 3) +
  geom_vline(xintercept = 0.90, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(limits = c(-5, 5), breaks = seq(-6, 6, by = 1)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
  scale_color_manual(values = c('red', 'black')) +
  annotate(geom = 'curve', xend = 0.89, y = 4, x = 0.7, yend = 2.5, curvature = 0.25, 
           arrow = arrow(length = unit(3.5, 'mm'))) +
  annotate('text', x = 0.65, y = 4.5, label = '0.90 cosine similarity', size = 5, colour = 'black', fontface = 2) +
  labs(title = paste0('SNV')
       ) +
  theme_bw() +
  theme_bigfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    legend.position = 'none',
    plot.margin = margin(t= 5, l = 5, b = 5, r = 5)) +
  # have to remove the axis labels here and not after applying ggMarginal (plot will not render somehow)
  rremove('xlab') +
  rremove('ylab')

snv_p2 <- ggExtra::ggMarginal(snv_p2, type = 'histogram', margins = 'y')

#########################################################################################
### DBS #################################################################################

# get the percentage of tumors that have a consine similarity over 0.90
dbs_percent_over_90 <- dbs_cos_tot_mut_combined %>%
  summarise(percentage_over_90 = mean(over_90, na.rm = TRUE)) %>% 
  pull()

# get the reverse value (percentage under 0.90)
dbs_n_below_90 <- floor((1 - dbs_percent_over_90) * dbs_n_samples)

# plot the mutational burden ratio against the cosine similarity of 96 dbs context profiles
dbs_p2 <- ggplot(dbs_cos_tot_mut_combined, aes(x = dbs_cosine_similarity, y = log2(mut_burden_ratio), color = over_90)) +
  geom_point(alpha = 0.4, shape = 16, size = 3) +
  geom_vline(xintercept = 0.90, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(limits = c(-2, 6), breaks = seq(-6, 6, by = 1)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
  scale_color_manual(values = c('red', 'black')) +
  labs(title = paste0('DBS')
       ) +
  theme_bw() +
  theme_bigfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    legend.position = 'none') +
  rremove('xlab') +
  rremove('ylab')

dbs_p2 <- ggExtra::ggMarginal(dbs_p2, type = 'histogram', margins = 'y')

#########################################################################################
### INDEL ###############################################################################

# get the percentage of tumors that have a consine similarity over 0.90
indel_percent_over_90 <- indel_cos_tot_mut_combined %>%
  summarise(percentage_over_90 = mean(over_90, na.rm = TRUE)) %>% pull()

# get the reverse value (percentage under 0.90)
indel_n_below_90 <- floor((1 - indel_percent_over_90) * indel_n_samples)

# plot the mutational burden ratio against the cosine similarity of 96 indel context profiles
indel_p2 <- ggplot(indel_cos_tot_mut_combined, aes(x = indel_cosine_similarity, y = log2(mut_burden_ratio), color = over_90)) +
  geom_point(alpha = 0.4, shape = 16, size = 3) +
  geom_vline(xintercept = 0.90, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(limits = c(-4, 5), breaks = seq(-6, 6, by = 1)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
  scale_color_manual(values = c('red', 'black')) +
  labs(title = paste0('INDEL')
       ) +
  theme_bw() +
  theme_bigfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    legend.position = 'none') +
  rremove('xlab') +
  rremove('ylab')

indel_p2 <- ggExtra::ggMarginal(indel_p2, type = 'histogram', margins = 'y')

#########################################################################################
### SV $#################################################################################

# get the percentage of tumors that have a cosine similarity over 0.90
sv_percent_over_90_len_cutoff <- sv_cos_tot_mut_combined %>%
  summarise(percentage_over_90_len_cutoff = mean(over_90_len_cutoff, na.rm = TRUE)) %>% pull()

# get the reverse value (percentage under 0.90)
sv_n_below_90_len_cutoff <- floor((1 - sv_percent_over_90_len_cutoff) * sv_n_samples)

# filter the high ratios out so the plot does not look like a pancake
sv_cos_tot_mut_combined_filtered_high_ratios_len_cutoff <- sv_cos_tot_mut_combined

# new sample number
sv_n_samples_len_cutoff <- nrow(sv_cos_tot_mut_combined_filtered_high_ratios_len_cutoff)

# get the percentage of tumors that have a cosine similarity over 0.90
sv_percent_over_90_len_cutoff <- sv_cos_tot_mut_combined_filtered_high_ratios_len_cutoff %>%
  summarise(percentage_over_90_len_cutoff = mean(over_90_len_cutoff, na.rm = TRUE)) %>% pull()

# get the reverse value (percentage under 0.90)
sv_n_below_90_len_cutoff <- floor((1 - sv_percent_over_90_len_cutoff) * sv_n_samples)

# plot the mutational burden ratio against the cosine similarity of 96 sv context profiles
sv_p2_len_cutoff <- ggplot(sv_cos_tot_mut_combined_filtered_high_ratios_len_cutoff, 
                           aes(x = sv_cosine_similarity_len_cutoff, 
                               y = log2(mut_burden_ratio_len_cutoff), 
                               color = over_90_len_cutoff)) +
  geom_point(alpha = 0.4, shape = 16, size = 3) +
  geom_vline(xintercept = 0.90, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(limits = c(-5, 5), breaks = seq(-6, 6, by = 1)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
  scale_color_manual(values = c('red', 'black')) +
  labs(title = paste0('SV length > ', sv_length_cutoff, ' bp')
       ) +
  theme_bw() +
  theme_bigfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    legend.position = 'none') +
  rremove('xlab') +
  rremove('ylab')

sv_p2_len_cutoff <- ggExtra::ggMarginal(sv_p2_len_cutoff, type = 'histogram', margins = 'y')

# stitch plots together
ratio_plots <- ggarrange(snv_p2,
                         dbs_p2, 
                         indel_p2, 
                         sv_p2_len_cutoff,
                         labels = NULL,
                         ncol = 1, nrow = 4,
                         align = 'hv', 
                         font.label = list(size = 10, color = 'black', face = 'bold', family = NULL, position = 'top'))

ratio_plots <- annotate_figure(ratio_plots, 
                               left = text_grob('Log2(TMB ratio Hartwig/PCAWG)', size = 18, rot = 90, vjust = 0.5),
                               bottom = text_grob('Cosine similarity', size = 18),
                               fig.lab = 'C', fig.lab.size = 24, fig.lab.face = 'bold')
