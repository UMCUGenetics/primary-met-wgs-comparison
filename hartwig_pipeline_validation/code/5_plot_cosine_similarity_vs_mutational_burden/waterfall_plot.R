### WATERFALL PLOT: plot cos sim by Hartwig tmb

#########################################################################################
### SNV $################################################################################

snv_p3 <- ggplot(snv_cos_tot_mut_combined, aes(x = log10(snv_hartwig_total_mutational_burden), y = snv_cosine_similarity)) +
  geom_point(alpha = 0.3) +
  scale_x_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 1)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(title = paste0('SNV'),
       x = 'Log10(Hartwig Total mutational burden)',
       y = 'Cosine similarity') +
  theme_bw() +
  theme_bigfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.margin = margin(t= 5, l = 5, b = 5, r = 15))

#########################################################################################
### DBS #################################################################################

dbs_p3 <- ggplot(dbs_cos_tot_mut_combined, aes(x = log10(dbs_hartwig_total_mutational_burden), y = dbs_cosine_similarity)) +
  geom_point(alpha = 0.3) +
  scale_x_continuous(limits = c(0, 4.5), breaks = seq(0, 6, by = 1)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(title = paste0('DBS'),
       x = 'Log10(Hartwig Total mutational burden)',
       y = 'Cosine similarity') +
  theme_bw() +
  theme_bigfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.margin = margin(t= 5, l = 5, b = 5, r = 15))

#########################################################################################
### INDEL ###############################################################################

indel_p3 <- ggplot(indel_cos_tot_mut_combined, aes(x = log10(indel_hartwig_total_mutational_burden), y = indel_cosine_similarity)) +
  geom_point(alpha = 0.3) +
  scale_x_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 1)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(title = paste0('INDEL'),
       x = 'Log10(Hartwig Total mutational burden)',
       y = 'Cosine similarity') +
  theme_bw() +
  theme_bigfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.margin = margin(t= 5, l = 5, b = 5, r = 15))

#########################################################################################
### SV ##################################################################################

sv_p3_len_cutoff <- ggplot(sv_cos_tot_mut_combined, aes(x = log10(sv_hartwig_total_mutational_burden_len_cutoff), y = sv_cosine_similarity_len_cutoff)) +
  geom_point(alpha = 0.3) +
  scale_x_continuous(limits = c(0, 4), breaks = seq(0, 6, by = 1)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(title = paste0('SV length >', sv_length_cutoff, 'bp'),
       x = 'Log10(Hartwig Total mutational burden)',
       y = 'Cosine similarity') +
  theme_bw() +
  theme_bigfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.margin = margin(t= 5, l = 5, b = 5, r = 15))

# stitch plots together
waterfall_plots <- ggarrange(snv_p3 + rremove('ylab') + rremove('xlab'),
                             dbs_p3 + rremove('ylab') + rremove('xlab'), 
                             indel_p3 + rremove('ylab') + rremove('xlab'), 
                             sv_p3_len_cutoff+ rremove('ylab') + rremove('xlab'),
                    labels = NULL,
                    ncol = 1, nrow = 4,
                    align = 'hv', 
                    font.label = list(size = 10, color = 'black', face = 'bold', family = NULL, position = 'top'))

waterfall_plots <- annotate_figure(waterfall_plots, 
                                   left = text_grob('Cosine similarity', size = 18, rot = 90, vjust = 0.5),
                                   bottom = text_grob('Log10(Hartwig TMB)', size = 18),
                                   fig.lab = 'B', fig.lab.size = 24, fig.lab.face = 'bold')
                