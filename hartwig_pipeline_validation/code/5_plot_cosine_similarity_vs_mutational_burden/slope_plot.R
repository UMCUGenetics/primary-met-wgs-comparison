### SLOPE PLOTS: Plot the TMB of PCAWG vs Hartwig plus the slope

#########################################################################################
### SNV #################################################################################

# get the number of rows
snv_n_samples <- nrow(snv_cos_tot_mut_combined)

# build linear model Hartwig total burden explained by pcawg total mutational burden
snv_model <- lm(snv_hartwig_total_mutational_burden ~ snv_pcawg_total_mutational_burden, data = snv_cos_tot_mut_combined)

# get coefs of this model
snv_slope <- snv_model$coefficients[2]

snv_y_limit <- 1000000

# plot this relationship and the regression line
snv_p1 <- ggplot(snv_cos_tot_mut_combined, aes(x = snv_pcawg_total_mutational_burden, y = snv_hartwig_total_mutational_burden)) +
  geom_point(color = '#C0C0C0') +
  geom_smooth(method = 'lm', se = FALSE, color = '#000000') +
  scale_x_continuous(labels = unit_format(unit = 'M', scale = 1e-6), limits = c(0, snv_y_limit)) +
  scale_y_continuous(labels = unit_format(unit = 'M', scale = 1e-6), limits = c(0, snv_y_limit)) +
  annotate('text', x = snv_y_limit * 0.25, y = snv_y_limit * 0.75, label = paste0('Slope: ', round(snv_slope, 2)),
           size = 5, colour = 'black', fontface = 2) +
  labs(title = 'SNV',
       x = 'PCAWG SNV TMB',
       y = 'Hartwig SNV TMB') +
  theme_bw() +
  theme_bigfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.margin = margin(t= 5, l = 5, b = 5, r = 15))

#########################################################################################
### DBS #################################################################################

# get the number of rows
dbs_n_samples <- nrow(dbs_cos_tot_mut_combined)

# build linear model Hartwig total burden explained by pcawg total mutational burden
dbs_model <- lm(dbs_hartwig_total_mutational_burden ~ dbs_pcawg_total_mutational_burden, data = dbs_cos_tot_mut_combined)

# get coefs of this model
dbs_slope <- dbs_model$coefficients[2]

dbs_y_limit <- 20000

# plot this relationship and the regression line
dbs_p1 <- ggplot(dbs_cos_tot_mut_combined, aes(x = dbs_pcawg_total_mutational_burden, y = dbs_hartwig_total_mutational_burden)) +
  geom_point(color = '#C0C0C0') +
  geom_smooth(method = 'lm', se = FALSE, color = '#000000') +
  scale_x_continuous(labels = unit_format(unit = 'K', scale = 1e-3),
                     breaks = seq(0, dbs_y_limit, by = 4000),
                     limits = c(0, dbs_y_limit)) +
  scale_y_continuous(labels = unit_format(unit = 'K', scale = 1e-3),
                     breaks = seq(0, dbs_y_limit, by = 4000),
                     limits = c(0, dbs_y_limit)) +
  annotate('text', x = dbs_y_limit * 0.25, y = dbs_y_limit * 0.75, label = paste0('Slope: ', round(dbs_slope, 2)),
           size = 5, colour = 'black', fontface = 2) +
  labs(title = 'DBS',
       x = 'PCAWG DBS TMB',
       y = 'Hartwig DBS TMB') +
  theme_bw() +
  theme_bigfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.margin = margin(t= 5, l = 5, b = 5, r = 15))

#########################################################################################
### Indel ###############################################################################

# get the number of rows
indel_n_samples <- nrow(indel_cos_tot_mut_combined)

# build linear model Hartwig total burden explained by pcawg total mutational burden
indel_model <- lm(indel_hartwig_total_mutational_burden ~ indel_pcawg_total_mutational_burden, data = indel_cos_tot_mut_combined)

# get coefs of this model
indel_slope <- indel_model$coefficients[2]

indel_y_limit <- 15000

# plot this relationship and the regression line
indel_p1 <- ggplot(indel_cos_tot_mut_combined, aes(x = indel_pcawg_total_mutational_burden, y = indel_hartwig_total_mutational_burden)) +
  geom_point(color = '#C0C0C0') +
  geom_smooth(method = 'lm', se = FALSE, color = '#000000') +
  scale_x_continuous(labels = unit_format(unit = 'K', scale = 1e-3),
                     breaks = seq(0, indel_y_limit, by = 5000),
                     limits = c(0, indel_y_limit)) +
  scale_y_continuous(labels = unit_format(unit = 'K', scale = 1e-3),
                     breaks = seq(0, indel_y_limit, by = 5000),
                     limits = c(0, indel_y_limit)) +
  annotate('text', x = indel_y_limit * 0.25, y = indel_y_limit * 0.75, label = paste0('Slope: ', round(indel_slope, 2)),
           size = 5, colour = 'black', fontface = 2) +
  labs(title = 'INDEL',
       x = 'PCAWG Indel TMB',
       y = 'Hartwig Indel TMB') +
  theme_bw() +
  theme_bigfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.margin = margin(t= 5, l = 5, b = 5, r = 15))

#########################################################################################
### SV ##################################################################################

# get the number of rows
sv_n_samples <- nrow(sv_cos_tot_mut_combined)

# ------------------------------------ Model with length cutoff

# build linear model Hartwig total burden explained by pcawg total mutational burden
sv_model_len_cutoff <- lm(sv_hartwig_total_mutational_burden_len_cutoff ~ sv_pcawg_total_mutational_burden_len_cutoff, data = sv_cos_tot_mut_combined)

# get coefs of this model
sv_slope_len_cutoff <- sv_model_len_cutoff$coefficients[2]

sv_y_limit <- 3000

# plot this relationship and the regression line
sv_p1_len_cutoff <- ggplot(sv_cos_tot_mut_combined, aes(x = sv_pcawg_total_mutational_burden_len_cutoff, y = sv_hartwig_total_mutational_burden_len_cutoff)) +
  geom_point(color = '#C0C0C0') +
  geom_smooth(method = 'lm', se = FALSE, color = '#000000') +
  scale_x_continuous(labels = unit_format(unit = 'K', scale = 1e-3), 
                     limits = c(0, sv_y_limit)) +
  scale_y_continuous(labels = unit_format(unit = 'K', scale = 1e-3),
                     limits = c(0, sv_y_limit)) +
  annotate('text', x = sv_y_limit * 0.25, y = sv_y_limit * 0.75, label = paste0('Slope: ', round(sv_slope_len_cutoff, 2)),
           size = 5, colour = 'black', fontface = 2) +
  labs(title = paste0('SV length >', sv_length_cutoff, 'bp'),
       x = 'PCAWG SV TMB',
       y = 'Hartwig SV TMB') +
  theme_bw() +
  theme_bigfont +
  theme(
    axis.text = element_text(color = '#000000', family = 'Helvetica'),
    axis.line = element_line(color = '#000000'),
    plot.margin = margin(t= 5, l = 5, b = 5, r = 15))

# stitch plots together
slope_plots <- ggarrange(snv_p1,
                         dbs_p1,
                         indel_p1,
                         sv_p1_len_cutoff,
                         labels = NULL,
                         ncol = 4, nrow = 1,
                         align = "hv",
                         font.label = list(size = 8, color = "black", face = "bold", family = NULL, position = "top"))

slope_plots <- annotate_figure(slope_plots,
                               fig.lab = 'A', fig.lab.size = 24, fig.lab.face = 'bold')
