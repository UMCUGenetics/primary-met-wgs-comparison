############### Clean actionable variants ###############
# author: remy (sascha)
# date: 12/01/2022
# last updated: 26/01/2022

### Description
# input: actionable variants table that was parsed from VCF files by parse_actionable_variants.R script.
# This script cleans the actionable variants table. It uses the amino acid code table to change the amino
# acid residue change from three letter code to single letter code to match the actionability table from Hartwig.

# libraries
library(tidyverse) # data manipulation and plotting
library(naniar) # replace with NA function

#========= Path prefixes =========#
base_dir <- list(
  hpc='/base/path',
  mnt='/base/path',
  umc='/base/path'
)
for(i in base_dir){
  if(dir.exists(i)){
    base_dir <- i
    break
  }
}

# define dates of the input files and other variables
actionable_variants_results_date <- '2022_01_26' # old: 2022_01_20, new: 2022_01_26
linx_drivers_resutls_date <- '2021_12_03'
fusion_results_date <- '2021_12_03'
purity_results_date <- '2021_12_03'
exclude_hypermutators <- FALSE

# --------------------------------------------- ACTIONABILITY

actionability <- read_tsv(file = paste0(base_dir, '/path/to/hartwig_actionability.tsv')) %>%
  # remove variants for which we dont have data for (expression and hypermethylation)
  filter(!mutated_variant %in% c('over exp', 'dec exp', 'hypermethylation'))

n_per_mutated_variant <- actionability %>%
  group_by(mutated_variant) %>%
  count() %>%
  arrange(desc(n))

# --------------------------------------------- METADATA (MSI high/low marker)

metadata <- read_tsv(paste0(base_dir, '/path/to/metadata.tsv')) %>%
  # select relevant columns
  select(sample_id, cohort, cancer_type, msi_status) %>%
  # change the msi_status column to match the actionability table
  mutate(msi_status = if_else(msi_status == 'MSI', 'MSI high', 'MSI neg'))

# --------------------------------------------- PURITY (TMB high/low marker)

purity <- read_tsv(paste0(base_dir, '/path/to/',
                          purity_results_date,'_purple.purity.txt.gz')) %>%
  select(sample, tmbPerMb) %>%
  mutate(tmb_status = case_when(tmbPerMb < 1 ~ 'TMB low',
                                tmbPerMb >= 1 & tmbPerMb < 10 ~ 'TMB normal',
                                tmbPerMb >= 10 ~ 'TMB high')) %>%
  select(sample, tmb_status)

# join this data to the metadata
metadata <- metadata %>%
  inner_join(purity, by = c('sample_id' = 'sample'))
rm(purity)

# --------------------------------------------- PARSED ACTONABLE VARIANTS

# read in linx drivers
actionable_variants <- read_tsv(paste0(base_dir, '/path/to/',
                                       actionable_variants_results_date ,'/actionable_variants_prefilter.tsv.gz'),
                         col_names = TRUE)

# clean the intron_exon_rank column
actionable_variants <- actionable_variants %>%
  mutate(intron_exon_rank = case_when(str_detect(var_type, pattern = 'insertion') ~ str_replace_all(intron_exon_rank,
                                                                                                    pattern = '([0-9]+)\\/[0-9]+',
                                                                                                    replacement = 'exon \\1 ins'),
                                      str_detect(var_type, pattern = 'deletion') ~ str_replace_all(intron_exon_rank,
                                                                                                   pattern = '([0-9]+)\\/[0-9]+',
                                                                                                   replacement = 'exon \\1 del'),
                                      TRUE ~ str_replace_all(intron_exon_rank,
                                                             pattern = '([0-9]+)\\/[0-9]+',
                                                             replacement = 'exon \\1 mut')))

# ----------------------------- LINX DRIVER

linx_drivers <- read_tsv(file = paste0(base_dir, '/path/to/',
                                       linx_drivers_resutls_date ,'/linx_drivers.tsv'))

# modify the driver column to match the actionability table
linx_drivers <- linx_drivers %>%
  mutate(driver = str_to_lower(driver),
         driver = if_else(driver == 'mutation' & category == 'TSG' & biallelic, 'del', driver),
         driver = if_else(driver == 'mutation', 'mut', driver)) %>%
  drop_na(driver)

# ----------------------------- FUSIONS

linx_fusions_high <- read_tsv(paste0(base_dir, '/path/to/',
                                     fusion_results_date ,'/linx_fusions.tsv'),
                         col_names = TRUE)

# left join the metadata table, rename 'name' column to 'gene'
linx_fusions_high <- linx_fusions_high %>%
  semi_join(metadata, by = 'sample_id') %>%
  rename(gene = 'name') %>%
  select(sample_id, gene)

########################## WILD-TYPE BIOMARKERS

actionability_wild_type <- actionability %>%
  filter(mutated_variant == 'wild_type') %>%
  distinct(gene, mutated_variant)

# create a complete wild-type markers df with all the samples
actionable_wild_type_marker <- expand_grid(sample_id = metadata$sample_id, gene = actionability_wild_type$gene)

# get the samples with a mutated wild-type marker gene
samples_with_mutated_wild_type_markers <- actionable_variants %>%
  semi_join(actionability_wild_type, by = 'gene') %>%
  distinct(sample_id, gene)

# now anti join this mutated wild type marker table to the complete table to get the sample which have a wild-type marker gene
actionable_wild_type_marker <- actionable_wild_type_marker %>%
  anti_join(samples_with_mutated_wild_type_markers, by = c('sample_id', 'gene')) %>%
  mutate(mutated_variant = 'wild_type')

########################## FUSION

actionability_fusions <- actionability %>%
  filter(mutated_variant == 'fusion') %>%
  filter(!str_detect(gene, pattern = '\\_')) %>%
  distinct(gene) %>% pull(gene)

start_list <- list()

for (i in 1:length(actionability_fusions)) {
  start_list[[i]] <- if_else(startsWith(linx_fusions_high$gene, actionability_fusions[i]),
                            str_replace(linx_fusions_high$gene, pattern = regex('^([A-Z0-9]+)\\_.*'), replacement = '\\1'),
                            'drop')
}

start_df <- as.data.frame(do.call(cbind, start_list))

end_list <- list()

for (i in 1:length(actionability_fusions)) {
  end_list[[i]] <- if_else(endsWith(linx_fusions_high$gene, actionability_fusions[i]),
                            str_replace(linx_fusions_high$gene, pattern = regex('.*\\_([A-Z0-9]+)$'), replacement = '\\1'),
                            'drop')
}

end_df <- as.data.frame(do.call(cbind, end_list))

linx_fusions_hightt <- cbind(linx_fusions_high, start_df, end_df)
colnames(linx_fusions_hightt) <- c('sample_id', 'gene', str_c('start_', 1:length(actionability_fusions)), str_c('end_', 1:length(actionability_fusions)))

linx_fusions_hightt <- linx_fusions_hightt %>%
  pivot_longer(., cols = 3:ncol(linx_fusions_hightt), names_to = 'names', values_to = 'gene2') %>%
  select(-names, -gene) %>%
  rename(gene = 'gene2') %>%
  filter(gene != 'drop') %>%
  mutate(mutated_variant = 'fusion')

actionable_fusions <- linx_fusions_high %>%
  semi_join(actionability, by = 'gene') %>%
  mutate(mutated_variant = 'fusion')

actionable_fusions <- actionable_fusions %>%
  union(., linx_fusions_hightt)

######################### PROTEIN RESIDUE CHANGES

## ARBITRARY

# create a second table that contains the unique arbitrary protein residue changes (e.g. BRAF V600E -> BRAF V600X)
arbitrary_actionable_variants <- actionable_variants %>%
  mutate(arbitrary_protein_residue_change = str_replace_all(protein_residue_change, pattern = '([A-Z][0-9]+)[A-Z]', replacement = '\\1X')) %>%
  filter(str_detect(arbitrary_protein_residue_change, pattern = 'X$')) %>%
  select(sample_id, gene, arbitrary_protein_residue_change) %>%
  distinct()

arbitrary_actionable_protein_residue_changes <- arbitrary_actionable_variants %>%
  semi_join(actionability, by = c('gene', 'arbitrary_protein_residue_change' = 'mutated_variant')) %>%
  select(sample_id, gene, arbitrary_protein_residue_change) %>%
  rename(mutated_variant = 'arbitrary_protein_residue_change')

## NON-ARBITRARY

actionable_protein_residue_changes1 <- actionable_variants %>%
  semi_join(actionability, by = c('gene', 'protein_residue_change' = 'mutated_variant')) %>%
  select(sample_id, gene, protein_residue_change) %>%
  rename(mutated_variant = 'protein_residue_change')

actionable_protein_residue_changes2 <- actionable_variants %>%
  semi_join(actionability, by = c('gene', 'intron_exon_rank' = 'mutated_variant')) %>%
  select(sample_id, gene, intron_exon_rank) %>%
  rename(mutated_variant = 'intron_exon_rank')

actionable_protein_residue_changes <- actionable_protein_residue_changes1 %>%
  union(., actionable_protein_residue_changes2)

######################### TMB / MSI BIOMARKERS

actionable_tmb_msi <- metadata %>%
  pivot_longer(., cols = c('msi_status', 'tmb_status'), names_to = 'names', values_to = 'mutated_variant') %>%
  filter(mutated_variant != 'TMB normal') %>%
  mutate(gene = mutated_variant) %>%
  select(sample_id, gene, mutated_variant)

######################### OTHER CHANGES

actionable_linx_driver <- linx_drivers %>%
  semi_join(actionability, by = c('gene', 'driver' = 'mutated_variant')) %>%
  select(sample_id, gene, driver) %>%
  rename(mutated_variant = 'driver')

########################## Combine all the different variants
actionable_variants_clean <- actionable_wild_type_marker %>%
  union(., arbitrary_actionable_protein_residue_changes) %>%
  union(., actionable_protein_residue_changes) %>%
  union(., actionable_fusions) %>%
  union(., actionable_tmb_msi) %>%
  union(., actionable_linx_driver)

write_tsv(actionable_variants_clean, file = paste0(base_dir, '/output/path/',
                                                   actionable_variants_results_date ,'/actionable_variants.tsv'))