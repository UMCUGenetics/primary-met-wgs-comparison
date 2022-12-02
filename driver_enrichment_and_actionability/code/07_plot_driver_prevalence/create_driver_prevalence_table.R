############### Create driver prevalence table ###############
# author: remy (sascha)
# date: 23/11/2021
# last updated: 24/10/2022

### Description
# This script creates a table that includes all reported LINX drivers and fusions. Input files are combined driver and fusion tables.
## For linx drivers, the point mutation drivers are split up according to their subtype (missense, nonsense, splice, inframe, frameshift; more than one can apply).
# Additionally, a mutation is classified as 'multihit' if both alleles are affected or calssified as 'promoter' if none of the criteria apply.
# Finally, several subtypes are collapsed into one: amplification = 'AMP' & 'PARTIAL_AMP', indel = 'FRAMESHIFT' & 'INFRAME'.
## The fusions are split up as well. The following steps were applied to simplify fusion gene names:
# 1. if reportedType == 'EXON_DEL_DUP': same gene; fusion name is simplified to only the gene involved
# 2. if reportedType == 'IG_KNOWN_PAIR': immunoglobulin is always the enhancer; fusion is renamed to IGH fusion partner
# 3. if reportedType == 'KNOWN_PAIR': look up the activated cancer driver gene(s) at https://www.intogen.org/search.
# 3. Should both genes be annotated as drivers, then one fusion is split up into two separate rows. If there is only one then only the driver gene name is retained. If there is none then this fusion is filtered out.
# 4. if reportedType == 'PROMISCUOUS_3': gene at the 3' end (reported as: geneEnd) is the activated gene; use gene name in 'geneEnd'
# 5. if reportedType == 'PROMISCUOUS_5': gene at the 5' end (reported as: geneStart) is the activated gene; use gene name in 'geneStart'
# 6. if reportedType == 'PROMISCUOUS_BOTH': both are activated; split up into two rows, one for 3' and one for 5' fusion partner.

### Input
# metadata.tsv
# filtered LINX drivers tabel from step 2
# LINX fusions table from step 2
# fusion lookup table which contains the activated cancer driver gene in the fusion pair according to https://www.intogen.org/search
# filtered driver clonality table from step 3

### Output 
# a driver prevalence table among the cohort

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

# define dates of the input files
driver_results_date <- '2021_12_03' # prefix of your generated LINX driver table
fusion_results_date <- '2021_12_03' # prefix of your generated fusion driver table
driver_clonality_date <- '2021_12_03' # prefix of your generated LINX driver clonality table

# ------------------------------ METADATA

metadata <- read_tsv(paste0(base_dir, '/data/processed/metadata/metadata.tsv'),
                     col_types = 'icccccccccccccccllllcccdiiildclccllllllliiiiiiiiiii') %>%
  # select relevant columns
  select(patient_id, sample_id, cohort, cancer_type, cancer_type_code, 
         cancer_subtype, progression_status_code, is_blacklisted_subtype, metastatic_location)

# calculate sample size
sample_size <- metadata %>%
  group_by(cancer_type, cohort) %>%
  summarise(sample_size = n(), .groups = 'drop') %>%
  pivot_wider(names_from = 'cohort', values_from = 'sample_size') %>%
  mutate(cancer_type_label = paste0(cancer_type, ' (', PCAWG, ' vs. ', Hartwig, ')'))

# ------------------------------ DRIVERS

linx_drivers <- read_tsv(paste0(base_dir, '/results/02_combine_driver_catalogs/', driver_results_date ,'/driver/results/linx_drivers.tsv'),
                         col_names = TRUE)

# ----------------------------- FUSIONS

linx_fusions_high <- read_tsv(paste0(base_dir, '/results/02_combine_driver_catalogs/', driver_results_date ,'/driver/results/linx_fusions.tsv'),
                         col_names = TRUE) %>%
  filter(reported == TRUE)

# prepare the linx_fusions data for the stacked barplot: which one is the activated cancer driver gene?
linx_fusions_high_list <- linx_fusions_high %>%
  group_by(reportedType) %>%
  group_split()
# 1. if reportedType == 'EXON_DEL_DUP': same gene; delete everything after the first gene
linx_fusions_high_list[[1]] <- linx_fusions_high_list[[1]] %>%
  mutate(name = geneStart)
# 2. if reportedType == 'IG_KNOWN_PAIR': always immunoglobulin that is the enhancer; rename to IGH fusion partner
linx_fusions_high_list[[2]] <- linx_fusions_high_list[[2]] %>%
  mutate(name = str_replace(name, pattern = '^.+\\_', replacement = ''))
# 3. if reportedType == 'KNOWN_PAIR': look up the activated cancer driver gene at https://www.intogen.org/search
fusion_lookup_table <- read_tsv(file = paste0(base_dir, '/data/external/fusion_lookup_table.tsv'))

tmp_tbl <- linx_fusions_high_list[[3]] %>%
  inner_join(fusion_lookup_table, by = 'name')

# split the fusions where cancer_driver == 'both'
tmp_tbl_both <- tmp_tbl %>%
  filter(cancer_driver == 'both') %>%
  separate(col = 'name', into = c('name1', 'name2'), sep = '_') %>%
  pivot_longer(., cols = c('name1', 'name2'), names_to = 'name_drop', values_to = 'name') %>%
  select(-cancer_driver, -n, -name_drop)

tmp_tbl_single <- tmp_tbl %>%
  filter(!cancer_driver %in% c('both', 'none')) %>%
  mutate(name = cancer_driver) %>%
  select(-cancer_driver, -n)

# combine both cleaned tmp tbls
linx_fusions_high_list[[3]] <- rbind(tmp_tbl_both, tmp_tbl_single)
rm(tmp_tbl, tmp_tbl_both, tmp_tbl_single)
# 4. if reportedType == 'PROMISCUOUS_3': gene at the geneEnd is the activated gene; use gene name in 'geneEnd'
linx_fusions_high_list[[4]] <- linx_fusions_high_list[[4]] %>%
  mutate(name = geneEnd)
# 5. if reportedType == 'PROMISCUOUS_5': gene at the geneStart is the activated gene; use gene name in 'geneStart'
linx_fusions_high_list[[5]] <- linx_fusions_high_list[[5]] %>%
  mutate(name = geneStart)
# 6. if reportedType == 'PROMISCUOUS_BOTH': both are activated; use both
linx_fusions_high_list[[6]] <- linx_fusions_high_list[[6]] %>%
  separate(col = 'name', into = c('5_prime_gene', '3_prime_gene'), sep = '_') %>%
  pivot_longer(., cols = c('5_prime_gene', '3_prime_gene'), names_to = 'prime', values_to = 'name') %>%
  select(-prime)

# combine the list into a df again
linx_fusions_hightt <- do.call(rbind, linx_fusions_high_list)

# add category == 'ONCO' for these fusions
linx_fusions_hightt <- linx_fusions_hightt %>%
  mutate(category = 'ONCO')

# --------------------------- DRIVER CLONALITY (point mutations only)

driver_clonality <- read_tsv(paste0(base_dir, '/processed/linx_driver/03_parse_driver_clonality/results/', driver_clonality_date, '/driver_clonality_clean.tsv')) %>%
  distinct()

# -------------------------- PREPARE DATA

# left join metadata to get additional information (necessary for next step)
linx_drivers_high <- metadata %>%
  inner_join(linx_drivers_high, by = 'sample_id')

# add a new column 'driver2' which basically splits up the MUTATION drivers into mutation type drivers (e.g. MISSENSE, NONSENSE etc.)
linx_drivers_hightt <- linx_drivers_high %>%
  # reduce excess in certain mutations (e.g. missense) to either present / not present (1 / 0) 
  mutate(across(missense:frameshift, ~if_else(.x > 0, 1, 0))) %>%
  # add up all the different point and indel mutations to create total_presence (helper column)
  mutate(total_presence = missense + nonsense + splice + inframe + frameshift) %>%
  # then make the table longer for easier manipulation
  pivot_longer(., cols = missense:frameshift, names_to = 'mutation_drivers', values_to = 'presence') %>%
  # make mutation_drivers upper case
  mutate(mutation_drivers = toupper(mutation_drivers)) %>%
  # filter out those rows where total_presence > 0 & presence == 0 
  # (meaning there is another mutation type active, but not this one so it can be dropped)
  filter(!(total_presence > 0 & presence == 0)) %>% ### NOTE: DONT FORGET TO WRAP EXPRESSION IN BRACKETS ()
  # finally add the new 'driver' column which shows MUTATION driver subtypes in place of MUTATION
  mutate(driver = case_when(driver != 'MUTATION' ~ driver,
                             driver == 'MUTATION' & presence == 1 ~ mutation_drivers,
                             driver == 'MUTATION' & presence == 0 & biallelic ~ 'MULTIHIT',
                             driver == 'MUTATION' & presence == 0 & biallelic == FALSE ~ 'PROMOTER',
                             TRUE ~ 'nix')) %>%
  # drop the columns used for the manipulation and then drop duplicate rows, then update
  select(-mutation_drivers, -presence, -total_presence) %>%
  distinct()

# --- add the fusions to the drivers

# rename 'name' column to 'gene' in clean fusion table
linx_fusions_hightt <- linx_fusions_hightt %>%
  rename(gene = 'name') %>%
  mutate(driver = 'FUSION')

# get the names of columns that are the same in linx drivers and fusions table
common_column_names <- names(linx_drivers_hightt)[names(linx_drivers_hightt) %in% names(linx_fusions_hightt)]

# filter linx drivers table based on the character vector above
linx_drivers_hightt <- linx_drivers_hightt %>%
  select(common_column_names)
linx_fusions_high_plot <- linx_fusions_hightt %>%
  select(common_column_names)

# combine fusions and the other drivers
linx_drivers_hightt <- linx_drivers_hightt %>%
  union(., linx_fusions_high_plot)

# get sample size per cancer type in long format (necessary for joining later)
sample_size <- sample_size %>%
  pivot_longer(., cols = c('PCAWG', 'Hartwig'), names_to = 'cohort', values_to = 'sample_size')

# get cohort sample size
cohort_size <- sample_size %>%
  group_by(cohort) %>%
  summarise(cohort_size = sum(sample_size), .groups = 'drop')

# join metadata, sample size and cohort size to the drivers table
linx_drivers_hightt <- linx_drivers_hightt %>%
  right_join(metadata, by = 'sample_id') %>%
  inner_join(sample_size, by = c('cancer_type', 'cohort')) %>%
  inner_join(cohort_size, by = 'cohort')

# create output dir
dir.create(path = paste0(base_dir, '/results/07_plot_driver_prevalence'), recursive = TRUE)

# write table to disk
write_tsv(linx_drivers_hightt, file = paste0(base_dir, '/results/07_plot_driver_prevalence/driver_prevalence.tsv'))
