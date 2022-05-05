############### Clean driver clonality table ###############
# author: remy (sascha)
# date: 09/11/2021
# last updated: 09/11/2021

### Description
# This script reads the parsed driver clonality variants table in. the problem with this parsed table is that it still includes
# multiple variants per gene. this script finds the variant that was added to the LINX driver catalog. 
# This clean table is then written to disk.

# libraries
library(tidyverse) # data manipulation and plotting

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

# input variables
driver_clonality_date <- '2021_12_03'

# --------------------------------------------- DRIVER CLONALITY

# read in driver clonality
driver_clonality <- read_tsv(paste0(base_dir, '/path/to/', driver_clonality_date, '/driver_clonality.tsv')) %>%
  distinct()

# add new column tier2 that is PANEL if any variant in a driver gene shows tier == PANEL, 
# HOTSPOT for the same reason, or REST if the none of the variants fall into PANEL or HOTSPOT
driver_clonality <- driver_clonality %>%
  group_by(sample_id, gene) %>%
  mutate(tier2 = case_when(any(tier == 'PANEL') ~ 'PANEL',
                           any(tier == 'HOTSPOT') ~ 'HOTSPOT',
                           TRUE ~ 'REST')) %>%
  ungroup()

# filter the drivers where the two tier columns match, or they belong to the 'REST'
driver_clonality <- driver_clonality %>%
  filter(tier == tier2 | tier2 == 'REST')

# if all variants of a sample + gene combo are subclonal, then create a new column and label them all 'subclonal'
# if this is not the case then it is a 'clonal' variant
driver_clonality <- driver_clonality %>%
  group_by(sample_id, gene) %>%
  mutate(clonality2 = if_else(all(clonality == 'subclonal'), 'subclonal', 'clonal')) %>%
  ungroup()

# get rid of duplicate rows + add the type of driver
driver_clonality <- driver_clonality %>%
  distinct(sample_id, gene, clonality2) %>%
  mutate(driver = 'MUTATION')

# write clean driver clonality table to disk
write_tsv(driver_clonality, 
          file = paste0(base_dir, 
                        '/path/to/output/', 
                        driver_clonality_date, '/driver_clonality_clean.tsv'))
