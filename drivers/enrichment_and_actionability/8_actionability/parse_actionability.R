############### Investigate LINX driver output ###############
# author: remy (sascha)
# date: 16/12/2021
# last updated: 17/01/2021

### Description
# This script parses the table of actionable tumor variants and extracts (i) the name of the 
# cancer driver gene from the 'tumorProfile' column, then makes the table long format so each gene is one row.
# Wild type genes are extracted as well and moved to a different column called 'wild_type'.

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

# read in actionable variants
actionability <- read_tsv(paste0(base_dir, '/path/to/actionable_variants.tsv.gz'))

# get the cancer types names used in the meatadata table
source(paste0(base_dir, '/path/to/r_objects/cancer_type_order.R'))

## test del: what cancer types are there?
# actionability %>%
#   filter(str_detect(cancerType, pattern = '.*laryn.*')) %>% pull(cancerType) %>% unique()

## relabeling cancer types

# split the hierarchically higher cancer type categories like 'sarcoma' into 
# cancer types we use in this study (e.g. Leiomyosarcoma, Liposarcoma)
split_cancerType <- data.frame(cancerType = c(rep('adenocarcinoma', 16),
                                              rep('sarcoma', 2), 
                                              rep(c('gastroesophageal junction adenocarcinoma',
                                                    'gastroesophageal adenocarcinoma',
                                                    'gastroesophageal cancer'), 2)),
                               cancer_type = c('Breast cancer', 'Colorectum carcinoma', 'Esophagus cancer', 'Stomach cancer', 
                                               'Cholangiocarcinoma', 'Hepatocellular carcinoma', 'Pancreas carcinoma', 
                                               'Cervix carcinoma', 'Ovarian cancer', 'Uterus carcinoma', 
                                               'Kidney clear cell carcinoma', 'Non small cell lung cancer',
                                               'Prostate carcinoma', 'Skin melanoma', 'Thyroid cancer', 'Urothelial cancer',
                                               'Leiomyosarcoma', 'Liposarcoma', 
                                               rep('Stomach cancer', 3), rep('Esophagus cancer', 3)))

# complete the table
actionability <- actionability %>%
  left_join(split_cancerType, by = 'cancerType') %>%
  mutate(cancer_type = as.character(cancer_type)) %>%
  mutate(cancerType = if_else(is.na(cancer_type), cancerType, cancer_type)) %>%
  select(-cancer_type)

# relabel cancer types to match the metadata labels
actionability <- actionability %>%
  mutate(cancerType = fct_collapse(cancerType,
                                   'Breast cancer' = c('breast cancer',
                                                       'triple-receptor negative breast cancer',
                                                       'Her2-receptor positive breast cancer',
                                                       'estrogen-receptor positive breast cancer',
                                                       'Her2-receptor negative breast cancer',
                                                       'breast ductal carcinoma',
                                                       'estrogen-receptor negative breast cancer',
                                                       'breast carcinoma',
                                                       'progesterone-receptor positive breast cancer',
                                                       'breast lobular carcinoma',
                                                       'breast adenocarcinoma',
                                                       'inflammatory breast carcinoma',
                                                       'breast papillary carcinoma',
                                                       'breast metaplastic carcinoma'), 
                                   'Glioblastoma multiforme' = c('glioblastoma',
                                                                 'brain glioblastoma multiforme',
                                                                 'glioblastoma proneural subtype',
                                                                 'glioblastoma mesenchymal subtype'), 
                                   'Colorectum carcinoma' = c('colorectal cancer',
                                                              'colorectal adenocarcinoma',
                                                              'colon cancer',
                                                              'colon carcinoma',
                                                              'colon adenoma',
                                                              'colorectal carcinoma',
                                                              'colon adenocarcinoma',
                                                              'colorectal adenoma',
                                                              'large intestine cancer',
                                                              'rectum cancer',
                                                              'rectosigmoid cancer',
                                                              'rectum adenocarcinoma',
                                                              'anal canal cancer'), 
                                   'Esophagus cancer' = c('esophagus adenocarcinoma',
                                                          'esophagus squamous cell carcinoma',
                                                          'esophageal cancer',
                                                          'esophageal carcinoma'),
                                   'Stomach cancer' = c('gastric adenocarcinoma',
                                                        'gastric signet ring cell adenocarcinoma',
                                                        'stomach cancer',
                                                        'stomach carcinoma'), 
                                   'Cholangiocarcinoma' = c('intrahepatic cholangiocarcinoma',
                                                            'cholangiocarcinoma',
                                                            'biliary tract cancer',
                                                            'extrahepatic bile duct carcinoma',
                                                            'bile duct cancer'), 
                                   'Hepatocellular carcinoma' = c('hepatocellular carcinoma',
                                                                  'hepatocellular adenoma',
                                                                  'liver cancer',
                                                                  'liver carcinoma'), 
                                   'Pancreas carcinoma' = c('pancreatic cancer',
                                                            'pancreatic ductal adenocarcinoma',
                                                            'pancreatic adenocarcinoma',
                                                            'pancreatic carcinoma',
                                                            'pancreatic acinar cell adenocarcinoma',
                                                            'pancreatic ductal carcinoma',
                                                            'pancreatic adenosquamous carcinoma'),
                                   'Pancreas neuroendocrine' = c('pancreatic endocrine carcinoma',
                                                                 'neuroendocrine tumor',
                                                                 'neuroendocrine carcinoma'), 
                                   'Cervix carcinoma' = c('cervix carcinoma',
                                                          'cervical cancer',
                                                          'endocervical carcinoma',
                                                          'cervical adenocarcinoma',
                                                          'cervical squamous cell carcinoma'), 
                                   'Ovarian cancer' = c('ovarian cancer',
                                                        'ovary serous adenocarcinoma',
                                                        'ovarian clear cell adenocarcinoma',
                                                        'ovarian serous carcinoma',
                                                        'ovary epithelial cancer',
                                                        'ovarian carcinoma',
                                                        'ovarian clear cell carcinoma',
                                                        'ovarian mucinous neoplasm',
                                                        'endometrioid ovary carcinoma',
                                                        'endometrioid ovary carcinoma',
                                                        'ovary adenocarcinoma',
                                                        'small-cell carcinoma of the ovary of hypercalcemic type'),
                                   'Uterus carcinoma' = c('uterine cancer',
                                                          'endometrial cancer',
                                                          'endometrial adenocarcinoma',
                                                          'endometrial carcinoma',
                                                          'endometrial serous adenocarcinoma',
                                                          'endometrial mixed adenocarcinoma',
                                                          'endometrial clear cell adenocarcinoma',
                                                          'endometrioid ovary carcinoma'),
                                   'Upper respiratory tract cancer' = c('nasopharynx carcinoma',
                                                                        'hypopharynx cancer',
                                                                        'pharynx squamous cell carcinoma',
                                                                        'tongue squamous cell carcinoma',
                                                                        'nasal cavity cancer',
                                                                        'head and neck squamous cell carcinoma',
                                                                        'head and neck cancer',
                                                                        'head and neck carcinoma',
                                                                        'larynx cancer',
                                                                        'laryngeal carcinoma',
                                                                        'oral cavity cancer',
                                                                        'oral squamous cell carcinoma'), 
                                   'Kidney clear cell carcinoma' = c('kidney cancer',
                                                                     'renal cell carcinoma',
                                                                     'clear cell renal cell carcinoma',
                                                                     'renal carcinoma'), 
                                   'Non small cell lung cancer' = c('lung cancer',
                                                                    'lung non-small cell carcinoma',
                                                                    'lung adenocarcinoma',
                                                                    'lung squamous cell carcinoma',
                                                                    'lung papillary adenocarcinoma',
                                                                    'lung non-squamous non-small cell carcinoma',
                                                                    'lung carcinoma',
                                                                    'lung large cell carcinoma',
                                                                    'adenosquamous lung carcinoma',
                                                                    'mucinous lung adenocarcinoma'),
                                   'Diffuse large B-cell lymphoma' = c('diffuse large B-cell lymphoma',
                                                                       'marginal zone B-cell lymphoma',
                                                                       'B-cell lymphoma',
                                                                       'primary mediastinal B-cell lymphoma',
                                                                       'lymphoma'), 
                                   'Prostate carcinoma' = c('prostate cancer',
                                                            'prostate adenocarcinoma',
                                                            'prostate carcinoma',
                                                            'castration-resistant prostate carcinoma'), 
                                   'Skin melanoma' = c('skin cancer',
                                                       'skin melanoma',
                                                       'skin carcinoma',
                                                       'melanoma'), 
                                   'Leiomyosarcoma' = c('uterus leiomyosarcoma',
                                                        'leiomyosarcoma',
                                                        'uterine corpus myxoid leiomyosarcoma'),
                                   'Liposarcoma' = c('liposarcoma',
                                                     'myxoid liposarcoma',
                                                     'dedifferentiated liposarcoma',
                                                     'adult liposarcoma',
                                                     'pleomorphic liposarcoma'), 
                                   'Thyroid cancer' = c('thyroid gland medullary carcinoma',
                                                        'thyroid gland papillary carcinoma',
                                                        'thyroid gland follicular carcinoma',
                                                        'thyroid gland carcinoma',
                                                        'thyroid gland cancer',
                                                        'thyroid gland anaplastic carcinoma',
                                                        'thyroid gland Hurthle cell carcinoma',
                                                        'differentiated thyroid gland carcinoma'), 
                                   'Urothelial cancer' = c('urinary bladder cancer',
                                                           'urinary system cancer',
                                                           'bladder carcinoma',
                                                           'bladder urothelial carcinoma',
                                                           'bladder papillary transitional cell neoplasm',
                                                           'invasive bladder transitional cell carcinoma',
                                                           'bladder transitional cell papilloma',
                                                           'bladder carcinoma in situ')))

## filtering
actionability <- actionability %>%
  # exclude 'Advanced Solid Tumor' entries
  filter(cancerType != 'Advanced Solid Tumor') %>%
  # only 'actionable' drug target combinations
  filter(evidenceType == 'Actionable') %>%
  # only 'sensitive' and 'resistant' actionable targets
  filter(responseType %in% c('sensitive', 'predicted - sensitive', 'resistant', 'resistand - predicted')) %>%
  # only include Tier A (PDA approved) and B (late clinical)
  filter(ampCapAscoEvidenceLevel %in% c('A', 'B'))

######################################### VARIANT LEVEL ANALYSIS ######################################### 

## handling wild-type genes

# add a new column wild_type
actionabilitytt <- actionability %>%
  # extract all the values in 'genes' that have the suffix 'wild-type' at least once
  mutate(wild_type = str_extract_all(tumorProfile, pattern = regex('[A-Z0-9]+\\swild-type.*')),
         .after = 'tumorProfile') %>%
  # cast list column as character
  mutate(wild_type = as.character(wild_type)) %>%
  # replace empty character strings with NA
  mutate(wild_type = if_else(wild_type == 'character(0)', NA_character_, wild_type))

# parsing of the 'wild_type' column
actionabilitytt <- actionabilitytt %>%
  # remove all non-wild-type genes at the end of the string (e.g. 'BRAF wild-type KRAS' -> 'BRAF wild-type')
  mutate(wild_type = str_remove(wild_type, pattern = regex('\\s[A-Z0-9]+(\\s[A-Z0-9]+)?$'))) %>%
  # remove the 'wild-type' label after each gene (e.g. 'BRCA1 wild-type' -> 'BRCA1')
  mutate(wild_type = str_replace_all(wild_type, pattern = regex('.*?(\\s?[A-Z0-9]+)\\swild-type(\\s([A-Z0-9]+)\\swild-type)?.*?'),
                                           replacement = '\\1 \\3')) %>%
  # replace double white space with single white space
  mutate(wild_type = str_replace_all(wild_type, pattern = regex('(\\s)\\s'),
                                           replacement = '\\1')) %>%
  # remove the non-wild-type genes at the end of the string (e.g. 'BRCA1 RAD50 loss' -> 'BRCA1')
  mutate(wild_type = str_remove_all(wild_type, pattern = regex('\\s[A-Z0-9]+\\s[a-z]+(\\s[a-z]+)?$'))) %>%
  # trim white spaces around genes
  mutate(wild_type = str_trim(wild_type)) %>%
  # replace the NAs in the 'wild_type' column (prevents them from being deleted later)
  replace_na(., replace = list(wild_type = 'none')) %>%
  # split up the 'wild_type' into 4 columns
  separate(., col = 'wild_type', into = str_c('wild_type_', c(1:4))) %>%
  # then make the table longer based on the newly created 4 columns
  pivot_longer(., cols = str_c('wild_type_', c(1:4)),names_to = 'names', values_to = 'wild_type') %>%
  # drop 'names' column
  select(1:4, wild_type, everything(), -names) %>%
  # eliminate NA rows in long 'wild_type' column
  filter(!is.na(wild_type))

## handling mutated genes

# add new column 'mutated'
actionabilitytt2 <- actionabilitytt %>%
  # ------------------------------------ Amino acid residue changes (sub/del/ins)
  # remove the simple and complex protein residue substitution strings (e.g. A750G, but also E23Vfs*17 & V600E/K)
  mutate(mutated = str_replace_all(tumorProfile, pattern = regex('\\s([A-Z][0-9]+[A-Z])(\\/[A-Z])?'), replacement = ':\\1'),
         .after = tumorProfile) %>%
  mutate(mutated = str_replace_all(mutated, pattern = regex('\\s([A-Z][0-9]+([A-Z])?fs(\\*[0-9]+)?)'), replacement = ':\\1')) %>%
  # remove protein residue deletion, duplication and insertion strings 
  # (e.g. A750_P755del, S622del, Y772_A775dup, H773_V774insGHPH, G776delinsIC or E709_T710delinsD)
  mutate(mutated = str_replace_all(mutated, pattern = regex('\\s(([A-Z0-9]+\\_)?[A-Z0-9]+del([a-zA-Z]+)?)'), replacement = ':\\1')) %>%
  mutate(mutated = str_replace_all(mutated, pattern = regex('\\s([A-Z0-9]+\\_.*(ins|dup|splice)([a-zA-Z]+)?)'), replacement = ':\\1')) %>%
  # ------------------------------------ fusions
  # fuse fusions together: adopt LINX annotation (e.g. EGFR - ERG -> EGFR_ERG)
  mutate(mutated = str_replace_all(mutated, pattern = regex('\\s\\-\\s'), replacement = '_')) %>%
  # handling self-fusions & rearrangements: adopt LINX annotation 
  # (e.g. 'FGFR3 fusion' -> 'FGFR3_FGFR3' or 'FGFR3 rearrange' -> 'FGFR3_FGFR3')
  mutate(mutated = str_replace_all(mutated, pattern = regex('rearrange'), replacement = 'fusion')) %>%
  # ------------------------------------ wild-type
  # extract all the values in 'genes' that have the suffix 'wild-type' at least once
  mutate(mutated = str_remove_all(mutated, pattern = regex('[A-Z0-9]+\\swild-type'))) %>%
  # get rid of excess white space
  mutate(mutated = str_squish(mutated)) %>%
  # mark empty cells in 'mutated' as 'none'
  mutate(mutated = if_else(mutated == '', 'none', mutated)) %>%
  # ------------------------------------ Exon alterations
  #  handling one exon typo (ckbEntryId: 75485)
  mutate(mutated = str_replace_all(mutated, pattern = regex('(.*?)del exon14(.*?)'),
                                         replacement = '\\1exon 14 del\\2')) %>%
  # handling cases where exon number is fused with exon (e.g. 'exon11' -> 'exon 11 mut')
  mutate(mutated = str_replace_all(mutated, pattern = '(exon)([0-9]+)',
                                   replacement = '\\1 \\2 mut')) %>%
  # handling exon annotation: remove the annotation (e.g. 'FLT3 exon 14 ins' -> FLT3)
  mutate(mutated = str_replace_all(mutated, pattern = regex('([A-Z0-9]+)\\s(exon)\\s([0-9]+)\\s(del|ins|mut)'),
                                         replacement = '\\1:\\2:\\3:\\4')) %>%
  # ------------------------------------ all the other cancer driver gene alterations (neg/pos/mut/del/over_exp etc.)
  # remove the remaining cancer gene alteration annotations to keep only gene names 
  # (e.g. 'PTEN negative' -> 'PTEN', 'FGFR1 act mut' -> 'FGFR1', 'CD274 over exp KRAS mut' -> 'CD274 KRAS')
  mutate(mutated = str_replace_all(mutated, pattern = regex('(\\s?[A-Z0-9]+)\\s([a-z0-9]+)(\\s([a-z0-9]+)?)?'),
                                         replacement = '\\1:\\2\\3')) %>%
  mutate(mutated = str_replace_all(mutated, pattern = regex('(\\:[a-z]+)\\s([a-z]+)'),
                                         replacement = '\\1:\\2')) %>%
  # split up the 'mutated' into 4 columns
  separate(., col = 'mutated', into = str_c('mutated_', c(1:4)), sep = ' ') %>%
  # then make the table longer based on the newly created 4 columns
  pivot_longer(., cols = str_c('mutated_', c(1:4)), names_to = 'names', values_to = 'mutated') %>%
  # drop 'names' column
  select(1:4, mutated, everything(), -names) %>%
  # eliminate NA rows in long 'mutated' column
  filter(!is.na(mutated)) %>%
  distinct()

# handle TMB and MSI cases
actionabilitytt3 <- actionabilitytt2 %>%
  # unify the labels (e.g. 'MSI:negative' -> 'MSI:neg')
  mutate(mutated = str_replace_all(mutated, pattern = regex('^MSI:negative|^MSI:low'),
                                   replacement = 'MSI:neg')) %>%
  # replace the ':' in MSI and TMB strings, so they wont get separated downstream like the other variants
  mutate(mutated = if_else(str_detect(mutated, pattern = regex('^MSI|^TMB')), 
                           str_replace(mutated, pattern = ':', replacement = ' '), mutated))

# handle NAs and harmonize some labels
actionabilitytt4 <- actionabilitytt3 %>%
  # replace the first dash ':' with a plus '+'
  mutate(mutated = str_replace(mutated, pattern = regex(':'),
                                         replacement = '+')) %>%
  # split 'mutated' into 'mutated' & 'mutated_variant' by using the + as a separator
  separate(., col = 'mutated', into = c('mutated', 'mutated_variant'), sep = '\\+') %>%
  # replace all ':' in mutated_variant with ' '
  mutate(mutated_variant = str_replace_all(mutated_variant, pattern = regex(':'),
                                           replacement = ' ')) %>%
  # deleted duplicate rows based on the following columns:
  distinct(., ckbEntryId, tumorProfile, mutated, mutated_variant, wild_type, treatment, cancerType, ampCapAscoEvidenceLevel,
           .keep_all =  TRUE) %>%
  # fill the NA values in 'mutated_variant' with 'wild_type' if mutated is 'none', 'MSI ...'/'TMB ...' or else 'fusion'
  mutate(mutated_variant = case_when(mutated == 'none' ~ 'wild_type',
                                     str_detect(mutated, pattern = regex('^MSI|^TMB')) ~ mutated,
                                     is.na(mutated_variant) ~ 'fusion',
                                     TRUE ~ mutated_variant)) %>%
  # include wild_type variant genes in 'gene' column
  mutate(mutated = if_else(mutated == 'none', wild_type, mutated)) %>%
  # add mutated genes that are marked as 'pos' or 'positive' to the wild type column
  mutate(wild_type = if_else(str_detect(mutated_variant, pattern = regex('^pos')), 'wild_type', wild_type)) %>%
  mutate(mutated_variant = if_else(str_detect(mutated_variant, pattern = regex('^pos')), 'wild_type', mutated_variant)) %>%
  # collapse similar factors into 'del' and 'mut'
  mutate(mutated_variant = fct_collapse(mutated_variant,
                                        'mut' = c('mut', 'mutant', 'inact mut', 'act mut'),
                                        'del' = c('del', 'neg', 'negative', 'loss'))) %>%
  # rename column 'mutated' to 'gene'
  rename(gene = 'mutated') %>%
  select(1:4, gene, mutated_variant, everything(), -wild_type) %>%
  distinct()

# write this clean table to disk
write_tsv(actionabilitytt4,
          file = paste0(base_dir, '/output/path/hartwig_actionability.tsv'))

##################################################################################################