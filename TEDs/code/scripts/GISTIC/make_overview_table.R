#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(reshape2)
library(glue)
library(ggplot2)
library(scales)
library(ggpubr)
library(rjson)
library(readxl)
library(GenomicRanges)
library(VariantAnnotation)

library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
chromosomes <- seqnames(get(ref_genome))[1:23]
library(data.table)

colorpalette=c("Carboplatin"="#000000",
               "Cisplatin"="#FFFF00",
               "Gemcitabine"="#1CE6FF",
               "Platinum"="#FF34FF",
               "Pyrimidine_antagonist"="#FF4A46",
               "untreated"="#008941",
               "Alkaloid"="#006FA6",
               "Alkylating"="#A30059",
               "Anastrozole"="#FFDBE5",
               "Anthracycline"="#7A4900",
               "Anti_HER2"="#0000A6",
               "Anti_VEGF"="#63FFAC",
               "Aromatase_inhibitor"="#B79762",
               "Bevacizumab"="#004D43",
               "Capecitabine"="#8FB0FF",
               "CDK4__6_inhibitor"="#997D87",
               "Cyclophosphamide"="#5A0007",
               "Docetaxel"="#809693",
               "Doxorubicin"="#FEFFE6",
               "Epirubicin"="#1B4400",
               "Eribulin"="#4FC601",
               "Everolimus"="#3B5DFF",
               "Exemestane"="#4A3B53",
               "Fluorouracil"="#FF2F80",
               "Folate_antagonist"="#61615A",
               "Fulvestrant"="#BA0900",
               "GnRH_antagonist"="#6B7900",
               "Goserelin"="#00C2A0",
               "Letrozole"="#FFAA92",
               "Leuprorelin"="#FF90C9",
               "Megestrol"="#B903AA",
               "Methotrexate"="#D16100",
               "Microtubule_inhibitor"="#DDEFFF",
               "mTOR_inhibitor"="#000035",
               "Paclitaxel"="#7B4F4B",
               "Palbociclib"="#A1C299",
               "Pertuzumab"="#300018",
               "Selective_ER_modulator"="#0AA6D8",
               "Tamoxifen"="#013349",
               "Taxane"="#00846F",
               "Trastuzumab"="#372101",
               "Vinorelbine"="#FFB500",
               "Anti_EGFR"="#C2FFED",
               "Irinotecan"="#A079BF",
               "Leucovorin"="#CC0744",
               "Oxaliplatin"="#C0B9B2",
               "Panitumumab"="#C2FF99",
               "Topoisomerase_inhibitor"="#001E09",
               "Multikinase_nhibitor"="#00489C",
               "Sunitinib"="#6F0062",
               "Anti_PD_1"="#0CBD66",
               "Erlotinib"="#EEC3FF",
               "Gefitinib"="#456D75",
               "Immunotherapy"="#B77B68",
               "Osimertinib"="#7A87A1",
               "Pemetrexed"="#788D66",
               "Doxorubicin_liposomal"="#885578",
               "Abiraterone"="#FAD09F",
               "Anti_AR"="#FF8A9A",
               "Bicalutamide"="#D157A0",
               "Cabazitaxel"="#BEC459",
               "Enzalutamide"="#456648",
               "Luteinizing_hormone_releasing_hormone_LHRH"="#0086ED",
               "Anti_CTLA_4"="#886F4C",
               "Ipilimumab"="#34362D",
               "Nivolumab"="#B4A8BD")

cancercolors=c("BLCA"='#8dd3c7',
               "BRCA"='#ffffb3',
               "CESC"='#bebada',
               "CHOL"='#fb8072',
               "COREAD"='#80b1d3',
               "DLBCL"='#fdb462',
               "ESCA"='#b3de69',
               "GBM"='#fccde5',
               "HNSC"='#d9d9d9',
               "KIRC"='#ff1417',
               "LIHC"='#ff6611',
               "NSCLC"='#c4ff00',
               "OV"='#ff8844',
               "PAAD"='#ffee55',
               "PRAD"='#ffff99',
               "SKCM"='#78FA37',
               "STAD"='#aacc22',
               "UCEC"='#bbdd77')


mechanisms <- c('Alkaloid',
                'Alkylating',
                'Anthracycline',
                'Anti_AR__GnRH',
                'Anti_CTLA_4',
                'Anti_EGFR',
                'Anti_HER2',
                'Anti_PD_1',
                'Anti_VEGF',
                'Aromatase_inhibitor',
                'CDK4__6_inhibitor',
                'Folate_antagonist',
                'GnRH_antagonist',
                'Immunotherapy',
                'Microtubule_inhibitor',
                'Multikinase_inhibitor',
                'Platinum',
                'Pyrimidine_antagonist',
                'Selective_ER_modulator',
                'Taxane',
                'Topoisomerase_inhibitor',
                'mTOR_inhibitor')


drugs <- c('Abiraterone',
  'Anastrozole',
  'Bevacizumab',
  'Bicalutamide',
  'Cabazitaxel',
  'Capecitabine',
  'Carboplatin',
  'Cisplatin',
  'Cyclophosphamide',
  'Degarelix',
  'Docetaxel',
  'Doxorubicin',
  'Doxorubicin_liposomal',
  'Enzalutamide',
  'Epirubicin',
  'Eribulin',
  'Erlotinib',
  'Etoposide',
  'Everolimus',
  'Exemestane',
  'Fluorouracil',
  'Fulvestrant',
  'Gefitinib',
  'Gemcitabine',
  'Goserelin',
  'Ipilimumab',
  'Irinotecan',
  'Letrozole',
  'Leucovorin',
  'Leuprorelin',
  'Luteinizing_hormone_releasing_hormone_LHRH',
  'Megestrol',
  'Methotrexate',
  'Nilutamide',
  'Nivolumab',
  'Osimertinib',
  'Oxaliplatin',
  'Paclitaxel',
  'Palbociclib',
  'Panitumumab',
  'Pembrolizumab',
  'Pemetrexed',
  'Pertuzumab',
  'Sunitinib',
  'Tamoxifen',
  'Trastuzumab',
  'Vinorelbine')

GISTIC_runs <- c('/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BLCA/GISTIC_BLCA_Alkaloid/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BLCA/GISTIC_BLCA_Platinum/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BLCA/GISTIC_BLCA_Pyrimidine_antagonist/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BLCA/GISTIC_BLCA_untreated/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-ERpos/GISTIC_BRCA-ERpos_Alkaloid/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-ERpos/GISTIC_BRCA-ERpos_Alkylating/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-ERpos/GISTIC_BRCA-ERpos_Anthracycline/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-ERpos/GISTIC_BRCA-ERpos_Anti_AR__GnRH/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-ERpos/GISTIC_BRCA-ERpos_Anti_HER2/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-ERpos/GISTIC_BRCA-ERpos_Anti_VEGF/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-ERpos/GISTIC_BRCA-ERpos_Aromatase_inhibitor/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-ERpos/GISTIC_BRCA-ERpos_CDK4__6_inhibitor/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-ERpos/GISTIC_BRCA-ERpos_Folate_antagonist/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-ERpos/GISTIC_BRCA-ERpos_Microtubule_inhibitor/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-ERpos/GISTIC_BRCA-ERpos_mTOR_inhibitor/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-ERpos/GISTIC_BRCA-ERpos_Platinum/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-ERpos/GISTIC_BRCA-ERpos_Pyrimidine_antagonist/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-ERpos/GISTIC_BRCA-ERpos_Selective_ER_modulator/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-ERpos/GISTIC_BRCA-ERpos_Taxane/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-ERpos/GISTIC_BRCA-ERpos_untreated/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA/GISTIC_BRCA_Alkaloid/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA/GISTIC_BRCA_Alkylating/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA/GISTIC_BRCA_Anthracycline/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA/GISTIC_BRCA_Anti_AR__GnRH/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA/GISTIC_BRCA_Anti_HER2/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA/GISTIC_BRCA_Anti_VEGF/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA/GISTIC_BRCA_Aromatase_inhibitor/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA/GISTIC_BRCA_CDK4__6_inhibitor/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA/GISTIC_BRCA_Folate_antagonist/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA/GISTIC_BRCA_Microtubule_inhibitor/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA/GISTIC_BRCA_mTOR_inhibitor/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA/GISTIC_BRCA_Platinum/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA/GISTIC_BRCA_Pyrimidine_antagonist/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA/GISTIC_BRCA_Selective_ER_modulator/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA/GISTIC_BRCA_Taxane/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA/GISTIC_BRCA_untreated/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-HER2pos/GISTIC_BRCA-HER2pos_Alkylating/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-HER2pos/GISTIC_BRCA-HER2pos_Anthracycline/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-HER2pos/GISTIC_BRCA-HER2pos_Anti_AR__GnRH/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-HER2pos/GISTIC_BRCA-HER2pos_Anti_HER2/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-HER2pos/GISTIC_BRCA-HER2pos_Aromatase_inhibitor/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-HER2pos/GISTIC_BRCA-HER2pos_Pyrimidine_antagonist/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-HER2pos/GISTIC_BRCA-HER2pos_Selective_ER_modulator/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-HER2pos/GISTIC_BRCA-HER2pos_Taxane/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-HER2pos/GISTIC_BRCA-HER2pos_untreated/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-TNB/GISTIC_BRCA-TNB_Alkylating/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-TNB/GISTIC_BRCA-TNB_Anthracycline/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-TNB/GISTIC_BRCA-TNB_Platinum/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-TNB/GISTIC_BRCA-TNB_Pyrimidine_antagonist/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-TNB/GISTIC_BRCA-TNB_Taxane/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/BRCA-TNB/GISTIC_BRCA-TNB_untreated/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/CESC/GISTIC_CESC_Platinum/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/CESC/GISTIC_CESC_Taxane/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/CESC/GISTIC_CESC_untreated/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/CHOL/GISTIC_CHOL_Platinum/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/CHOL/GISTIC_CHOL_Pyrimidine_antagonist/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/CHOL/GISTIC_CHOL_untreated/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/COREAD/GISTIC_COREAD_Anti_EGFR/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/COREAD/GISTIC_COREAD_Anti_VEGF/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/COREAD/GISTIC_COREAD_Platinum/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/COREAD/GISTIC_COREAD_Pyrimidine_antagonist/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/COREAD/GISTIC_COREAD_Topoisomerase_inhibitor/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/COREAD/GISTIC_COREAD_untreated/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/COREAD-MSI-POLE/GISTIC_COREAD-MSI-POLE_untreated/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/COREAD-MSS/GISTIC_COREAD-MSS_Anti_EGFR/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/COREAD-MSS/GISTIC_COREAD-MSS_Anti_VEGF/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/COREAD-MSS/GISTIC_COREAD-MSS_Platinum/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/COREAD-MSS/GISTIC_COREAD-MSS_Pyrimidine_antagonist/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/COREAD-MSS/GISTIC_COREAD-MSS_Topoisomerase_inhibitor/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/COREAD-MSS/GISTIC_COREAD-MSS_untreated/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/DLBCL/GISTIC_DLBCL_untreated/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/ESCA/GISTIC_ESCA_Platinum/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/ESCA/GISTIC_ESCA_Pyrimidine_antagonist/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/ESCA/GISTIC_ESCA_Taxane/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/ESCA/GISTIC_ESCA_untreated/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/GBM/GISTIC_GBM_Alkylating/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/GBM/GISTIC_GBM_untreated/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/HNSC/GISTIC_HNSC_Platinum/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/HNSC/GISTIC_HNSC_untreated/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/KIRC/GISTIC_KIRC_Multikinase_inhibitor/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/KIRC/GISTIC_KIRC_untreated/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/LIHC/GISTIC_LIHC_untreated/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/LMS/GISTIC_LMS_Anthracycline/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/LMS/GISTIC_LMS_untreated/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/LPS/GISTIC_LPS_untreated/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/LUAD/GISTIC_LUAD_Anti_EGFR/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/LUAD/GISTIC_LUAD_Folate_antagonist/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/LUAD/GISTIC_LUAD_Immunotherapy/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/LUAD/GISTIC_LUAD_Platinum/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/LUAD/GISTIC_LUAD_Taxane/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/LUAD/GISTIC_LUAD_untreated/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/LUSC/GISTIC_LUSC_Platinum/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/LUSC/GISTIC_LUSC_Pyrimidine_antagonist/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/LUSC/GISTIC_LUSC_untreated/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/OV/GISTIC_OV_Anthracycline/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/OV/GISTIC_OV_Anti_VEGF/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/OV/GISTIC_OV_Platinum/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/OV/GISTIC_OV_Pyrimidine_antagonist/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/OV/GISTIC_OV_Selective_ER_modulator/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/OV/GISTIC_OV_Taxane/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/OV/GISTIC_OV_untreated/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/PAAD/GISTIC_PAAD_Platinum/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/PAAD/GISTIC_PAAD_Pyrimidine_antagonist/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/PAAD/GISTIC_PAAD_Topoisomerase_inhibitor/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/PAAD/GISTIC_PAAD_untreated/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/PANET/GISTIC_PANET_untreated/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/PRAD/GISTIC_PRAD_Anti_AR__GnRH/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/PRAD/GISTIC_PRAD_Immunotherapy/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/PRAD/GISTIC_PRAD_Taxane/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/PRAD/GISTIC_PRAD_untreated/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/SKCM/GISTIC_SKCM_BRAF_inhibitor/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/SKCM/GISTIC_SKCM_Immunotherapy/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/SKCM/GISTIC_SKCM_MEK_inhibitor/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/SKCM/GISTIC_SKCM_untreated/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/STAD/GISTIC_STAD_Platinum/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/STAD/GISTIC_STAD_Pyrimidine_antagonist/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/STAD/GISTIC_STAD_untreated/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/THCA/GISTIC_THCA_untreated/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/UCEC/GISTIC_UCEC_Platinum/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/UCEC/GISTIC_UCEC_Taxane/',
                 '/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/UCEC/GISTIC_UCEC_untreated/')



google_clinical <- as.data.frame(readRDS(file = "/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/clinical_overview.rds"))
CN_genes_df <- as.data.frame(t(as.data.frame(readRDS(file = "/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/CN_gene_database/gene_CN_strict25_final.rds"))))
overview_all <- as.data.frame(read_excel("/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/Drivers_per_sample/ploidy_plus_2.5/filtered_peak_overview_ploidy_plus_2.5_qvalue05_final_NOV2022.xlsx"))
outdir <- "/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/Drivers_per_sample/ploidy_plus_2.5/Driver_per_sample_qvalue05_NOV2022/"


process_gistic <- function(GISTIC_repo="/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/PRAD/GISTIC_PRAD_Anti_AR__GnRH/"){
  
  Sampleoverview_final <- NULL
  
  ## Checks --------------------------------
  if(!dir.exists(GISTIC_repo)){
    stop("provide correct GISTIC_repo path")
  }else{
    #locate all purple and linx files
    GISTIC_files <- list.files(GISTIC_repo, recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
  }

  groupfolder = basename(GISTIC_repo) %>% gsub(pattern = "\\..*$",replacement =  "")
  stringsplitlength <- length(strsplit(groupfolder, '_')[[1]])
  
  
  if(stringsplitlength==3){
    Drug <- strsplit(groupfolder, '_')[[1]][3]
    cancerType <- strsplit(groupfolder, '_')[[1]][2]
  } else if(stringsplitlength==4){
    Drug <- sprintf("%s_%s",strsplit(groupfolder, '_')[[1]][3],strsplit(groupfolder, '_')[[1]][4])
    cancerType <- strsplit(groupfolder, '_')[[1]][2]
  }else if(stringsplitlength==5){
    Drug <- sprintf("%s_%s_%s",strsplit(groupfolder, '_')[[1]][3],strsplit(groupfolder, '_')[[1]][4],strsplit(groupfolder, '_')[[1]][5])
    cancerType <- strsplit(groupfolder, '_')[[1]][2]
  }else if(stringsplitlength==6){
    Drug <- sprintf("%s_%s_%s_%s",strsplit(groupfolder, '_')[[1]][3],strsplit(groupfolder, '_')[[1]][4],strsplit(groupfolder, '_')[[1]][5],strsplit(groupfolder, '_')[[1]][6])
    cancerType <- strsplit(groupfolder, '_')[[1]][2]
  }else if(stringsplitlength==7){
    Drug <- sprintf("%s_%s_%s_%s_%s",strsplit(groupfolder, '_')[[1]][3],strsplit(groupfolder, '_')[[1]][4],strsplit(groupfolder, '_')[[1]][5],strsplit(groupfolder, '_')[[1]][6],strsplit(groupfolder, '_')[[1]][7])
    cancerType <- strsplit(groupfolder, '_')[[1]][2]
  }
  

  pretreated_samples_hmf <- fromJSON(file = sprintf("/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/json_files/%s/%s.json",cancerType,Drug))
  nontreated_samples_PCAWG <- fromJSON(file = sprintf("/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/json_files/%s/untreated.json",cancerType))
  samples <- c(pretreated_samples_hmf,nontreated_samples_PCAWG)
  
  #proces AMP
  CN_genes_df_amp = CN_genes_df %>% tibble::rownames_to_column(var = "sample") %>%  tidyr::separate(sample, c("sampleId", "eventype"))
  CN_genes_df_amp = CN_genes_df_amp %>% dplyr::filter(eventype=="AMP")
  CN_genes_df_del = CN_genes_df %>% tibble::rownames_to_column(var = "sample") %>%  tidyr::separate(sample, c("sampleId", "eventype"))
  CN_genes_df_del = CN_genes_df_del %>% dplyr::filter(eventype=="DEL")
  head(CN_genes_df_del[1:3,1:3])
  
  df_sample_amp = setNames(data.frame(matrix(ncol = 1, nrow = length(samples))),c("sampleId"))
  df_sample_amp$sampleId <- samples

  #1) add known
  SampleoverviewknownAMP <- df_sample_amp
  SampleoverviewunknownAMP <- df_sample_amp
  SampleoverviewknownDEL <- df_sample_amp
  SampleoverviewunknownDEL <- df_sample_amp
  
  
  overview_all_drug_tumortype_known_amp <- overview_all %>% dplyr::filter(cancertype==cancerType,drug==Drug,q.value_TT_fdr<0.05,SV_type=="Amp",Known=="Yes")
  overview_all_drug_tumortype_unknown_amp <- overview_all %>% dplyr::filter(cancertype==cancerType,drug==Drug,q.value_TT_fdr<0.05,SV_type=="Amp",is.na(Known))
  overview_all_drug_tumortype_known_del <- overview_all %>% dplyr::filter(cancertype==cancerType,drug==Drug,q.value_TT_fdr<0.05,SV_type=="Del",Known=="Yes")
  overview_all_drug_tumortype_unknown_del <- overview_all %>% dplyr::filter(cancertype==cancerType,drug==Drug,q.value_TT_fdr<0.05,SV_type=="Del",is.na(Known)) 
  
  if(nrow(overview_all_drug_tumortype_known_amp)>0){
    CN_genes_df_amp_known <- CN_genes_df_amp %>% dplyr::select(sampleId,overview_all_drug_tumortype_known_amp$gene) 
    colnames(CN_genes_df_amp_known)[1] <- "sampleId"
    CN_genes_df_amp_known[CN_genes_df_amp_known==2] <- "AMP_known"
    SampleoverviewknownAMP <- left_join(df_sample_amp,CN_genes_df_amp_known,by="sampleId")
  }
  if(nrow(overview_all_drug_tumortype_unknown_amp)>0){
    CN_genes_df_amp_unknown <- CN_genes_df_amp %>% dplyr::select(sampleId,overview_all_drug_tumortype_unknown_amp$gene) 
    colnames(CN_genes_df_amp_unknown)[1] <- "sampleId"
    CN_genes_df_amp_unknown[CN_genes_df_amp_unknown==2] <- "AMP_unknown"
    SampleoverviewunknownAMP <- left_join(df_sample_amp,CN_genes_df_amp_unknown,by="sampleId")
  }
  SampleoverviewAMP <- left_join(SampleoverviewknownAMP,SampleoverviewunknownAMP,by="sampleId")

  if(nrow(overview_all_drug_tumortype_known_del)>0){
    CN_genes_df_del_known <- CN_genes_df_del %>% dplyr::select(sampleId,overview_all_drug_tumortype_known_del$gene) 
    colnames(CN_genes_df_del_known)[1] <- "sampleId"
    CN_genes_df_del_known[CN_genes_df_del_known==2] <- "DEL_known"
    SampleoverviewknownDEL <- left_join(df_sample_amp,CN_genes_df_del_known,by="sampleId")
  }
  if(nrow(overview_all_drug_tumortype_unknown_del)>0){
    CN_genes_df_del_unknown <- CN_genes_df_del %>% dplyr::select(sampleId,overview_all_drug_tumortype_unknown_del$gene) 
    colnames(CN_genes_df_del_unknown)[1] <- "sampleId"
    CN_genes_df_del_unknown[CN_genes_df_del_unknown==2] <- "DEL_unknown"
    SampleoverviewunknownDEL <- left_join(df_sample_amp,CN_genes_df_del_unknown,by="sampleId")
  }
  SampleoverviewDEL <- left_join(SampleoverviewknownDEL,SampleoverviewunknownDEL,by="sampleId")
  
  Sampleoverview <- left_join(SampleoverviewAMP,SampleoverviewDEL,by="sampleId")
  
  
  if(ncol(Sampleoverview)>1) {
    if(nrow(Sampleoverview)!=length(samples)) {
      print("CHECK")}
    
    df_sample_amp_subs = setNames(data.frame(matrix(ncol = 3, nrow = length(samples))),c("sampleId","drug","cancertype"))
    df_sample_amp_subs$sampleId <- samples
    df_sample_amp_subs <- df_sample_amp_subs %>% dplyr::mutate(treated = ifelse(grepl('^DO', sampleId),TRUE,FALSE))
    df_sample_amp_subs$drug <- Drug
    df_sample_amp_subs$cancertype <- cancerType
    Sampleoverview_final <- left_join(df_sample_amp_subs,Sampleoverview,by="sampleId")
  }else{
    Sampleoverview_final <- NULL
  }
  return(Sampleoverview_final)
}

for(GISTIC in GISTIC_runs){
  out_dir_cancertype_drug = basename(GISTIC) %>% gsub(pattern = "\\..*$",replacement =  "")
  is.not.null <- function(x) !is.null(x)
  print(GISTIC)
  out <- process_gistic(GISTIC_repo=GISTIC)
  print(out)
  if(!is.null(out)){
    write.table(out, file=sprintf("%s/Drivers_per_sample_%s.txt",outdir,out_dir_cancertype_drug),sep = "\t", col.names = TRUE,qmethod = "double", quote = FALSE,row.names = FALSE)
  }
}
  
