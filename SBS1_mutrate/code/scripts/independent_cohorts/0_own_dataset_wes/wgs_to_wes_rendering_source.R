
## Loading basic packages --------------------------------

library(readxl)
library(tidyr)
library(ggplot2)
library(stringr)
library(plyr)
library(dplyr)
library(magrittr)
library(mutSigExtractor)


## Define functions --------------------------------

`%notin%` <- Negate(`%in%`)


## Assign values --------------------------------

# Cancer type order for fig1 layout
cancer_type_order_fig1 <- c("Pan-cancer", "Glioblastoma multiforme", "Upper respiratory tract carcinoma", "Thyroid carcinoma",
                            "Lung adenocarcinoma", "Lung squamous cell carcinoma","Diffuse large B-cell lymphoma", "Breast carcinoma",
                            "Cholangiocarcinoma", "Hepatocellular carcinoma", "Pancreas carcinoma",
                            "Pancreas neuroendocrine", "Colorectal carcinoma", "Esophageal carcinoma", "Stomach carcinoma", 
                            "Kidney renal clear cell carcinoma", "Cervical carcinoma", "Ovarian serous adenocarcinoma", "Uterus carcinoma", 
                            "Bladder urothelial carcinoma", "Prostate carcinoma", "Leiomyosarcoma", "Liposarcoma", "Skin melanoma")

cancer_type_code_order_fig1 <- c("PAN", "GBM", "HNSC", "THCA", "LUAD", "LUSC", "DLBCL", "BRCA", "CHOL", "LIHC", "PAAD", "PANET", "COREAD", "ESCA",
                                 "STAD", "KIRC", "CESC", "OV", "UCEC", "BLCA", "PRAD", "LMS", "LPS", "SKCM")




# Cancer type order for the rest of the analysis + color palette
cancer_type_order <- c("Breast carcinoma", "Glioblastoma multiforme", "Colorectal carcinoma", "Esophageal carcinoma",
                       "Stomach carcinoma", "Cholangiocarcinoma", "Hepatocellular carcinoma", "Pancreas carcinoma",
                       "Pancreas neuroendocrine", "Cervical carcinoma", "Ovarian serous adenocarcinoma", "Uterus carcinoma",
                       "Upper respiratory tract carcinoma", "Kidney renal clear cell carcinoma", "Lung adenocarcinoma",
                       "Lung squamous cell carcinoma", "Diffuse large B-cell lymphoma", "Prostate carcinoma", "Skin melanoma",
                       "Leiomyosarcoma", "Liposarcoma", "Thyroid carcinoma", "Bladder urothelial carcinoma")

cancer_type_code_order <- c("BRCA", "GBM", "COREAD", "ESCA", "STAD", "CHOL", "LIHC", "PAAD", "PANET", "CESC", "OV",
                            "UCEC", "HNSC", "KIRC", "LUAD", "LUSC", "DLBCL", "PRAD", "SKCM", "LMS", "LPS", "THCA", "BLCA")


# Color palette for cancer types
cancer_type_palette <- c("#702963", "#808080", "#B4C424", "#98FB98", "#228B22", "#008080", "#000080", "#40E0D0", "#87CEEB", "#FFB6C1", "#FF00FF",
                         "#F33A6A", "#36454F", "#6F8FAF", "#A95C68", "#C38E96", "#8A9A5B", "#8B0000", "#FFC000", "#80461B", "#E49B0F", "#C2B280", "#EE4B2B")


# Tissue group order
included_tissue_group <- c("Breast", "CNS", "GI_core", "GI_dev", "Gyn", "Head_and_neck", "Kidney", "Lung", "Lymphoid",
                           "Prostate", "Skin", "Soft_tissue", "Thyroid", "Urothelial")



# Cohort order + cohort color palette
cohort_order <- c("PCAWG", "Hartwig")
cohort_color_palette <- c("#f58134", "#9966CC")

