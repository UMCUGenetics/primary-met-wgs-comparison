############### Load libraries ######################
### General
if (!require("devtools")) install.packages("devtools"); library(devtools) # installing from other sources
if (!require("googlesheets4")) install.packages("googlesheets4"); library(googlesheets4) # interact with google sheets
### Data manipulation
if (!require("naniar")) install.packages("naniar"); library(naniar) # cleaning data
if (!require("tidyr")) install.packages("tidyr"); library(tidyr) # tidy model output
if (!require("broom")) install.packages("broom"); library(broom) # handling NAs
if (!require("mutSigExtractor")) devtools::install_github('https://github.com/UMCUGenetics/mutSigExtractor/'); library(mutSigExtractor) # reading and manipulating VCF files
if (!require("matrixTests")) install.packages("matrixTests"); library(matrixTests) # needed for mutSigExtractor
if (!require("tidyverse")) install.packages("tidyverse"); library(tidyverse) # data manipulation and plotting
### Plotting
if (!require("ggthemes")) install.packages("ggthemes"); library(ggthemes) # plot theme
if (!require("ggrepel")) install.packages("ggrepel"); library(ggrepel) # handling labels
if (!require("ggpubr")) install.packages("ggpubr"); library(ggpubr) # plot layouts
if (!require("patchwork")) install.packages("patchwork"); library(patchwork) # plot layouts
if (!require("ggnewscale")) install.packages("ggnewscale"); library(ggnewscale) # multiple color scales
if (!require("gt")) install.packages("gt"); library(gt) # pretty tables
if (!require("sysfonts")) install.packages("sysfonts"); library(sysfonts) # google fonts
