library("readr")
library(devtools)
library(tidyr)
library(dplyr)
library(rjson)


args = commandArgs(trailingOnly = TRUE)
input_csv_file <- fromJSON(file = args[1])
output_csv_file <- args[2]
samples <- input_csv_file
CN_files <-read.table(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/CopyNumber_links.txt",sep = "\t",header = TRUE)
CN_files <- CN_files[CN_files$SampleID %in% samples, ]
if(nrow(CN_files)!=length(samples)){
  stop("nrow(CN_files)!=length(samples)")
}

print("test")
google_clinical <- as.data.frame(readRDS(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/clinical_overview.rds"))
print("test")
metadata_hmf <- read.csv2("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/cancer_types_HMF_PCAWG_2_metadata.csv",header = TRUE,sep = ",")
print("test")

head(CN_files$path)

datalist = lapply(as.character(CN_files$path), function(x){read.table(file=x,sep="\t", header=T)})
names(datalist) <- CN_files$SampleID
driver_data_2 = do.call(rbind, datalist)
driver_data_2 <- driver_data_2 %>% tibble::rownames_to_column(var = "sample_id") %>% tidyr::separate(sample_id, c("sample_id", "mut_nr") )
driver_data_2 <- left_join(driver_data_2,metadata_hmf[,c("sample_id","gender")],by = "sample_id")
regions <- driver_data_2 %>% dplyr::select(sample_id,chromosome,copyNumber,start,end,baf,gender)

# Fix absolute CN below zero.
regions$copyNumber <- base::ifelse(regions$copyNumber < 0, 0, regions$copyNumber)

unique(regions$gender)
# Add a copy-number to the sex-chromosomes (so a log2(1) -1 equals 0 instead of -1 for the X and Y chromosomes)
regions$copyNumber[regions$gender == 'MALE' & (regions$chromosome == "X" | regions$chromosome == "Y")] <- regions$copyNumber[regions$gender == 'MALE' & (regions$chromosome == "X" | regions$chromosome == "Y")] + 1

# Convert to GISTIC format.
gistic.regions <- data.frame(sample = regions$sample_id, chr = regions$chromosome, start = regions$start, end = regions$end, bafCount = regions$baf, copyNumber.log2 = log2(regions$copyNumber) - 1)

gistic.regions$copyNumber.log2 <- base::ifelse(gistic.regions$copyNumber.log2 <= -15, -10, gistic.regions$copyNumber.log2)


#print(gistic.regions)
write_tsv(gistic.regions, output_csv_file)