
## In this script we retrieve the mutation and patient information (their positions and reference/alternative bases) of enriched non-coding drivers of interests


## Init --------------------------------

expanding_path <- '/home/ali313/Documents/studies/master/umc-project'

path_prefix <- c(
  hpc  ='/hpc/',
  local = paste0(expanding_path, '/hpc/')
)
path_prefix <- path_prefix[min(which(dir.exists(path_prefix)))]
path_prefix <- sub('/hpc/','',path_prefix)

base_dir <- paste0(path_prefix,'/hpc/cuppen/projects/P0025_PCAWG_HMF/')


## Read in the significantly mutated elements that we are interested in retrieving patient and mutation info for  --------------------------------
treatment_sig_hits <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external-files/significant_hits_noncoding.tsv"),
                               header =T, stringsAsFactors = F, sep = "\t")


# determine the corresponding mutation file and element ranges of the significantly mutated elements of interest
for (i in 1:nrow(treatment_sig_hits)){
  
  element_id <- treatment_sig_hits[i,"id"]
  cancer_type <- treatment_sig_hits[i,"cancer_type_code"]
  treatment <- treatment_sig_hits[i,"treatment"]
  element_type <- treatment_sig_hits[i,"genomic_element"]
  
  
  # read in the element info
  element_input_path <- paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external-files/genomic-elements/", element_type , ".tsv.gz")
  element_input <- read.csv(file = element_input_path, stringsAsFactors = F, header = T, sep = "\t")
  
  # make sure the rows corresponding to chromosome M are excluded
  element_input <- element_input[element_input$chr != "chrM",]
  
  # select the element of interest
  element_input_of_interest <- element_input[element_input$id == element_id,]
  
  # save the data
  write.table(element_input_of_interest[,-4], file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/r-objects/non-coding-drivers/version-202211/mutation-info-retrieve/20221123-request/", i, ".element.bed"),
              col.names = T, quote = F, sep = "\t", row.names = F)
  
  # read in the mutation info
  mutation_input_path <- paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external-files/input-som-muts-ActiveDriverWGS/version-202211/", cancer_type, "/", cancer_type, ".", treatment, ".dndscv.activedriverwgs.input.tsv.gz")
  mutation_input <- read.csv(file = mutation_input_path, stringsAsFactors = F, header = T, sep = "\t")
  
  # add the additional metadata of significantly mutated elements to mutations for later reference
  mutation_input[,colnames(treatment_sig_hits)] <- treatment_sig_hits[i,]
  
  # save the data
  write.table(mutation_input, file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/r-objects/non-coding-drivers/version-202211/mutation-info-retrieve/20221123-request/", i, ".mut-input.bed"),
              col.names = T, quote = F, sep = "\t", row.names = F)
}


# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

## Bash code to retrieve the intersecting mutations

## Intersect significantly mutated element ranges with mutations
# for i in {1..nr_files}; do echo ${i} & ~/bedtools intersect -a ${i}.mut-input.bed -b ${i}.element.bed > ./${i}.intersect.tsv; done

## Concatanate the files
# for i in {1..nr_files}; do cat ${i}.intersect.tsv >> ../treatment-enriched-elements-full-info.tsv; done

## Assigning the correct field names using sed
# sed -i '1i chr\pos1\pos2\ref\alt\patient_id\element_id\_colnames(treatment_sig_hits)[-1]_' ../treatment-enriched-elements-full-info.tsv

# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------------------


# read in the data
sig_hit_info <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/r-objects/non-coding-drivers/version-202211/mutation-info-retrieve/treatment-enriched-elements-full-info.tsv"),
                         header = F, stringsAsFactors = F, sep = "\t")

# remove hypermutated samples
sig_hit_info <- sig_hit_info[sig_hit_info$patient_id %notin% sample_metadata$sample_id[sample_metadata$smnv_load > 90000],]



# save the data
write.table(sig_hit_info, file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/r-objects/non-coding-drivers/version-202211/mutation-info-retrieve/20221123-request/treatment-enriched-elements-full-info.tsv"),
            row.names = F, quote = F, sep = "\t")