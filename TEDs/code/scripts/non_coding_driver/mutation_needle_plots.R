

## In this script we make needle plots for genes of interest in different cancer type and treatment groups

## Init --------------------------------

expanding_path <- '/home/ali313/Documents/studies/master/umc-project'

path_prefix <- c(
  hpc  ='/hpc/',
  local = paste0(expanding_path, '/hpc/')
)
path_prefix <- path_prefix[min(which(dir.exists(path_prefix)))]
path_prefix <- sub('/hpc/','',path_prefix)

base_dir <- paste0(path_prefix,'/hpc/cuppen/projects/P0025_PCAWG_HMF/')


## Select the treatment group, cancer type, and genes you want to depict the mutation in (gene lengths were obtained from UniProtKB)  --------------------------------

treatments <- c("Anti_EGFR", "untreated",
                "Anti_AR__GnRH", "untreated",
                "Aromatase_inhibitor", "untreated",
                "Aromatase_inhibitor", "untreated")

cancer_type_codes <- c("LUAD", "LUAD",
                       "PRAD", "PRAD",
                       "BRCA", "BRCA",
                       "BRCA_ERpos", "BRCA_ERpos")

genes_of_interest <- c("EGFR", "EGFR",
                       "AR", "AR",
                       "ESR1", "ESR1",
                       "ESR1", "ESR1")

protein_length <- c(1210, 1210,
                    920, 920,
                    595, 595,
                    595, 595)



# find the corresponding mutations from the mutation file, assign their mutation type, and calculate the frequencies of each mutation type
for (i in 1:8){
  
  print(i)
  
  treatment <- treatments[i]
  cancer_type_code <- cancer_type_codes[i]
  gene_of_interest <- genes_of_interest[i]
  
  # read in the mutation info
  mutation_input <- read.csv(file = paste0(base_dir, "drivers/processed/resistance/dndscv/", cancer_type_code, "/", treatment, ".dndscv.annotmuts.tsv.gz"),
                             header = T, stringsAsFactors = F, sep = "\t")
  
  # select for mutations that happen in the gene of interest
  mutations_of_interet <- mutation_input[mutation_input$gene == gene_of_interest,]
  
  
  if (nrow(mutations_of_interet) > 0){
    
    
    # assign the correct mutation type to mutations
    mutations_of_interet$impact[nchar(mutations_of_interet$ref) > nchar(mutations_of_interet$mut) & str_detect(mutations_of_interet$ntchange, pattern = "inf")] <- "Inframe Deletion"
    mutations_of_interet$impact[nchar(mutations_of_interet$ref) > nchar(mutations_of_interet$mut) & !str_detect(mutations_of_interet$ntchange, pattern = "inf")] <- "Frameshift Deletion"
    mutations_of_interet$impact[nchar(mutations_of_interet$ref) < nchar(mutations_of_interet$mut) & str_detect(mutations_of_interet$ntchange, pattern = "inf")] <- "Inframe Insertion"
    mutations_of_interet$impact[nchar(mutations_of_interet$ref) < nchar(mutations_of_interet$mut) & !str_detect(mutations_of_interet$ntchange, pattern = "inf")] <- "Frameshift Insertion"
    mutations_of_interet$impact[nchar(mutations_of_interet$ref) == nchar(mutations_of_interet$mut) & nchar(mutations_of_interet$mut) > 1] <- "MNV"
    
    mutations_of_interet$aachange[mutations_of_interet$aachange == "."] <- mutations_of_interet$impact[mutations_of_interet$aachange == "."]
    
    # reshape the data frame
    mutations_of_interet <- mutations_of_interet[,c(3, 15, 17, 1:2, 4:14, 16)]
    colnames(mutations_of_interet)[1:2] <- c("coord", "category")
    
    mutations_of_interet$coord[mutations_of_interet$codonsub != "."] <- str_sub(mutations_of_interet$aachange[mutations_of_interet$codonsub != "."],2,-2)
    mutations_of_interet$coord[mutations_of_interet$codonsub == "."] <- ceiling(as.numeric(lapply(str_split(mutations_of_interet$ntchange[mutations_of_interet$codonsub == "."], pattern = "-"), "[[", 1))/3)
    
    
    dups_coord <- mutations_of_interet$coord[duplicated(paste(mutations_of_interet$coord, mutations_of_interet$category))]
    
    
    # the value is the number of instances of each mutation type
    mutations_of_interet$value <- 1
    
    for (i in 1:length(unique(dups_coord))){
      dups_aachange <- mutations_of_interet$category[mutations_of_interet$coord == unique(dups_coord)[i]]
      for (j in 1:length(unique(dups_aachange)))
        mutations_of_interet$value[mutations_of_interet$coord == unique(dups_coord)[i] & mutations_of_interet$category == unique(dups_aachange)[j]] <- nrow(mutations_of_interet[mutations_of_interet$coord == unique(dups_coord)[i] & mutations_of_interet$category == unique(dups_aachange)[j],])
    }
    
    
    # remove duplicate mutations
    mutations_of_interet <- mutations_of_interet[!duplicated(paste(mutations_of_interet$coord, mutations_of_interet$category)),]
    
    
  } else {
    mutations_of_interet$value <- numeric(0)
    
    mutations_of_interet <- mutations_of_interet[,c(3, 15, 17, 1:2, 4:14, 16)]
    
    colnames(mutations_of_interet)[1:2] <- c("coord", "category")
    
    mutations_of_interet$coord <- as.character(mutations_of_interet$coord)
    
  }
  
  # save the mutation list
  write.table(mutations_of_interet, file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/r-objects/lollipops/positions/20221123-request/", gene_of_interest, ".", treatment, ".", cancer_type_code, ".tsv"),
              row.names = F, quote = F, sep = "\t")
}









## Make the needle plots  --------------------------------

# load necessary packages
library("org.Hs.eg.db")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(mutsneedle)
library(shiny)
library(stringr)



# assign the index of the groups to plot
i <- 1
treatment <- treatments[i]
cancer_type_code <- cancer_type_codes[i]
gene_of_interest <- genes_of_interest[i]
length_of_protein <- protein_length[i]

# get the gene ids
genes_id <- as.character(unlist(mget(x=genes_of_interest,envir=org.Hs.egALIAS2EG)))
gene_id <- genes_id[i]

# read in the mutation list to be plotted
mutations_of_interet <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/r-objects/lollipops/positions/20221123-request/", gene_of_interest, ".", treatment, ".", cancer_type_code, ".tsv"),
                                 header = T, stringsAsFactors = F, sep = "\t")
mutations_of_interet$coord <- as.character(mutations_of_interet$coord)


# remove Essential_Splice type of mutations
mutations_of_interet <- mutations_of_interet[mutations_of_interet$category != "Essential_Splice",]




## Prepare the input data frame of shinyApp  --------------------------------
regiondata <- exampleRegionData()

regiondata[1,] <- c(paste0("0-", length_of_protein), "Gene")
regiondata <- regiondata[1,]


# assign the gene elements for each gene (poistions obtained from EMBL_EBI Database)
if (i %in% c(1,2)){
  
  regiondata[2,] <- c("57-168", "Recep_L_domain")
  regiondata[3,] <- c("177-338", "Furin-like")
  regiondata[4,] <- c("361-481", "Recep_L_domain")
  regiondata[5,] <- c("505-637", "GF_recep_IV")
  regiondata[6,] <- c("712-968", "PK_Tyr_Ser-Thr")
} else if (i %in% c(3,4)){
  
  regiondata[2,] <- c("6-449", "Androgen_recep")
  regiondata[3,] <- c("558-627", "zf-C4")
  regiondata[4,] <- c("690-881", "Hormone_recep")
  
} else if (i %in% 5:8){
  
  regiondata[2,] <- c("42-181", "Oest_recep")
  regiondata[3,] <- c("183-252", "zf_C4")
  regiondata[4,] <- c("332-531", "Hormone_recep")
  regiondata[5,] <- c("552-595", "ESR1_C")
  
  
} 


# define the category_vec of interest
category_vec <- c("Nonsense", "Missense", "Synonymous", "Frameshift Deletion", "Inframe Deletion", "Frameshift Insertion", "Inframe Insertion",
                  "MNV", "Stop_loss")


# this is only for representation purposes
for (k in 1:length(category_vec)){
  if (category_vec[k] %notin% mutations_of_interet$category){
    mutations_of_interet <- rbind(mutations_of_interet, c(5000, category_vec[k], 0))
  }
}


# set the order of the data frame
index <- order(match(mutations_of_interet$category, category_vec, nomatch = F))
index <- index[index != 0]
mutations_of_interet <- mutations_of_interet[index,]



# this is only for representation purposes (trying to have the y-axis of treated and untreated group in the same scale)
if (i %in% c(2)){
  while (j < 19){
    j <- j + 1
    mutations_dummy <- mutations_of_interet[2,]
    mutations_dummy$coord <- 4000
    mutations_of_interet <- rbind(mutations_of_interet, mutations_dummy)
  }
} else if (i %in% c(4)){
  while (j < 20){
    j <- j + 1
    mutations_dummy <- mutations_of_interet[2,]
    mutations_dummy$coord <- 4000
    mutations_of_interet <- rbind(mutations_of_interet, mutations_dummy)
  }
} else if (i %in% c(6)){
  while (j < 46){
    j <- j + 1
    mutations_dummy <- mutations_of_interet[2,]
    mutations_dummy$coord <- 4000
    mutations_of_interet <- rbind(mutations_of_interet, mutations_dummy)
  }
} else if (i %in% c(8)){
  while (j < 40){
    j <- j + 1
    mutations_dummy <- mutations_of_interet[6,]
    mutations_dummy$coord <- 4000
    mutations_of_interet <- rbind(mutations_of_interet, mutations_dummy)
  }
} 


# set the row names
row.names(mutations_of_interet) <- 1:nrow(mutations_of_interet)



## Plot the mutations using shinyApp and then save  --------------------------------

mm <- shinyApp(
  ui = mutsneedleOutput("id",width=800,height=500),
  server = function(input, output) {
    output$id <- renderMutsneedle(
      mutsneedle(mutdata=mutations_of_interet,domains=regiondata,
                 gene = gene_of_interest,
                 # colorMap = ccc,
                 xlab="Protein Location",
                 ylab="Mutation Frequency",
                 maxlength = length_of_protein)
    )
  }
)


