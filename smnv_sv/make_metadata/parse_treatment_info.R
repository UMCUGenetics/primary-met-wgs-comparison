## Init ================================
options(stringsAsFactors=F)

WRITE_OUTPUT <- FALSE

## Paths --------------------------------
path_prefix <- c(
   hpc  ='/hpc/',
   local=path.expand('~/hpc/')
)
path_prefix <- path_prefix[ min(which(dir.exists(path_prefix))) ]
path_prefix <- sub('/hpc/','/',path_prefix)

base_dir <- paste0(path_prefix,'/hpc/cuppen/projects/P0025_PCAWG_HMF/')

## Load data ================================
treatment_info <- read.delim(paste0(base_dir,'/passengers/processed/metadata/treatment/info_treatment.tsv'))
treatment_info$dummy_bool <- TRUE

treatment_info <- subset(treatment_info, mechanism!='') ## Remove Isolated Limb Perfusion (ILP) right leg

## Drug type ================================
##  Preliminary one hot encoded treatment data
m_treatment_pre <- reshape2::acast(
   data=treatment_info, 
   formula=sample_id~mechanism, 
   value.var='dummy_bool', fill=NA_real_
)
m_treatment_pre[is.na(m_treatment_pre)] <- 0
m_treatment_pre[m_treatment_pre>=1] <- 1 ## Some samples can have e.g. more than one platinum therapy
m_treatment_pre <- !!m_treatment_pre

## --------------------------------
## Some columns represent multiple treatments
## Reassign these TRUE values to the correct column with a single treatment.
#treatment_names <- tolower(colnames(m_treatment_pre))
treatments_lower <- tolower(colnames(m_treatment_pre))
treatments_lower_groups <- lapply(treatments_lower, function(i){
   strsplit(i,', ')[[1]]
})
names(treatments_lower_groups) <- colnames(m_treatment_pre)
treatments_lower_uniq <- sort(unique(unlist(treatments_lower_groups, use.names=F)))

## Get treatment name, preserving capitals
treatments_comma_no_yes <- c(
   grep(',', colnames(m_treatment_pre), invert=T, value=T),
   grep(',', colnames(m_treatment_pre), invert=F, value=T)
)
treatments_comma_no_yes <- paste0('__',treatments_comma_no_yes,'__')
treatments_comma_no_yes <- gsub(', ','__, __',treatments_comma_no_yes)

names(treatments_lower_uniq) <- sapply(treatments_lower_uniq, function(i){
   #i=treatments_lower_uniq[70]
   sel_index <- grep(
      paste0('__',i,'__'),
      tolower(treatments_comma_no_yes),
      fixed=T
   )
   treatments_comma_no_yes[ sel_index[1] ]
})

names(treatments_lower_uniq) <- gsub('__','',names(treatments_lower_uniq))

## --------------------------------
## Find the columns that need to be merged
merge_colnames <- lapply(treatments_lower_uniq, function(i){
   #i='anaplastic lymphoma kinase inhibitor'
   bool <- sapply(treatments_lower_groups, function(j){
      i %in% j
   })
   names(bool)[bool]
})
names(merge_colnames) <- names(treatments_lower_uniq)
merge_colnames <- merge_colnames[sort(names(merge_colnames))]
merge_colnames <- merge_colnames[!duplicated(merge_colnames)]

## Merge TRUE values
m_treatment <- do.call(cbind, lapply(merge_colnames, function(i){
   i_m <- m_treatment_pre[,i,drop=F]
   rowSums(i_m)
}))
m_treatment <- m_treatment>0
m_treatment <- m_treatment[,order(colnames(m_treatment))]

## Rename treatments
renameCols <- function(m, old.name, new.name){
   colnames(m)[colnames(m) %in% old.name] <- new.name
   #return(colnames(m))
   return(m)
}

m_treatment <- renameCols(
   m_treatment,
   old.name=c('Other','Unknown'),
   new.name=c('Other treatment(s)','Unknown treatment(s)')
)

## Export
write.table(
   m_treatment,
   gzfile(paste0(base_dir,'/passengers/processed/metadata/treatment/m_treatment.txt.gz')),
   sep='\t',quote=F
)

## Treatment type ================================
m_treatment_type <- reshape2::acast(
   data=treatment_info, 
   formula=sample_id~type, 
   value.var='dummy_bool', fill=NA_real_
)
m_treatment_type[is.na(m_treatment_type)] <- 0
m_treatment_type[m_treatment_type>=1] <- 1 ## Some samples can have e.g. more than one platinum therapy
m_treatment_type <- !!m_treatment_type

## Export
write.table(
   m_treatment_type,
   gzfile(paste0(base_dir,'/passengers/processed/metadata/treatment/m_treatment_type.txt.gz')),
   sep='\t',quote=F
)








