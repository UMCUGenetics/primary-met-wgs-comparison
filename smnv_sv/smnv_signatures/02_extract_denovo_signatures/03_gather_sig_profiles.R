options(stringsAsFactors=F)

## Paths ================================
path_prefix <- c(
   hpc  ='/hpc/',
   local=path.expand('~/hpc/')
)
path_prefix <- path_prefix[ min(which(dir.exists(path_prefix))) ]
path_prefix <- sub('/hpc/','/',path_prefix)

base_dir <- paste0(path_prefix,'/hpc/cuppen/projects/P0025_PCAWG_HMF/')
wd <- paste0(base_dir,'/passengers/processed/sigs_denovo/extractions/12_fixed_seeds/')

devtools::load_all(paste0(path_prefix,'/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORD/processed/scripts_main/mutSigExtractor'))

## Load data ================================
##
selected_ranks <- read.delim(paste0(wd,'/metadata/selected_ranks.txt'))
context_translations <- read.delim(paste0(base_dir,'/passengers/processed/sigs_denovo/dep/mut_context_translations.txt'))

##
nmf_metadata <- (function(){
   df <- data.frame(
      dir=list.dirs(paste0(wd,'/output/'), recursive=F)
   )
   df$extraction_name <- basename(df$dir)
   df$tissue_type <- sub('[.].+$','',df$extraction_name)
   
   df$mut_type <- sub('^.+[.]','',df$extraction_name)
   df$mut_type <- c(snv='SBS',dbs='DBS',indel='ID')[df$mut_type]
   
   df$rank <- selected_ranks$optimum_rank[ match(df$extraction_name, selected_ranks$extraction_name) ]
   #df$extraction_name <- NULL
   return(df)
})()

## Get matrices ================================
getMatrixPath <- function(dir, rank, what='signatures'){
   #dir='/Users/lnguyen//hpc/cuppen/projects/P0025_PCAWG_HMF//passengers/processed/sigs_denovo/extractions/06_DR104-update4_pcawgSoftFilt//output//Biliary.dbs'
   #what='signatures'
   #rank=8
   
   dir.mut_type <- grep(
      '(SBS|DBS|ID)\\d+$',
      list.dirs(dir, recursive=F),
      value=T
   )
   if(length(dir.mut_type)>1){
      dir.mut_type <- dir.mut_type[1]
      warning('Multiple mut type dirs were found selecting the first: ', dir.mut_type)
   }
   
   ## Get paths to matrices
   dir.solutions <- paste0(dir.mut_type,'/All_Solutions/')
   if(what=='signatures'){
      paths <- list.files(dir.solutions,'*_Signatures.txt$', recursive=T, full.names=T)
   } else if(what=='activities'){
      paths <- list.files(dir.solutions,'*_Activities.txt$', recursive=T, full.names=T)
   } else {
      stop("`what` must be 'signatures' or 'activities'")
   }
   
   ## Get NMF rank
   rank_name <- sapply( strsplit(basename(paths),'_'), `[[`, 2 )
   rank_name <- sub('^S','',rank_name)
   
   paths[ as.integer(rank_name)==rank ]
}

PROFILES.REF <- list(
   SBS=mutSigExtractor::SBS_SIGNATURE_PROFILES_V3,
   DBS=mutSigExtractor::DBS_SIGNATURE_PROFILES,
   ID=mutSigExtractor::INDEL_SIGNATURE_PROFILES
)

sig_profiles <- lapply(1:nrow(nmf_metadata), function(i){
   #i=28
   
   ## Args
   dir <- nmf_metadata$dir[i]
   tissue_type <- nmf_metadata$tissue_type[i]
   mut_type <- nmf_metadata$mut_type[i]
   rank <- nmf_metadata$rank[i]
   message(i,': ',tissue_type,' ',mut_type,', rank=',rank)
   
   if(rank==0){ 
      ## Make dummy matrix for tissue types without any signatures
      context_names <- rownames(PROFILES.REF[[mut_type]])
      
      m <- matrix(
         data=0, 
         nrow=length(context_names), 
         dimnames=list(context_names, paste0(mut_type,'_0'))
      )
      
      class(m) <- c(class(m), 'dummy')
      
   } else {
      ## Read sig profile
      path <- getMatrixPath(dir, rank)
      m <- read.delim(path)
      
      ## Convert to mutSigExtractor context names
      m$MutationType <- context_translations$name.mutSigExtractor[
         match(m$MutationType, context_translations$name.sigProfiler)
      ]
      
      ## Convert to matrix
      rownames(m) <- m$MutationType
      m$MutationType <- NULL
      m <- as.matrix(m)
      
      ## Fix denovo sig names
      colnames(m) <- paste0(mut_type,'_',LETTERS[1:ncol(m)])  #sub('\\d+','_',colnames(m))
   }
   
   class(m) <- c(class(m), mut_type)
   return(m)
})
names(sig_profiles) <- with(nmf_metadata, paste0(tissue_type,'.',mut_type,'.',rank))
#names(sig_profiles) <- with(nmf_metadata, paste0(tissue_type,'.',mut_type))

## Add ref sig to denovo sig name ================================
sig_metadata <- read.delim(mutSigExtractor:::SIG_METADATA_PATH)

matchDenovoToRefSigs <- function(profiles.denovo, output='profiles'){
   if(F){
      profiles.denovo=sig_profiles$Lymphoid_DBS_0
      profiles.denovo=sig_profiles$Biliary_DBS_3
      output='profiles'
   }

   ## Init
   output <- match.arg(output, c('cossims','profiles'))

   ##
   if(inherits(profiles.denovo,'dummy')){
      
      max_sim_ref_sig <- 'NA'
      sig_etiology <- 'NA'
      max_sim <- 'NA'
      
   } else {
      ## Get relevant ref signature profiles
      mut_type <- grep(c('SBS|DBS|ID'), class(profiles.denovo), value=T)
      profiles.ref <- PROFILES.REF[[mut_type]]
      
      ## Cos sim function
      calcProfileCosSim <- function(profiles.denovo, profiles.ref){
         cosSim <- function(x, y) { sum(x*y)/sqrt(sum(x^2)*sum(y^2)) }
         
         out <- apply(profiles.denovo, 2, function(i){
            #i=profiles.de_novo[,1]
            apply(profiles.ref, 2, function(j){
               cosSim(i, j)
            })
         })
         
         t(out)
      }
      
      ## Match each denovo sigs to a ref sig
      cos_sims <- calcProfileCosSim(profiles.denovo, profiles.ref)
      
      max_sim <- round(apply(cos_sims,1,max),2) * 100
      max_sim_ref_sig <- colnames(cos_sims)[ max.col(cos_sims) ]
      
      sig_etiology <- sig_metadata$sig_etiology[ match(max_sim_ref_sig, sig_metadata$sig_name) ]
   }
   
   ## Add info to denovo sig names
   sig_names <- paste0(
      colnames(profiles.denovo),'.',
      max_sim_ref_sig,'.',
      sig_etiology,'.',
      max_sim
   )
   colnames(profiles.denovo) <- sig_names
   
   return(profiles.denovo)
}

sig_profiles_2 <- lapply(sig_profiles, matchDenovoToRefSigs)

## Export ================================
sig_profiles_dir <- paste0(wd,'/sig_profiles/')
dir.create(sig_profiles_dir, showWarnings=F)

for(i in names(sig_profiles_2)){
   #i=names(sig_profiles_2)[1]
   out_path <- paste0(sig_profiles_dir,'/',i,'.txt')
   write.table(sig_profiles_2[[i]], out_path, sep='\t',quote=F)
}
