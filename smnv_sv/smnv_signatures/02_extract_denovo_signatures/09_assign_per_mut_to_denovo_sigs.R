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

## Load inputs ================================
## Args --------------------------------
args <- commandArgs(trailingOnly=T)

sample_name <- args[1]
som_vcf_path <- args[2]
sig_path_SBS <- args[3]
sig_path_DBS <- args[4]
sig_path_ID <- args[5]

out_dir <- args[6]
if(is.na(out_dir)){ out_dir <- paste0(wd,'/sig_contrib/') }

## Misc data --------------------------------
fit_metadata <- read.delim(paste0(wd,'/metadata/fit_metadata.txt.gz'))
sig_metadata <- read.delim(paste0(wd,'/sig_contrib/fit_lsq.post_processed/sig_metadata.post_processed.txt.gz'))

if(F){
   fit_metadata_ss <- fit_metadata[fit_metadata$sample=='WIDE01011147T',]
   
   ##
   sample_name <- fit_metadata_ss$sample
   som_vcf_path <- paste0(path_prefix, fit_metadata_ss$sample_dir, fit_metadata_ss$som_vcf_path)
   sig_path_SBS <- paste0(path_prefix, fit_metadata_ss$sig_profile_dir, fit_metadata_ss$sig_path_SBS)
   sig_path_DBS <- paste0(path_prefix, fit_metadata_ss$sig_profile_dir, fit_metadata_ss$sig_path_DBS)
   sig_path_ID <- paste0(path_prefix, fit_metadata_ss$sig_profile_dir, fit_metadata_ss$sig_path_ID)
}

## Main ================================
message('\n## Loading variants')
variants <- variantsFromVcf(som_vcf_path, vcf.filter='PASS', verbose=T)

message('\n## Extracting contexts')
contexts <- list(
   SBS=extractSigsSnv(df=variants, output='contexts')[,1],
   DBS=extractSigsDbs(df=variants, output='contexts')[,1],
   ID=extractSigsIndel(df=variants, output='contexts', method='PCAWG')[,1]
)

message('## Loading denovo sig profiles')
sig_profiles <- list(
   SBS=read.delim(sig_path_SBS),
   DBS=read.delim(sig_path_DBS),
   ID=read.delim(sig_path_ID)
)

## Assignment ================================
message('\n## Performing per mutation signature assignment')
sig_assignments <- lapply(names(sig_profiles), function(i){
   #i=names(sig_profiles)[2]

   df <- assignSigPerMut(
      df=variants, fit.method='lsq',
      mode=i, signature.profiles=sig_profiles[[i]],
      verbose=T
   )
   
   if(nrow(df)>=1){
      df$mut_type <- i
      
      ## Make unique sig name
      df$denovo_name_uniq <- sub('[.].+$','',df$assigned_sig)
      
      tissue_type_group <- fit_metadata[
         match(sample_name, fit_metadata$sample),
         'tissue_type_group'
      ]
      df$denovo_name_uniq <- paste0(tissue_type_group,'.',df$denovo_name_uniq)
      
      ## Add metdata
      index <- match(df$denovo_name_uniq, sig_metadata$denovo_name_uniq)
      df$sig_name <- sig_metadata$sig_name[ index ]
      df$etiology <- sig_metadata$etiology[ index ]
      df$matched_ref_cos_sim <- sig_metadata$matched_ref_cos_sim[ index ]
      df$mean_within_clust_cos_sim <- sig_metadata$mean_within_clust_cos_sim[ index ]
      
      ## Clean up columns
      df$assigned_sig <- NULL
      df$denovo_name_uniq[is.na(df$sig_name)] <- NA
   } else {
      df$mut_type <-
         df$denovo_name_uniq <-
         df$sig_name <-
         df$etiology <-
         df$matched_ref_cos_sim <-
         df$mean_within_clust_cos_sim <-
         character()
   }
   
   return(df)
})
sig_assignments <- do.call(rbind, sig_assignments)

##
muts_assigned_dir <- paste0(out_dir,'/muts_assigned/')
dir.create(muts_assigned_dir, showWarnings=F)
write.table(
   sig_assignments,
   gzfile(paste0(muts_assigned_dir,'/',sample_name,'.txt.gz')),
   sep='\t',quote=F, row.names=F
)
