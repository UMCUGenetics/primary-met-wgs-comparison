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

## Args ================================
args <- commandArgs(trailingOnly=T)

sample_name <- args[1]
som_vcf_path <- args[2]
sig_path_SBS <- args[3]
sig_path_DBS <- args[4]
sig_path_ID <- args[5]

out_dir <- args[6]
if(is.na(out_dir)){ out_dir <- paste0(wd,'/sig_contrib/') }

if(F){
   ##
   fit_metadata <- read.delim(paste0(wd,'/metadata/fit_metadata.txt.gz'))
   
   fit_metadata_ss <- fit_metadata[fit_metadata$tissue_type_group=='Lymphoid',]
   fit_metadata_ss <- fit_metadata_ss[1,]
   
   fit_metadata_ss <- fit_metadata[fit_metadata$sample=='DO36119',]
   
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

## lsq fit ================================
message('## Performing least squares fitting')
sig_contribs <- lapply(names(sig_profiles), function(i){
   #i=names(sig_profiles)[1]
   fitToSignatures(contexts[[i]], sig_profiles[[i]], scale.contrib=T)
})
names(sig_contribs) <- names(sig_profiles)

##
contexts_dir <- paste0(out_dir,'/contexts/')
dir.create(contexts_dir, showWarnings=F)
contexts_ex <- t(do.call(c, unname(contexts)))
write.table(
   contexts_ex, 
   paste0(contexts_dir,'/',sample_name,'.txt'),
   sep='\t',quote=F, row.names=F
)

##
fit_lsq_dir <- paste0(out_dir,'/fit_lsq/')
dir.create(fit_lsq_dir, showWarnings=F)
sig_contribs_ex <- t(do.call(c, unname(sig_contribs)))
write.table(
   sig_contribs_ex, 
   paste0(fit_lsq_dir,'/',sample_name,'.txt'),
   sep='\t',quote=F, row.names=F
)

# ## Assignment ================================
# message('\n## Performing per mutation signature assignment')
# sig_assignments <- lapply(names(sig_profiles), function(i){
#    #i=names(sig_profiles)[2]
# 
#    out <- assignSigPerMut(
#       df=variants, fit.method='lsq',
#       mode=i, signature.profiles=sig_profiles[[i]],
#       verbose=T
#    )
# 
#    if(nrow(out)>=1){
#       out$mut_type <- i
#    } else {
#       out$mut_type <- character()
#    }
# 
#    return(out)
# })
# sig_assignments <- do.call(rbind, sig_assignments)
# 
# ##
# muts_assigned_dir <- paste0(out_dir,'/muts_assigned/')
# dir.create(muts_assigned_dir, showWarnings=F)
# write.table(
#    sig_assignments,
#    gzfile(paste0(muts_assigned_dir,'/',sample_name,'.txt.gz')),
#    sep='\t',quote=F, row.names=F
# )
