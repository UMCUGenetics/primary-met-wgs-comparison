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
wd <- paste0(base_dir,'/passengers/processed/metadata/gene_status/')

## Dependencies --------------------------------
source(paste0(base_dir,'/passengers/analysis/common_functions/common_functions.R'))

## Load data ================================
##
linx_drivers <- cacheAndReadData(
   paste0(base_dir,'/passengers/processed/linx_tables/linx_drivers.merged.txt.gz')
)

fragile_sites <- read.csv(paste0(wd,'/fragile_sites_hmf.csv'))

## Main ================================
## Subset linx drivers --------------------------------
df <- linx_drivers
df <- df[df$sample %in% sample_metadata$sample_id,]
df <- df[df$driverLikelihood>0.5,]

## Cast to matrix form --------------------------------
m <- reshape2::acast(data=df, formula=sample~gene, value.var='sample', fun.aggregate=length, fill=NA_real_)
m <- as.data.frame(m)

## Remove fragile site genes
m <- m[,!(colnames(m) %in% fragile_sites$Gene)]

## Add missing samples
m <- m[sample_metadata$sample_id,]; rownames(m) <- sample_metadata$sample_id 

## Fill NAs with 0
m[is.na(m)] <- 0 

## Convert to bool
m <- m>=1 

## Remove infrequent genes
m <- m[,colSums(m)>=15]

## Export --------------------------------
## HPC
write.table(
   m,
   gzfile(paste0(base_dir,'/passengers/processed/metadata/gene_status/has_gene_mut.txt.gz')),
   sep='\t',row.names=T,quote=F
)