options(stringsAsFactors=F)

## Path prefixes ================================
path_prefix <- c(
   hpc  ='/hpc/',
   local=path.expand('~/hpc/')
)
path_prefix <- path_prefix[ min(which(dir.exists(path_prefix))) ]
path_prefix <- sub('/hpc/','/',path_prefix)

base_dir <- paste0(path_prefix,'/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/')

devtools::load_all(paste0(base_dir,'/CHORD/processed/scripts_main/mutSigExtractor/'))

## Main ================================
write.tsv <- function(...){ write.table(..., sep='\t', quote=F) }

args <- commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
out_dir <- args[2]
vcf_som <- args[3]
vcf_sv <- args[4]

if(F){
  sample_name <- 'test'
  vcf_som <- '/Users/lnguyen/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/DO51542-from-jar//purplesoft3.2/DO51542T.purple.somatic.vcf.gz'
}

#out_dir <- paste0(base_dir,'/datasets/processed/HMF_DR104_update3/matrices/')
#mut_types <- c('snv','dbs','indel_chord','indel')
mut_types <- c('snv','dbs','indel_chord','indel','sv')
out_paths <- paste0(out_dir,'/',mut_types,'/',sample_name,'_',mut_types,'.txt') ## subdirs were already created (manually)
names(out_paths) <- mut_types
for(i in paste0(out_dir,'/',mut_types,'/')){ dir.create(i, showWarnings=F) }

main <- function(sample_name, vcf_som, vcf_sv){

   contexts <- list()

   message('## Extracting SNV contexts...')
   if(!file.exists(out_paths['snv'])){
      contexts$snv <- extractSigsSnv(vcf_som, output='contexts', vcf.filter='PASS', sample.name=sample_name)
      write.tsv(contexts$snv, out_paths['snv'])
   }

   message('## Extracting indel contexts (CHORD)...')
   if(!file.exists(out_paths['indel_chord'])){
      contexts$indel_chord <- extractSigsIndel(vcf_som, vcf.filter='PASS', sample.name=sample_name, method='CHORD')
      write.tsv(contexts$indel_chord, out_paths['indel_chord'])
   }
   
   message('## Extracting indel contexts (PCAWG)...')
   if(!file.exists(out_paths['indel'])){
      contexts$indel <- extractSigsIndel(vcf_som, vcf.filter='PASS', sample.name=sample_name, method='PCAWG')
      write.tsv(contexts$indel, out_paths['indel'])
   }

   message('## Extracting DBS contexts...')
   if(!file.exists(out_paths['dbs'])){
      contexts$dbs <- extractSigsDbs(vcf_som, output='contexts', vcf.filter='PASS', sample.name=sample_name)
      write.tsv(contexts$dbs, out_paths['dbs'])
   }
   
   message('## Extracting SV contexts...')
   if(!file.exists(out_paths['sv'])){
      contexts$sv <- extractSigsSv(vcf_sv, output='contexts', vcf.filter='PASS', sample.name=sample_name)
      write.tsv(contexts$sv, out_paths['sv'])
   }
}

main(sample_name, vcf_som, vcf_sv)

# ## Local execution ================================
# vcf_paths <- read.delim(paste0(base_dir,'/datasets/processed/HMF_DR104_update3/metadata/vcf_paths.txt.gz'))
# for(i in 1:nrow(vcf_paths)){
#    message('Processing [',i,']: ',vcf_paths$sample[i])
#    main(
#       vcf_paths$sample[i],
#       paste0('/Users/lnguyen/',vcf_paths$som[i]),
#       paste0('/Users/lnguyen/',vcf_paths$sv[i])
#    )
# }
