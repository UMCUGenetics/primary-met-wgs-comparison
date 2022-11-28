options(stringsAsFactors=F)

## Paths ================================
path_prefix <- c(
   hpc='/hpc/',
   local=path.expand('~/hpc/')
)
path_prefix <- path_prefix[ min(which(dir.exists(path_prefix))) ]
path_prefix <- sub('/hpc/','/',path_prefix)

base_dir <- paste0(path_prefix,'/hpc/cuppen/projects/P0025_PCAWG_HMF/')
wd <- paste0(base_dir,'/passengers/processed/sigs_denovo/extractions/12_fixed_seeds/')

## Load data ================================
## Select samples --------------------------------
message('Loading metadata...')
manifest <- read.delim(paste0(base_dir,'/passengers/processed/manifest/05_20220216/manifest_HMF_PCAWG.txt.gz'))
metadata <- read.delim(paste0(wd,'/metadata/SuppTable1_sample_metadata.txt.gz'), na.strings=c('NA','#N/A'))

# ## Blacklist
# metadata <- subset(metadata, !is_blacklisted)

## Tissue type groups
large_tt_groups <- with(subset(metadata, !is_blacklisted),{
   tt_counts <- sort(table(tissue_type), decreasing=T)
   large_tt_groups <- names(tt_counts)[ tt_counts>=30 ]
   return(large_tt_groups)
})

metadata$tissue_type_group <- metadata$tissue_type
metadata$tissue_type_group[ !(metadata$tissue_type_group %in% large_tt_groups) ] <- 'Other'
metadata$tissue_type_group[ metadata$tissue_type_group=='Unknown' ] <- 'Other'

## Export
write.table(
   metadata, paste0(wd,'/metadata/nmf_metadata.txt'),
   sep='\t',row.names=F,quote=F
)

## Contexts --------------------------------
message('Loading contexts...')
contexts <- readRDS(paste0(base_dir,'/passengers/processed/mut_contexts/04_fixed_smnvs/matrices/contexts_merged.rds'))

mut_types <- c('snv','dbs','indel')
contexts <- contexts[mut_types]

## sigProfiler context names (only needed for indels)
context_translations <- read.delim(paste0(base_dir,'/passengers/processed/sigs_denovo/dep/mut_context_translations.txt'))
colnames(contexts$indel) <- context_translations$name.sigProfiler[
   match(colnames(contexts$indel), context_translations$name.mutSigExtractor)
]

## Use sigProfiler style matrix
contexts <- lapply(contexts, function(i){
   #i=contexts$snv
   df <- t(i)
   df <- data.frame(
      MutationType=rownames(df),
      #df[,1:10,drop=F],
      df,
      row.names=NULL
   )
   return(df)
})
#lapply(contexts,function(i){ i[,1:10] })

## Make param grid ================================
message('Spawning jobs...')

param_grid <- expand.grid(
   mut_type=mut_types,
   tissue_type=unique(metadata$tissue_type_group)
)

param_grid$max_rank <- c(snv=21, dbs=8, indel=10)[ param_grid$mut_type ]

## Exceptions
param_grid$max_rank[param_grid$mut_type=='indel' & param_grid$tissue_type=='Lung'] <- 15

##
param_grid$index <- 1:nrow(param_grid)
param_grid$index <- formatC(param_grid$index, format='d', flag='0', width=nchar(max(param_grid$index)))

param_grid$context_type_arg <- c(snv='96',dbs='DINUC',indel='ID')[param_grid$mut_type]

##
sample_groups <- with(
   subset(metadata, !is_blacklisted),
   split(sample_id, tissue_type_group)
)
sample_groups$Other <- subset(
   metadata, 
   tissue_type_group=='Other' & ((tissue_type=='Unknown' & is_selected_biopsy) | (tissue_type!='Unknown' & !is_blacklisted_sample)),
   sample_id, drop=T
)

master_out_dir <- paste0(wd,'/output/')
dir.create(master_out_dir, showWarnings=F)

## Spawn jobs ================================
"#!/bin/bash
#SBATCH --job-name=${tissue_type}_${context_type_arg}
#SBATCH --output=${out_dir}/slurm.out
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --ntasks-per-node=20

if [[ ! -f ${out_dir}/job.done ]]; then

source /hpc/cuppen/projects/P0025_PCAWG_HMF/passengers/processed/sigs_denovo/dep/env_sigProfiler/bin/activate
python ${wd}/scripts/01_run_SigProfiler.py \\
--input_data ${input_data_path} \\
--context_type ${context_type_arg} \\
--output ${out_dir} \\
--maximum_signatures ${max_rank} && touch ${out_dir}/job.done

else

echo Done file exists: ${out_dir}/job.done

fi

" -> job_template

writeString <- function(string, path){
   fileConn <- file(path)
   writeLines(string, fileConn)
   close(fileConn)
}

for(i in 1:nrow(param_grid)){
   #i=1
   mut_type <- as.character(param_grid$mut_type[i])
   context_type_arg <- as.character(param_grid$context_type_arg[i])
   tissue_type <- as.character(param_grid$tissue_type[i])
   max_rank <- param_grid$max_rank[i]
   #index <- param_grid$index[i]
   
   ## Make output dir
   out_dir <- paste0(master_out_dir,'/',tissue_type,'.',mut_type,'/')
   dir.create(out_dir, showWarnings=F)
   
   if(file.exists(paste0(out_dir,'/job.done'))){
      message('## NMF already complete @: ', basename(out_dir))
   } else {
      message('## Spawning at: ', basename(out_dir))
      
      ## Write input contexts matrix
      sel_samples <- sample_groups[[tissue_type]]
      input_data <- contexts[[mut_type]]
      input_data <- input_data[,c('MutationType',sel_samples)]
      input_data_path <- paste0(out_dir,'/contexts.txt')
      write.table(input_data, input_data_path, sep='\t', quote=F, row.names=F)
      
      ##
      job_string <- stringr::str_interp(
         job_template, 
         env=list(
            tissue_type=tissue_type,
            context_type_arg=context_type_arg,
            out_dir=out_dir,
            wd=wd,
            input_data_path=input_data_path,
            max_rank=max_rank
         )
      )
      job_string <- gsub('/Users/lnguyen/','/',job_string)
      job_string <- gsub('/+','/',job_string)
      
      writeString(job_string, paste0(out_dir,'/job.sh'))
   }
}

## Job submit script ================================
"#!/bin/bash
for i in ${master_out_dir}/*/; do 
   if [[ ! -f $i/job.done ]]; then
      sbatch $i/job.sh
   else
      echo Done file exists: $i/job.done
   fi
done
" -> submit_template

submit_template <- stringr::str_interp(
   submit_template, 
   env=list(master_out_dir=master_out_dir)
)
submit_template <- gsub('/Users/lnguyen/','/',submit_template)
submit_template <- gsub('/+','/',submit_template)

writeString(submit_template, paste0(master_out_dir,'/submit_jobs.sh'))
