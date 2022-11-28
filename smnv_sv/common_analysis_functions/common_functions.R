## Init ================================
path_prefix <- c(
   hpc  ='/hpc/',
   local=path.expand('~/hpc/')
)
path_prefix <- path_prefix[ min(which(dir.exists(path_prefix))) ]
path_prefix <- sub('/hpc/','/',path_prefix)

base_dir <- paste0(path_prefix,'/hpc/cuppen/projects/P0025_PCAWG_HMF/')

## Dependencies --------------------------------
devtools::load_all(paste0(base_dir,'/passengers/analysis/common_functions/statsExtra'))

## Sample metadata ================================
sample_metadata <- read.delim(
   paste0(base_dir,'/passengers/processed/metadata/main/metadata_20221111_1251.txt.gz'),
   na.strings=c('NA','#N/A')
)

## Make factor but retain order in data
factor2 <- function(x, ...){ factor(x=x, levels=unique(x), ...) } 

## Force order: tissue group > cancer type > cancer subtype --------------------------------
sample_metadata <- sample_metadata[
  with(sample_metadata, order(tissue_group, cancer_type, cancer_subtype)) 
,]

sample_metadata$tissue_type <- factor2(sample_metadata$tissue_type)
sample_metadata$tissue_group <- factor2(sample_metadata$tissue_group)
sample_metadata$cancer_type <- factor2(sample_metadata$cancer_type)
sample_metadata$cancer_subtype <- factor2(sample_metadata$cancer_subtype)
sample_metadata$cancer_type_code <- factor2(sample_metadata$cancer_type_code)

## Force order for other columns --------------------------------
sample_metadata$cohort <- as.factor(sample_metadata$cohort)
sample_metadata$cancer_stage <- as.factor(ifelse(sample_metadata$cohort=='Hartwig', 'Metastatic','Primary'))

sample_metadata$metastatic_location <- factor(sample_metadata$metastatic_location, c('Distant','Lymph','Local','Primary'))
sample_metadata$progression_status <- as.factor(sample_metadata$progression_status)
sample_metadata$progression_status_code <- as.factor(sample_metadata$progression_status_code)

## Formatting of other columns --------------------------------
## Primary/metastatic
sample_metadata$cohort[sample_metadata$cohort=='HMF'] <- 'Hartwig'
sample_metadata$is_hmf_sample <- sample_metadata$cohort=='Hartwig'

## HRD
sample_metadata$is_hrd <- sample_metadata$hr_status=='HR_deficient'
sample_metadata$is_hrd[sample_metadata$hr_status=='cannot_be_determined'] <- NA

## MSI
sample_metadata$is_msi <- sample_metadata$msi_status=='MSI'
sample_metadata$is_msi[is.na(sample_metadata$msi_status)] <- NA

## Various group splits --------------------------------
sample_metadata$cancer_type_subtype <- (function(){
   v <- paste0(
      sample_metadata$cancer_type,', ',
      sample_metadata$cancer_subtype
   )
   
   v <- sub('(, NA)|(, None)', '',v)
   v <- factor2(v)
   return(v)
})()

sample_metadata$cancer_type_code_subtype <- (function(){
   v <- paste0(
      sample_metadata$cancer_type_code,', ',
      sample_metadata$cancer_subtype
   )
   
   v <- sub('(, NA)|(, None)', '',v)
   v <- factor2(v)
   return(v)
})()

sample_metadata$cancer_type_code_supertype_subtype <- (function(){
   v <- paste0(
      sample_metadata$cancer_type_code,', ',
      sample_metadata$cancer_subtype
   )
   
   #v <- sub(', NA', '',v)
   v[grepl(', NA', v)] <- NA
   v <- sub(', None', '',v)
   v <- factor2(v)
   return(v)
})()

sample_metadata$cancer_type_code_subtype_subset <- (function(){
   v <- paste0(
      sample_metadata$cancer_type_code,', ',
      sample_metadata$cancer_subtype
   )
   
   v[grep('(, NA)|(, None)',v)] <- NA
   v <- factor2(v)
   return(v)
})()

sample_metadata$cancer_type_code_subtype_subset_hr_status <- (function(){
   v <- paste0(
      sample_metadata$cancer_type_code_subtype_subset, ', ',
      ifelse(sample_metadata$is_hrd,'HRD','non-HRD')
   )
   v[grep('NA',v)] <- NA
   v <- factor(v, sort(unique(v)))
   return(v)
})()

sample_metadata$cancer_type_and_code <- with(sample_metadata,{
   factor2( paste0(cancer_type,'\n(',cancer_type_code,')') )
})

## Gene status --------------------------------
m_gene_status <- read.delim(paste0(base_dir,'/passengers/processed/metadata/gene_status/has_gene_mut.txt.gz'), check.names=F)

## PURPLE data --------------------------------
purple_purity <- read.delim(paste0(base_dir,'/passengers/processed/metadata/purity/05_20220216/purple.purity.txt.gz'), check.names=F)
sample_metadata$ploidy <- purple_purity$ploidy[ match(sample_metadata$sample_id, purple_purity$sample) ]

## Treatment data --------------------------------
m_treatment <- read.delim(
   paste0(base_dir,'/passengers/processed/metadata/treatment/m_treatment.txt.gz'),
   check.names=F
)

## Sample retrieval functions ================================
## Retrieve metadata per sample
getSampleMetadata <- function(sample.names, keys=NULL, drop.dimensions=T, drop.factor.levels=T){
   #sample.names=c('CPCT02030502','CPCT02180026')
   #sample.names=pd$sample
   
   ## If keys not specified, return all columns
   if(is.null(keys)){ 
     keys <- colnames(sample_metadata) 
   }
   
   ## Get metadata for samples
   out <- sample_metadata[match(sample.names, sample_metadata$sample_id), keys, drop=drop.dimensions]
   
   ## Remove unused factor levels
   if(drop.factor.levels & (is.data.frame(out) | is.factor(out))){ 
     out <- droplevels(out) 
   }
   
   return(out)
}

## Subset samples by based on `sample_metadata` column conditions
subsetSamples <- function(x, ...){
   ## Sample selection
   sel_samples <- subset(
      sample_metadata, 
      ...,
      #is_blacklisted==F & is_blacklisted_subtype==F, 
      sample_id, 
      drop=T
   )
   
   ## Select rows in data
   main <- function(xx){
      out <- xx[match(sel_samples, rownames(xx)),]
      if(length(sel_samples)!=nrow(out)){
         warning('No. of selected samples (',length(sel_samples),') does not match the number of rows in the output (',nrow(out),')')
      }
      return(out)
   }
   
   ## Output
   if(is.data.frame(x) | is.matrix(x)){
      return( main(x) )
   }
   
   if(is.list(x)){
      return( lapply(x, main) )
   }
   
   stop('`x` must be a dataframe, matrix, or list')
}

## Long form data frame containing sample and sample group info
getSampleGroups <- function(
   label.group.var='cancer_type',
   comparison.group.var='cancer_stage', 
   sample.metadata=sample_metadata, 
   switch.comparison.group.var=F, rm.msi.subtypes=T,
   output='group_info', min.group.size=NULL, show.label.group.size=T,
   verbose=2
){
   
   if(F){
      ## Args
      sample.metadata=sample_metadata
      label.group.var='cancer_type_code_supertype_subtype'
      #label.group.var='cancer_type_code_subtype_subset'
      
      #comparison.group.var='metastatic_location'
      comparison.group.var='cancer_stage'
      #comparison.group.var='metastatic_progression_code'
      
      switch.comparison.group.var=F
      rm.msi.subtypes=T
      
      output='group_info'
      min.group.size=NULL
      show.label.group.size=T
      verbose=T
   }
   
   ## Init --------------------------------
   require(reshape2)
   
   ## Initialize new metadata table for editing
   sample_metadata_ss <- sample.metadata
   
   ## Check group vars
   if(!(output %in% c('group_info','group_counts'))){
      stop('`output` must be group_info or group_counts')
   }
   
   if(!is.factor(sample_metadata_ss[[comparison.group.var]])){
      stop('`',comparison.group.var,'` must be a factor')
   }
   
   if(!is.factor(sample_metadata_ss[[label.group.var]])){
      stop('`',label.group.var,'` must be a factor')
   }
   
   ## Remove blacklisted samples
   sample_metadata_ss <- subset(sample_metadata_ss, !is_blacklisted)
   
   ## Other variables
   MIN_GROUP_SIZE <- if(!is.null(min.group.size)){ min.group.size } else { 0 }
   group_info <- NULL
   
   ## Announce sample group name
   if(verbose){
      message(
         '\nSelecting samples for: ',comparison.group.var,' x ',label.group.var
      )
   }
   
   ## Reverse comparison.group.var case and ctrl levels --------------------------------
   if(switch.comparison.group.var){
      new_comparison_var_levels <- (function(){
         v <- sample_metadata_ss[[comparison.group.var]]
         v <- c(
            levels(v)[2:length(levels(v))],
            levels(v)[1]
         )
         return(v)
      })()
      
      sample_metadata_ss[[comparison.group.var]] <- factor(
         sample_metadata_ss[[comparison.group.var]], 
         new_comparison_var_levels
      )
   }
   
   ## Helper functions --------------------------------
   ## Split samples into groups
   splitSamples <- function(output='data.frame', min.group.size=0){
      
      if(!(output %in% c('data.frame','list'))){
         stop('`output` must be data.frame or list')
      }
      
      ## Get sample ids
      l1 <- split(sample_metadata_ss, sample_metadata_ss[[label.group.var]]) ## samples with NA `label.group.var` are automatically removed upon split
      l2 <- lapply(l1, function(i){
         #i=l1[[1]]
         l_i <- split(i, i[[comparison.group.var]])
         groups <- lapply(l_i,`[[`,'sample_id')
         
         ## Assign NULL to groups with too few samples
         group_sizes <- sapply(groups, length)
         if(any(group_sizes<min.group.size)){ return(NULL) }
         
         return(groups)
      })
      
      ## Remove comparisons with too few samples
      l2 <- l2[!sapply(l2, is.null)]
      
      ## Rm zero length lists
      l2 <- l2[ sapply(l2, length)>0 ]
      
      ## Output nested list of sample ids
      if(output=='list'){ return(l2) }
      
      ## Output long form data frame
      group_info <- reshape2::melt(l2)
      colnames(group_info) <- c('sample_id','comparison_group','label_group')
      
      return(group_info)
   }
   
   ## Split samples into e.g.: Primary vs metastatic sites [Distant, Local, Lymph]
   splitSamplesOneVsRest <- function(ctrl.comparison.group.name='Primary', min.group.size=0){

      ## Get rest group names, e.g. [Distant, Local, Lymph]
      rest_comparison_group_names <- levels(sample_metadata_ss[[comparison.group.var]])
      rest_comparison_group_names <- rest_comparison_group_names[rest_comparison_group_names!=ctrl.comparison.group.name]
      
      nested_sample_list <- splitSamples(output='list')
      
      group_info <- lapply(names(nested_sample_list), function(i){
         #i=names(nested_sample_list)[[1]]
         i_l <- lapply(rest_comparison_group_names, function(j){
            #j=rest_comparison_group_names[[1]]
            pair <- c(
               nested_sample_list[[i]][j], 
               nested_sample_list[[i]][ctrl.comparison.group.name]
            )
            
            names(pair) <- c(comparison.group.var, ctrl.comparison.group.name)
            
            ##
            if(switch.comparison.group.var){
               pair <- pair[c(ctrl.comparison.group.name, comparison.group.var)]
            }
            
            ## Assign NULL to groups with too few samples
            group_sizes <- sapply(pair, length)
            if(any(group_sizes<min.group.size)){ return(NA) }
            
            return(pair)
         })
         
         if(!switch.comparison.group.var){
            names(i_l) <- paste0(i,', ',rest_comparison_group_names, ' vs ', ctrl.comparison.group.name)
         } else {
            names(i_l) <- paste0(i,', ',ctrl.comparison.group.name, ' vs ', rest_comparison_group_names)
         }
         
         return(i_l)
      })
      names(group_info) <- names(nested_sample_list)
      
      group_info <- reshape2::melt(group_info)
      
      ## Assign colnames
      colnames_new <- c(value='sample_id', L1='label_group_raw', L2='label_group', L3='comparison_group')
      colnames(group_info) <- colnames_new[ colnames(group_info) ]
      group_info$label_group_raw <- NULL
      
      ## Remove comparisons with few samples
      group_info <- group_info[!is.na(group_info$comparison_group),]
      
      return(group_info)
   }
   
   ## Caner type labels --------------------------------
   valid_label_group_vars <- c(
      'cancer_type','cancer_type_code',
      'cancer_type_subtype','cancer_type_code_subtype','cancer_type_code_subtype_subset',
      'cancer_type_code_supertype_subtype', 'cancer_type_and_code', 
      'cancer_type_code_subtype_subset_hr_status'
   )
   
   ## Check for valid cancer types
   if(!(label.group.var %in% valid_label_group_vars)){
      stop("`label.group.var` must be one of the following: ",paste(valid_label_group_vars, collapse=', '))
   }
   
   ## Remove MSI subtypes
   if(rm.msi.subtypes & grepl('subtype', label.group.var)){
      sample_metadata_ss <- sample_metadata_ss[
         !grepl('MSI', sample_metadata_ss$cancer_type_code_subtype)
      ,]
   }
   
   ## Groups --------------------------------
   ## cohort x cancer type
   if(label.group.var %in% valid_label_group_vars & comparison.group.var %in% c('cohort','cancer_stage')){
      #if(is.null(min.group.size)){ MIN_GROUP_SIZE <- 15 }
      group_info <- splitSamples(min.group.size=MIN_GROUP_SIZE)
      
   # } else if(label.group.var %in% valid_label_group_vars & comparison.group.var %in% c('cohort','cancer_stage')){
   #    ## Rm samples with unknown subtype
   #    sample_metadata_ss <- subset(sample_metadata_ss, !is.na(cancer_subtype))
   # 
   #    ##
   #    #if(is.null(min.group.size)){ MIN_GROUP_SIZE <- 5 }
   #    group_info <- splitSamples(min.group.size=MIN_GROUP_SIZE)
   
   ## metastatic_location x cancer type
   } else if(label.group.var %in% valid_label_group_vars & comparison.group.var=='metastatic_location'){
      ## Rm samples with unknown subtype
      sample_metadata_ss <- subset(sample_metadata_ss, !is.na(metastatic_location)) 
      
      ##
      #if(is.null(min.group.size)){ MIN_GROUP_SIZE <- 5 }
      group_info <- splitSamplesOneVsRest(
         ctrl.comparison.group.name='Primary', 
         min.group.size=MIN_GROUP_SIZE
      )
      
   ## progression status x cancer type
   } else if(label.group.var %in% valid_label_group_vars & comparison.group.var %in% c('progression_status','progression_status_code')){
   
      ## Rm samples with unknown progression status
      if(comparison.group.var=='progression_status'){
         ctrl_comparison_group_name <- 'metastatic'
         sample_metadata_ss <- subset(sample_metadata_ss, !is.na(progression_status)) 
      } else if(comparison.group.var=='progression_status_code'){
         ctrl_comparison_group_name <- 'MET'
         sample_metadata_ss <- subset(sample_metadata_ss, !is.na(progression_status_code)) 
      }
      
      #if(is.null(min.group.size)){ MIN_GROUP_SIZE <- 5 }
      group_info <- splitSamplesOneVsRest(
         ctrl.comparison.group.name=ctrl_comparison_group_name, 
         min.group.size=MIN_GROUP_SIZE
      )
   }
   
   if(is.null(group_info)){
      stop('Invalid `label.group.var` and `comparison.group.var` combinations')
   }
   
   ## Output --------------------------------
   ## Tabulate group counts
   group_counts <- unclass(table(group_info$label_group, group_info$comparison_group))
   
   ## Add e.g. (10 | 21) to the end of label group
   if(show.label.group.size){
      new_label_groups <- structure(
         paste0(rownames(group_counts),' (',group_counts[,1],' | ',group_counts[,2],')'),
         names=rownames(group_counts)
      )
      group_info$label_group <- new_label_groups[ as.character(group_info$label_group) ]
   }
   
   ## Force (original) factor levels
   group_info$label_group <- factor2(group_info$label_group)
   group_info$comparison_group <- factor2(group_info$comparison_group)
   
   if(output=='group_counts'){ 
      return(group_counts) 
   }
   
   ## Info
   if(verbose==1){ 
      message(
         'Total samples: ',nrow(sample_metadata_ss),'\n',
         'Comparisons: ',nrow(group_counts),'\n',
         paste(colnames(group_counts), collapse=' vs ')
      )
   }
   if(verbose==2){ print(group_counts) }
   
   return(group_info)
}

## Add sample values from a matrix or data frame to sample group info data frame
getSampleData <- function(group.info, data, out.format='long'){
   if(F){
      group.info=sample_groups$cancer_type_code.metastatic_location
      data=mut_type_load
   }
   
   ## Checks --------------------------------
   if(!is.matrix(data) & !is.data.frame(data)){
      stop('`data` must be a matrix or data.frame')
   }
   
   if(!(all(group.info$sample_id %in% rownames(data)))){ 
      stop('Some sample ids in `group.info` are not in `data`') 
   }
   
   if(!(out.format %in% c('long','wide'))){
      stop("`out.format` must be 'long' or 'wide'")
   }
   
   ## Match sample order of `data` to `group.info` --------------------------------
   data_matched <- as.data.frame(data, check.names=F)
   data_matched <- cbind(sample_id=rownames(data_matched), data_matched, row.names=NULL)
   data_matched <- data_matched[match(group.info$sample_id, data_matched$sample_id),]
   
   ## Append feature values to `group.info` --------------------------------
   feature_names <- colnames(data)
   
   if(out.format=='wide'){
      return(
         cbind(group.info, data_matched[,feature_names] )
      )
   }
   
   feature_names <- colnames(data)
   l <- lapply(feature_names, function(i){
      cbind(group.info, feature=i, value=data_matched[[i]])
   })
   do.call(rbind, l)
}

## Sample groups ================================
sample_groups <- list()
MIN_GROUP_SIZE_LARGE <- 15
MIN_GROUP_SIZE_SMALL <- 5

## By cancer types --------------------------------
sample_groups$cancer_stage.cancer_type <- getSampleGroups(
   label.group.var='cancer_type', 
   comparison.group.var='cancer_stage',
   min.group.size=MIN_GROUP_SIZE_LARGE
)

# if(F){
#    tab <- getSampleGroups(
#       label.group.var='cancer_type', 
#       comparison.group.var='cancer_stage',
#       min.group.size=MIN_GROUP_SIZE_LARGE,
#       output='group_counts'
#    )
#    sum(tab)
#    colSums(tab)
# }

sample_groups$cancer_stage.cancer_type_and_code <- (function(){
   df <- getSampleGroups(
      label.group.var='cancer_type_and_code', 
      comparison.group.var='cancer_stage',
      min.group.size=MIN_GROUP_SIZE_LARGE
   )
   df$label_group <- sub('\\) \\(', '; ', df$label_group)
   df$label_group <- factor2(df$label_group)
   return(df)
})()

sample_groups$cancer_stage.cancer_type_code_subtype <- getSampleGroups(
   label.group.var='cancer_type_code_subtype', 
   comparison.group.var='cancer_stage',
   min.group.size=MIN_GROUP_SIZE_SMALL
)

sample_groups$cancer_stage.cancer_type_code_subtype_subset <- getSampleGroups(
   label.group.var='cancer_type_code_subtype_subset', 
   comparison.group.var='cancer_stage',
   min.group.size=MIN_GROUP_SIZE_SMALL
)

sample_groups$cancer_stage.cancer_type_code_supertype_subtype <- getSampleGroups(
   label.group.var='cancer_type_code_supertype_subtype', 
   comparison.group.var='cancer_stage',
   min.group.size=0 ## No need to provide threshold. All cancer supertypes with <15 samples already blacklisted in metadata
)

## By metastatic location --------------------------------
sample_groups$metastatic_location.cancer_type_code <- getSampleGroups(
   label.group.var='cancer_type_code', 
   comparison.group.var='metastatic_location',
   min.group.size=MIN_GROUP_SIZE_SMALL
)

sample_groups$metastatic_location.cancer_type_code_subtype_subset <- getSampleGroups(
   label.group.var='cancer_type_code_subtype_subset', 
   comparison.group.var='metastatic_location',
   min.group.size=MIN_GROUP_SIZE_SMALL
)

## By primary progression --------------------------------
sample_groups$progression_status_code.cancer_type_code <- getSampleGroups(
   label.group.var='cancer_type_code', 
   comparison.group.var='progression_status_code', 
   switch.comparison.group.var=T,
   min.group.size=MIN_GROUP_SIZE_SMALL
)

## Misc --------------------------------
sample_groups$cancer_stage.cancer_type_code_subtype_subset_hr_status <- getSampleGroups(
   label.group.var='cancer_type_code_subtype_subset_hr_status', 
   comparison.group.var='cancer_stage',
   min.group.size=MIN_GROUP_SIZE_SMALL
)

sample_groups$PRAD_MET.cohort <- getSampleGroups(
   label.group.var='cancer_type',
   comparison.group.var='cohort',
   sample.metadata=subset(sample_metadata, cancer_type_code=='PRAD' & progression_status_code=='MET'),
   min.group.size=MIN_GROUP_SIZE_SMALL
)

## Misc functions ================================
## Convert to HMF ids for Hartwig samples
anonymizeSampleIds <- function(v){
   #v <- rownames(sig_contribs_raw)
   is_pcawg_sample <- grepl('^DO',v)
   sample_id_2 <- getSampleMetadata(v, 'sample_id_2')
   v[!is_pcawg_sample] <- sample_id_2[!is_pcawg_sample]
   return(v)
}

## Save tables locally
cacheAndReadData <- function(
   remote.path, local.path=NULL,
   cache.dir=path.expand('~/Documents/R_cache/'),
   overwrite=F
){
   
   #remote.path=paste0(base_dir,'/misc/processed/Chromatin_modifiers/scripts/annotate_svs/vis_sv_data_compact.txt.gz')
   
   ## Init --------------------------------
   if(!dir.exists(cache.dir)){
      dir.create(cache.dir, recursive=T)
   }
   
   ext <- c('.rds','.txt','.txt.gz','.csv','.csv.gz')
   regex <- gsub('[.]','[.]',ext)
   regex <- paste0(regex,'$')
   
   is_valid_ext <- sapply(regex, function(i){ grepl(i, remote.path) })
   if(all(!is_valid_ext)){
      stop('File extension must be one of the following:\n  ', paste(ext, collapse=', '))
   }
   
   ext <- ext[ which(is_valid_ext) ]
   regex <- regex[ which(is_valid_ext) ]
   
   set.seed(nchar(remote.path))
   
   ## Copy --------------------------------
   if(dir.exists('/hpc/')){
      local.path <- remote.path
   } else {
      if(is.null(local.path)){
         local.path <- paste0(
            cache.dir,
            sub(regex,'',basename(remote.path)),
            '.',paste(sample(letters, 8), collapse=''),ext
         )
      }
      
      local.path <- gsub('/+', '/', local.path)
      if(!file.exists(local.path) | overwrite){
         if(!file.exists(local.path)){
            message('Making local copy: ', local.path)
         } else {
            message('Updating local copy: ', local.path)
         }
         
         system(sprintf(
            'rsync -a %s %s',
            remote.path,
            local.path
         ))
      }
   }
   
   ## Read --------------------------------
   message('Reading local copy: ', local.path)
   if(ext=='.rds'){
      readRDS(local.path)
   } else if(grepl('txt',ext)){
      read.delim(local.path, check.names=F)
   } else {
      read.csv(local.path, check.names=F)
   }
}

## Relative contribution
relContrib <- function(m, row.sums=NULL){
   #m <- all_svs
   
   if(is.null(row.sums)){ row.sums <- rowSums(m) }
   m <- m/row.sums
   m[is.na(m)] <- 0
   return(m)
}

## Enrichment functions ================================
enrTrueVsFalse <- function(
   m, bool=NULL, bool.var.name='is_hmf_sample', alternative='two.sided', 
   include.pancancer.group=FALSE, p.adjust.per='none', #p.adjust.per.feature=FALSE,
   show.sample.size=T, order.by.pvalue=F, avg.numeric.func='median', show.conting=FALSE
){
   if(F){
      m=all_svs
      m=m_viral_ins
      m=sig_contribs.renamed.rel
      alternative='greater';
      bool=NULL
      bool.var.name='is_hmf_sample'
      show.sample.size=T; order.by.pvalue=F; avg.numeric.func='median'; show.conting=F; include.pancancer.group=F; p.adjust.per.feature=F
   }
   
   if(is.numeric(as.matrix(m))){
      data_type <- 'numeric'
   } else if(is.logical(as.matrix(m))) {
      data_type <- 'logical'
   } else {
      stop('m must contain only numeric or logical data')
   }
   
   ## Add metadata --------------------------------
   feature_names <- colnames(m)
   m <- data.frame(m, check.names=F)
   
   m$tissue_cancer_type <- getSampleMetadata(rownames(m), 'tissue_cancer_type')
   
   if(!is.null(bool)){
      m$bool <- bool
   } else if(!is.null(bool.var.name)){
      m$bool <- getSampleMetadata(rownames(m), bool.var.name)
   }
   
   if(anyNA(m$bool)){
      warning('Boolean variable contains NA values. Samples with NA boolean values will be removed')
      m <- m[!is.na(m$bool),]
   }
   
   m$tissue_cancer_type <- as.factor(m$tissue_cancer_type)
   if(include.pancancer.group){
      m <- rbind(
         m,
         within(m, { tissue_cancer_type <- 'Pancancer' })
      )
   }
   
   ## Calculate stats --------------------------------
   m_split <- split(m, m$tissue_cancer_type)
   
   enr <- lapply(names(m_split), function(i){
      #i=names(m_split)[1]
      #i='BRCA'
      #print(i)
      m_ss <- m_split[[i]]
      out <- univarFeatSel.default(
         x=m_ss[,feature_names,drop=F],
         y=m_ss[,'bool'],
         alternative=alternative,
         show.sample.size=show.sample.size, order.by.pvalue=F, avg.numeric.func=avg.numeric.func,
         show.conting=show.conting
      )
      cbind(tissue_cancer_type=i, out)
   })
   enr <- do.call(rbind, enr)

   NULL -> enr$is_pass_feature -> enr$is_keep_feature
   
   ##
   if(!is.null(p.adjust.per)){
      if(p.adjust.per=='none'){
         enr$qvalue <- p.adjust(enr$pvalue, method='bonferroni')
      } else if(p.adjust.per=='feature'){
         enr$index <- 1:nrow(enr)
         enr <- do.call(rbind, lapply(split(enr, enr$feature), function(i){
            i$qvalue <- p.adjust(i$pvalue, method='bonferroni')
            return(i)
         }))
         enr <- enr[order(enr$index),]
         enr$index <- NULL
      } else if(p.adjust.per=='cancer_type') {
         enr$index <- 1:nrow(enr)
         enr <- do.call(rbind, lapply(split(enr, enr$tissue_cancer_type), function(i){
            i$qvalue <- p.adjust(i$pvalue, method='bonferroni')
            return(i)
         }))
         enr <- enr[order(enr$index),]
         enr$index <- NULL
      }
   }
   
   ##
   if(order.by.pvalue){ enr <- enr[order(enr$pvalue),] }
   
   return(enr)
}

## --------------------------------
plotEnr <- function(
   ## Inputs
   enr=NULL, m=NULL,
   
   ## Enrichment args
   bool=NULL, bool.var.name='is_hmf_sample', alternative='two.sided', 
   include.pancancer.group=FALSE, p.adjust.per='none',
   
   ## Plot args
   thres.pvalue=0.01, thres.qvalue=0.05, thres.eff.size=0.1, show.p.signif=FALSE, rm.non.signif=FALSE,
   value.name='Average', false.name='primary',true.name='metastatic',
   feat.type.parse.func=NULL, plot.type='moon', facet.scales='free_x', bool.colors=c('FALSE'='#F58134', 'TRUE'='#6B469A')
){
   if(F){
      m <- cbind(
         SV_LOAD=sv_load,
         DEL=rowSums(dels),
         DUP=rowSums(dups),
         COMPLEX=rowSums(complex_svs),
         LINE=all_svs[,'LINE']
      )
      
      m <- log10(m+1)
      
      #enr <- enrTrueVsFalse(log10(m+1), bool.var.name='is_hmf_sample')
      bool.var.name='is_hmf_sample'; bool.colors=c('FALSE'='#F58134', 'TRUE'='#6B469A')
      alternative='greater'
      thres.pvalue=0.01; thres.qvalue=0.01; thres.eff.size=0.1; show.p.signif=F
      value.name='Average'; false.name='primary';true.name='metastatic';
      feat.type.parse.func=NULL; plot.type='moon'; facet.scales='free'; include.pancancer.group=F;
      rm.non.signif=F
   }
   
   require(ggplot2)
   
   if(!(plot.type %in% c('moon','bar','violin','volcano'))){
      stop("Plot type must be 'moon', 'bar', 'violin', or 'volcano")
   }
   
   ## Enrichment --------------------------------
   if(is.null(enr) & is.null(m)){ stop('Either `enr` or `m` must be specified') }
   if(is.null(enr)){
      enr <- suppressWarnings({
         enrTrueVsFalse(
            m, bool=bool, bool.var.name=bool.var.name, 
            alternative=alternative, 
            include.pancancer.group=include.pancancer.group, 
            p.adjust.per.feature=p.adjust.per
         )
      }) 
   }
   
   ## Format enr table to long form for plotting --------------------------------
   pd <- (function(){
      template <- enr[,c('tissue_cancer_type','feature','pvalue','qvalue','eff_size')]
      
      case <- structure(
         enr[,c('avg_case','n_case')],
         names=c('avg','n')
      )
      
      ctrl <- structure(
         enr[,c('avg_ctrl','n_ctrl')],
         names=c('avg','n')
      )
      
      rbind(
         cbind(bool=TRUE, template, case),
         cbind(bool=FALSE,template, ctrl)
      )
   })()
   
   ##  --------------------------------
   data_type <- enr$feature_type[1]
   if(data_type=='numeric'){ 
      test_name <- 'Wilcox test'
      effsize_name <- "Cliff's delta"
      #effsize_letter <- "D"
   } else {
      test_name <- 'Fisher test'
      effsize_name <- "Cramer's V"
      #effsize_letter <- "V"
   }
   
   ## Significance --------------------------------
   ##
   signif_types <- c(
      'ns',
      
      sprintf('%s>%s; p<%s; effsize>=%s',true.name,false.name,thres.pvalue,thres.eff.size),
      sprintf('%s>%s; q<%s; effsize>=%s',true.name,false.name,thres.qvalue,thres.eff.size),
      
      sprintf('%s>%s; p<%s; effsize>=%s',false.name,true.name,thres.pvalue,thres.eff.size),
      sprintf('%s>%s; q<%s; effsize>=%s',false.name,true.name,thres.qvalue,thres.eff.size)
   )
   
   pd$signif_type <- signif_types[1]
   
   if(show.p.signif){
      pd$signif_type[pd$pvalue<thres.pvalue & pd$eff_size>=thres.eff.size] <- signif_types[2]
      pd$signif_type[pd$pvalue<thres.pvalue & pd$eff_size<=thres.eff.size] <- signif_types[4]
   }
   
   pd$signif_type[pd$qvalue<thres.qvalue & pd$eff_size>=thres.eff.size] <- signif_types[3]
   pd$signif_type[pd$qvalue<thres.qvalue & pd$eff_size<=thres.eff.size] <- signif_types[5]
   
   pd$signif_type <- droplevels(factor(pd$signif_type, signif_types))
   
   ##
   pd$signif_label <- (function(){
      df <- list()
      df$p_label <- sprintf(
         'p=%s\nES=%s',
         signif(pd$pvalue,2),
         round(pd$eff_size,2)
      )
      
      df$q_label <- sprintf(
         'q=%s\nES=%s',
         signif(pd$qvalue,2),
         round(pd$eff_size,2)
      )
      df <- as.data.frame(df)
      
      df$label <- ''
      if(show.p.signif){
         p_sel <- pd$pvalue<thres.pvalue & abs(pd$eff_size)>=thres.eff.size
         df$label[ p_sel ] <- df$p_label[ p_sel ]
      }
      q_sel <- pd$qvalue<thres.qvalue & abs(pd$eff_size)>=thres.eff.size
      df$label[ q_sel ] <- df$q_label[ q_sel ]
      
      return(df$label)
   })()
   
   ## --------------------------------
   pd$feature <- factor(pd$feature, unique(pd$feature))
   pd$tissue_group <- sample_metadata$tissue_group[ match(pd$cancer_type_code, sample_metadata$cancer_type_code) ]
   
   if(include.pancancer.group){
      pd$tissue_group[pd$cancer_type_code=='Pancancer'] <- 'Pancancer'
      pd$tissue_group <- as.factor(pd$tissue_group)
      pd$tissue_group <- relevel(pd$tissue_group,'Pancancer')
   }
   
   if(rm.non.signif){
      enriched_features <- unique(pd[pd$signif_type!='ns','feature'])
      pd <- pd[pd$feature %in% enriched_features,]
   }
   
   if(!is.null(feat.type.parse.func)){
      pd$feature_type <- feat.type.parse.func(pd$feature)
      pd$feature_type <- factor(pd$feature_type, unique(pd$feature_type))
   }
   
   ## Plot params --------------------------------
   ##
   xlabels <- (function(){
      df <- unique( pd[,c('bool', 'cancer_type_code', 'n')] )
      df <- reshape2::dcast(df, cancer_type_code~bool, value.var='n')
      colnames(df) <- c('cancer_type_code','false','true')
      v <- with(df, paste0(cancer_type_code,'\n(',false,' vs ',true,')'))
      names(v) <- df$cancer_type_code
      return(v)
   })()
   
   xtitle <- sprintf('Cancer type\n(%s vs %s)',false.name,true.name)
   
   ## Plot --------------------------------
   ## Moon plot
   if(plot.type=='moon'){
      pd$eff_size[is.na(pd$eff_size)] <- 0
      pd$pvalue[is.na(pd$pvalue)] <- 1
      pd$pvalue[is.na(pd$qvalue)] <- 1
      pd$avg[is.na(pd$avg)] <- 0
      
      p <- ggplot(subset(pd, tissue_group=='GI_dev'), aes(x=cancer_type_code, y=feature))
      p <- ggplot(pd, aes(x=cancer_type_code, y=feature))
      
      if(!is.null(feat.type.parse.func)){
         p <- p + facet_grid(feature_type~tissue_group, scales='free',space='free', switch='y')
      } else {
         p <- p + facet_grid(~tissue_group, scales='free',space='free', switch='y')
      }
      
      ##
      aes_values <- data.frame(signif_type=signif_types)
      aes_values$color <- c('darkgrey','red','red','blue','blue')
      aes_values$linetype <- c('solid','twodash','solid','twodash','solid')
      aes_values <- aes_values[aes_values$signif_type %in% levels(pd$signif_type),]
      
      ##
      p <- p + 
         gggibbous::geom_moon( 
            aes(ratio=0.5, right=bool, fill=avg, size=abs(eff_size), color=signif_type, linetype=signif_type),
            stroke=0.3
         ) +
         #geom_point( aes(size=abs(eff_size), color=signif_type), shape=21, fill=NA, stroke=0.3) +
         
         ##
         scale_fill_distiller(
            name=value.name, 
            palette='Spectral'
         ) +
         scale_size_continuous(
            name=paste0('Eff size (',effsize_name,')'), 
            range=c(0.1,7)
         ) +
         scale_color_manual(
            name=paste0('Signif (',test_name,'; ',effsize_name,')'), 
            values=structure(aes_values$color, names=aes_values$signif_type),
            drop=TRUE
         ) +
         scale_linetype_manual(
            name=paste0('Signif (',test_name,'; ',effsize_name,')'), 
            values=structure(aes_values$linetype, names=aes_values$signif_type),
            drop=TRUE
         ) +
         guides(
            fill=guide_colorbar(order=1, barheight=4, frame.colour='black'),
            size=guide_legend(order=2),
            color=guide_legend(order=3, override.aes=list(size=5)),
            linetype=guide_legend(order=3)
         ) +
         #ggtitle('test') +
         
         ## Axes
         scale_y_discrete(limits=rev, position='right', name='Feature') +
         scale_x_discrete(labels=xlabels, name=sprintf('Cancer type\n(%s vs %s)',false.name,true.name)) +
         
         theme_bw() +
         theme(
            panel.grid.minor=element_blank(),
            panel.spacing=unit(0,'pt'),
            strip.text.x=element_text(angle=90, hjust=0, vjust=0.5),
            axis.text.x.top=element_text(angle=90, hjust=0, vjust=0.5),
            axis.text.x.bottom=element_text(angle=90, hjust=1, vjust=0.5)
         )
      
      return(p)
   }
   
   if(plot.type=='volcano'){
      ## Prep data
      signif_thres <- if(show.p.signif){ thres.pvalue } else { thres.qvalue }
      #es_thres <- thres.eff.size
      pd <- subset(pd, bool==TRUE)
      pd$point_label <- paste0(pd$cancer_type_code,', ',pd$feature)
      pd$point_label[pd$signif_type=='ns'] <- ''
      
      ## Colors
      gg_color_hue <- function(n) {
         hues = seq(15, 375, length = n + 1)
         hcl(h = hues, l = 65, c = 100)[1:n]
      }
      
      enriched_features <- unique(subset(pd, signif_type!='ns', feature, drop=T))
      not_enriched_features <- droplevels(unique(pd$feature[!(pd$feature %in% enriched_features)]))
      
      color_pal <- c(
         structure(gg_color_hue(length(enriched_features)), names=as.character(enriched_features)),
         structure(rep('grey', length(not_enriched_features)), names=levels(not_enriched_features))
      )
      
      ## Plot
      if(show.p.signif){
         p <- ggplot(pd, aes(x=eff_size, y=-log10(pvalue)))
      } else {
         p <- ggplot(pd, aes(x=eff_size, y=-log10(qvalue)))
      }
      p <- p +
         geom_hline(yintercept=-log10(signif_thres), linetype='dotted') +
         geom_vline(xintercept=c(-thres.eff.size, thres.eff.size), linetype='dotted') +
         geom_point(aes(color=feature)) +
         scale_color_manual(values=color_pal) +
         ggrepel::geom_text_repel(aes(label=point_label, color=feature), size=2.5, max.overlaps=25) +
         guides(color='none') +
         xlab(effsize_name) +
         ggtitle(bool.var.name) +
         theme_bw()
      
      return(p)
   }
   
   if(plot.type=='bar'){
      ## Bar plot
      p <- ggplot(pd, aes(x=cancer_type_code, y=avg))
      
      if(length(unique(pd$feature))>1){
         p <- p + facet_grid(feature~tissue_group, scales=facet.scales, space='free_x')
      } else {
         p <- p + facet_grid(.~tissue_group, scales=facet.scales, space='free')
      }
      
      p <- p +
         geom_bar(aes(fill=bool), stat='identity', position='dodge', size=0.2, color='black') +
         scale_fill_manual(values=c('FALSE'='#f9b385', 'TRUE'='#a690c2')) +
         
         geom_text(
            data=subset(pd, bool==TRUE), 
            mapping=aes(label=signif_label, y=Inf), 
            angle=90, hjust=1, vjust=0.5, size=2.5
         ) +
         
         scale_x_discrete(labels=xlabels, name=xtitle) +
         ylab(value.name) +
         xlab(xtitle) +
         
         theme_bw() +
         theme(
            panel.grid.minor=element_blank(),
            panel.spacing.x=unit(0,'pt'),
            strip.text.x=element_text(angle=90, hjust=0, vjust=0.5),
            strip.text.y=element_text(angle=0, hjust=0, vjust=0.5),
            axis.text.x.top=element_text(angle=90, hjust=0, vjust=0.5),
            axis.text.x.bottom=element_text(angle=90, hjust=1, vjust=0.5),
            legend.position='none'
         )
      
      return(p)
   }
   
   if(plot.type=='violin'){
      ## Matrix to long dataframe
      pd_raw <- reshape2::melt(as.matrix(m))
      colnames(pd_raw) <- c('sample','feature','value')
      
      ## Get metadata
      pd_raw <- cbind(
         getSampleMetadata(pd_raw$sample, c('cancer_type_code','tissue_group',bool.var.name)),
         pd_raw
      )
      colnames(pd_raw)[colnames(pd_raw)==bool.var.name] <- 'bool'
      
      pd_raw <- pd_raw[!is.na(pd_raw$bool),]
      
      pd_raw$cancer_type_code <- as.factor(pd_raw$cancer_type_code)
      pd_raw$tissue_group <- as.factor(pd_raw$tissue_group)
      
      if(include.pancancer.group){
         pd_raw <- rbind(
            pd_raw,
            within(pd_raw,{
               cancer_type_code <- 'Pancancer'
               tissue_group <- 'Pancancer'
            })
         )
         
         pd_raw$cancer_type_code <- relevel(as.factor(pd_raw$cancer_type_code),'Pancancer')
         pd_raw$tissue_group <- relevel(as.factor(pd_raw$tissue_group),'Pancancer')
      }
      
      pd_raw <- pd_raw[order(pd_raw$tissue_group, pd_raw$cancer_type_code),]
      
      #gghalves::geom_half_violin(data=pd_l, side='l', size=0.3, width=0.8) +
      #gghalves::geom_half_violin(data=pd_r, side='r', size=0.3, width=0.8) +
      
      ## Plot
      p <- ggplot(pd_raw, aes(x=cancer_type_code, y=value, fill=bool))
      if(length(unique(pd$feature))>1){
         p <- p + facet_grid(feature~tissue_group, scales=facet.scales, space='free_x')
      } else {
         p <- p + facet_grid(.~tissue_group, scales=facet.scales, space='free')
      }
      
      p <- p +
         gghalves::geom_half_violin(data=subset(pd_raw, bool==FALSE), side='l', size=0.3, width=0.8, draw_quantiles=0.5) +
         gghalves::geom_half_violin(data=subset(pd_raw, bool==TRUE ), side='r', size=0.3, width=0.8, draw_quantiles=0.5) +
         #ggbeeswarm::geom_beeswarm(dodge.width=0.8, corral='wrap', corral.width=0.2, size=0.1, color='grey') +
         
         #ggbeeswarm::geom_beeswarm(dodge.width=0.8, corral='wrap', corral.width=0.2, size=0.1, color='grey') +
         #geom_boxplot(outlier.shape=NA, alpha=0.5) +
         
         scale_x_discrete(labels=xlabels, name=sprintf('Cancer type\n(%s vs %s)',false.name,true.name)) +
         scale_fill_manual(values=c('FALSE'='#f9b385', 'TRUE'='#a690c2')) +
         ylab(value.name) +
         
         geom_text(
            data=subset(pd, bool==TRUE), 
            mapping=aes( label=signif_label, y=Inf, color=as.factor(sign(eff_size)) ), 
            vjust=1.05, #angle=90, hjust=1, vjust=0.5, 
            size=2.5 #, fontface='bold', color='black'
         ) +
         scale_color_manual(values=c('-1'='blue','0'='black','1'='red')) +
         
         theme_bw() +
         theme(
            panel.grid.minor=element_blank(),
            panel.spacing.x=unit(0,'pt'),
            strip.text.x=element_text(angle=90, hjust=0, vjust=0.5),
            strip.text.y=element_text(angle=0, hjust=0, vjust=0.5),
            axis.text.x.top=element_text(angle=90, hjust=0, vjust=0.5),
            axis.text.x.bottom=element_text(angle=90, hjust=1, vjust=0.5),
            legend.position='none'
         )
      
      return(p)
   }
}


# 
# ####
# if(F){
#    df <- subset(
#       sample_metadata, 
#       !is_blacklisted_sample,
#       c(sample_id, patient_id, cancer_type_code, sbs_load, sbs_load.clonal)
#    )
#    df$sbs_clonal_fraction <- df$sbs_load.clonal / df$sbs_load
#    patient_ids.multi_biopsy <- df$patient_id[duplicated(df$patient_id)]
#    df <- subset(df, patient_id %in% patient_ids.multi_biopsy)
#    df_split <- split(df, df$patient_id)
#    l <- lapply(df_split, function(i){
#       #i=df_split$HMF000809
#       outer(i$sbs_clonal_fraction, i$sbs_clonal_fraction, '-') |> as.dist() |> as.vector() |> abs()
#    })
#    
#    mean(unlist(l))
#    sd(unlist(l))
# }

