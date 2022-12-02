#!/usr/bin/env Rscript
library(GenomicRanges)
library(VariantAnnotation)
library(reshape2)
library(dplyr)
library(stringr)
library("ggplot2")
library(rjson)
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
chromosomes <- seqnames(get(ref_genome))[1:23]
options(stringsAsFactors = F)

myCols <- c("darkgray","steelblue3","red3","lightgray")
names(myCols) <- c("neutral1","gain","loss","neutral2")



cramerV.default <- function(case.true, case.false, ctrl.true, ctrl.false, show.sign=T){
  
  ## Vectorized version
  observed <- cbind(case.true, case.false, ctrl.true, ctrl.false)
  n <- rowSums(observed)
  
  rowsums1 <- case.true  + ctrl.true
  rowsums2 <- case.false + ctrl.false
  colsums1 <- case.true  + case.false
  colsums2 <- ctrl.true  + ctrl.false
  
  expected <- cbind(
    case.true  = rowsums1*colsums1,
    case.false = rowsums2*colsums1,
    ctrl.true  = rowsums1*colsums2,
    ctrl.false = rowsums2*colsums2
  ) / n
  
  ## Chi-squared statistic
  chi2 <- rowSums(
    (observed - expected)^2 / expected
  )
  
  ## Calculate Cramers V --------------------------------
  ## Based on this tutorial: https://www.real-statistics.com/chi-square-and-f-distributions/effect-size-chi-square/
  ## V = sqrt(chi2/(n*df))
  ## where df* = min(r – 1, c – 1) and r = the number of rows and c = the number of columns in the contingency table.
  V <- sqrt(chi2/n)
  V[is.na(V)] <- 0
  
  if(show.sign){
    ## Add sign based on log(odds ratio)
    observed_1 <- observed + 1
    ratio_case <- observed_1[,'case.true'] / observed_1[,'case.false']
    ratio_ctrl <- observed_1[,'ctrl.true'] / observed_1[,'ctrl.false']
    log_odds <- log2(ratio_case / ratio_ctrl)
    V[log_odds<0] <- -V[log_odds<0]
  }
  
  return(V)
}

get_stats <- function(df=samples_treatment_gene_gistic,drug,cancertype) {
  genes <- names(df)[7:length(names(df))]
  df_sample= setNames(data.frame(matrix(ncol = 12, nrow = 0)),c("drug","cancertype","gene","HMF_treat_2","HMF_treat_0","HMF_notreat_2","HMF_notreat_0","PCAWG_notreat_2","PCAWG_notreat_0","Fisher","Cramer_Est","Fisher_Est"))
  for(gene in genes){
    df_sub <- cbind(df[,1:6],dplyr::select(df,gene))
    df_sub <- df_sub %>% group_by_at(c("treated","cohort",sprintf("%s",gene))) %>% tally() %>% as.data.frame()
    is.integer0 <- function(x){is.integer(x) && length(x) == 0L}
  
    if(is.integer0(dplyr::filter(df_sub, cohort=="HMF" , treated == "Yes" ) %>% filter_(sprintf("`%s` == 2",gene))%>% dplyr::pull(n))){HMFtreat2=0}else{HMFtreat2=dplyr::filter(df_sub, cohort=="HMF" , treated == "Yes" ) %>% filter_(sprintf("`%s` == 2",gene))%>% dplyr::pull(n)}
    if(is.integer0(dplyr::filter(df_sub, cohort=="HMF" , treated == "Yes" ) %>% filter_(sprintf("`%s` == 0",gene))%>% dplyr::pull(n))){HMFtreat0=0}else{HMFtreat0=dplyr::filter(df_sub, cohort=="HMF" , treated == "Yes" ) %>% filter_(sprintf("`%s` == 0",gene))%>% dplyr::pull(n)}
    if(is.integer0(dplyr::filter(df_sub, cohort=="HMF" , treated == "No" ) %>% filter_(sprintf("`%s` == 2",gene))%>% dplyr::pull(n))){HMFnotreat2=0}else{HMFnotreat2=dplyr::filter(df_sub, cohort=="HMF" , treated == "No" ) %>% filter_(sprintf("`%s` == 2",gene))%>% dplyr::pull(n)}
    if(is.integer0(dplyr::filter(df_sub, cohort=="HMF" , treated == "No" ) %>% filter_(sprintf("`%s` == 0",gene))%>% dplyr::pull(n))){HMFnotreat0=0}else{HMFnotreat0=dplyr::filter(df_sub, cohort=="HMF" , treated == "No" ) %>% filter_(sprintf("`%s` == 0",gene))%>% dplyr::pull(n)}
    if(is.integer0(dplyr::filter(df_sub, cohort=="PCAWG" , treated == "No" ) %>% filter_(sprintf("`%s` == 2",gene))%>% dplyr::pull(n))){PCAWGnotreat2=0}else{PCAWGnotreat2=dplyr::filter(df_sub, cohort=="PCAWG" , treated == "No" ) %>% filter_(sprintf("`%s` == 2",gene))%>% dplyr::pull(n)}
    if(is.integer0(dplyr::filter(df_sub, cohort=="PCAWG" , treated == "No" ) %>% filter_(sprintf("`%s` == 0",gene))%>% dplyr::pull(n))){PCAWGnotreat0=0}else{PCAWGnotreat0=dplyr::filter(df_sub, cohort=="PCAWG" , treated == "No" ) %>% filter_(sprintf("`%s` == 0",gene))%>% dplyr::pull(n)}
    
    stats <- data.frame(drug=drug,
                        cancertype=cancertype,
                        gene = gene, 
               HMF_treat_2 =     HMFtreat2,
               HMF_treat_0 =     HMFtreat0,
               HMF_notreat_2 =   HMFnotreat2,
               HMF_notreat_0 =   HMFnotreat0,
               PCAWG_notreat_2 = PCAWGnotreat2,
               PCAWG_notreat_0 = PCAWGnotreat0)
    if(stats$PCAWG_notreat_2==0){
      stats$PCAWG_notreat_2 <- 1
    }
    stats$Fisher <- fisher.test(
      rbind(
        c(stats$HMF_treat_2,stats$HMF_treat_0),
        c(stats$PCAWG_notreat_2,stats$PCAWG_notreat_0)
      ))$p.value
    stats$Cramer_Est <- cramerV.default(stats$HMF_treat_2,stats$HMF_treat_0,stats$PCAWG_notreat_2,stats$PCAWG_notreat_0)
    stats$Fisher_Est <- fisher.test(
      rbind(
        c(stats$HMF_treat_2,stats$HMF_treat_0),
        c(stats$PCAWG_notreat_2,stats$PCAWG_notreat_0)
      ))$estimate
    df_sample <- rbind(df_sample,stats)
  }
  return(df_sample)
  df_sample <- NULL
}

rank_genes<- function(df=gr_genes_scores,df2=add_ampl_context){
  df_sample= setNames(data.frame(matrix(ncol = 3, nrow = 0)),c("gene","dist","rank"))
  for(peak in unique(sort(df$peakID))){
    
    dfsub <- subset(df,peakID==peak)
    dfsub2 <- subset(df2,Unique.Name==peak)
    
    gr_dfsub <- with(dfsub, GRanges(seqnames, IRanges(start, end)))
    seqlevelsStyle(gr_dfsub) = "UCSC"
    gr_dfsub <- gr_dfsub[ seqnames(gr_dfsub)  %in% chromosomes,]
    
    gr_dfsub2 <- with(dfsub2, GRanges(Chromosome, IRanges(TOPpeak_start, TOPpeak_end), peakID = Unique.Name,peakStart=TOPpeak_start , peakEnd=TOPpeak_end))
    seqlevelsStyle(gr_dfsub2) = "UCSC"
    gr_dfsub2 <- gr_dfsub2[ seqnames(gr_dfsub2)  %in% chromosomes,]
    
    distance(gr_dfsub,gr_dfsub2)
    dfsub$dist <- distance(gr_dfsub,gr_dfsub2)
    order.scores<-order(dfsub$dist)
    dfsub$rank <- NA
    dfsub$rank[order.scores] <- 1:nrow(dfsub)
    dfsub <- dfsub %>% dplyr::select(gene,dist,rank)
    df_sample <- rbind(df_sample,dfsub)
  }
  return(df_sample)
}

process_gistic <- function(GISTIC_repo="/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/LUAD/GISTIC_LUAD_Anti_EGFR/results/",outdir,hpc){
  outdir <- GISTIC_repo
  print(GISTIC_repo)
  print(outdir)
  print(hpc)
  
  ## Checks --------------------------------

  if(!dir.exists(GISTIC_repo)){
    stop("provide correct GISTIC_repo path")
  }else{
    #locate all purple and linx files
    GISTIC_files <- list.files(GISTIC_repo, recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
  }
  
  if(hpc=="Yes"){
    google_clinical <- as.data.frame(readRDS(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/clinical_overview.rds"))
    hg19 <- readRDS(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/CN_gene_database/hg19.rds")
    genes <- data.frame(read.csv("/hpc/cuppen/shared_resources/HMF_data/resources/External_Resources/HMFTools-Resources/ensembl_db/ensembl_gene_data.csv",sep = "," ,na=""))
    CN_genes_df = as.data.frame(t(as.data.frame(readRDS(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/CN_gene_database/gene_CN_strict25_final.rds"))))
  }else{
    google_clinical <- as.data.frame(readRDS(file = "/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/clinical_overview.rds"))
    hg19 <- readRDS(file = "/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/hg19.rds")
    genes <- data.frame(read.csv("/Users/avanhoeck/hpc/cuppen/shared_resources/HMF_data/resources/External_Resources/HMFTools-Resources/ensembl_db/ensembl_gene_data.csv",sep = "," ,na=""))
    CN_genes_df = as.data.frame(t(as.data.frame(readRDS(file = "/Users/avanhoeck/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/CN_gene_database/gene_CN_strict25_final.rds"))))
  }
  
  genes <- genes %>% dplyr::select(GeneId,GeneName,Chromosome,GeneStart,GeneEnd)
  genes <- genes[!duplicated(genes$GeneName),]
  print(length(genes$GeneId))
  gr_genes <- with(genes, GRanges(Chromosome, IRanges(GeneStart, GeneEnd), gene = GeneName, genelenght = GeneEnd-GeneStart, geneStart = GeneStart, geneEnd = GeneEnd))
  seqlevels(gr_genes) <- paste("chr", seqlevels(gr_genes), sep="")
  seqlevelsStyle(gr_genes) = "UCSC"
  gr_genes <- gr_genes[ seqnames(gr_genes) %in% chromosomes,]
  
  groupfolder = basename(GISTIC_repo) %>% gsub(pattern = "\\..*$",replacement =  "")
  stringsplitlength <- length(strsplit(groupfolder, '_')[[1]])
  
  
  if(stringsplitlength==3){
      drug <- strsplit(groupfolder, '_')[[1]][3]
      cancertype <- strsplit(groupfolder, '_')[[1]][2]
      } else if(stringsplitlength==4){
        drug <- sprintf("%s_%s",strsplit(groupfolder, '_')[[1]][3],strsplit(groupfolder, '_')[[1]][4])
        cancertype <- strsplit(groupfolder, '_')[[1]][2]
        }else if(stringsplitlength==5){
          drug <- sprintf("%s_%s_%s",strsplit(groupfolder, '_')[[1]][3],strsplit(groupfolder, '_')[[1]][4],strsplit(groupfolder, '_')[[1]][5])
          cancertype <- strsplit(groupfolder, '_')[[1]][2]
          }else if(stringsplitlength==6){
            drug <- sprintf("%s_%s_%s_%s",strsplit(groupfolder, '_')[[1]][3],strsplit(groupfolder, '_')[[1]][4],strsplit(groupfolder, '_')[[1]][5],strsplit(groupfolder, '_')[[1]][6])
            cancertype <- strsplit(groupfolder, '_')[[1]][2]
          }else if(stringsplitlength==7){
            drug <- sprintf("%s_%s_%s_%s_%s",strsplit(groupfolder, '_')[[1]][3],strsplit(groupfolder, '_')[[1]][4],strsplit(groupfolder, '_')[[1]][5],strsplit(groupfolder, '_')[[1]][6],strsplit(groupfolder, '_')[[1]][7])
            cancertype <- strsplit(groupfolder, '_')[[1]][2]
          }
  
  dir.create(sprintf("%s/results_percancer_strict25/",outdir))
  tempoutdir <- sprintf("%s/results_percancer_strict25/",outdir)
  pretreated_samples_hmf <- fromJSON(file = sprintf("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/json_files/%s/%s.json",cancertype,drug))
  nontreated_samples <- fromJSON(file = sprintf("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/json_files/%s/untreated.json",cancertype))
    
  GISTICscores <- GISTIC_files[grepl("scores.gistic",GISTIC_files)]
  GISTICscoresdf <- as.data.frame(read.table(GISTICscores,header = T,sep = "\t",colClasses=c("character","numeric","numeric","numeric","numeric","numeric","numeric")))
  colnames(GISTICscoresdf)[5] <- "log10.q.value"
  GISTICscoresdf[GISTICscoresdf$Chromosome==23,]$Chromosome <- "X"
  GISTICscoresdf_amp <- GISTICscoresdf %>% dplyr::filter(Type=="Amp")
  gr_GISTICscoresdf_amp <- with(GISTICscoresdf_amp, GRanges(Chromosome, IRanges(Start, End), log10 = log10.q.value,scoreStart=Start , scoreEnd=End))
  seqlevels(gr_GISTICscoresdf_amp) <- paste("chr", seqlevels(gr_GISTICscoresdf_amp), sep="")
  seqlevelsStyle(gr_GISTICscoresdf_amp) = "UCSC"
  gr_GISTICscoresdf_amp <- gr_GISTICscoresdf_amp[ seqnames(gr_GISTICscoresdf_amp)  %in% chromosomes,]
  
  GISTICscoresdf_del <- GISTICscoresdf %>% dplyr::filter(Type=="Del")
  gr_GISTICscoresdf_del <- with(GISTICscoresdf_del, GRanges(Chromosome, IRanges(Start, End), log10 = log10.q.value,scoreStart=Start , scoreEnd=End))
  seqlevels(gr_GISTICscoresdf_del) <- paste("chr", seqlevels(gr_GISTICscoresdf_del), sep="")
  seqlevelsStyle(gr_GISTICscoresdf_del) = "UCSC"
  gr_GISTICscoresdf_del <- gr_GISTICscoresdf_del[ seqnames(gr_GISTICscoresdf_del)  %in% chromosomes,]
  
  
  
  lesionfiles <- GISTIC_files[grepl("all_lesions.conf_95",GISTIC_files)]
  all_lesions.conf_95 <- read.table(lesionfiles[which(nchar(lesionfiles)==min(nchar(lesionfiles)))],header = T,sep = "\t")
  
  #get the peaks
  all_peaks <- all_lesions.conf_95 %>% dplyr::select(Unique.Name,Wide.Peak.Limits) %>% dplyr::filter(!str_detect(Unique.Name, "CN values") )
  amp_peaks <- all_peaks %>% dplyr::select(Unique.Name,Wide.Peak.Limits) %>% dplyr::filter(str_detect(Unique.Name, "Amplification") )
  amp_peaks <- amp_peaks %>% tidyr::separate(Wide.Peak.Limits,sep = "\\(",c("Wide.Peak.Limits",NA))
  amp_peaks <- amp_peaks %>% tidyr::separate(Wide.Peak.Limits,sep = "\\:",c("Chromosome","region"))
  amp_peaks <- amp_peaks %>% tidyr::separate(region,sep = "\\-",c("start","end"))
  amp_peaks$Unique.Name <- gsub(amp_peaks$Unique.Name,pattern = " ",replacement =  ".")
  amp_peaks$start <- as.numeric(amp_peaks$start)
  amp_peaks$end <- as.numeric(amp_peaks$end)
  amp_peaks$peaklenght = amp_peaks$end-amp_peaks$start
  amp_peaks <- amp_peaks  %>% dplyr::filter(peaklenght<1000000)
  amp_peaks$peaklenght <- NULL
  
  del_peaks <- all_peaks %>% dplyr::select(Unique.Name,Wide.Peak.Limits) %>% dplyr::filter(str_detect(Unique.Name, "Deletion") )
  del_peaks <- del_peaks %>% tidyr::separate(Wide.Peak.Limits,sep = "\\(",c("Wide.Peak.Limits",NA))
  del_peaks <- del_peaks %>% tidyr::separate(Wide.Peak.Limits,sep = "\\:",c("Chromosome","region"))
  del_peaks <- del_peaks %>% tidyr::separate(region,sep = "\\-",c("start","end"))
  del_peaks$Unique.Name <- gsub(del_peaks$Unique.Name,pattern = " ",replacement =  ".")
  del_peaks$start <- as.numeric(del_peaks$start)
  del_peaks$end <- as.numeric(del_peaks$end)
  del_peaks$peaklenght = del_peaks$end-del_peaks$start
  del_peaks <- del_peaks  %>% dplyr::filter(peaklenght<1000000)
  del_peaks$peaklenght <- NULL
  

  pretreated_samples_hmf_df= data.frame(sampleId =pretreated_samples_hmf, 
                                        ttype = cancertype, 
                                        drugtype =drug,
                                        cohort = "HMF",
                                        treated = "Yes") 
  
  

  nontreated_samples_pcawg_df= data.frame(sampleId = nontreated_samples, 
                                          ttype = cancertype, 
                                          drugtype =drug,
                                          cohort = "PCAWG",
                                          treated = "No")
  
  samples_treatment <- rbind(pretreated_samples_hmf_df,nontreated_samples_pcawg_df)
  
  sample_overview = data.frame(ttype = cancertype, 
                               drugtype =drug,
                               type="amplification",
                               pretreated_samples_hmf = nrow(pretreated_samples_hmf_df),
                               nontreated_samples_pcawg = if(is.na(nontreated_samples_pcawg_df$sampleId)){nontreated_samples_pcawg = 0}else {nontreated_samples_pcawg = nrow(nontreated_samples_pcawg_df)},
                               number_amp_wide_peaks = nrow(amp_peaks),
                               number_del_wide_peaks = nrow(del_peaks))
  
  write.table(sample_overview, file=sprintf("%s/sample_overview_%s.txt",tempoutdir,groupfolder),sep = "\t", col.names = TRUE,qmethod = "double", quote = FALSE,row.names = FALSE)
  
  
  if(nrow(amp_peaks)>0){
    gr_amp_peaks <- with(amp_peaks, GRanges(Chromosome, IRanges(start, end), peakID = Unique.Name,peakStart=start , peakEnd=end))
    seqlevelsStyle(gr_amp_peaks) = "UCSC"
    gr_amp_peaks <- gr_amp_peaks[ seqnames(gr_amp_peaks)  %in% chromosomes,]
    olaps <- GenomicRanges::findOverlaps(gr_GISTICscoresdf_amp, gr_amp_peaks, ignore.strand=TRUE)
    qh_amp <- S4Vectors::queryHits(olaps)
    sh_amp <- S4Vectors::subjectHits(olaps)
    GISTICscores_peak_amp <- gr_GISTICscoresdf_amp[qh_amp]
    GISTICscores_peak_amp$peakID <- gr_amp_peaks$peakID[sh_amp]
    GISTICscores_peak_amp <- as.data.frame(GISTICscores_peak_amp)
    GISTICscores_peak_amp$Cols <- "gain"
  

  #process amplification
  
    olaps <- GenomicRanges::findOverlaps(gr_genes,gr_amp_peaks, ignore.strand=TRUE)
    qh <- S4Vectors::queryHits(olaps)
    sh <- S4Vectors::subjectHits(olaps)
    
    gr_genes_scores <- gr_genes[qh]
    gr_genes_scores$peakID <- gr_amp_peaks$peakID[sh]
    gr_genes_scores$peakStart <- gr_amp_peaks$peakStart[sh]
    gr_genes_scores$peakEnd <- gr_amp_peaks$peakEnd[sh]
    gr_genes_scores <- as.data.frame(gr_genes_scores)
  
    add_context <- all_lesions.conf_95 %>% dplyr::select(Unique.Name,Peak.Limits,q.values) %>% dplyr::filter(!str_detect(Unique.Name, "CN values") )
    add_ampl_context <- add_context %>% dplyr::select(Unique.Name,Peak.Limits,q.values) %>% dplyr::filter(str_detect(Unique.Name, "Amplification") )
    add_ampl_context <- add_ampl_context %>% tidyr::separate(Peak.Limits,sep = "\\(",c("Peak.Limits",NA))
    add_ampl_context <- add_ampl_context %>% tidyr::separate(Peak.Limits,sep = "\\:",c("Chromosome","region"))
    add_ampl_context <- add_ampl_context %>% tidyr::separate(region,sep = "\\-",c("TOPpeak_start","TOPpeak_end"))
    add_ampl_context$Unique.Name <- gsub(add_ampl_context$Unique.Name,pattern = " ",replacement =  ".")
    add_ampl_context$TOPpeak_start <- as.numeric(add_ampl_context$TOPpeak_start)
    add_ampl_context$TOPpeak_end <- as.numeric(add_ampl_context$TOPpeak_end)
    ranked_genes <- rank_genes(df=gr_genes_scores,df2=add_ampl_context)
    add_ampl_context$Chromosome <- NULL
    add_ampl_context <- add_ampl_context %>% dplyr::rename(peakID=Unique.Name)
    gr_genes_scores <- dplyr::left_join(gr_genes_scores,add_ampl_context,by="peakID") 
    gr_genes_scores <- dplyr::left_join(gr_genes_scores,ranked_genes,by="gene") 
    gr_genes_scores <- unique(gr_genes_scores)
    
    CN_genes_df_amp = CN_genes_df %>% tibble::rownames_to_column(var = "sample") %>%  tidyr::separate(sample, c("sampleId", "eventype"))
    CN_genes_df_amp = CN_genes_df_amp %>% dplyr::filter(eventype=="AMP")
    
    samples_treatment_gene <- left_join(samples_treatment,CN_genes_df_amp,by="sampleId")
    samples_treatment_gene <- samples_treatment_gene[complete.cases(samples_treatment_gene), ]
    rm(CN_genes_df_amp)
    
    samples_treatment_gene_gistic = cbind(samples_treatment_gene[,1:6],samples_treatment_gene[,names(samples_treatment_gene) %in% gr_genes_scores$gene])
    samples_treatment_gene_gistic <- samples_treatment_gene_gistic[complete.cases(samples_treatment_gene_gistic), ]
    
    stat_genes <- get_stats(df=samples_treatment_gene_gistic,drug,cancertype)
    gr_genes_scores <- left_join(gr_genes_scores,stat_genes,by="gene")
    gr_genes_scores <- gr_genes_scores %>% arrange(Fisher) %>% dplyr::select(drug,cancertype,gene,seqnames,peakID,dist,HMF_treat_2,HMF_treat_0,PCAWG_notreat_2,PCAWG_notreat_0,q.values,
                                                                             Fisher,Cramer_Est,Fisher_Est,dist,rank)
    
    
    colnames(amp_peaks)[1] <- "peakID"
    gr_genes_scores <- left_join(gr_genes_scores,amp_peaks,by="peakID")
    gr_genes_scores <- gr_genes_scores %>% arrange(Fisher) %>% dplyr::select(drug,cancertype,gene,seqnames,peakID,start,end,dist,HMF_treat_2,HMF_treat_0,PCAWG_notreat_2,PCAWG_notreat_0,q.values,
                                                                             Fisher,Cramer_Est,Fisher_Est,dist,rank)

    
    
    gz1 <- gzfile(sprintf("%s/GISTIC_genes_results_amp_%s.csv.gz",tempoutdir,groupfolder), "w")
    write.csv(gr_genes_scores, gz1,quote = FALSE,sep="\t",col.names = TRUE,row.names = FALSE,qmethod = "double")
    close(gz1)
  }
  
  
  
  
  #process deletions
  if(nrow(del_peaks)>0){
    gr_del_peaks <- with(del_peaks, GRanges(Chromosome, IRanges(start, end), peakID = Unique.Name,peakStart=start , peakEnd=end))
    seqlevelsStyle(gr_del_peaks) = "UCSC"
    gr_del_peaks <- gr_del_peaks[ seqnames(gr_del_peaks)  %in% chromosomes,]
    olaps <- GenomicRanges::findOverlaps(gr_GISTICscoresdf_del, gr_del_peaks, ignore.strand=TRUE)
    qh_del <- S4Vectors::queryHits(olaps)
    sh_del <- S4Vectors::subjectHits(olaps)
    GISTICscores_peak_del <- gr_GISTICscoresdf_del[qh_del]
    GISTICscores_peak_del$peakID <- gr_del_peaks$peakID[sh_del]
    GISTICscores_peak_del <- as.data.frame(GISTICscores_peak_del)
    GISTICscores_peak_del$Cols <- "loss"
    
    
    
    samples_treatment_gene_gistic <- NULL
    gr_genes_scores <- NULL
    
    olaps <- GenomicRanges::findOverlaps(gr_genes,gr_del_peaks, ignore.strand=TRUE)
    qh <- S4Vectors::queryHits(olaps)
    sh <- S4Vectors::subjectHits(olaps)
    
    gr_genes_scores <- gr_genes[qh]
    gr_genes_scores$peakID <- gr_del_peaks$peakID[sh]
    gr_genes_scores$peakStart <- gr_del_peaks$peakStart[sh]
    gr_genes_scores$peakEnd <- gr_del_peaks$peakEnd[sh]
    gr_genes_scores <- as.data.frame(gr_genes_scores)
    
    add_context <- all_lesions.conf_95 %>% dplyr::select(Unique.Name,Peak.Limits,q.values) %>% dplyr::filter(!str_detect(Unique.Name, "CN values") )
    add_del_context <- add_context %>% dplyr::select(Unique.Name,Peak.Limits,q.values) %>% dplyr::filter(str_detect(Unique.Name, "Deletion") )
    add_del_context <- add_del_context %>% tidyr::separate(Peak.Limits,sep = "\\(",c("Peak.Limits",NA))
    add_del_context <- add_del_context %>% tidyr::separate(Peak.Limits,sep = "\\:",c("Chromosome","region"))
    add_del_context <- add_del_context %>% tidyr::separate(region,sep = "\\-",c("TOPpeak_start","TOPpeak_end"))
    add_del_context$Unique.Name <- gsub(add_del_context$Unique.Name,pattern = " ",replacement =  ".")
    add_del_context$TOPpeak_start <- as.numeric(add_del_context$TOPpeak_start)
    add_del_context$TOPpeak_end <- as.numeric(add_del_context$TOPpeak_end)
    ranked_genes <- rank_genes(df=gr_genes_scores,df2=add_del_context)
    add_del_context$Chromosome <- NULL
    add_del_context <- add_del_context %>% dplyr::rename(peakID=Unique.Name)
    gr_genes_scores <- dplyr::left_join(gr_genes_scores,add_del_context,by="peakID") 
    gr_genes_scores <- dplyr::left_join(gr_genes_scores,ranked_genes,by="gene") 
    gr_genes_scores <- unique(gr_genes_scores)
    
    CN_genes_df_del = CN_genes_df %>% tibble::rownames_to_column(var = "sample") %>%  tidyr::separate(sample, c("sampleId", "eventype"))
    CN_genes_df_del = CN_genes_df_del %>% dplyr::filter(eventype=="DEL")
    
    
    samples_treatment_gene <- left_join(samples_treatment,CN_genes_df_del,by="sampleId")
    samples_treatment_gene <- samples_treatment_gene[complete.cases(samples_treatment_gene), ]
    rm(CN_genes_df_del)
    
    samples_treatment_gene_gistic = cbind(samples_treatment_gene[,1:6],samples_treatment_gene[,names(samples_treatment_gene) %in% gr_genes_scores$gene])
    samples_treatment_gene_gistic <- samples_treatment_gene_gistic[complete.cases(samples_treatment_gene_gistic), ]
  
    stat_genes <- get_stats(df=samples_treatment_gene_gistic,drug,cancertype)
    gr_genes_scores <- left_join(gr_genes_scores,stat_genes,by="gene")
    gr_genes_scores <-  gr_genes_scores %>% arrange(Fisher) %>% dplyr::select(drug,cancertype,gene,seqnames,peakID,dist,HMF_treat_2,HMF_treat_0,PCAWG_notreat_2,PCAWG_notreat_0,q.values,
                                                                          Fisher,Cramer_Est,Fisher_Est,dist,rank)
    
    
    colnames(del_peaks)[1] <- "peakID"
    gr_genes_scores <- left_join(gr_genes_scores,del_peaks,by="peakID")
    gr_genes_scores <- gr_genes_scores %>% arrange(Fisher) %>% dplyr::select(drug,cancertype,gene,seqnames,peakID,start,end,dist,HMF_treat_2,HMF_treat_0,PCAWG_notreat_2,PCAWG_notreat_0,q.values,
                                                                             Fisher,Cramer_Est,Fisher_Est,dist,rank)

    gz1 <- gzfile(sprintf("%s/GISTIC_genes_results_del_%s.csv.gz",tempoutdir,groupfolder), "w")
    write.csv(gr_genes_scores, gz1,quote = FALSE,sep="\t",col.names = TRUE,row.names = FALSE,qmethod = "double")
    close(gz1)
    }
  
  
  
  
  df <- rbind(amp_peaks,del_peaks)
  
  if(nrow(df)==0){
    df <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("ttype","drugtype","Unique.Name","Chromosome","start","end"))
    write.table(as.data.frame(df), file=sprintf("%s/wide_peak_overview_%s.txt",tempoutdir,groupfolder),sep = "\t", col.names = TRUE,qmethod = "double", quote = FALSE,row.names = FALSE)
  }else{
    df <- rbind(amp_peaks,del_peaks)
    df <- cbind(drugtype = drug, df)
    df <- cbind(ttype = cancertype, df)
    write.table(df, file=sprintf("%s/wide_peak_overview_%s.txt",tempoutdir,groupfolder),sep = "\t", col.names = TRUE,qmethod = "double", quote = FALSE,row.names = FALSE)
  }
  
  
  print("process done")
  
}


#####

args <- commandArgs(trailingOnly = TRUE)

process_gistic(GISTIC_repo=args[1], outdir=args[2],hpc=args[3])
