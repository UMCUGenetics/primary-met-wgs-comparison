#!/usr/bin/env Rscript
library(GenomicRanges)
library(VariantAnnotation)
#library(ggplot2)
library(reshape2)
library(dplyr)
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
chromosomes <- seqnames(get(ref_genome))[1:23]
options(stringsAsFactors = F)


#========= Query =========#
args = commandArgs(trailingOnly=TRUE)


sampleid <- args[1]
sprintf("START PROGRAM FOR %s", sampleid)
google_clinical <- readRDS(file = "/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/clinical_overview.rds")

#genes coordinates
genes <- data.frame(read.csv("/hpc/cuppen/shared_resources/HMF_data/resources/External_Resources/HMFTools-Resources/ensembl_db/ensembl_gene_data.csv",sep = "," ,na=""))
genes <- genes %>% dplyr::select(GeneId,GeneName,Chromosome,GeneStart,GeneEnd)
genes <- genes[!duplicated(genes$GeneName),]
gr_genes <- with(genes, GRanges(Chromosome, IRanges(GeneStart, GeneEnd), gene = GeneName, genelenght = GeneEnd-GeneStart, geneStart = GeneStart, geneEnd = GeneEnd))
seqlevels(gr_genes) <- paste("chr", seqlevels(gr_genes), sep="")
seqlevelsStyle(gr_genes) = "UCSC"
gr_genes <- gr_genes[ seqnames(gr_genes) %in% chromosomes,]

cohort <- google_clinical %>% dplyr::filter(sampleId==sampleid) %>% dplyr::pull(cohort)
if(cohort == "HMF"){
  setname <- google_clinical %>% dplyr::filter(sampleId==sampleid) %>% dplyr::pull(setName)
  CNV <- data.frame(read.csv(sprintf("/hpc/cuppen/shared_resources/HMF_data/DR-104-update5/somatics/%s/purple/%s.purple.cnv.somatic.tsv",sampleid,sampleid),sep = "\t"))
  purple <- data.frame(read.csv(sprintf("/hpc/cuppen/shared_resources/HMF_data/DR-104-update5/somatics/%s/purple/%s.purple.purity.tsv",sampleid,sampleid),sep = "\t"))
  genomeploidy <- purple %>% dplyr::pull(ploidy)
}
if(cohort == "PCAWG"){
  CNV <- data.frame(read.csv(sprintf("/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/%s-from-jar/purplesoft3.3/%sT.purple.cnv.somatic.tsv",sampleid,sampleid),sep = "\t"))
  purple <- data.frame(read.csv(sprintf("/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/%s-from-jar/purplesoft3.3/%sT.purple.purity.tsv",sampleid,sampleid),sep = "\t"))
  genomeploidy <- purple %>% dplyr::pull(ploidy)
}

gr_CNV <- with(CNV, GRanges(chromosome, IRanges(start, end), CN = copyNumber,CNStart=start , CNEnd=end,minorAlleleCN=minorAlleleCopyNumber,majorAlleleCN=majorAlleleCopyNumber))
seqlevels(gr_CNV) <- paste("chr", seqlevels(gr_CNV), sep="")
seqlevelsStyle(gr_CNV) = "UCSC"
gr_CNV <- gr_CNV[ seqnames(gr_CNV)  %in% chromosomes,]

olaps <- GenomicRanges::findOverlaps(gr_genes, gr_CNV, ignore.strand=TRUE)
qh <- S4Vectors::queryHits(olaps)
sh <- S4Vectors::subjectHits(olaps)

cgenes_withBPingenes <- gr_genes[qh]
cgenes_withBPingenes$CN <- gr_CNV$CN[sh]
cgenes_withBPingenes$CNStart <- gr_CNV$CNStart[sh]
cgenes_withBPingenes$CNEnd <- gr_CNV$CNEnd[sh]
cgenes_withBPingenes$minorAlleleCN <- gr_CNV$minorAlleleCN[sh]
cgenes_withBPingenes$majorAlleleCN <- gr_CNV$majorAlleleCN[sh]
cgenes_withBPingenes_df <- as.data.frame(cgenes_withBPingenes)

# start processing each gene
df_all <- NULL
df_all= setNames(data.frame(matrix(ncol = 9, nrow = 0)),c("seqnames","start","end","width","gene","genelenght","CN","minorAlleleCN","majorAlleleCN"))
for(geneID in unique(cgenes_withBPingenes_df$gene)){
  df <- cgenes_withBPingenes_df%>% dplyr::filter(gene==geneID)
  if(nrow(df)==1){
    df_temp <- df %>% dplyr::select(seqnames,start,end,width,gene,genelenght,CN,minorAlleleCN,majorAlleleCN)
    }
  if(nrow(df)>1){
    lowest_CN <- min(df$CN)
    df_temp <- df %>% dplyr::filter(CN==lowest_CN)%>% dplyr::select(seqnames,start,end,width,gene,genelenght,CN,minorAlleleCN,majorAlleleCN)
    df_temp <- df_temp[1,]
  }
  df_all=rbind(df_all,df_temp)
}

#assign GISTIC bins
df_all <- df_all %>% dplyr::mutate(GIST_AMP=ifelse(df_all$CN > sum(2.5,genomeploidy),2,0)) #sum(0.9,genomeploidy),2,0) 2.5*genomeploidy
df_all <- df_all %>% dplyr::mutate(GIST_DEL=ifelse(df_all$CN < 0.3,2,0)) #sum(genomeploidy+(-2))
names(df_all)[names(df_all) == "GIST_AMP"] <- sprintf("%s_AMP",sampleid)
names(df_all)[names(df_all) == "GIST_DEL"] <- sprintf("%s_DEL",sampleid)

df_all <- df_all %>% dplyr::select(gene,sprintf("%s_AMP",sampleid),sprintf("%s_DEL",sampleid))
rownames(df_all) <- df_all$gene
df_all$gene <- NULL
write.table(df_all, file = sprintf("/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/GISTIC_rebuttal/scripts/CN_gene_database/%s_CNgene_strict25.txt",sampleid),sep = "\t", col.names = NA,qmethod = "double", quote = FALSE)
