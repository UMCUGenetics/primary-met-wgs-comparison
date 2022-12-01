library(dndscv)
args = commandArgs(trailingOnly=TRUE)
mut_file <- args[1]
output_path <- args[2]
treatment_id <- args[3]
muts <- read.csv(gzfile(mut_file), as.is = TRUE, sep="\t")

# gene specific 


# compute analysis

out <- dndscv(muts, outmats=T)
outd<-paste(output_path,"/",treatment_id,".dndscv.annotmuts.tsv.gz",sep="")
outr<-paste(output_path,"/",treatment_id,".dndscv.results.tsv.gz",sep="")
outci<-paste(output_path,"/",treatment_id,".dndscv.ci.tsv.gz",sep="")
ci = geneci(out)

# save results

write.table(out$annotmuts,file=gzfile(outd),sep="\t",row.names=FALSE, quote = FALSE)
write.table(out$sel_cv,file=gzfile(outr),sep="\t",row.names=FALSE, quote = FALSE)
write.table(ci,file=gzfile(outci),sep="\t",row.names=FALSE, quote = FALSE)




 

