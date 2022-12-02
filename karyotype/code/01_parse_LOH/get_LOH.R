############### Get LOH ############### 
# author: arne

### Description
# This script parses the Loss of Heterozygosity (LOH) proportion from the .purple.cnv.somatic.tsv files (PURPLE output).

### Input
# individual .purple.cnv.somatic.tsv files

### Output
# LOH proportion per sample table

# global options
options(stringsAsFactors = FALSE, scipen = 999)

# libs
here::set_here(path='../../')
source(paste0(here::here(), '/code/r_objects/libs.R'))

#========= Path prefixes =========#
## Hartwig
# dirs<-Sys.glob("/path/to/Hartwig/data_request/somatics/*/purple/")
## PCAWG
# dirs<-Sys.glob("/path/to/PCAWG/samples")

# list all files
files=list.files(dirs, pattern = ".purple.cnv.somatic.tsv$",full.names=TRUE, recursive = TRUE)

# custom functions
'%!in%' <- function(x,y)!('%in%'(x,y))
doubt_df <- data.frame(file = character(),chromArm = character(), highest= integer(), second= integer(), tp = character(),
                       stringsAsFactors = FALSE)

readFile <- function(filename){
  first_line  <- readLines(filename, n=1)
  sample_data <- read.table(filename, skip = 1, header=FALSE, sep= "\t")
  names(sample_data) <- as.vector(unlist(strsplit(first_line[length(first_line)],"\t")))
  return(sample_data)
}

addMajorMinorAllele <- function(sample_data){
  sample_data$copyNumber[sample_data$copyNumber < 0] <- 0
  sample_data["unboundMinorAllele"] <- (1 - sample_data$baf) * sample_data$copyNumber
  sample_data$unboundMinorAllele[sample_data$unboundMinorAllele < 0.25 ] <- 0
  sample_data["majorAllele"] <- sample_data$copyNumber - sample_data$unboundMinorAllele
  sample_data$majorAllele[sample_data$majorAllele < 0 ] <- 0
  return(sample_data)
}

chromArmPloidy <- function(sample_data, sample_name){
  majorPloidy_list <- c()
  minorPloidy_list <- c()
  LOH_list <- c()
  
  lapply(c(1:22,"X"), function(c){
    if(c %!in% c(13,14,15,21,22)){
      subset_data <- subset(sample_data, sample_data$`chromosome` == c)
      loc_centromeer <- which(subset_data$segmentEndSupport == "CENTROMERE")
      # Calculates everything for P arm
      subset_p  <- subset_data[1:loc_centromeer,]
      df_subset <- ploidyWeight(subset_p)

      ploidyFrameP <- createPloidyFrame(df_subset, c, "p")
      return_listP <- aneuPloidy(ploidyFrameP, c, "p", sample_name)
      tail(ploidyFrameP)
      LOHarmP <- ploidyFrameP %>% dplyr::filter(lohStatus == TRUE) %>% mutate(segment_size = as.integer(end) - as.integer(start)) %>% summarise(Total = sum(segment_size))  %>% dplyr::pull(Total)

      majorPloidy_list[length(majorPloidy_list)+1] <<- return_listP[1]
      minorPloidy_list[length(minorPloidy_list)+1] <<- return_listP[2]
      LOH_list[length(LOH_list)+1] <<- as.list(LOHarmP)
      
      # Calculates everything for Q arm
      subset_q   <- subset_data[(loc_centromeer+1):length(subset_data$start),]
      df_subset2 <<- ploidyWeight(subset_q)

      ploidyFrameQ <<- createPloidyFrame(df_subset2,c, "q")
      return_listQ <- aneuPloidy(ploidyFrameQ, c, "q", sample_name)
      LOHarmQ <- ploidyFrameQ %>% dplyr::filter(lohStatus == TRUE) %>% mutate(segment_size = as.integer(end) - as.integer(start)) %>% summarise(Total = sum(segment_size))  %>% dplyr::pull(Total)
      

      majorPloidy_list[length(majorPloidy_list)+1] <<- return_listQ[1]
      minorPloidy_list[length(minorPloidy_list)+1] <<- return_listQ[2]
      LOH_list[length(LOH_list)+1] <<- as.list(LOHarmQ)

    }else{
      
      # Calculates everything for the autosomes that are treated as one armed
      subset_data <- subset(sample_data, sample_data$`chromosome` == c)
      df_subset   <- ploidyWeight(subset_data)
      ploidyFrame <- createPloidyFrame(df_subset, c, "q")
      return_list <- aneuPloidy(ploidyFrame, c, "q", sample_name)
      LOHarm <- ploidyFrame %>% dplyr::filter(lohStatus == TRUE) %>% mutate(segment_size = as.integer(end) - as.integer(start)) %>% summarise(Total = sum(segment_size))  %>% dplyr::pull(Total)
      
      majorPloidy_list[length(majorPloidy_list)+1] <<- return_list[1]
      minorPloidy_list[length(minorPloidy_list)+1] <<- return_list[2]
      LOH_list[length(LOH_list)+1] <<- as.list(LOHarm)
      
    }
  })
  return(list(majorPloidy_list, minorPloidy_list, LOH_list))
}

ploidyWeight <- function(data_frame){
  data_frame["majorGroup"] <- cut(data_frame$majorAllele,
                                  breaks= c(-1.5,-0.5,0.5,c(1*(1.5:20))), right = FALSE,
                                  include.lowest = TRUE, labels = c(-1:19))

  data_frame["group"] <- rleid(data_frame$majorGroup)
  data_frame$majorAllele[data_frame$majorAllele <= 0 ] <- NA
  data_frame["mapWeight"] <- rep(1, length(data_frame[,1]))
  data_frame$mapWeight[is.na(data_frame$mapWeight)] <- 1

  return(data_frame)
}

createPloidyFrame <- function(arm_subset, chromosome, arm){
  grouplist <- unique(arm_subset$group)

  tmp_list  <- lapply(1:length(unique(arm_subset$group)), function(i){
    group_set <- subset(arm_subset, arm_subset$group == grouplist[i])
    start <- group_set$start[1]
    end   <- group_set$end[length(group_set$end)]

    segment_size <- (group_set$end - group_set$start) + 1
    weight <-  segment_size * group_set$mapWeight

    averageMAP <- group_set$majorAllele * weight
    averageMAP <- sum(averageMAP) / sum(weight)

    averageMinor <- group_set$unboundMinorAllele * weight
    averageMinor <- sum(averageMinor) / sum(weight)

    averageCN <- group_set$copyNumber * weight
    averageCN <- sum(averageCN) / sum(weight)

    if(!is.na(averageMinor) && averageMinor <= 0.25 && !is.na(averageMAP) && averageMAP > 0.8){
      loh <- TRUE
    }else{
      loh <- FALSE
    }
    tmp_vector <- c(-1:19)
    return(c(paste0(chromosome,arm), as.integer(start), as.integer(end), as.numeric(averageMAP), as.numeric(averageMinor),
             as.numeric(averageCN), as.integer(tmp_vector[group_set$majorGroup[1]]), as.logical(loh)))
  })

  ploidyFrame <- as.data.frame(do.call(rbind, tmp_list), stringsAsFactors = FALSE)
  colnames(ploidyFrame) <- c("chromArm", "start", "end", "averageMajorAP", "averageMinorAP", "averageCN", "majorAPcategory", "lohStatus")
  ploidyFrame["minorAPcategory"] <- round(as.numeric(ploidyFrame$averageMinorAP))
  return(ploidyFrame)
}


aneuPloidy <- function(ploidyFrame, chromosome, arm, filename){
  ploidyFrame["segmentSize"] <- (as.numeric(ploidyFrame$end) - as.numeric(ploidyFrame$start))

  major_list <-  aggregate(segmentSize~majorAPcategory, data = ploidyFrame, FUN = sum)
  major_list["percentage"]   <- (major_list$segmentSize / sum(major_list$segmentSize)) * 100

  minor_list <-  aggregate(segmentSize~minorAPcategory, data = ploidyFrame, FUN = sum)
  minor_list["percentage"]   <- (minor_list$segmentSize / sum(minor_list$segmentSize)) * 100

  tmpList  <- list(major_list, minor_list)
  category <- c("major", "minor")

  counter <- 1
  ploidyList <- lapply(tmpList, function(subList){
    column <- paste0(category[counter], "APcategory")
    if(length(order(subList$percentage)) > 1){
      if(max(subList$percentage) >= 50  && (subList[which.max(subList$percentage), "percentage"] - sort(subList$percentage, decreasing = TRUE)[2]) > 10){
        ploidyLevel <- subList[which.max(subList$percentage), column]
      }else{
        ploidyLevel <- NA
        doubt_df[nrow(doubt_df)+1,] <<- c(filename, paste0(chromosome, arm),  subList[which.max(subList$percentage), column],
                                          subList[which(subList$percentage == (sort(subList$percentage,decreasing=TRUE)[2])), column], category[counter])
      }
      ploidyLevel <- subList[which.max(subList$percentage), "majorAPcategory" ] 
    }else{
      ploidyLevel <- subList[which.max(subList$percentage), column]
    }
    counter <<- counter + 1
    return(ploidyLevel)
  })

  return(list(as.numeric(unlist(ploidyList[1])), as.numeric(as.character(unlist(ploidyList[2])))))
}

getChromArms <- function(){
  chromArms = c()
  for(i in c(1:22,"X")){
    if(i %!in% c(13,14,15,21,22)){
      chromArms <- c(chromArms, c(paste0(i,"p"), paste0(i,"q")))
    }else{
      chromArms <- c(chromArms, c(paste0(i,"q")))
    }
  }
  chromArms <- append(chromArms, "genome_level")
  return(chromArms)
}

ploidyEstimator <- function(files){
  chromArms <- getChromArms()
  df_LOH  <- data.frame(chromArms[0:41])
  names(df_LOH) <- "chromArms"
  
  counter <<- 0
  return_list <- lapply(files, function(file){
    counter <<- counter + 1
    sample_data <- readFile(file)
    sample_name <- unlist(strsplit(file,"/"))[length(unlist(strsplit(files[1],"/")))]
    sample_name <- unlist(strsplit(sample_name,"\\."))[1]
    sample_data <- addMajorMinorAllele(sample_data)
    return_list <- chromArmPloidy(sample_data, sample_name)
    df_LOH[sample_name]  <<- round(as.numeric(unlist(return_list[3])))
    

  })
  write.table(df_LOH , file = "LOH.tsv" , sep="\t", col.names = TRUE, row.names = FALSE)
  
}

ploidyEstimator(files)