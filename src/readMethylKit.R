library(methylKit)
library(ggplot2)
library(reshape2)
library(ggpubr)
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #Make pretty colours
locs<-as.list(list.files("output/methylation_calling/",full.names = T)[grep("CX_report.txt$",list.files("output/methylation_calling/",full.names = T))]) #reads the data
indNamesVector<-sub("_trimmed.*$","",sub("^.*/","",locs))
indNames<-as.list(indNamesVector)


####### Run the first time (This loads all methylation files and normalizes)

methylationFiles<-methRead(location = locs,indNames,assembly = "arabidopsisRef",pipeline = "bismarkCytosineReport",treatment = c(rep(1,times=length(indNames))),mincov = 10,dbdir = "methylKitdb")
methylDB<-makeMethylDB(methylationFiles,dbdir = "methylKitdb")

normalized<-normalizeCoverage(methylDB,dbdir="methylKitdb")
 


#This calculates number of sites covered for each individual (with 10x coverage)
#Only runs on linux!!! 
#TODO Make this more elegant
nLociData<-data.frame(nloci=as.numeric(system("for file in methylKitdb/*normed.txt.bgz; do zcat $file | wc -l; done",TRUE)),file=system("for file in methylKitdb/*normed.txt.bgz; do echo $file; done",TRUE))

#Add individual names to the data (Individuals should be named following the pop_ind format)
for(i in 1:length(nLociData[,1])){
  nLociData$ind[i]<-sub("^.*/","",sub("_normed.txt.bgz","",nLociData$file[i]))
  nLociData$indNumber[i]<-sub("^.*_","",sub("_normed.txt.bgz","",nLociData$file[i]))
  nLociData$loc[i]<-grep(paste0(nLociData$indName[i],"_"),locs,value = T)
}
max(nLociData$nloci)

write.table(nLociData,"results/nLociPerInd.tsv")