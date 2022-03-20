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

methylationFilesCG<-methRead(location = locs,indNames,assembly = "arabidopsisRef",pipeline = "bismarkCytosineReport",treatment = c(rep(1,times=length(indNames))),mincov = 10,dbdir = "filtCG")
methylationFilesCHG<-methRead(location = locs,indNames,assembly = "arabidopsisRef",pipeline = "bismarkCytosineReport",treatment = c(rep(1,times=length(indNames))),mincov = 10,dbdir = "filtCHG",context="CHG")
methylationFilesCHH<-methRead(location = locs,indNames,assembly = "arabidopsisRef",pipeline = "bismarkCytosineReport",treatment = c(rep(1,times=length(indNames))),mincov = 10,dbdir = "filtCHH",context="CHH")
nLociData<-data.frame()
for(i in 1:length(methylationFilesCG)){
  nSitesCG<-length(methylationFilesCG@.Data[[i]]$chr)
  nSitesCHG<-length(methylationFilesCHG@.Data[[i]]$chr)
  nSitesCHH<-length(methylationFilesCHH@.Data[[i]]$chr)
  nSites<-data.frame(sampleID=methylationFilesCHG@.Data[[i]]@sample.id,nSitesCG,nSitesCHG,nSitesCHH)
  nLociData<-rbind(nLociData,nSites)
}
nLociData$tot<-rowSums(nLociData[,2:4])
write.table(nLociData,"results/nLociPerInd.tsv")
