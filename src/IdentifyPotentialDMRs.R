rm(list=ls())
library(ggplot2)
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #Make pretty colours
annot<-read.csv("results/mergedAnnot.csv") #Load annotation data
annot$chr<-sub("chr","",annot$X) #Change the format of the contig names
siteDSS<-readRDS("results/siteDSS.RData") #Load the DSS output
siteDSS<-siteDSS[!is.na(siteDSS$fdrs),] #Remove NAs from DSS output
annot<-annot[annot$chr%in%siteDSS$chr,] #Remove non covered contigs from the analysis
siteDSS$qval<-p.adjust(siteDSS$pvals,method="fdr") #Redo the multiple testing correction to 

dataDMCs<-siteDSS[siteDSS$qval<0.1,]#Get all DMCs
dataDMCs<-dataDMCs[order(dataDMCs$chr,dataDMCs$pos),] #Order based on positions and chromosome
dataDMCs<-dataDMCs[duplicated(dataDMCs)==F,] #Remove duplications?????
write.table(annot,"results/coveredAnnotation.tsv",row.names=F)

#NOTE: THIS ONLY WORKS IF YOU HAVE A treated-nonTreated experimental setup
#otherwise you have to identify the DMRs for each specific treatment!

windowSize<-500 #This is the number of bases you look to (forwards)
dataOut<-data.frame()
i<-1

pb <- txtProgressBar(min = 0, max = nrow(dataDMCs), style = 3)
total<-0
  
for(i in 1:nrow(dataDMCs)){
  if(total==0){
    startRegionPos<-dataDMCs$pos[i]
    startRegionChr<-dataDMCs$chr[i]
    
  }  
  counts<-nrow(dataDMCs[which(dataDMCs$chr==dataDMCs$chr[i]&
                                dataDMCs$pos>dataDMCs$pos[i]&
                                dataDMCs$pos<(dataDMCs$pos[i]+windowSize)),])
  if(counts==0){
    dataOut<-rbind(dataOut,data.frame(chr=startRegionChr,startPos=startRegionPos,endPos=dataDMCs$pos[i],count=total+1))
    total<-0}else{
      total<-total+1
    }
  setTxtProgressBar(pb, i)
}


DMRs<-dataOut[dataOut$count>10,]


annotatedDMRs<-merge(DMRs,annot,by="chr")
write.table(annotatedDMRs,"results/DMRs.tsv")



