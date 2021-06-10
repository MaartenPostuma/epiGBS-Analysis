library(methylKit)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(DSS)
library(bsseq)
library(GenomicRanges)

colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #Make pretty colours
#################

Design<-read.table("results/Design_filtered.tsv") #Read design table


################
#Load filtered methylKit data to import into DSS (Using the same filters)
positions<-data.frame() 
for (context in c("CG","CHG","CHH")){
  i<-i+1
  united<-readMethylDB(list.files(paste0("Filt",context),
                                  pattern = paste0("methylBase_",context,".txt.bgz$")
                                  ,full.names=T)
  )
  positions<-rbind(positions,data.frame(getData(united),context=context))
}


#Make a dataframe with order of methylKit samples (Methylkit orders differently i.e. sample_1 goes after sample_10 while R tends to do it the otherway around)
methylKitOrder<-data.frame(sample=getSampleID(united),num=1:length(getSampleID(united)))

#Make a dataframe that has the context per cytosine
contexts<-data.frame(chr=positions$chr,start=positions$start,context=positions$context,loci=paste(positions$chr,positions$start,sep="_"))

#####################
#Load epiGBS data into epiGBS based on the filtering performed in methylKit

sampleList<-list()#Make a list which will contain all samples
i<-1
for(sample in Design$loc){ #Load all bismark data
  sampleID<-Design$indName[Design$loc==sample]#Get the sample ID from the design file
  covColumn<-grep(paste0("coverage",methylKitOrder$num[methylKitOrder$sample==sampleID],"$"),
                  colnames(positions)) #get the coverage column from methylKit
  
  
  #Create GenomicRanges object containing the loci that passed the filters we employed in methylKit
    #I.e. only sites that have min cov > 10 (is not NA in the dataset) and occur in 80% of individuals
  lociFilt<-GRanges(positions$chr[is.na(positions[,covColumn])==F],
                    ranges=IRanges(start=positions$start[is.na(positions[,covColumn])==F],
                                   end=positions$end[is.na(positions[,covColumn])==F]))
  
  sampleList[[i]]<-read.bismark(sample,loci=lociFilt) #Read them using bsseq
  i<-i+1  
}

myDSS<-combineList(sampleList) #Combine them into a BSseq dataset
##########################################


#Fit the model you want to fit!

myFit <- DMLfit.multiFactor(myDSS, Design,~site+treatment.y)
test<-DMLtest.multiFactor(myFit,term="site")
test$loci<-paste(test$chr,test$pos,sep="_")

testContext<-merge(test,contexts,by="loci") #Merge the test results with the contexts data.frame 

#plot pvalue distribution for different contexts
ggplot(testContext,aes(x=pvals,fill=context))+
  geom_histogram(binwidth = 0.01)+
  facet_grid(context~.,scales="free")+
  scale_fill_manual(values = colorBlindBlack8)


