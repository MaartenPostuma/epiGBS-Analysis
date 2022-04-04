library(methylKit)
library(ggplot2)
library(reshape2)
library(ggpubr)
#########################################
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #Make pretty colours

#Takes output from the readMethylKit + a design table:
#The Design table should contain a column called sampleID with the names of the samples (same as the names used in epiGBS) and can contain all other treatments etc. These will be combined with nLoci data after filtering to create filtered design. Which you can use to determine what filtering steps leave enough samples per treatment group
#
######## Run the first time (This loads all methylation files) ###########


locs<-as.list(list.files("output/methylation_calling/",full.names = T)[grep("CX_report.txt$",list.files("output/methylation_calling/",full.names = T))]) #reads the data
indNamesVector<-sub("_trimmed.*$","",sub("^.*/","",locs))
indNames<-as.list(indNamesVector)

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

##############################################
nLociData<-read.table("results/nLociPerInd.tsv") #Load nLoci per sample generated previously
Design<-read.table("data/Design.tsv") #Load design file

#############################################

max(nLociData)
filterCutOff<-mean(nLociData$tot)*0.25 #Filtering cut off. Can be anything from a number of loci (in this case)
                                       #or based on the max / mean / min etc. 
                                       #I Haven't figured out what works best yet.
ind_miss<-0.8       #Fraction of individuals that should have 10x coverage at a certain site for it to be included

#####################
locs<-as.list(list.files("output/methylation_calling/",full.names = T)[grep("CX_report.txt$",list.files("output/methylation_calling/",full.names = T))]) #reads the data
indNamesVector<-sub("_trimmed.*$","",sub("^.*/","",locs))
indNames<-as.list(indNamesVector)




#Plots where the filtering cut off occurs and how many individuals would get filtered
#Check output location of the plot
ggplot(nLociData,aes(x=tot))+geom_histogram(binwidth = 50000,fill="white",col="black")+
  theme_classic()+
  geom_vline(xintercept = filterCutOff,col="red",size=2)+xlab("Number of Loci")+
  ggtitle(paste("n Individuals filtered =",
                length(nLociData$tot[nLociData$tot<filterCutOff]),"out of",
                length(nLociData[,1])))
  ggsave("results/nLociDistrib.png")
#################################



nLociDataFiltFinal<-nLociData[nLociData$tot>filterCutOff,] #Filters the samples that don't meet the                                                               filtering criteria
nLociDataFiltFinal$indNumber<-sub("^.*_","",nLociDataFiltFinal$sampleID)
#You can add other filtering steps here. In this experiment every sample was twice. So we only kept samples which passed filters twice
nLociDataFiltFinal<-nLociDataFiltFinal[nLociDataFiltFinal$indNumber%in%names(which(table(nLociDataFiltFinal$indNumber) > 1)),]





##########################################################
merged<-merge(nLociDataFiltFinal,Design,by="sampleID") #Filters the design



aggregate(tot~site+treatment,
          merged,length) #Creates a table with the number of individuals                                                        passing filters per experimental treatment,                                                           need to change the treatment into the column                                                          names in your design file

write.table(merged,"results/Design_filtered.tsv") #writes the filtered Design

#################
#Create the methylKit databases for each different context

finalMethylationFilesCG<-reorganize(methylationFilesCG, sample.ids=nLociDataFiltFinal$sampleID, 
                                    treatment = c(1:nrow(nLociDataFiltFinal))) #Keep individuals that pass filtering
filteredCG<-filterByCoverage(finalMethylationFilesCG, lo.count = 10,  #Filter sites for coverage
                             lo.perc = NULL, hi.count = NULL, context = "CpG", hi.perc = 99.9)
unitedCG<-unite(filteredCG,min.per.group=as.integer(nrow(nLociDataFiltFinal)*ind_miss),suffix="CG")#Unite


finalMethylationFilesCHG<-reorganize(methylationFilesCHG, sample.ids=nLociDataFiltFinal$sampleID, 
                                     treatment = c(1:nrow(nLociDataFiltFinal)))
filteredCHG<-filterByCoverage(finalMethylationFilesCHG, lo.count = 10, 
                              lo.perc = NULL, hi.count = NULL, context = "CpG", hi.perc = 99.9)
unitedCHG<-unite(filteredCHG,min.per.group=as.integer(nrow(nLociDataFiltFinal)*ind_miss),suffix="CHG")



finalMethylationFilesCHH<-reorganize(methylationFilesCHH, sample.ids=nLociDataFiltFinal$sampleID, 
                                     treatment = c(1:nrow(nLociDataFiltFinal)))
filteredCHH<-filterByCoverage(finalMethylationFilesCHH, lo.count = 10, 
                              lo.perc = NULL, hi.count = NULL, context = "CpG", hi.perc = 99.9)
unitedCHH<-unite(filteredCHH,min.per.group=as.integer(nrow(nLociDataFiltFinal)*ind_miss),suffix="CHH")


saveRDS(unitedCG, file = "results/united.filtered.CG.RData")
saveRDS(unitedCHG, file = "results/united.filtered.CHG.RData")
saveRDS(unitedCHH, file = "results/united.filtered.CHH.RData")
