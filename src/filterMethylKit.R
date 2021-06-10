library(methylKit)
library(ggplot2)
library(reshape2)
library(ggpubr)
#########################################

#Takes output from the readMethylKit + a design table:
#The Design table should contain a column called sampleID with the names of the samples (same as the names used in epiGBS) and can contain all other treatments etc. These will be combined with nLoci data after filtering to create filtered design. Which you can use to determine what filtering steps leave enough samples per treatment group
#


##############################################
nLociData<-read.table("results/nLociPerInd.tsv") #Load nLoci per sample generated previously
Design<-read.table("data/Design.tsv") #Load design file

#############################################


filterCutOff<-50000 #Filtering cut off. Can be anything from a number of loci (in this case)
                                       #or based on the max / mean / min etc. 
                                       #I Haven't figured out what works best yet.
ind_miss<-0.8            #Fraction of individuals a 


#####################

#Plots where the filtering cut off occurs and how many individuals would get filtered
#Check output location of the plot
ggplot(nLociData,aes(x=nloci))+geom_histogram(binwidth = 5000,fill="white",col="black")+
  theme_classic()+
  geom_vline(xintercept = filterCutOff,col="red",size=2)+
  ggtitle(paste("n Individuals filtered =",
                length(nLociData$nloci[nLociData$nloci<filterCutOff]),"out of",
                length(nLociData[,1])))+
  ggsave("results/nLociDistrib.png")
#################################



nLociDataFiltFinal<-nLociData[nLociData$nloci>filterCutOff,] #Filters the samples that don't meet the                                                               filtering criteria

#You can add other filtering steps here. In this experiment every sample was twice. So we only kept samples which passed filters twice
nLociDataFiltFinal<-nLociDataFiltFinal[nLociDataFiltFinal$indNumber%in%names(which(table(nLociDataFiltFinal$indNumber) > 1)),]





##########################################################
merged<-merge(nLociDataFiltFinal,Design,by.x="indName",by.y="sampleID") #Filters the design



aggregate(nloci~site+treatment.y,
          merged,length) #Creates a table with the number of individuals                                                        passing filters per experimental treatment,                                                           need to change the treatment into the column                                                          names in your design file

write.table(merged,"results/Design_filtered.tsv") #writes the filtered Design

#################
#Create the methylKit databases for each different context
methylationFilesCG<-methRead(location = as.list(nLociDataFiltFinal$loc),
                             as.list(nLociDataFiltFinal$indName),
                             assembly = "arabidopsisRef",
                             pipeline = "bismarkCytosineReport",
                             treatment = rep(1,times=length(nLociDataFiltFinal$nloci)),
                             mincov = 10,
                             dbdir = "FiltCG")
methylDBCG<-makeMethylDB(methylationFilesCG,dbdir = "FiltCG")

normalizedCG<-normalizeCoverage(methylDBCG,dbdir="FiltCG")
unitedCG<-unite(normalizedCG,min.per.group=as.integer(length(nLociDataFiltFinal$loc)*ind_miss),suffix="CG")



methylationFilesCHG<-methRead(location = as.list(nLociDataFiltFinal$loc),
                              as.list(nLociDataFiltFinal$indName),
                              assembly = "arabidopsisRef",
                              pipeline = "bismarkCytosineReport",
                              treatment = rep(1,times=length(nLociDataFiltFinal$nloci)),
                              mincov = 10,context="CHG",
                              dbdir = "FiltCHG")
methylDBCHG<-makeMethylDB(methylationFilesCHG,dbdir = "FiltCHG")

normalizedCHG<-normalizeCoverage(methylDBCHG,dbdir="FiltCHG")
unitedCHG<-unite(normalizedCHG,
                 min.per.group=as.integer(length(nLociDataFiltFinal$loc)*ind_miss),suffix="CHG")



# unitedCHG<-readMethylDB("/scratch2/maarten/scabiosa/FiltCHG/methylBase_2311ef4affd682.txt.bgz")


# 
methylationFilesCHH<-methRead(location = as.list(nLociDataFiltFinal$loc),
                              as.list(nLociDataFiltFinal$indName),
                              assembly = "arabidopsisRef",
                              pipeline = "bismarkCytosineReport",
                              treatment = rep(1,times=length(nLociDataFiltFinal$nloci)),
                              mincov = 10,context="CHH",
                              dbdir = "FiltCHH")
methylDBCHH<-makeMethylDB(methylationFilesCHH,dbdir = "FiltCHH")

normalizedCHH<-normalizeCoverage(methylDBCHH,dbdir="FiltCHH")
unitedCHH<-unite(normalizedCHH,min.per.group=as.integer(length(nLociDataFiltFinal$loc)*ind_miss),suffix="CHH")

