rm(list=ls())
library(methylKit)
library(reshape2)

#Adjusted from seepers et al. (2021?)
Design<-read.table("results/Design_filtered.tsv") #Read design table to eventually combine with methylation data

for(context in c("CG","CHH","CHG")){ #Make the data for all 3 contexts
united<-readMethylDB(paste0("Filt",context,"/methylBase_",context,".txt.bgz")) #Read the methylBase objects made when filtering

methylation<-getData(united) #Load data


columnsC<- grep(pattern="^numCs", colnames(methylation)) #Get num C's
columnsT<-grep(pattern="^numTs", colnames(methylation)) #Get num T's
methylation.Pos<-methylation[,c(1,3,columnsC,columnsT)] #Make data.frame with CHR_POS and and Cs and Ts columns


firstC<-which(colnames(methylation.Pos)=="numCs1") #Determine which column is the first C colummn
lastC<-which(colnames(methylation.Pos)==paste0("numCs",length(columnsC))) #Determine the last C column
methPatternLong<-melt(methylation.Pos[c(1,2,c(firstC:lastC))], id.vars = c(1,2), 
                                variable.name = "ID", value.name = "numCs")#Make the data long format



firstT<-which(colnames(methylation.Pos)=="numTs1") #Determine which column is the first C colummn
lastT<-which(colnames(methylation.Pos)==paste0("numTs",length(columnsT))) #Determine the last C column
methPatternLongT<-melt(methylation.Pos[c(1,2,c(firstT:lastT))],id.vars=c(1,2),
                       variable.name = "ID", value.name = "numTs") #Make the data long format

methPatternLong$numTs<-methPatternLongT$numTs #Add the Ts column to the C columns

methPatternLong$ID<-gsub(methPatternLong$ID, pattern = "numCs", #Change IDs to something more like IDs 
                         replacement = "ID_")           



IDFrame<-data.frame(indName=united@sample.ids,ID=paste0("ID_",1:length(united@sample.ids))) #Get the real IDs
DesignWithIDs<-merge(IDFrame,Design,by="indName") #Merge real IDs with IDs created above

methPatternLongCorrectID<-merge(methPatternLong,DesignWithIDs,by="ID") #Combine Design with methylation data

write.table(methPatternLongCorrectID,paste0("Filt",context,"/longFormat",context,".tsv")) #Write the methylation data
}

