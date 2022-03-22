rm(list=ls())
library(methylKit)
library(ggplot2)
library(reshape2)
library(ggpubr)
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #Make pretty colours
#################
#performs basic percentage methylation and PCA analysis
#Make sure to change variables specifying colour and shape in plots to match the variables
#In your own Design table
#################
unitedCG<-readRDS("results/united.filtered.CG.RData")
unitedCHG<-readRDS("results/united.filtered.CHG.RData")
unitedCHH<-readRDS("results/united.filtered.CHH.RData")
#Loads the filtered design

Design<-read.table("results/Design_filtered.tsv")


i<-0
percMethData<-pcaPlot<-list()

#Reads the united data for each different context
for (context in c("CG","CHG","CHH")){
  i<-i+1
  # #########################################
  # # PCA's for each context
  # pca<-PCASamples(get(paste0("united",context)),scale =T,obj.return = T,screeplot = T) #create the PCA
  # 
  # pcaPerc<-round(pca$sdev/sum(pca$sdev)*100,1) #Calculate percentage of variante per PC
  # pcaData<-data.frame(PC1=pca$x[,1],PC2=pca$x[,2],PC3=pca$x[,3],indName=row.names(pca$x))
  # pcaPlotData<-merge(pcaData,Design,by="indName")
  # pcaPlot[[i]]<-ggplot(pcaPlotData,aes(x=PC1,y=PC2,col=site,shape=treatment.y))+
  #   geom_point(size=3)+theme_light()+
  #   scale_colour_manual(values = colorBlindBlack8)+ggtitle(context)+
  #   xlab(paste0("PC1 (",pcaPerc[1],"%)"))+
  #   ylab(paste0("PC2 (",pcaPerc[2],"%)"))
  # ##########################################
  # Percentage methylation dataframes
  percMethData[[i]]<-as.data.frame(t(percMethylation(get(paste0("united",context)))))
  
  
}

# g <- ggarrange(pcaPlot[[1]],pcaPlot[[2]],pcaPlot[[3]], ncol=3, nrow=1, widths=c(6,6), align="h", labels="AUTO", common.legend = TRUE, legend="right", legend.grob=get_legend(pcaPlot[[1]]))
# g
# ggsave("results/pca.png",g,height=4.5,width=9,dpi="retina")


#Combine the percentage methylation with the Design.
percMethDataOut<-data.frame(sampleID=row.names(percMethData[[1]]),
                            CG=rowMeans(percMethData[[1]],na.rm=T),
                            CHG=rowMeans(percMethData[[2]],na.rm=T),
                            CHH=rowMeans(percMethData[[3]],na.rm=T))
percMethDataOut$num<-sub("^.*_","",percMethDataOut$sampleID)

meltPlot<-melt(percMethDataOut,id.vars=c("sampleID","num"))
percMethPlot<-merge(meltPlot,Design,by="sampleID")

ggplot(percMethPlot,aes(x=treatment,y=value,col=site))+geom_boxplot()+facet_grid(.~variable,scales="free",space = "free")+ylab("percentage methylation")+
  scale_colour_manual(values = colorBlindBlack8)+theme_light()
  ggsave("results/percMeth.png")



  


ggplot(percMethPlot,aes(x=site,y=value,col=treatment,group=num))+geom_point()+geom_line()+facet_grid(variable~treatment,scales="free")+theme_light()+
  scale_colour_manual(values=colorBlindBlack8)
ggsave("results/perPairDiffMeth.png",width=14,height=9)
  