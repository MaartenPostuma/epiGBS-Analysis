rm(list=ls())
library(lme4);library(lmerTest);library(optimx);library(base);
library(parallel);library(stringr);library(reshape2)


#The models contain numCs,numTs in stead of flat percentages to account in variation in coverage
nullModel<-as.formula("cbind(numCs, numTs)~1 + (1|indNumber)") #Write your null model
fullModel<-as.formula("cbind(numCs, numTs)~1 + site + (1|indNumber)") #Write your full model


for(context in c("CG","CHG","CHH")){
function_glmm <- function(methIn,nullModel,fullModel) {tryCatch({
  
  meth_temp_final <- as.data.frame(methIn) #load data
  CHR_POS <- as.character(meth_temp_final[1,2]) #Get the position
  
  null <- glmer(nullModel, meth_temp_final, #Create the null model
                family="binomial")
  model <- glmer(fullModel, meth_temp_final, #Create the full mdoel
                 family="binomial")
  
  lrt <- anova(null, model) #Run anova between null and full model
  pvalue=lrt$"Pr(>Chisq)"[2] #Get pvalue
  deviance <- lrt$deviance[2] #Get deviance
  Df <- lrt$Df[2] #get degrees of freedom
  Chisq <- lrt$Chisq[2] #Get ChiSquared
  
  message<-data.frame(model@optinfo$conv$lme4)  #Save all messages output by the glmer
  mess1<-str_c(" ",message$messages[1]) 
  mess2<-str_c(" ",message$messages[2]) 
  mess3<-str_c(" ",message$messages[3])
  
  mess1.1<-str_c(" ",message[1,1])
  mess2.1<-str_c(" ",message[1,2])
  mess3.1<-str_c(" ",message[1,3])
  mess4.1<-str_c(" ",message[1,4])
  
  message1.2 <- data.frame(model@optinfo$message)
  mess1.2 <-str_c(" ",message1.2[1,1])
  
  rdf <- df.residual(model) #Save residuals
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2) #calculate Pearson.chisq
  ratio <- Pearson.chisq/rdf #calculate ratio
  pval_disp <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE) #calculate pval dispersal?
  
  
  sum<-summary(model) #Get model summary
  
  coefs<-t(melt(data.frame(sum$coefficients))) #Save the coefficients of the model
  coefsOut<-data.frame(matrix(coefs[2,],ncol=length(coefs[2,])))
  colnames(coefsOut)<-paste(row.names(sum$coefficients),coefs[1,],sep="_")
  
  random<-t(melt(data.frame(sum$varcor))) #Save the information in the random effects of the model
  randomOut<-data.frame(matrix(random[5,],ncol=length(random[5,])))
  colnames(randomOut)<-paste(random[1,],random[4,],sep="_")
  
  return(data.frame(CHR_POS=CHR_POS, df=Df, Chisq=Chisq, pval_lrt=as.numeric(pvalue), df_residual=rdf, 
                    disp_ratio=ratio, message1=mess1, pval_disp=pval_disp,ratio=ratio,Pearson.chisq=Pearson.chisq,
                    message2=mess2, message3=mess3, message1.1=mess1.1, message2.1=mess2.1, 
                    message3.1=mess3.1, message4.1=mess4.1, message1.2=mess1.2,coefsOut,randomOut))
}, error = function(e){cat("ERROR :", conditionMessage(e), "\n");print(CHR_POS)})
}




methylationData<-read.table(paste0("Filt",context,"/longFormat",context,".tsv")) #Loads the data created using makeGLMMData.R
CHR_POS<-paste(methylationData$chr,methylationData$end)#Creates a vector with all of the loci name
glmm_in<-split(methylationData,f = as.factor(CHR_POS)) #Split the large data.frame in to a list based on loci

glmm_out <- mclapply(glmm_in, function_glmm,nullModel,fullModel, mc.cores=15) 
data_glmm<-do.call("rbind",glmm_out)
data_glmm$qval <- p.adjust(data_glmm$pval_lrt, method="fdr", n = nrow(data_glmm))
write.table(data_glmm,paste0("output/glmer",context,"tsv"))
}

