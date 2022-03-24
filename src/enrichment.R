rm(list=ls())

annotatedDMRs<-read.table("results/DMRs.tsv")
annot<-read.table("results/coveredAnnotation.tsv",h=T)


m<-sum(annot$repMasVerySimple=="transposon")
n<-sum(!annot$repMasVerySimple=="transposon")
k<-nrow(annotatedDMRs)
x<-1:k
diffMethTransposons<-sum(annotatedDMRs$repMasVerySimple=="transposon")

probabilities <- dhyper(x, m, n, k, log = FALSE)
plot(probabilities)
lines(x=rep(diffMethTransposons,2),y=0:1,col="red")
pvalue <- sum(probabilities[diffMethTransposons:length(probabilities)])
pvalue
