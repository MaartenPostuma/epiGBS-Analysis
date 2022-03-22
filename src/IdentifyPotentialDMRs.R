data<-read.table("results/subDSS.tsv")#This is the data.frame with the significant DMCs


windowSize<-500 #This is the number of bases you look to (forwards)
dataOut<-data.frame()
total<-0



for(i in 1:nrow(data)){
if(total==0){
  startRegionChr<-data$chr[i]
  startRegionPos<-data$start[i]
}  
counts<-nrow(data[data$chr==data$chr[i]&data$start>data$start[i]&data$start<data$start[i]+windowSize,])
if(counts==0){
    dataOut<-rbind(dataOut,data.frame(startChr=startRegionChr,startPos=startRegionPos,endChr=data$chr[i],endPos=data$start[i],count=total+1))
    total<-0}else{
          total<-total+1
}
  
}
