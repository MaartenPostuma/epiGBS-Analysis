data<-read.table("yourDMCs.tsv") #This should be the file with your DMCs

windowSize<-500 #This is the number of bases you look to (forwards)
dataOut<-data.frame()
total<-0


pb <- txtProgressBar(min = 0, max = nrow(data), style = 3)

for(i in 1:nrow(data)){
  if(total==0){
    posRegionChr<-data$chr[i]
    posRegionPos<-data$pos[i]
  }  
  counts<-nrow(data[which(data$chr==data$chr[i]&data$pos>data$pos[i]&data$pos<(data$pos[i]+windowSize)),])
  if(counts==0){
    dataOut<-rbind(dataOut,data.frame(posChr=posRegionChr,posPos=posRegionPos,endChr=data$chr[i],endPos=data$pos[i],count=total+1))
    total<-0}else{
      total<-total+1
    }
  setTxtProgressBar(pb, i)
}



