library(dplyr)
library(data.table)

readResult<-function(resultFile,dirName){
  dt<-fread(paste0(dirName,"//",resultFile),sep=",",header=TRUE)
  dt[,name:=resultFile]
  browser()
  writedt<-dt[.(iters,cost,accuracy,ari,elapsed,InitName,FCM,DATASET)];
  return(writedt)
}

directory<-dirname(sys.frame(1)$ofile)
resultDir<-paste0(directory,"//results")
files<-list.files(resultDir)

resultList<-lapply(files,
                   function(x,dirName) fread(paste0(dirName,"//",x),
                                             sep=",",header=TRUE),
                   dirName=resultDir)
#browser()
finalRes<-rbindlist(resultList)

finalRes<-finalRes[,which(duplicated(colnames(finalRes))):=NULL]
finalRes<-finalRes[,.(iters.mean=mean(iters),iters.std=sd(iters),
                      cost.mean=mean(cost),cost.std=sd(cost),
                      accuracy.mean=mean(accuracy),accuracy.std=sd(accuracy),
                      ari.mean=mean(ari),ari.std=sd(ari),
                      elapsed.mean=mean(elapsed),elapsed.std=sd(elapsed)
                      ),by=.(DATASET,FCM,InitName)]
#finalRes<-finalRes[,.(iters,cost,accuracy,ari,elapsed,InitName),by=.(DATASET,FCM)]
#browser()
write.csv(finalRes,paste0(resultDir,"//finalResult.csv"))