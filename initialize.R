library("cluster")
library(xlsx)
directory<-dirname(sys.frame(1)$ofile)
#source("C://Users//wangc//workspace//EFCMPSO//FCM.R")
source(paste0(directory,"//FCM.R"))
dirPath<-paste0(directory,"//Dataset//")
#dirPath<-"C://Users//wangc//LaptopBackUp//Wangjin//UmassMed//Code//Dataset//"
dataFileList<-list();
id<-1
# # ##***************iris**********************
# dataFile<-"iris.data"
# dataName<-paste(dirPath,dataFile,sep="")
# #browser()
# dataFile<-list(name=dataName,fmt="txt",sep=",",label=5,cn=3,attrCol=-5,attrTypNum=1:4)
# 
# dataFileList[[id]]<-dataFile
# id<-id+1
##********************bank market******************
dataFile<-"bankMarket//bank-full.csv"
dataName<-paste0(dirPath,dataFile)
dataFile<-list(name=dataName,fmt="csv",sep=";",
               label=17,cn=2,attrCol=c(1,12,13,14,15)
               )

dataFileList[[id]]<-dataFile
id<-id+1

##*****************credit card client****************
dataFile<-"credit_card_client//default_of_credit_card_clients.csv"
dataName<-paste0(dirPath,dataFile)
dataFile<-list(name=dataName,fmt="csv",sep=",",
               label=25,cn=2,attrCol=c(2,6:24))

dataFileList[[id]]<-dataFile
id<-id+1
##****************Shuttle*************************
dataFile<-"shuttle//shuttle.trn"
dataName<-paste0(dirPath,dataFile)
dataFile<-list(name=dataName,fmt="csv",sep=" ",label=10,cn=7,attrCol=1:9)

dataFileList[[id]]<-dataFile
id<-id+1
##****************Occupancy*************************
dataFile<-"occupancy//datatraining.csv"
dataName<-paste0(dirPath,dataFile)
dataFile<-list(name=dataName,fmt="csv",sep=",",label=8,cn=2,attrCol=3:7)

dataFileList[[id]]<-dataFile
id<-id+1
##******************MAGIC gamma telescope ***********************
dataFile<-"MAGIC_Gamma_Telescope//magic04.data"
dataName<-paste0(dirPath,dataFile)
dataFile<-list(name=dataName,fmt="csv",sep=",",label=11,cn=2,attrCol=1:10)

dataFileList[[id]]<-dataFile
id<-id+1
##******************skin segmentation ***********************
dataFile<-"skin_segmation//Skin_NonSkin.csv"
dataName<-paste0(dirPath,dataFile)
dataFile<-list(name=dataName,fmt="csv",sep=",",label=4,cn=2,attrCol=1:3)

dataFileList[[id]]<-dataFile
id<-id+1
# ##****************BreastCancer*************************
# dataFile<-"breast-cancer-wisconsin//breast-cancer-wisconsin.data"
# dataName<-paste0(dirPath,dataFile)
# dataFile<-list(name=dataName,fmt="txt",sep=",",label=11,cn=2,attrCol=2:10)
# 
# dataFileList[[id]]<-dataFile
# id<-id+1
# ##****************pima-indians-diabetes*************************
# dataFile<-"pima-indians-diabetes//pima-indians-diabetes.data"
# dataName<-paste0(dirPath,dataFile)
# dataFile<-list(name=dataName,fmt="txt",sep=",",label=9,cn=2,attrCol=-9)
# 
# dataFileList[[id]]<-dataFile
# id<-id+1
##****************pageblock*************************
# dataFile<-"pageblock//page-blocks.data"
# dataName<-paste0(dirPath,dataFile)
# dataFile<-list(name=dataName,fmt="other",sep=" ",label=11,cn=5,attrCol=-11)
# 
# dataFileList[[id]]<-dataFile
# id<-id+1
##********************Simulated Data****************************
tf<-seq(1,10,2)*10000
for (i in tf){
  dataFile<-sprintf("SimulateData//%dx3.csv",i)
  dataName<-paste0(dirPath,dataFile)
  dataFile<-list(name=dataName,fmt="csv",sep=",",label=3,cn=5,attrCol=-3)
  
  dataFileList[[id]]<-dataFile
  id<-id+1  
}

##********************end************************
singleRun<-function(df,runType,initType,param){
  myData<-ImportData(df);
  ptm<-proc.time()
  initFunc<-match.fun(initType)
  regularized<-regularize(myData[,df$attrCol])
  #browser()
  particles<-initFunc(regularized,param)
  #particles<-initFunc(as.matrix(myData[,df$attrCol]),param)
  #particles<-randomParInit(as.matrix(myData[,1:4]),param,pn,randSeedNum)
  
  #particles<-denParInit(as.matrix(myData[,1:4]),param,pn,dc,krt)
  coreFunc<-match.fun(runType)
  #browser()
  result<-coreFunc(regularized,param,particles)
  #result<-coreFunc(as.matrix(myData[,df$attrCol]),param,particles)
  #result<-fpsofcm(as.matrix(myData[,-1*df$label]),param,particles)
  elapsed<-proc.time()-ptm;
  initLabel<-apply(result$f,1,function(x) match(max(x),x))
  ari<-calibrate(initLabel,myData[,df$label])
  return(list(iters=result$iters,cost=result$cost,accuracy=ari$acc,ari=ari$ari,elapsed=elapsed))
}

mcRun<-function(paramDataFile,mcNum,runType,initType,resultDir){
  param<-lapply(1:mcNum,function(i) list(m=2,e=1e-5,cn=paramDataFile$cn,dc=2,pn=10,krt=10,randSeed=i))
  result<-lapply(param,singleRun,df=paramDataFile,runType=runType,initType=initType)
  res<-data.frame(matrix(unlist(result),nrow=length(result),byrow=TRUE))
  colnames(res)<-c("iters","cost","accuracy","ari","elapsed")
  fileFg<-last(unlist(strsplit(paramDataFile$name,"//")))
  
  resFile<-paste0(resultDir,"//",initType,"_",runType,"_",fileFg,".csv")
  #browser()
  write.csv(res,resFile,sep=",",col.names=TRUE)
}
initName<-c("randomParInit","denParInit")
coreName<-c("fpsofcm","fpso","fcm")
#coreName<-c("fpso")
saveIn<-paste0(directory,"//results")
for (init in initName){
  for (cf in coreName){
    lapply(dataFileList,mcRun,mcNum=1,cf,init,saveIn)    
  }  
}
  


# myData<-ImportData(dataFile);
# param<-list(m=2,e=1e-5,cn=dataFile$cn)
# # result<-fcm(as.matrix(myData[,1:4]),param)
# # initLabel<-apply(result$f,1,function(x) match(max(x),x))
# # ari<-calibrate(initLabel,myData[,dataFile$label])
# pn<-10;
# particles<-randomParInit(as.matrix(myData[,1:4]),param,pn)
# dc<-2
# krt<-pn
# 
# particles<-denParInit(as.matrix(myData[,1:4]),param,pn,dc,krt)
# 
# result<-fpsofcm(as.matrix(myData[,1:4]),param,particles)
# initLabel<-apply(result$f,1,function(x) match(max(x),x))
# ari<-calibrate(initLabel,myData[,dataFile$label])

# test<-fanny(as.matrix(myData[,1:4]),3,diss=F)
#  initLabel<-apply(test$membership,1,function(x) match(max(x),x))
#  ari<-calibrate(initLabel,myData[,dataFile$label])
