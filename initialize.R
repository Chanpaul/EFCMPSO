library("cluster")
library(xlsx)
library(data.table)
library(pryr)
library(Rcpp)
Sys.setenv(SPARK_HOME="C://Users//wangj//Downloads//spark-2.0.2-bin-hadoop2.7")
.libPaths(c(file.path(Sys.getenv("SPARK_HOME"),"R","lib"),.libPaths()))

library(SparkR)
#Rcpp.package.skeleton("psofcm")

directory<-dirname(sys.frame(1)$ofile)
#source("C://Users//wangc//workspace//EFCMPSO//FCM.R")
source(paste0(directory,"//FCM.R"))
dirPath<-paste0(directory,"//Dataset//")
#dirPath<-"C://Users//wangc//LaptopBackUp//Wangjin//UmassMed//Code//Dataset//"
dataFileList<-list();
id<-1
#browser()

##********************QP******************
dataFile<-"ppp_zero_wide//QP_dur_clusterlabel.csv"
dataName<-paste0(dirPath,dataFile)
dataFile<-list(name=dataName,fmt="csv",sep=",",
               label=27,cn=4,attrCol=5:22
)
dataFileList[[id]]<-dataFile
id<-id+1
##*******************TDTA*******************
dataFile<-"TDTA//TDTA_mi10_10_dur_idx_single.csv"
dataName<-paste0(dirPath,dataFile)
dataFile<-list(name=dataName,fmt="csv",sep=",",
               label=13,cn=3,attrCol=4:12
)
dataFileList[[id]]<-dataFile
id<-id+1

##********************glass******************
dataFile<-"glass//glass.data"
dataName<-paste0(dirPath,dataFile)
dataFile<-list(name=dataName,fmt="csv",sep=",",
               label=11,cn=6,attrCol=2:10
)
dataFileList[[id]]<-dataFile
id<-id+1
##********************Image segmentation******************
# dataFile<-"imageSegmentation//segmentation.data"
# dataName<-paste0(dirPath,dataFile)
# dataFile<-list(name=dataName,fmt="csv",sep=",",
#                label=1,cn=7,attrCol=2:20
# )
# dataFileList[[id]]<-dataFile
# id<-id+1

# # ##***************iris**********************
dataFile<-"iris//iris.data"
dataName<-paste(dirPath,dataFile,sep="")
#browser()
dataFile<-list(name=dataName,fmt="txt",sep=",",
               label=5,cn=3,attrCol=-5)

dataFileList[[id]]<-dataFile
id<-id+1
#
##******************MAGIC gamma telescope ***********************
dataFile<-"MAGIC_Gamma_Telescope//magic04.data"
dataName<-paste0(dirPath,dataFile)
dataFile<-list(name=dataName,fmt="csv",sep=",",label=11,cn=2,attrCol=1:10)

dataFileList[[id]]<-dataFile
id<-id+1
# ##****************inonsphere*************************
dataFile<-"ionosphere//ionosphere.data"
dataName<-paste0(dirPath,dataFile)
dataFile<-list(name=dataName,fmt="txt",sep=",",label=35,cn=2,attrCol=c(-2,-35))

dataFileList[[id]]<-dataFile
id<-id+1
# ##****************Simulated inonsphere*************************
# dataFile<-"ionosphere//inon10000x35.csv"
# dataName<-paste0(dirPath,dataFile)
# dataFile<-list(name=dataName,fmt="csv",sep=",",label=35,cn=2,attrCol=c(-2,-35))
# 
# dataFileList[[id]]<-dataFile
# id<-id+1

# ##****************pima-indians-diabetes*************************
dataFile<-"pima-indians-diabetes//pima-indians-diabetes.data"
dataName<-paste0(dirPath,dataFile)
dataFile<-list(name=dataName,fmt="txt",sep=",",label=9,cn=2,attrCol=-9)

dataFileList[[id]]<-dataFile
id<-id+1
# ##****************wine*************************
# dataFile<-"wine//wine.data"
# dataName<-paste0(dirPath,dataFile)
# dataFile<-list(name=dataName,fmt="txt",sep=",",label=1,
#                cn=3,attrCol=2:13)
# 
# dataFileList[[id]]<-dataFile
# id<-id+1
##********************Simulated Data****************************
# simDir<-paste0(dirPath,"QPSim");
# simFile<-grep("\\d+x\\d+.csv",list.files(simDir),perl=TRUE,value=TRUE);
# for (i in simFile){
#   dataName<-paste0(simDir,"//",i);
#   dataFile<-list(name=dataName,fmt="csv",sep=",",label=19,cn=4,attrCol=-19);
#   dataFileList[[id]]<-dataFile
#   id<-id+1
# }

# ##****************BreastCancer*************************
dataFile<-"breast-cancer-wisconsin//breast-cancer-wisconsin.data"
dataName<-paste0(dirPath,dataFile)
dataFile<-list(name=dataName,fmt="txt",sep=",",label=11,cn=2,attrCol=2:10)

dataFileList[[id]]<-dataFile
id<-id+1
# 
# ##****************pageblock*************************
dataFile<-"pageblock//page-blocks.data"
dataName<-paste0(dirPath,dataFile)
dataFile<-list(name=dataName,fmt="others",sep=" ",label=11,cn=5,attrCol=-11)

dataFileList[[id]]<-dataFile
id<-id+1
##********************Simulated Pageblock Data****************************
# tf<-c(1,1.5,2)*10000
# for (i in tf){
#   dataFile<-sprintf("pageblock//%dx11.csv",i)
#   dataName<-paste0(dirPath,dataFile)
#   dataFile<-list(name=dataName,fmt="csv",sep=",",label=11,cn=5,attrCol=-11)
# 
#   dataFileList[[id]]<-dataFile
#   id<-id+1
# }
# browser()
##********************end************************
postProcess<-function(src,des){
  require(data.table)
  f<-list.files(path=src,full.names=TRUE);
  res<-lapply(f,
              function(fn) {temp<-fread(fn,sep=",",header=TRUE);
              temp[,"type"]=tail(unlist(strsplit(fn,"/")),n=1);
              return(temp)})
  fres<-res[[1]]
  for(i in 2:length(res)){
    fres<-rbind(fres,res[[i]])
  }
  #browser()
  write.csv(fres,des,sep=",")
}


singleRun<-function(df,runType,initType,param){
  
  myData<-na.omit(ImportData(df));
  #browser()
  myData<-batchCategoryToReal(myData);
  
  ptm<-proc.time()
  initFunc<-match.fun(initType)
  #browser()
  #regularized<-regularize(myData[,df$attrCol])
  regularized<-myData[,df$attrCol]
  #browser()
  particles<-initFunc(regularized,param)
  #browser()
  #particles<-initFunc(as.matrix(myData[,df$attrCol]),param)
  #particles<-randomParInit(as.m atrix(myData[,1:4]),param,pn,randSeedNum)
  
  #particles<-denParInit(as.matrix(myData[,1:4]),param,pn,dc,krt)
  coreFunc<-match.fun(runType)
  #browser()
  result<-coreFunc(regularized,param,particles)
  #result<-coreFunc(as.matrix(myData[,df$attrCol]),param,particles)
  #result<-fpsofcm(as.matrix(myData[,-1*df$label]),param,particles)
  elapsed<-proc.time()-ptm;
  initLabel<-apply(result$f,1,function(x) match(max(x),x))
  ari<-calibrate(initLabel,myData[,df$label])
  rm(myData)
  return(list(iters=result$iters,cost=result$cost,accuracy=ari$acc,ari=ari$ari,elapsed=elapsed))
}

mcRun<-function(paramDataFile,mcNum,runType,initType,resultDir){
  param<-lapply(1:mcNum,function(i) list(m=2,e=1e-5,cn=paramDataFile$cn,dc=6,pn=10,krt=10,randSeed=i))
  
  result<-lapply(param,singleRun,df=paramDataFile,runType=runType,initType=initType)
  #browser()
  res<-data.frame(matrix(unlist(result),nrow=length(result),byrow=TRUE))
  colnames(res)<-c("iters","cost","accuracy","ari","elapsed")
  fileFg<-last(unlist(strsplit(paramDataFile$name,"//")))
  
  resFile<-paste0(resultDir,"//",initType,"_",runType,"_",fileFg,".csv")
  #browser()
  resdatatable<-setDT(res)
  resdatatable[,c("InitName","FCM","DATASET"):=list(initType,runType,fileFg)]
  write.csv(resdatatable,resFile,sep=",",col.names=TRUE)
  
}
# initName<-c("denParInit")#,"randomParInit")
# coreName<-c("fpsofcm")#,"fpso","fcm")
initName<-c("randomParInit")
coreName<-c("fpsofcm","fcm")
#coreName<-c("fpso")
saveIn<-paste0(directory,"//results")
for (init in initName){
  for (cf in coreName){
    lapply(dataFileList,mcRun,mcNum=5,cf,init,saveIn)    
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
