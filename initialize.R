library("cluster")
source("D://UmassMed//Paper related//IEEE-TBD//Code//FCM.R")
dirPath<-"D://UmassMed//Code//Dataset//"
# ##***************iris**********************
dataFile<-"iris.data"
dataName<-paste(dirPath,dataFile,sep="")
#browser()
dataFile<-list(name=dataName,fmt="txt",sep=",",label=5,cn=3)
##********************Wine******************
# dataFile<-"wine//wine.data"
# dataName<-paste(dirPath,dataFile)
# dataFile<-list(name=dataName,fmt="txt",sep=",",label=1,cn=3)
##*****************Glass****************
# dataFile<-"glass//glass.data"
# dataName<-paste(dirPath,dataFile)
# dataFile<-list(name=dataName,fmt="txt",sep=",",label=11,cn=6,id=1)
##****************ionosphere*************************
# dataFile<-"ionosphere//ionosphere.data"
# dataName<-paste(dirPath,dataFile)
# dataFile<-list(name=dataName,fmt="txt",sep=",",label=35)
##****************BreastCancer*************************
# dataFile<-"breast-cancer-wisconsin//breast-cancer-wisconsin.data"
# dataName<-paste(dirPath,dataFile)
# dataFile<-list(name=dataName,fmt="txt",sep=",",label=11)
##****************pima-indians-diabetes*************************
# dataFile<-"pima-indians-diabetes//pima-indians-diabetes.data"
# dataName<-paste(dirPath,dataFile)
# dataFile<-list(name=dataName,fmt="txt",sep=",",label=9)
##****************pageblock*************************
# dataFile<-"pageblock//page-blocks.data"
# dataName<-paste(dirPath,dataFile)
# dataFile<-list(name=dataName,fmt="txt",sep=" ",label=11)
##****************TDTA*************************
# dataFile<-"01_12_2015//TDTA//TDTA_mi10_10_dur_idx.xls"
# dataName<-paste(dirPath,dataFile)
# dataFile<-list(name=dataName,fmt="excel",sep=" ",label=11)
##********************end************************
#dataPath<-"D://UmassMed//Paper related//Data//iris.data";
myData<-ImportData(dataFile);
#myData[,1:4]<-normalit(myData[,1:4],byrow=TRUE);
#browser();
param<-list(m=2,e=1e-5,cn=3)
# result<-fcm(as.matrix(myData[,1:4]),param)
# initLabel<-apply(result$f,1,function(x) match(max(x),x))
# ari<-calibrate(initLabel,myData[,dataFile$label])
result<-fpsofcm(as.matrix(myData[,1:4]),param)
initLabel<-apply(result$f,1,function(x) match(max(x),x))
ari<-calibrate(initLabel,myData[,dataFile$label])

# test<-fanny(as.matrix(myData[,1:4]),3,diss=F)
#  initLabel<-apply(test$membership,1,function(x) match(max(x),x))
#  ari<-calibrate(initLabel,myData[,dataFile$label])
