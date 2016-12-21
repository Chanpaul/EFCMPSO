##********************bank market******************
dataFile<-"bankMarket//bank-full.csv"
dataName<-paste0(dirPath,dataFile)
dataFile<-list(name=dataName,fmt="csv",sep=";",
               label=17,cn=2,attrCol=1:16
)
dataFileList[[id]]<-dataFile
id<-id+1
##*****************credit card client****************
dataFile<-"credit_card_client//default_of_credit_card_clients.csv"
dataName<-paste0(dirPath,dataFile)
dataFile<-list(name=dataName,fmt="csv",sep=",",
               label=25,cn=2,attrCol=2:24)

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
dataFile<-list(name=dataName,fmt="csv",sep=",",label=7,cn=2,attrCol=2:6)

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