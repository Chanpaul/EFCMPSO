library(MASS)
library(dplyr)
library(data.table)
library(xlsx)

singleGen<-function(dataNum,param){
  totalNum<-dataNum
  ratio<-param$classRatio
  cn<-param$cn
  tDir<-param$dataDir
  mu<-param$mu
  sigma<-param$sigma
  dt<-lapply(1:cn,function(i,p_mu,p_sigma,p_ratio) 
    as.data.table(mvrnorm(n=totalNum*ratio[i],
                          unlist(mu[i]),
                          sigma[[i]],
                          tol=1e-6,
                          empirical=FALSE))%>%mutate(className=paste0("C",i)),
             p_mu=mu,
             p_sigma=sigma,
             p_ratio=ratio);
  dt<-rbindlist(dt)
  dataFile<-paste0(tDir,nrow(dt),"x",ncol(dt),".csv");
  #browser()
  write.csv(dt,dataFile,sep=",",col.names=FALSE,row.names=FALSE)
}
#****************Method I based on specific dataset**************#
directory<-dirname(sys.frame(1)$ofile)
dtMiddleName<-"//Dataset//QPSim//"
dataDir<-paste0(directory,dtMiddleName);
dtFileName<-"fittedParam.csv"
dtName<-paste0(directory,dtMiddleName,dtFileName);
#gaussianParam<-read.xlsx(dtName,sheetIndex=1)
gaussianParam<-fread(dtName,sep=",",header=TRUE)
className<-c("C1","C2","C3","C4");
mu<-lapply(className,function(cname,dt) dt[METRIC=="avg",.(eval(parse(text=cname)))],dt=gaussianParam)
sigma<-lapply(className,function(cname,dt) diag(dt[METRIC=="std",eval(parse(text=cname))]), 
              dt=gaussianParam)
#browser()
classRatio<-c(0.27,0.14,0.37,0.22);
cn<-4;
params<-list(classRatio=classRatio,
             cn=cn,dataDir=dataDir,
             mu=mu,sigma=sigma)
lenOfDatasets<-seq(1,10,2)*1000
lapply(lenOfDatasets,singleGen,param=params)

#****************Method II based on specific parameters**************#
# mu<-list(list(-0.4004,1.1997),
#          list(-0.6385,0.9283),
#          
#          list(-0.1972,0.9012),
#          list(-0.5480,0.4899),
#          
#          list(-0.4078,0.7603),
#          list(-0.3544,0.6999),
#          list(-0.3085,0.6602))
# sigma<-list(matrix(c(0.8993,0.0320,0.0320,0.7424),nrow=2,byrow=TRUE)*1e-3,
#             matrix(c(0.0008,0.0007,0.0007,0.0019),nrow=2,byrow=TRUE),
#             
#             matrix(c(0.5803,-0.0151,-0.0151,0.5789),nrow=2,byrow=TRUE)*1e-3,
#             matrix(c(0.0023,-0.0023,-0.0023,0.0046),nrow=2,byrow=TRUE),
#             
#             matrix(c(0.0014,-0.0009,-0.0009,0.0015),nrow=2,byrow=TRUE),
#             matrix(c(0.0030,-0.0017,-0.0017,0.0027),nrow=2,byrow=TRUE),
#             matrix(c(0.0020,-0.0015,-0.0015,0.0020),nrow=2,byrow=TRUE))
# singleGen<-function(param){
#   totalNum<-param$num
#   ratio<-param$classRatio
#   cn<-param$cn
#   tDir<-param$dataDir
#   dt<-data.frame()
#   for (i in 1:length(mu)){
#     tempDt<-mvrnorm(n=totalNum*ratio[i],unlist(mu[i]),sigma[[i]],tol=1e-6,empirical=FALSE) 
#     tempDt<-data.frame(tempDt)%>%mutate(class=min(i,cn))
#     dt<-rbind(dt,tempDt)
#   }
#   dataFile<-paste0(tDir,nrow(dt),"x",ncol(dt),".csv");
#   #browser()
#   write.csv(dt,dataFile,sep=",",col.names=FALSE,row.names=FALSE)
# }
# 
# 
# dataDir<-"C://Users//wangc//LaptopBackUp//Wangjin//UmassMed//Code//Dataset//SimulateData//"
# classRatio<-c(0.06,0.1,0.18,0.13,0.26,0.1,0.17)
# cn<-5;
# param<-list(num=0,classRatio=classRatio,cn=cn,dataDir=dataDir)
# lenOfDatasets<-seq(1,10,2)*10000
# params<-lapply(lenOfDatasets,function(tdataDir,tclassRatio,tcn,tnum) list(num=tnum,
#                                                                           classRatio=tclassRatio,
#                                                                           cn=tcn,
#                                                                           dataDir=tdataDir),
#                tclassRatio=classRatio,
#                tcn=cn,
#                tdataDir=dataDir)
# 
# lapply(params,singleGen)


