library(base)
library(stats)
library(gtools)  # for permu
library(pdfCluster) #for adjusted rand index;

#@param particle = list(pbest=pbest,muX=muX,muV=muV,fitness=0)
#@particle$pbest    the best of the particle;
#@particle$muX      the current position of the particle;
#@particle$muV      the current speed of the particle;
#@particle$fitness  the current fitness of the particle;
#@ param=list(m=m,e=e,cn=cn)
#@m fuziness
#@e error bound
#@cn cluster number

# fcm<- function (dataset,param)  #param is a list
# {
#   X<-dataset;
#   initValidityXB<-10000000;
#   cn<-param$cn;#cn is cluster number
#   m<-param$m;
#   e<-param$e;
#   N<-nrow(X)
#   n<-ncol(X)
#   #****************Initialization randomly select centroids********
#   X1 <-c(rep(1,N));
#   set.seed(1000);
#   centroidIndex=sample(1:N,cn);   #Generate N unrepeatable integers, then we select c of them.
#   v<-c();
#   d<-c();
#   #browser();
#   for (j in 1:cn){
#     tempRow<-X[centroidIndex[j],];
#     v<-rbind(v,tempRow);
#     xv = X - tcrossprod(X1,v[j,])#X1%*%t(v[j,]);
#     tempDist<-apply(xv,1,function(x) sum(x^2));
#     d<-cbind(d,tempDist);
#   }
#   d0<-d;
#   #browser()
#   d <- (d+1e-10)^(-1/(m-1));
#   f0 <- d/tcrossprod(apply(d,1,sum),rep(1,cn));
#   iter <- 0;                       # iteration counter
#   initialJ<-10000000000000000;
#   
#   J0<-sum(diag(t(f0^m)%*%d0));
#   distout<-c();
#   J<-c();
#   while (abs(J0-initialJ)>e){
#     iter = iter + 1;
#     initialJ=J0;
#     # Calculate centers
#     fm <-f0^m;
#     sumf<-apply(fm,2,sum);
#     #sumf = sum(fm);
#     v<-crossprod(fm,X)/(sumf%*%t(rep(1,n)));
#     #v = (fm'*X)./(sumf'*ones(1,n)); 
#     d<-c();
#     for (j in 1:cn){
#       xv = X - tcrossprod(X1,v[j,]);
#       #xv = X - X1*v[j,];
#       tempDist<-apply(xv,1,function(x) sum(x^2));
#       d<-cbind(d,tempDist);
#       #d[,j] = (xv*eye(n).*xv)*(w(j,:).^2)';
#     } 
#     distout<-sqrt(d);
#     d0<-d;  
#     J[iter]<-sum(diag(t(f0^m)%*%d0));
#     d = (d+1e-10)^(-1/(m-1));  #L2 distance
#     f0 <- d/tcrossprod(apply(d,1,sum),rep(1,cn));
#     #f0 = (d/ (sum(d,2)*ones(1,c)));      #method I
#     J0=J[iter];
#   }
#   results<-list(f=f0,d=distout,v=v,iter=iter,cost=J);
# }
ImportData<-function(dataFile)
{
  dataPath=dataFile$name;
  fileType=dataFile$fmt;
  sepStr=dataFile$sep;
  #import excel file
  if (fileType=="excel")
  {
    library(gdata)
    mydata = read.xls(dataPath)    
  }
  else if (fileType=="txt")
  {
    
    mydata<-read.table(dataPath,sep=sepStr)
    
  }
  else if (fileType=="csv")
  {
    mydata=read.csv(dataPath,sep=sepStr,header=TRUE)
  } else {
    mydata<-strsplit(readLines(dataPath),"\n")
    regPat<-"\\d*\\.?\\d+"
    mydata<-lapply(mydata,function(line) as.numeric(grep(regPat,unlist(line),perl=TRUE,value=TRUE)))
    mydata<-data.frame(matrix(unlist(mydata),nrow=length(mydata),byrow=TRUE))
  }
  return(mydata)  
}

calibrate<-function(yte, true_label){
  # Function for calculating clustering accuray and matching found 
  # labels with true labels. Assumes yte and y both are Nx1 vectors with
  # clustering labels. Does not support fuzzy clustering.
  #
  # Algorithm is based on trying out all reorderings of cluster labels, 
  # e.g. if yte = [1 2 2], try [1 2 2] and [2 1 1] so see witch fit 
  # the truth vector the best. Since this approach makes use of perms(),
  # the code will not run for unique(yte) greater than 10, and it will slow
  # down significantly for number of clusters greater than 7.
  #
  # Input:
  #   yte - result from clustering (y-test)
  #   y   - truth vector
  #
  # Output:
  #   accuracy    -   Overall accuracy for entire clustering (OA). For
  #                   overall error, use OE = 1 - OA.
  #   true_labels -   Vector giving the label rearangement witch best 
  #                   match the truth vector (y).
  #   CM          -   Confusion matrix. If unique(yte) = 4, produce a
  #                   4x4 matrix of the number of different errors and  
  #                   correct clusterings done.
  N <- length(true_label);
  cluster_names<-unique(true_label);
  #cluster_names = unique([yte' y']);
  #cluster_names = unique(yte);
  accuracy <- 0;
  maxInd <- 1;
  #browser()
  c_n<-length(cluster_names);
  perm<-permutations(c_n,c_n,as.vector(cluster_names))
  #perm = perms(unique([yte' y']));
  #perm = perms(unique(y));
  pN<-nrow(perm);
  pM<-ncol(perm);
  #[pN pM] = size(perm);
  
  pred_labels <- true_label;
  
  for (i in 1:pN){
    flipped_labels <- rep(0,N) #zeros(1,N);
    for (cl in 1 : pM){
      flipped_labels[yte==cl] = perm[i,cl];
    }
    testAcc = sum(flipped_labels == true_label)/N;
    if (testAcc > accuracy){
      accuracy = testAcc;
      maxInd = i;
      pred_labels = flipped_labels;
      
    }
  }
  corLabel<-perm[maxInd,];  # the corresponding relations between true labels and the test results.
  CM <-matrix(rep(0,pM*pM),nrow=pM,byrow=T)# zeros(pM,pM);
  for (rc in 1 : pM) {
    for (cc in 1 : pM){
      CM[rc,cc] <- sum(((true_label==rc)*(pred_labels==cc)));
    }
  }
  ari<-adj.rand.index(pred_labels,true_label);
  res<-list(ari=ari,acc=accuracy, p_lab=pred_labels, CM=CM, corLab=corLabel);
}
randomParInit<-function(dataset,param){  #random initialization;
  X<-dataset;
  cn<-param$cn;  #cluster number
  N<-nrow(X)     #number of data
  n<-ncol(X)     #number ofattribute
  p<-param$pn;         #number of particle
  m<-param$m
  randSeedNum<-param$randSeed
  gfitness<-0;
  particle<-list();
  for (k in 1:p){
    set.seed(k+randSeedNum*100);
    pbest<-matrix(runif(N*cn,min=0,max=1),nrow=N,byrow=T)
    muX<-matrix(runif(N*cn,min=0,max=1),nrow=N,byrow=T)
    muV<-matrix(runif(N*cn,min=0,max=1),nrow=N,byrow=T)
    #browser()
    pbest<-pbest/tcrossprod(apply(pbest,1,sum),rep(1,cn));
    muX<-muX/tcrossprod(apply(muX,1,sum),rep(1,cn));
    muV<-muV/tcrossprod(apply(muV,1,sum),rep(1,cn));
    
    fm <- pbest^m;
    sumf<-apply(fm,2,sum);
    v<-crossprod(fm,X)/(sumf%*%t(rep(1,n)));
    d<-c()
    X1<-rep(1,N);
    for (j in 1 : cn){
      xv = X - tcrossprod(X1,v[j,]);
      tempDist<-apply(xv,1,function(x) sum(x^2));
      d<-cbind(d,tempDist);
    }
    d<-sqrt(d)
    fitness=sum(diag(t(fm)%*%d));
    particle[[k]]=list(pbest=pbest,muX=muX,muV=muV,fitness=fitness,curFitness=fitness);
  }
  return(particle)  #return particle;
}

# mixDist<-function(dt,){
#   
# }

denParInit<-function(dataset,param) {  #krt is number of groups
  X<-dataset;
  cn<-param$cn;  #cluster number
  N<-nrow(X)     #number of data
  n<-ncol(X)     #number of attribute
  p<-param$pn; #number of particle
  m<-param$m
  krt<-param$krt;
  dc<-param$dc
  distMx<-dist(X, method = "euclidean", p = 2)
  result<-initCentroid(as.matrix(distMx),dc);
  clustRt<-result$clustRt;
  rho<-result$rho;
  delta<-result$delta;
  gbest<-1000000000;
  gfitness<-0;
  X1 <-c(rep(1,N));
  v<-list();
  particles<-list();
  temp1<-as.integer(length(clustRt)/cn);   #number of groups that centroids could be div
  if (temp1<1){
    t_v<-c()
    t_v<-rbind(t_v,X[clustRt,]);
    tt1<-1:N;
    tt2<-sample(tt1[-1*clustRt],cn-length(clustRt));
    t_v<-rbind(t_v,X[tt2,]);
    v[[1]]<-t_v
  } else if (temp1>=1 &temp1<p){
    temp2<-min(temp1,krt);
    for (jj in 1:temp2){
      id<-((jj-1)*cn+1):(jj*cn)
      v[[jj]]<-rbind(X[clustRt[id],])
    }
  } else {
    for (jj in 1:krt){
      id<-((jj-1)*cn+1):(jj*cn)
      v[[jj]]<-rbind(X[clustRt[id],])
    }
  }
  if (length(v)<p){
    t_len<-length(v);
    diffl<-p-t_len;
    for (jj in 1:diffl){
      tt1<-1:N;
      tt2<-sample(tt1[-1*clustRt],cn);
      v[[t_len+jj]]<-rbind(X[tt2,]);
    }
  }
  
  for (k in 1:p){
    temp_v<-v[[k]];
    d<-c();
    for (j in 1:cn){
      xv = X - tcrossprod(X1,temp_v[j,]);
      tempDist<-apply(xv,1,function(x) sum(x^2));
      d<-cbind(d,tempDist);
    } 
    d0<-sqrt(d);
    d = (d+1e-10)^(-1/(m-1));  #L2 distance
    #browser()
    pbest <- d/tcrossprod(apply(d,1,sum),rep(1,cn));
    f0<-pbest
    J<-sum(diag(t(f0^m)%*%d0))
    
    #browser()
    particles[[k]]=list(pbest=pbest,muX=pbest,muV=pbest,fitness=J,curFitness=J);
  }
  return(particles)
}
fpsofcm<-function (dataset,param,particles){
  #browser();
  X<-dataset;
  wmax<-0.9;
  wmin<-0.1
  #w<-0.1;
  c1<-2.0;
  c2<-2.0;
  cn<-param$cn;
  r1<-runif(1,0,1);
  r2<-runif(1,0,1);
  m = param$m;
  e = param$e;
  N<-nrow(X)
  n<-ncol(X)
  
  gbest<-particles[[1]]$pbest;
  oldGbest<-gbest;
  loopC<-0
  while(loopC<=500){
    fpsoRes<-fpsoModule(X,particles,param)     #fpsoRes=list(particle,iterations)
    #browser()
    fcmRes<-fcmModule(fpsoRes$particles,X,param); #fcmRes=list(particle,iterations)
    loopC<-loopC+fpsoRes$iters+fcmRes$iters
    gbest<-fcmRes$gbest
    if (max(abs(gbest-oldGbest))<e){
      break;
    }
    oldGbest<-gbest;
  }
  fm <- gbest^m;
  sumf<-apply(fm,2,sum);
  v<-crossprod(fm,X)/(sumf%*%t(rep(1,n)));
  d<-c()
  X1<-rep(1,N);
  for (j in 1 : cn){
    xv = X - tcrossprod(X1,v[j,]);
    tempDist<-apply(xv,1,function(x) sum(x^2));
    d<-cbind(d,tempDist);
  }
  d<-sqrt(d);
  J=sum(diag(t(fm)%*%d));
  result<-list(f=gbest,d=d,v=v,iters = loopC,cost = J)
}

operationPerParticle<-function(particle,dt,param){
  #browser()
  N<-nrow(dt)
  n<-ncol(dt)
  m<-param$m
  cn<-param$cn
  newParticle<-particle;
  fm <- particle$muX^m;
  sumf<-apply(fm,2,sum);
  v<-crossprod(fm,dt)/(sumf%*%t(rep(1,n)));
  X1<-rep(1,N);
  d<-c()
  for (j in 1 : cn){
    xv = dt - tcrossprod(X1,v[j,]);
    tempDist<-apply(xv,1,function(x) sum(x^2));
    d<-cbind(d,tempDist);
  }
  distout<-sqrt(d);
  f0<-particle$muX;
  tempFitness<-sum(diag(t(f0^m)%*%d));
  if (particle$fitness>tempFitness){
    newParticle$fitness<-tempFitness;
    newParticle$pbest<-f0;
  }
  return(newParticle)
}
updatePerParticle<-function(particle,w,gbest){
  N<-nrow(gbest)
  cn<-ncol(gbest)
  set.seed(100)
  r1<-runif(1,0,1)
  r2<-runif(1,0,1)
  c1<-2.0
  c2<-2.0
  newParticle<-particle
  #browser();
  newParticle$muV=w*particle$muV+c1*r1*(particle$pbest-particle$muX)+c2*r2*(gbest-particle$muX);
  newParticle$muX=particle$muX+newParticle$muV;
  newParticle$muX[newParticle$muX<0]=0;
  for (s in 1:N){
    if (sum(newParticle$muX[s,])==0){
      set.seed(1000)
      newParticle$muX[s,]=runif(cn,min=0,max=1)#unifrnd(0,1,1,c);
    }
  }
  #browser()
  newParticle$muX=newParticle$muX/(apply(newParticle$muX,1,sum)%*%t(rep(1,cn)));
  return(newParticle)
}

fcmPerParticle<-function(dt,particle,fcmParam){
  newParticle<-particle
  N<-nrow(dt)
  n<-ncol(dt)
  cn<-fcmParam$cn
  m<-fcmParam$m
  newParticle<-particle;
  fm <- particle$muX^(fcmParam$m);
  sumf<-apply(fm,2,sum);
  v<-crossprod(fm,dt)/(sumf%*%t(rep(1,n)));
  X1<-rep(1,N);
  d<-c()
  for (j in 1 : cn){
    xv = dt - tcrossprod(X1,v[j,]);
    tempDist<-apply(xv,1,function(x) sum(x^2));
    d<-cbind(d,tempDist);
  }
  distout<-sqrt(d);
  f0<-t(apply(distout,1,function(y,m) y^(-2/(m-1))/sum(y^(-2/(m-1))),m=fcmParam$m))
  newParticle$muX<-f0;
  tempFitness<-sum(diag(t(f0^m)%*%d));
  newParticle$curFitness<-tempFitness
  if (particle$fitness>tempFitness){
    newParticle$fitness<-tempFitness;
    newParticle$pbest<-f0;
  }
  return(newParticle)
}

fcmModule<-function(particles,dataset,fcmParam){
  iters<-0;
  pn<-length(particles)
  gbest<-particles[[1]]$pbest;     #Just for initialization
  tempFitness<-particles[[1]]$fitness;
  while (iters<6){
    iters<-iters+1;
    newParticles<-lapply(particles,fcmPerParticle,dt=dataset,fcmParam=fcmParam)
    for (k in 2:length(newParticles)){
      if(newParticles[[k]]$fitness<tempFitness){
        gbest<-newParticles[[k]]$pbest;
        tempFitness<-newParticles[[k]]$fitness
      }
    }
  }
  return(list(particles=newParticles,iters=iters,gbest=gbest))
}

fpsoModule<-function (dataset,particles,param){
  #browser();
  X<-dataset;
  wmax<-0.9;
  wmin<-0.1
  c1<-2.0;
  c2<-2.0;
  cn<-param$cn;
  r1<-runif(1,0,1);
  r2<-runif(1,0,1);
  m = param$m;
  e = param$e;
  N<-nrow(X)
  n<-ncol(X)
  
  gFitness<-particles[[1]]$fitness;
  oldGFitness<-gFitness
  gbest<-particles[[1]]$pbest;
  oldGbest<-gbest;
  iter<-0;
  pn<-length(particles);
  maxLoops<-95;
  #maxConvergeLoops<-200;
  #**********************Initialize Particle*********************
  # initRes<-randomInit(dataset,param,pn);
  # gbest<-initRes$gbest;
  # particle<-initRes$particle;
  # p<-length(particle)
  #*****************end*****************************
  loopC<-0;
  loopC1<-0;
  distout<-c()
  while (loopC<maxLoops){
    loopC<-loopC+1;
    iter<-iter+1;
    #browser()
    newParticles<-lapply(particles,operationPerParticle,dt=X,param=param)
    for (k in 1:length(newParticles)){
      if(newParticles[[k]]$fitness<gFitness){
        gbest<-newParticles[[k]]$pbest;
        gFitness<-newParticles[[k]]$fitness
      }
    }
    w<-wmax-loopC*(wmax-wmin)/maxLoops;
    newParticles<-lapply(newParticles,updatePerParticle,w=w,gbest=gbest)
    
    if (max(abs(gFitness-oldGFitness))<=1e-5){ 
        break
      }
    
    oldGFitness<-gFitness;
  }
  return(list(particles=newParticles,gbest=gbest,iters = iter))
}

fcm<-function(dataset,fcmParam,particles){  #only 1 particle
  iters<-0;
  mu<-particles[[1]]$muX;     #Just for initialization
  oldMu<-mu+10;
  e<-fcmParam$e
  particle<-particles[[1]];
  while (max(abs(mu-oldMu))>e){
    iters<-iters+1;
    particle<-fcmPerParticle(dataset,particle,fcmParam)
    oldMu<-mu;
    mu<-particle$muX
  }
  return(list(f=particle$muX,d=NULL,v=NULL,iters = iters,cost = particle$curFitness))
}

fpso<-function (dataset,param,particles){
  #browser();
  X<-dataset;
  wmax<-0.9;
  wmin<-0.1
  c1<-2.0;
  c2<-2.0;
  cn<-param$cn;
  r1<-runif(1,0,1);
  r2<-runif(1,0,1);
  m = param$m;
  e = param$e;
  N<-nrow(X)
  n<-ncol(X)
  #browser()
  gFitness<-particles[[1]]$fitness;
  oldGFitness<-gFitness
  gbest<-particles[[1]]$pbest;
  oldGbest<-gbest;
  iter<-0;
  pn<-length(particles);
  maxLoops<-1000;
  convergeLoops<-200
  #browser()
  loopC<-0;
  loop1<-0;
  distout<-c()
  while (loopC<maxLoops){
    loopC<-loopC+1;
    iter<-iter+1;
    newParticles<-lapply(particles,operationPerParticle,dt=X,param=param)
    tempMinFitness<-newParticles[[1]]$curFitness;
    for (k in 1:length(newParticles)){
      tempMinFitness<-min(tempMinFitness,newParticles[[k]]$curFitness);
      if(newParticles[[k]]$fitness<gFitness){
        gbest<-newParticles[[k]]$pbest;
        gFitness<-newParticles[[k]]$fitness
      }
    }
    w<-wmax-loop1*(wmax-wmin)/convergeLoops;
    newParticles<-lapply(newParticles,updatePerParticle,w=w,gbest=gbest)
    if (max(abs(gbest-oldGbest))<=1e-5){
      loop1<-loop1+1;
      if (loop1==convergeLoops){
        break  
      }
    } else{
      loop1<-0;
    }
    oldGbest<-gbest;
  }
  return(list(f=gbest,d=NULL,v=NULL,iters = iter,cost = gFitness))
}



initCentroid<-function(distMx,dc){
  ND<-nrow(distMx)
  dThresh<-sort(as.vector(distMx),decreasing=TRUE)[as.integer(dc*0.01*ND^2)]
  
  rho<-rep(0,ND)
  delta<-rep(0,ND)
  nneigh<-rep(0,ND);
  isBorder<-rep(1,ND);    # the border points are not directed neighbors of any points.    
  initCenterAt<-rep(0,ND);
  finalCenterAt<-c();
  for (i in 1:(ND-1)){
    for (j in (i+1):ND){
      rho[i]=rho[i]+exp(-(distMx[i,j]/dThresh)*(distMx[i,j]/dThresh))+1e-8;
      rho[j]=rho[j]+exp(-(distMx[i,j]/dThresh)*(distMx[i,j]/dThresh))+1e-8; 
    }
  } 
  maxd<-max(distMx);
  s_d_rho<-sort(rho,decreasing=T,index.return=T)
  rho_sorted<-s_d_rho$x
  ordrho<-s_d_rho$ix
  #browser()
  #[rho_sorted,ordrho]=sort(rho,'descend');
  delta[ordrho[1]]<--1;
  #delta(ordrho(1))=-1.;
  nneigh[ordrho[1]]=ordrho[1]
  #nneigh(ordrho(1))=ordrho(1);      # 0
  initCenterAt[ordrho[1]]=ordrho[1];
  #initCenterAt(ordrho(1))=ordrho(1);
  potClustRt<-c(ordrho[1]);
  #potClustRt=[ordrho(1)];
  
  for (ii in 2:ND){
    delta[ordrho[ii]]<-maxd;
    for (jj in 1:(ii-1)){
      if(distMx[ordrho[ii],ordrho[jj]]<delta[ordrho[ii]]){
        delta[ordrho[ii]]=distMx[ordrho[ii],ordrho[jj]]+0.000000000000001;
      }
      nneigh[ordrho[ii]]=ordrho[jj]; 
    }
      
    isBorder[nneigh[ordrho[ii]]]=0;
    if (distMx[nneigh[ordrho[ii]],ordrho[ii]]>dc){
      potClustRt=c(potClustRt,ordrho[ii]); 
      initCenterAt[ordrho[ii]]=ordrho[ii];
    }
    else {
      initCenterAt[ordrho[ii]]=initCenterAt[nneigh[ordrho[ii]]];
    }
      
  }
  delta[ordrho[1]]=max(delta);
  plot(delta,rho,'o');
  #browser()
   # %*******************For BSN 2015******************************
  gama=(delta/max(delta))*(rho/max(rho));
  sorted_clustCent<-sort(gama[potClustRt],decreasing=T,index.return=T);
  sorted_clustCent_gama<-sorted_clustCent$x
  ordinx<-sorted_clustCent$ix
  fClustRt<-potClustRt[ordinx];
  #finalCenterAt=initCenterAt;
#*********************************The end********************************
  result<-list(rho=rho,delta=delta,clustRt=fClustRt,nneigh=nneigh);
}

regularize<-function(X){
  zz<-apply(X,2,function(m) (m - min(m))/(max(m)-min(m)))
  return(zz) 
}

