library(base)
library(stats)
library(gtools)  # for permu
library(pdfCluster) #for adjusted rand index;
fcm<- function (dataset,param)  #param is a list
{
  X<-dataset;
  initValidityXB<-10000000;
  cn<-param$cn;#cn is cluster number
  m<-param$m;
  e<-param$e;
  N<-nrow(X)
  n<-ncol(X)
  #****************Initialization randomly select centroids********
  X1 <-c(rep(1,N));
  set.seed(1000);
  centroidIndex=sample(1:N,cn);   #Generate N unrepeatable integers, then we select c of them.
  v<-c();
  d<-c();
  #browser();
  for (j in 1:cn){
    tempRow<-X[centroidIndex[j],];
    v<-rbind(v,tempRow);
    xv = X - tcrossprod(X1,v[j,])#X1%*%t(v[j,]);
    tempDist<-apply(xv,1,function(x) sum(x^2));
    d<-cbind(d,tempDist);
  }
  d0<-d;
  #browser()
  d <- (d+1e-10)^(-1/(m-1));
  f0 <- d/tcrossprod(apply(d,1,sum),rep(1,cn));
  iter <- 0;                       # iteration counter
  initialJ<-10000000000000000;
  
  J0<-sum(diag(t(f0^m)%*%d0));
  distout<-c();
  J<-c();
  while (abs(J0-initialJ)>e){
    iter = iter + 1;
    initialJ=J0;
    # Calculate centers
    fm <-f0^m;
    sumf<-apply(fm,2,sum);
    #sumf = sum(fm);
    v<-crossprod(fm,X)/(sumf%*%t(rep(1,n)));
    #v = (fm'*X)./(sumf'*ones(1,n)); 
    d<-c();
    for (j in 1:cn){
      xv = X - tcrossprod(X1,v[j,]);
      #xv = X - X1*v[j,];
      tempDist<-apply(xv,1,function(x) sum(x^2));
      d<-cbind(d,tempDist);
      #d[,j] = (xv*eye(n).*xv)*(w(j,:).^2)';
    } 
    distout<-sqrt(d);
    d0<-d;  
    J[iter]<-sum(diag(t(f0^m)%*%d0));
    d = (d+1e-10)^(-1/(m-1));  #L2 distance
    f0 <- d/tcrossprod(apply(d,1,sum),rep(1,cn));
    #f0 = (d/ (sum(d,2)*ones(1,c)));      #method I
    J0=J[iter];
  }
  #tEnd=toc(tStart);
  #results
  results<-list(f=f0,d=distout,v=v,iter=iter,cost=J);
}
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
    mydata=read.table(dataPath,sep=sepStr)    
  }
  else if (fileType=="csv")
  {
    mydata=read.csv(dataPath)
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
randomInit<-function(dataset,param){  #random initialization;
  X<-dataset;
  cn<-param$cn;  #cluster number
  N<-nrow(X)     #number of data
  n<-ncol(X)     #number ofattribute
  p<-10;         #number of particle
  gfitness<-0;
  particle<-list();
  for (k in 1:p){
    browser
    set.seed(k);
    pbest<-matrix(runif(N*cn,min=0,max=1),nrow=N,byrow=T)
    muX<-matrix(runif(N*cn,min=0,max=1),nrow=N,byrow=T)
    muV<-matrix(runif(N*cn,min=0,max=1),nrow=N,byrow=T)
    #browser()
    pbest<-pbest/tcrossprod(apply(pbest,1,sum),rep(1,cn));
    muX<-muX/tcrossprod(apply(muX,1,sum),rep(1,cn));
    muV<-muV/tcrossprod(apply(muV,1,sum),rep(1,cn));
    particle[[k]]=list(pbest=pbest,muX=muX,muV=muV,center=c(),fitness=0);
    
  }
  gbest<-particle[[k]]$pbest;
  result<-list(particle=particle,gbest=gbest)  #return particle and gbest;
}

denInit<-function(dataset,param,dc,krt) {  #krt is number of groups
  X<-dataset;
  cn<-param$cn;  #cluster number
  N<-nrow(X)     #number of data
  n<-ncol(X)     #number ofattribute
  p<-10; #number of particle
  distMx<-dist(X, method = "euclidean", p = 2)
  result<-initCentroid(as.matrix(distMx),dc);
  clustRt<-result$clustRt;
  rho<-result$rho;
  delta<-result$delta;
  gbest<-1000000000;
  gfitness<-0;
  X1 <-c(rep(1,N));
  v<-list();
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
    for (j in 1:cn){
      xv = X - tcrossprod(X1,temp_v[j,]);
      tempDist<-apply(xv,1,function(x) sum(x^2));
      d<-cbind(d,tempDist);
    } 
    d = (d+1e-10)^(-1/(m-1));  #L2 distance
    d0<-d;  
    J<-sum(diag(t(f0^m)%*%d0))
    fitness<-1/J
    pbest <- d/tcrossprod(apply(d,1,sum),rep(1,cn));
    particle[[k]]=list(pbest=pbest,muX=pbest,muV=pbest,center=c(),fitness=fitness);
    if (gfitness<fitness){
      gfitness<-fitness;
      gbest<-pbest;
    }
    result<-list(particle=particle,gbest=gbest)
    
  }
  #result<-list(rho=rho,delta=delta,clustRt=fClustRt,nneigh=nneigh);
}
fpsofcm<-function (dataset,param){
  browser();
  X<-dataset;
  wmax<-0.9;
  wmin<-0.1
  #w<-0.1;
  c1<-2.0;
  c2<-2.0;
  cn<-param$cn;
  #r1<-0.2;
  #r2<-0.2;
  m = param$m;
  e = param$e;
  N<-nrow(X)
  n<-ncol(X)
  #**********************Initialize Particle*********************
  initRes<-randomInit(dataset,param);
  browser()
  gbest<-initRes$gbest;
  particle<-initRes$particle;
  #*****************end*****************************
  old_gbest<-matrix(rep(1,N*cn),nrow=N);
  loopC<-0;
  gfitness<-10000000;
  p<-length(particle)
  J_old<-rep(0,p);
  J_new<-rep(0,p);
  while (loopC<500){
    fpsofitness<-gfitness;
    fpsoloop<-0;
    while(fpsoloop<95){
      fpsoloop<-fpsoloop+1;
      for (k in 1:p){
        fm <- particle[[k]]$muX^m;
        sumf<-apply(fm,2,sum);
        particle[[k]]$center<-crossprod(fm,X)/(sumf%*%t(rep(1,n)));
        X1<-rep(1,N);
        v<-particle[[k]]$center;
        d<-c();
        for (j in 1:cn){
          xv = X - tcrossprod(X1,v[j,]);
          tempDist<-apply(xv,1,function(x) sum(x^2));
          d<-cbind(d,tempDist);
        }
        f0=particle[[k]]$muX;
        tempFitness<-sum(diag(t(f0^m)%*%d));
        J_new[k]<-tempFitness;
        if (particle[[k]]$fitness<tempFitness){
          particle[[k]]$fitness<-tempFitness;
          particle[[k]]$pbest<-f0;
        }
        if (particle[[k]]$fitness<gfitness){
          gfitness<-particle[[k]]$fitness;
          gbest<-particle[[k]]$pbest;
        }
      }
      
      #*****************************************************
      #browser()
      w<-wmax-fpsoloop*(wmax-wmin)/95;
      for (k in 1:p){
        set.seed(k)
        r1<-runif(1)
        r2<-runif(1)
        #browser();
        particle[[k]]$muV=w*particle[[k]]$muV+c1*r1*(particle[[k]]$pbest-particle[[k]]$muX)+c2*r2*(gbest-particle[[k]]$muX);
        particle[[k]]$muX=particle[[k]]$muX+particle[[k]]$muV;
        particle[[k]]$muX[particle[[k]]$muX<0]=0;
        for (s in 1:N){
          if (sum(particle[[k]]$muX[s,])==0){
            set.seed(1000)
            particle[[k]]$muX[s,]=runif(cn,min=0,max=1)#unifrnd(0,1,1,c);
          }
        }
        #browser()
        particle[[k]]$muX=particle[[k]]$muX/(apply(particle[[k]]$muX,1,sum)%*%t(rep(1,cn)));
      }
      #browser()
      if (min(abs(J_old-J_new))<e){
        break;
      }
      J_old<-J_new;
      #browser();
    }
    browser();
    fcmfitness<-gfitness;
    fcmloop<-0;
    browser()
    while (fcmloop<6){
      fcmloop<-fcmloop+1;
      for (k in 1:p){
        fm <- particle[[k]]$muX^m;
        sumf<-apply(fm,2,sum);
        particle[[k]]$center<-crossprod(fm,X)/(sumf%*%t(rep(1,n)));
        X1<-rep(1,N);
        v<-particle[[k]]$center
        d<-c()
        for (j in 1 : cn){
          xv = X - tcrossprod(X1,v[j,]);
          tempDist<-apply(xv,1,function(x) sum(x^2));
          d<-cbind(d,tempDist);
        }
        distout<-sqrt(d);
        f0<-particle[[k]]$muX;
        tempFitness<-sum(diag(t(f0^m)%*%d));
        J_new[k]<-tempFitness;
        #tempFitness=1/sum(sum(f0.^2.*d));
        if (particle[[k]]$fitness<tempFitness){
          particle[[k]]$fitness<-tempFitness;
          particle[[k]]$pbest<-f0;
        }
        if (particle[[k]]$fitness<gfitness){
          gfitness=particle[[k]]$fitness;
          gbest=particle[[k]]$pbest;
        }
      }
      if (abs(min(J_old-J_new))<e){
        break;
      }
      J_old<-J_new;
    }
    loopC<-fcmloop+fpsoloop;
    if (max(abs(gbest-old_gbest))<e){
      break;
    }
    old_gbest<-gbest;
  }
  browser()
  fm <- gbest^m;
  sumf<-apply(fm,2,sum);
  v<-crossprod(fm,X)/(sumf%*%t(rep(1,n)));
  J=gfitness;
  result<-list(f=gbest,d=distout,v=v,iter = loopC,cost = J)
}

fpso<-function (dataset,param){
  #browser();
  X<-dataset;
  w<-0.1;
  c1<-2.0;
  c2<-2.0;
  cn<-param$cn;
  r1<-0.2;
  r2<-0.2;
  m = param$m;
  e = param$e;
  N<-nrow(X)
  n<-ncol(X)
  gfitness<-100000000;
  iter<-0;
  #**********************Initialize Particle*********************
  initRes<-randomInit(dataset,param);
  gbest<-initRes$gbest;
  particle<-initRes$particle;
  p<-length(particle)
  #*****************end*****************************
  old_gbest<-matrix(rep(1,N*cn),nrow=N);
  loopC<-0;
  distout<-c()
  while (loopC<1000){
    iter<-iter+1;
    #*******************particle update*****************
    for (k in 1:p){
      fm <- particle[[k]]$muX^m;
      sumf<-apply(fm,2,sum);
      particle[[k]]$center<-crossprod(fm,X)/(sumf%*%t(rep(1,n)));
      X1<-rep(1,N);
      v<-particle[[k]]$center
      d<-c()
      for (j in 1:cn){
        xv = X - tcrossprod(X1,v[j,]);
        tempDist<-apply(xv,1,function(x) sum(x^2));
        d<-cbind(d,tempDist);
      }
      f0=particle[[k]]$muX;
      tempFitness<-sum(diag(t(f0^m)%*%d));
      if (particle[[k]]$fitness<tempFitness){
        particle[[k]]$fitness<-tempFitness;
        particle[[k]]$pbest<-f0;
      }
      if (particle[[k]]$fitness<gfitness){
        gfitness<-particle[[k]]$fitness;
        gbest<-particle[[k]]$pbest;
        distout<-d
      }
    }
    #*****************************************************
    for (k in 1:p){
      set.seed(k)
      r1<-runif(1);
      r2<-runif(1);
      particle[[k]]$muV=w*particle[[k]]$muV+c1*r1*(particle[[k]]$pbest-particle[[k]]$muX)+c2*r2*(gbest-particle[[k]]$muX);
      particle[[k]]$muX=particle[[k]]$muX+particle[[k]]$muV;
      particle[[k]]$muX[particle[[k]]$muX<0]=0;
      for (s in 1:N){
        if (sum(particle[[k]]$muX[s,])==0){
          set.seed(1000)
          particle[[k]]$muX[s,]=runif(cn,min=0,max=1)#unifrnd(0,1,1,c);
        }
      }
      #browser()
      particle[[k]]$muX=particle[[k]]$muX/(apply(particle[[k]]$muX,1,sum)%*%t(rep(1,cn)));
    }
    #*****************************************************
    if (max(abs(gbest-old_gbest))>1e-3){ 
      loopC<-0;
    } else {
      loopC<-loopC+1;
    }
    old_gbest<-gbest;
  }
  browser()
  fm <- gbest^m;
  sumf<-apply(fm,2,sum);
  v<-crossprod(fm,X)/(sumf%*%t(rep(1,n)));
  result<-list(f=gbest,d=distout,v=v,iter = iter,cost = gfitness)
}

initCentroid<-function(distMx,dc){
  ND<-nrow(distMx)
  rho<-rep(0,ND)
  delta<-rep(0,ND)
  nneigh<-rep(0,ND);
  isBorder<-rep(1,ND);    # the border points are not directed neighbors of any points.    
  initCenterAt<-rep(0,ND);
  finalCenterAt<-c();
  for (i in 1:(ND-1)){
    for (j in (i+1):ND){
      rho[i]=rho[i]+exp(-(distMx[i,j]/dc)*(distMx[i,j]/dc))+1e-8;
      rho[j]=rho[j]+exp(-(distMx[i,j]/dc)*(distMx[i,j]/dc))+1e-8; 
    }
  } 
  maxd<-max(distMx);
  s_d_rho<-sort(rho,decreasing=T,index.return=T)
  rho_sorted<-s_d_rho$x
  ordrho<-s_d_rho$ix
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
        delta[ordrho[ii]]=dist[ordrho[ii],ordrho[jj]]+0.000000000000001;
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
   # %*******************For BSN 2015******************************
  gama=(delta/max(delta))*(rho/max(rho));
  sorted_clustCent<-sort(gama(potClustRt),decreasing=T,index.return=T);
  sorted_clustCent_gama<-sorted_clustCent$x
  ordinx<-sorted_clustCent$ix
  fClustRt<-potClustRt(ordinx);
  #finalCenterAt=initCenterAt;
#*********************************The end********************************
  result<-list(rho=rho,delta=delta,clustRt=fClustRt,nneigh=nneigh);
}

normalit<-function(m){
  return (m - min(m))/(max(m)-min(m))
}

