
jin2=read.csv( "C:/Users/YDC/Desktop/上传数据/growdata_mic0.csv",header = T,row.names = 1)
jin3=read.csv( "C:/Users/YDC/Desktop/上传数据/growdata_mic2.csv",header = T,row.names = 1)
jin4=read.csv( "C:/Users/YDC/Desktop/上传数据/growdata_mic4.csv",header = T,row.names = 1)
jin5=read.csv( "C:/Users/YDC/Desktop/上传数据/growdata_mic6.csv",header = T,row.names = 1)
SNP=read.csv( "C:/Users/YDC/Desktop/上传数据/snpdata.csv",header = T,row.names = 1)
k8str <- read.csv("C:/Users/YDC/Desktop/上传数据/FASTSTRUCTURE.csv", header = F)

t <- jin2[,1]

#_________________________________Correction of population structure______________________________
get_strjin <- function(jin,si){
  
  
  BB <- apply(si,1,function(x){which(x==max(x))})
  
  BBt <- table(BB)
  BBn <- as.numeric(names(BBt))
  rjin2 <- c()
  for (i in 1:nrow(jin)) {
    
    
    pdata <-   jin[i,]
    subpop <- list()
    almea <- c()
    for( j in 1:length(BBt)){
      
      index1 <- which(BB==BBn[j])
      subpop[[j]] <- pdata[index1]
      almea[j] <- mean(as.numeric(pdata[index1]))
      
      
    }
    for (j in 1:length(BBt)) {
      chazhi=almea[j]-mean(almea)
      pdata[which(BB==BBn[j])] <- pdata[which(BB==BBn[j])]-chazhi
      
      
      
      
      
    }
    rjin2 <- rbind(rjin2,pdata)
    
    
  }
  
  
  rjin2 <- cbind(jin2[,1],rjin2)
  
  return(rjin2)}
rjin2 <- get_strjin(jin2[,-1],k8str)
rjin3 <- get_strjin(jin3[,-1],k8str)
rjin4 <- get_strjin(jin4[,-1],k8str)
rjin5 <- get_strjin(jin5[,-1],k8str)
#__________________________________Composite functional mapping_________________________________________________
library(mvtnorm)

library(pbapply)
get_miu3 =function(B,t){B[1]/(1+exp((4*B[2]*(B[3]-t)/B[1])+2))}
SAD1_get_matrix = function(par, times = t, options=list()) {   
  n <- ifelse (is.vector(times), length(times), NCOL(times) )   
  phi<- par[1]   
  v2 <- par[2]   
  tmp <- (1-phi^2)   
  sigma <- array(1, dim=c(n,n))   
  for(i in 1:n)   
  {     
    sigma[i,i:n] <- phi^( c(i:n) - i ) * (1-phi^(2*i))/tmp     
    sigma[i:n,i] <- sigma[i,i:n]   
  }   
  sigma <- sigma * abs(v2)   
  return(sigma); 
}
get_u = function(A,t){get_miu3(A[1:3],t)-get_miu3(A[4:6],t)}
get_sig= function(par,times = t, options=list()){SAD1_get_matrix(par = c(par[1],par[2]),times = t) + SAD1_get_matrix(par = c(par[3],par[4]),times = t) 
}
get_initial_par <- function(pheno,t){
  mean0 <- apply(pheno[,],2,mean)  
  c(max(mean0),
    max((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])),
    t[which.max(((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])))]-mean0[which.max(((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])))]/max((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])))
}
H0 = function(yt,t,par){
  miu=get_u(par[1:6],t)
  sigma=get_sig(par = par[7:10],times=t)
  L0 = c()
  L0 = sum(dmvnorm(t(yt),miu,sigma,log = TRUE))
  return(-L0)
}
H1 = function(yt,t,m0,m1,par){
  
  
  
  
  
  p1.0 <- yt[,c(m0)]
  p1.1 <- yt[,c(m1)]
  
  miu0 = get_u(par[1:6],t)
  miu1 = get_u(par[7:12],t)
  sigma = get_sig(par = par[13:16],times=t)
  
  L1.1 = sum(dmvnorm(t(p1.0),miu0,sigma,log = TRUE))
  L1.0=sum(dmvnorm(t(p1.1),miu1,sigma,log = TRUE))
  L1= L1.1+L1.0
  return(-L1)
}
parl0 <- optim(par = c(get_initial_par(t(rjin2[,-1]),jin2[,1]),get_initial_par(t(rjin5[,-1]),jin2[,1]),0.8,0.01,0.8,0.01),H0,yt = (rjin2[,-1]-rjin5[,-1]),t = jin2[,1])

parl0 <- optim(par = parl0$par,H0,yt = (rjin2[,-1]-rjin5[,-1]), t = jin2[,1])
parl0 <- optim(par = parl0$par,H0,yt = (rjin2[,-1]-rjin5[,-1]), t = jin2[,1])
parl0 <- optim(par = parl0$par,H0,yt = (rjin2[,-1]-rjin5[,-1]), t = jin2[,1])
parl0 <- optim(par = parl0$par,H0,yt = (rjin2[,-1]-rjin5[,-1]), t = jin2[,1])
parl0 <- optim(par = parl0$par,H0,yt = (rjin2[,-1]-rjin5[,-1]), t = jin2[,1])
parl0 <- optim(par = parl0$par,H0,yt = (rjin2[,-1]-rjin5[,-1]), t = jin2[,1])
parl0 <- optim(par = parl0$par,H0,yt = (rjin2[,-1]-rjin5[,-1]), t = jin2[,1])
parl0 <- optim(par = parl0$par,H0,yt = (rjin2[,-1]-rjin5[,-1]), t = jin2[,1])
parl0 <- optim(par = parl0$par,H0,yt = (rjin2[,-1]-rjin5[,-1]), t = jin2[,1])

optim_diff <- function(pheno_ck0,pheno_salt0,pheno_ck1,pheno_salt1,pheno_diff,t,m0,m1,xx){
  
  itime <- 100
  itimes <- 1
  par0 <-as.numeric(c(parl0$par[1:6],parl0$par))*xx
  repeat{
    a <- optim(par = par0 ,H1,yt = pheno_diff, t = t,m0=m0,m1=m1)
    
    b <- optim(a$par,H1,yt = pheno_diff, t = t,m0=m0,m1=m1)
    
    
    # cat("Logistic_diff",itimes,b$value,'\n')
    
    itimes <- itimes + 1
    
    if(all( abs(a$value-b$value) < 1 )||itimes == itime){ #itimesԽ?߾???Խ??
      break
    }else{
      par0 <- b$par
    }
  }
  b
  
}

get_lr <- function(SNP,control,stress,xx){
  vnp=SNP
  m0 <- which( vnp == 0)
  m1 <- which( vnp == 1)
  pheno_ck <- control[,-1]
  pheno_salt <- stress[,-1]
  pheno_ck0 <- pheno_ck[,m0]
  pheno_ck1 <- pheno_ck[,m1]
  pheno_salt0 <- pheno_salt[,m0]
  pheno_salt1 <- pheno_salt[,m0]
  
  parl1 <- optim_diff(pheno_ck0,pheno_salt0,pheno_ck1,pheno_salt1,(control[,-1]-stress[,-1]),t,m0,m1,xx)
  LR_plastic3 <- c(-2*(-parl0$value+parl1$value),parl1$par)
}

LR_plastic <- pbapply(SNP[,-1],1,get_lr,control=rjin2,stress=rjin5,xx=1)
lr <- c()
for (i in 1:1000) {
  sa <- sample(1:99,99)
 
    rjin2 <- rjin2[,sa]
rjin5 <- rjin5[,sa]
  
LR_plasticsa <- pbapply(SNP[,-1],1,get_lr,control=rjin2,stress=rjin5,xx=1)
lr[i] <- max(LR_plasticsa[1,])
  
}
Threshold_value <- sort(lr,decreasing = T)[50]
rjin2 <- get_strjin(jin2[,-1],k8str)
rjin3 <- get_strjin(jin3[,-1],k8str)
rjin4 <- get_strjin(jin4[,-1],k8str)
rjin5 <- get_strjin(jin5[,-1],k8str)

#____________________________Calculate genetic effect____________________________________________________________
get_VG <- function(jin2,jin5,marker_data,t,LR){
  
  diff_vg <- c() 
  
  
  
  
  
  for (a in 1:dim(marker_data)[1]) {
    
    AA <- as.numeric(which(marker_data[a,]==1))
    aa <- as.numeric(which(marker_data[a,]==0))
    
    NAA <- length(AA)
    Naa <- length(aa)
    
    p1 <- (NAA*2)/((NAA+Naa)*2) #A????Ƶ??
    p0 <- (Naa*2)/((NAA+Naa)*2) #a????Ƶ??
    
    mean_AA <- get_u(LR[2:7,a],t)
    mean_aa <- get_u(LR[8:13,a],t)
    AE <- (mean_AA - mean_aa)/2 
    
    Vg <- 2*p1*p0*(AE^2)  
    diff_vg <- rbind(diff_vg,Vg)
    cat(a,"finished","\n")
    
  }
  
  
  
  return(sqrt(diff_vg)) 
}

varh2tx <- get_VG(rjin2,rjin5,SNP[,-1],c(1:max(t)),LR_plastic)


#______________________________________Functional cluster____________________________________
requiredPackages = c("mvtnorm","reshape2","pbapply",'parallel','orthopolynom')
for(packages in requiredPackages){
  if(!require(packages,character.only = TRUE)) install.packages(packages)
  require(packages,character.only = TRUE)
}

age <- c(1:max(t))
get_init_par <- function(data,k){
  
  f3  <- function(y) {
    
    lop <- legendre.polynomials(n=13,normalized = F)
    lop_matrix <- as.matrix(as.data.frame(polynomial.values(polynomials = lop,x=scaleX(1:48,u=-1,v=1))))
    colnames(lop_matrix) <- paste0("L",0:13)
    lop_fit <- lm(y~lop_matrix[,2:14]) 
    return(as.numeric(lop_fit$coefficients))
  }
  
  
  
  
  
  #????k-means???ó?ʼ????--------------------------------------
  
  init_cluster <- kmeans(data,centers = k,iter.max = 100)
  pro <- table(init_cluster$cluster)/nrow(data) #????ģ?͸??ʸ???
  
  cuM <- init_cluster$centers
  cusd <- diag(cov(df))
  
  init_curve_para <-  t(apply(cuM, 1, f3))
  init_sd_para <- c(mean(cusd),0.4) #SAD1
  init_pro <- pro
  
  #????mclust
  #init_cluster <- mclust::densityMclust(df,k)
  #cuM <- t(init_cluster$parameters$mean)
  #init_pro <- init_cluster$parameters$pro
  #init_sd_para <- c(mean(cusd),0.5)
  
  return_object <- list(init_sd_para,init_curve_para,init_pro)
  names(return_object)<-c("init_sd_par","init_curve","init_pro")
  return(return_object)
}
pos <- which(as.numeric(LR_plastic[1,]) > Threshold_value)

get_cluster <- function(data,k,input){
  
  requiredPackages = c("mvtnorm","reshape2")
  for(packages in requiredPackages){
    if(!require(packages,character.only = TRUE)) install.packages(packages)
    require(packages,character.only = TRUE)
  }
  Delta <- 100; iter <- 0; itermax <- 2000;
  SAD1_get_matrix <- function(par,data){
    p <-  ncol(data)
    v2 <- par[1]
    phi <- par[2]
    tmp <- (1-phi^2)
    sigma <- array(dim=c(p,p))
    for(i in 1:p){
      sigma[i,i:p] <- phi^( c(i:p) - i ) * (1-phi^(2*i ))/tmp
      sigma[i:p,i] <- sigma[i,i:p]}
    sigma <- sigma*abs(v2)
    return(sigma)
  } 
  AR1_get_matrix <- function(par,data){
    n <- ncol(data)
    rho <- par[2]
    exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                      (1:n - 1))
    return(par[1]*rho^exponent)
  }
  Legendre.model <-function( t, mu, tmin=NULL, tmax=NULL )
  {
    u <- -1;
    v <- 1;
    if (is.null(tmin)) tmin<-min(t);
    if (is.null(tmax)) tmax<-max(t);
    ti    <- u + ((v-u)*(t-tmin))/(tmax - tmin);
    np.order <- length(mu)-1;
    L <- mu[1] + ti*mu[2];
    if (np.order>=2)
      L <- L + 0.5*(3*ti*ti-1)* mu[3] ;
    if (np.order>=3)
      L <- L + 0.5*(5*ti^3-3*ti)*mu[4] ;
    if (np.order>=4)
      L <- L + 0.125*(35*ti^4-30*ti^2+3)* mu[5];
    if (np.order>=5)
      L <- L + 0.125*(63*ti^5-70*ti^3+15*ti)*mu[6];
    if (np.order>=6)
      L <- L + (1/16)*(231*ti^6-315*ti^4+105*ti^2-5)* mu[7];
    if (np.order>=7)
      L <- L + (1/16)*(429*ti^7-693*ti^5+315*ti^3-35*ti)* mu[8];
    if (np.order>=8)
      L <- L + (1/128)*(6435*ti^8-12012*ti^6+6930*ti^4-1260*ti^2+35)* mu[9];
    if (np.order>=9)
      L <- L + (1/128)*(12155*ti^9-25740*ti^7+18018*ti^5-4620*ti^3+315*ti)* mu[10];
    if (np.order>=10)
      L <- L + (1/256)*(46189*ti^10-109395*ti^8+90090*ti^6-30030*ti^4+3465*ti^2-63)* mu[11];
    if (np.order>=11)
    {
      for(r in 11:(np.order))
      {
        kk <- ifelse(r%%2==0, r/2, (r-1)/2);
        for (k in c(0:kk) )
        {
          L <- L + (-1)^k*factorial(2*r-2*k)/factorial(k)/factorial(r-k)/factorial(r-2*k)/(2^r)*ti^(r-2*k)*mu[r+1];
        }
      }
    }
    return(L);
  }
  mle <- function(par,data,prob){
    par1 <- par[1:2]
    par2 <- matrix(par[-c(1:2)],nrow = k,ncol = 2)
    temp_S <- sapply(1:k, function(c) dmvnorm(data,
                                              Legendre.model(t=c(1:48),mu=par2[c,]),
                                              SAD1_get_matrix(par1,data))*prob[c] )
    LL <- sum(-log(rowSums(temp_S)))
    return(LL)
  }
  while ( Delta > 0.01 && iter <= itermax ) {
    # initiation
    if(iter == 0){
      init_sd_para <- input[[1]]
      init_curve_para <- input[[2]]
      pro <- input[[3]]
    }
    #E step, calculate the posterior probability
    old_par <- c(init_sd_para,init_curve_para)
    LL_mem <- mle(old_par,data,pro)
    mvn.c <- sapply(1:k, function(c) dmvnorm(data,
                                             Legendre.model(t=c(1:48),mu=init_curve_para[c,]),
                                             SAD1_get_matrix(init_sd_para,data))*pro[c] )
    omega <- mvn.c/rowSums(mvn.c)
    #M step, calculate parameters
    pro <- colSums(omega)/sum(omega)
    
    new_par <- try(optim(old_par, mle, data=data, prob=pro, method = "Nelder-Mead"))
    if ('try-error' %in% class(new_par))
      break
    L_Value <- new_par$value
    init_sd_para <- new_par$par[1:2]
    init_curve_para <- matrix(new_par$par[-c(1:2)],nrow = k)
    Delta <- abs(L_Value-LL_mem)
    if (Delta > 500)
      break
    cat('\n',"iter=",iter,"LL=",L_Value,'\n')
    iter <- iter+1; LL_mem <- L_Value
  } 
  
  BIC <- 2*(L_Value)+log(nrow(data))*length(old_par)
  
  cluster <- apply(omega,1,which.max)
  #clustered_df <- data.frame(cbind(row.names(data),df_fitted,cluster))
  clustered_df <- data.frame(cbind(row.names(data),df,cluster))
  colnames(clustered_df) <- c("row.names(data)",age,"cluster")
  long_df <- melt(clustered_df,id.vars=c("row.names(data)","cluster"))
  long_df[,4] <- as.numeric(long_df[,4])
  long_df[,3] <- as.numeric(long_df[,3])
  colnames(long_df) <- c("gene","cluster","time","fpkm")
  
  clustered_df <- clustered_df[,-1]
  return_object <- list(init_sd_para,init_curve_para,pro,LL_mem,BIC,clustered_df)
  names(return_object)<-c("sd_par", "curve_par", "pro", "LL", "BIC", "clustered_data")
  
  return(return_object)
  
}
df <- varh2tx
input=get_init_par(df,14)

a <- get_cluster(df,14,input) 
fitall <- a$clustered_data[,49]
get_kcl <- function(data,clust){
  
  k=length( table(clust))
  kcl <- list()
  for (i in 1:k) { posi<- which(clust==i)
  kcl[[i]] <- data[posi,]
  }
  
  return(kcl)
  
  
  
}
fclall <- get_kcl(varh2tx*100,fitall)
whipos <- fitall[as.numeric(pos)]
whipos <- table(whipos)
whipos <- as.numeric(names(whipos))

#_____________________________________Construction of gene network____________________________________
library(ggplot2)
library(reshape2)
library(patchwork)
library(scales)
library(splines)
library(orthogonalsplinebasis)
library(MASS)
library(grplasso)
library(parallel)














dat_F <- function(marker,pheno,pi1=1:8,pi2=9:16){
  
  ind <- pheno[,1]
  n_ind <- length(ind)
  npheno <- pheno[,-1]
  ind_t <- colnames(marker)[-(1:3)]
  nmarker <- marker[,-c(1:3)]
  ind_t1 <- c()
  for(i in 1:length(ind_t)){
    ind_t1 <- c(ind_t1,as.numeric(strsplit(ind_t[i],"X")[[1]][2]))
  }
  nm_i <- c()
  for(i in 1:n_ind){
    nm_i <- c(nm_i,which(ind[i]==ind_t1))
  }
  nmarker1 <- nmarker[,nm_i]
  pheno1 <- npheno[,pi1];pheno2 <- npheno[,pi2]
  
  return(list(Marker=nmarker1,H=pheno1,D=pheno2))
}



varf <- function(p1,p2,marker){
  
  nm <- dim(marker)[1]
  eff_h <- c()
  eff_d <- c()
  for(i in 1:nm){
    sm <- marker[i,]
    if(any(sm=="--")){
      sm[which(sm=="--")] <- NA
    }
    ni <- names(table(as.character(as.matrix(sm))))
    dat1 <- c();dat2 <- c();
    for(j in 1:length(ni)){
      dat1 <- rbind(dat1,colMeans(p1[which(sm==ni[j]),]))
      dat2 <- rbind(dat2,colMeans(p2[which(sm==ni[j]),]))
    }
    h <- apply(dat1,2,var)
    d <- apply(dat2,2,var)
    eff_h <- cbind(eff_h,h)
    eff_d <- cbind(eff_d,d)
  }
  
  return(list(eh=eff_h,ed=eff_d)) 
}


LMall <- function(NX,nt,nstep=30,order){
  
  stp <- (max(nt)-min(nt))/nstep
  res <- c()
  for(j in 1:nstep){
    
    tg1 <- Legendre.model11((j-1)*stp+1,np.order=order-1,tmin=min(nt), tmax=max(nt))
    tg2 <- Legendre.model11(j*stp/2+1,np.order=order-1,tmin=min(nt), tmax=max(nt))
    tg3 <- Legendre.model11(j*stp/2+1,np.order=order-1,tmin=min(nt), tmax=max(nt))
    tg4 <- Legendre.model11(j*stp+1,np.order=order-1,tmin=min(nt), tmax=max(nt))
    tmp1 <- rbind(tg1,tg2,tg3,tg4)
    res <- rbind(res,tmp1)
  }
  res
}

fitPKM <- function(para,NG,self,nconnect,nt,order,nstep,LL){
  
  odes <- ode.sovle.ind(NG,para,nconnect,nt,order,nstep,LL)
  sum((NG[,self]-(rowSums(odes)+NG[1,self]))^2
     # +5*sum(para^2)
  )
}
ode.sovle.ind <- function(NG,fitpar,nconnect,nt,order,nstep,LL){
  
  stp <- (max(nt)-min(nt))/nstep
  index <- which(nconnect==1)
  
  ind.par <- matrix(fitpar[1:(length(index)*(order-1))],ncol=order-1,byrow=T)
  allrep <- matrix(rep(0,length(index)),nrow=1)
  nn <- 1
  for(j in 1:nstep){
    tg1 <- (rowSums(t(apply(ind.par,1,"*",LL[nn,])))*NG[j,index])
    tg2 <- (rowSums(t(apply(ind.par,1,"*",LL[nn+1,])))*NG[j,index])
    tg3 <- (rowSums(t(apply(ind.par,1,"*",LL[nn+2,])))*NG[j,index])
    tg4 <- (rowSums(t(apply(ind.par,1,"*",LL[nn+3,])))*NG[j,index])
    tmp <- allrep[j,] +stp*(tg1+2*tg2+2*tg3+tg4)/6
    allrep <- rbind(allrep,tmp)
    nn <- nn + 4
  }
  allrep
}





optim.parallel <- function(connect,effect,n.cores,proc,order,times,nstep){
  
  diag(connect) <- 1
  nt1 <- min(times)
  nt2 <- max(times)
  
  LL <- LMall(NX=1,nt=seq(nt1,nt2,(nt2-nt1)/nstep),nstep=nstep,order=order)
  
  nx <- dim(effect)[2]
  
  grp <- floor(nx/n.cores)
  grp.i <- c()
  if(n.cores==1){
    grp.i <- c(grp.i,rep(1,nx))
  }else{
    for(ii in 1:n.cores){
      if(ii==n.cores){
        grp.i <- c(grp.i,rep(ii,nx-grp*(ii-1)))
      }else{
        grp.i <- c(grp.i,rep(ii,grp))
      }
    }
  }
  
  grp.ii <- unique(grp.i)
  
  res.list <- mclapply(grp.ii, function(i)
  {
    y.c <- 	which(grp.i==i)
    A <- sapply(y.c, proc, connect=connect,effect=effect,LL=LL,nstep=nstep,order=order,times=times);
    return (unlist(A));
  }, mc.cores=n.cores )
  
  res1 <- do.call("c", res.list)
  res2 <- parallel.data.optim(res1,connect,times)
  return(res2)
}

parallel.data.optim <- function(rd,nm,ntt){
  
  nrd <- matrix(rd,nrow=length(ntt))
  nn <- dim(nm)[1]
  ki <- 0
  allist <- list()
  for(i in 1:nn){
    iii <- (which(nm[i,]==1))
    iiil <- length(iii)
    tmp.d <- nrd[,(ki+1):(ki+iiil)]
    if(is.matrix(tmp.d)){
      colnames(tmp.d) <- iii
    }else{
      names(tmp.d) <- iii
    }
    
    allist[[i]] <- tmp.d
    ki <- ki + iiil
  }
  
  return(allist)
}


ode.optim <- function(y.c,connect,effect,LL,nstep,order,times){
  
  indexx <- which(connect[y.c,]==1)
  para <- rep(0.0001,length(indexx)*(order-1))
  res <- optim(para,fitPKM,NG=(effect),self=y.c,nconnect=connect[y.c,],nt=times,order=order,nstep=nstep,
               LL=LL,method="BFGS",control=list(maxit=2000,trace=T))
  cat("Gene=",y.c," ",res$value,"\n")
  A <- ode.sovle.ind(NG=(effect),res$par,nconnect=connect[y.c,],nt=times,order=order,nstep=nstep,LL=LL)
  return(A)
}


interType <- function(con,alle,sme){
  
  diag(con) <- 0
  nn <- dim(con)[1]
  connfest <- matrix(0,nrow=nn,ncol=nn)
  indp <- c()
  inter <- list()
  for(i in 1:nn){
    al <- alle[[i]]
    index <- which(as.numeric(colnames(al))==i)
    if(is.matrix(al)){
      lindp <- al[,index]
      linter <- al[,-index]
      indp <- cbind(indp,lindp)
      inter[[i]] <- linter
      rcor <- cor(sme[i,],linter)
    }else{
      indp <- cbind(indp,al)
      inter[[i]] <- 0
      rcor <- 0
    }
    
    
    connfest[i,which(con[i,]==1)] <- as.numeric(rcor)
  }
  
  return(list(connfest=connfest,connect=con,indp=indp,inter=inter))
  
}



Legendre.model11 <- function(t, np.order,tmin = NULL, tmax = NULL)
{
  u <- -1;
  v <- 1;
  if (is.null(tmin))
    tmin <- min(t);
  if (is.null(tmax))
    tmax <- max(t);
  ti    <- u + ((v - u) * (t - tmin)) / (tmax - tmin);
  L <- rep(NA,np.order)
  L[1] <- 1;
  if (np.order >= 2)
    L[2] <- 0.5 * (6 * ti) 
  if (np.order >= 3)
    L[3] <- 0.5 * (15 * ti ^ 2 - 3) 
  if (np.order >= 4)
    L[4] <-  0.125 * (35 * 4 * ti ^ 3 - 60 * ti) 
  if (np.order >= 5)
    L[5] <-  0.125 * (63 * 5 * ti ^ 4 - 210 * ti ^ 2 + 15)
  if (np.order >= 6)
    L[6] <-(1 / 16) * (231 * 6 * ti ^ 5 - 315 * 4 * ti ^ 3 + 105 * 2 * ti) 
  if (np.order >= 7)
    L[7] <- (1 / 16) * (429 * 7 * ti ^ 6 - 693 * 5 * ti ^ 4 + 315 * 3 *ti ^ 2 - 35)
  if (np.order>=8)
    L[8] <-  (1/128) * (6435 * 8 * ti ^ 7 - 12012 * 6 * ti ^ 5 + 6930 * 4 * ti ^ 3 - 1260 * 2 * ti)
  if (np.order>=9)
    L[9] <-  (1/128) * (12155 * 9 * ti ^ 8 - 25740 * 7 * ti ^ 6 + 18018 * 5 * ti ^ 4 - 4620 * 3 * ti ^ 2 + 315)
  if (np.order>=10)
    L[10] <-  (1/256) * (46189 * 10 * ti ^ 9 - 109395 * 8 * ti ^ 7 + 90090 * 6 * ti ^ 5 - 30030 * 4 * ti ^ 3 + 3465 * 2 * ti)
  if (np.order>=11)
    L[11] <-  (1/512) * (176358 * 11 * ti ^ 10 - 461890 * 9 * ti ^ 8 + 437580 * 7 * ti ^ 6 - 180180 * 5 * ti ^ 4 + 30030 * 3 * ti^ 2 - 1386)
  if (np.order>=12)
    L[12] <-  660.1943 * 12 * ti ^ 11 -1894.471 * 10 * ti ^ 9 + 2029.79 * 8 * ti ^ 7  - 997.0898 * 6 * ti ^ 5 + 219.9463 * 4 * ti^ 3 - 17.5957* 2 * ti
  if (np.order>=13)
    L[13] <-  1269.604 * 13 * ti ^ 12 -3961.166 * 11 * ti ^ 10 + 4736.177 * 9 * ti ^ 8  - 2706.387 * 7 * ti ^ 6 + 747.8174 * 5 * ti^ 4 - 87.97852* 3 * ti^ 2 + 2.932617
  if (np.order>=14)
    L[14] <-  2448.523 * 14 * ti ^ 13 -8252.429 * 12 * ti ^ 11 + 10893.21 * 10 * ti ^ 9  - 7104.265 * 8 * ti ^ 7 + 2368.088 * 6 * ti^ 5 - 373.9087* 4 * ti^ 3 + 21.99463* 2 * ti 
  if (np.order>=15)
    L[15] <-  4733.811 * 15 * ti ^ 14 -17139.66 * 13 * ti ^ 12 + 24757.29 * 11 * ti ^ 10  - 18155.34 * 9 * ti ^ 8 + 7104.265 * 7 * ti^ 6 - 1420.853* 5 * ti^ 4 + 124.6362* 3 * ti ^ 2 - 3.14209
  if (np.order>=16)
    L[16] <-  9171.759 * 16 * ti ^ 15 -35503.58 * 14 * ti ^ 13 + 55703.9 * 12 * ti ^ 11  - 45388.36 * 10 * ti ^ 9 + 20424.76 * 8 * ti^ 7 - 4972.986* 6 * ti^ 5 + 592.0221* 4 * ti ^ 3 - 26.70776* 2 * ti 
  if (np.order>=17)
    L[17] <-  17804 * 17 * ti ^ 16 -73374.07 * 15 * ti ^ 14 + 124262.5 * 13 * ti ^ 12  - 111407.8 * 11 * ti ^ 10 + 56735.45 * 9 * ti^ 8 - 16339.81* 7 * ti^ 6 + 2486.493* 5 * ti ^ 4 - 169.1492* 3 * ti ^ 2 +3.33847
  if (np.order>=18)
    L[18] <-  34618.89 * 18 * ti ^ 17 -151334 * 16 * ti ^ 15 + 275152.8 * 14 * ti ^ 13  - 269235.5 * 12 * ti ^ 11 + 153185.7 * 10 * ti^ 9 - 51061.91* 8 * ti^ 7 + 9531.556* 6 * ti ^ 5 - 888.0331* 4 * ti ^ 3 +31.71547* 2 * ti 
  if (np.order>=19)
    L[19] <-  67415.74 * 19 * ti ^ 18 -311570 * 17 * ti ^ 16 + 605336.1 * 15 * ti ^ 14  - 642023.1 * 13 * ti ^ 12 + 403853.3 * 11 * ti^ 10 - 153185.7* 9 * ti^ 8 + 34041.27* 7 * ti ^ 6 - 4084.952* 5 * ti ^ 4 +222.0083* 3 * ti ^ 2 -3.523941
  if (np.order>=20)
    L[20] <-  131460.7 * 20 * ti ^ 19 -640449.5 * 18 * ti ^ 17 + 1324173 * 16 * ti ^ 15  - 1513340 * 14 * ti ^ 13 + 1043288 * 12 * ti^ 11 - 444238.6* 10 * ti^ 9 + 114889.3* 8 * ti ^ 7 - 17020.64* 6 * ti ^ 5 +1276.548* 4 * ti ^ 3 -37.00138* 2 * ti 
  if (np.order>=21)
    L[21] <-  256661.4 * 21 * ti ^ 20 -1314607 * 19 * ti ^ 18 + 2882023 * 17 * ti ^ 16  - 3531127 * 15 * ti ^ 14 + 2648345 * 13 * ti^ 12 - 1251945* 11 * ti^ 10 + 370198.8* 9 * ti ^ 8 - 65651.02* 7 * ti ^ 6 +6382.738* 5 * ti ^ 4 -283.6773* 3 * ti^ 2 + 3.700138 
  if (np.order>=22)
    L[22] <-  501656.3 * 22 * ti ^ 21 -2694944 * 20 * ti ^ 19 + 6244383 * 18 * ti ^ 17  - 8165732 * 16 * ti ^ 15 + 6620863 * 14 * ti^ 13 - 3442849* 12 * ti^ 11 + 1147616* 10 * ti ^ 9 - 237985* 8 * ti ^ 7 +28722.32* 6 * ti ^ 5 -1772.983* 4 * ti^ 3 + 42.55159* 2 * ti
  if (np.order>=23)
    L[23] <-  981501.4 * 23 * ti ^ 22 -5518219 * 21 * ti ^ 20 + 13474721 * 19 * ti ^ 18  - 18733149 * 17 * ti ^ 16 + 16331463 * 15 * ti^ 14 - 9269209* 13 * ti^ 12 + 3442849* 11 * ti ^ 10 - 819725.9* 9 * ti ^ 8 +118992.5* 7 * ti ^ 6 -9574.107* 5 * ti^ 4 + 354.5966* 3 * ti^ 2 -3.868326
  if (np.order>=24)
    L[24] <-  1922107 * 24 * ti ^ 23 -11287266 * 22 * ti ^ 21 + 28970650 * 20 * ti ^ 19  - 42669950 * 18 * ti ^ 17 + 39807941 * 16 * ti^ 15 - 24497195* 14 * ti^ 13 + 10041643* 12 * ti ^ 11 - 2705096* 10 * ti ^ 9 +461095.8* 8 * ti ^ 7 -46274.85* 6 * ti^ 5 + 2393.527* 4 * ti^ 3 -48.35408* 2 * ti
  if (np.order>=25)
    L[25] <-  3767330 * 25 * ti ^ 24 -23065284 * 23 * ti ^ 22 + 62079965 * 21 * ti ^ 20  - 96568835 * 19 * ti ^ 18 + 96007388 * 17 * ti^ 16 - 63692706* 15 * ti^ 14 + 28580061* 13 * ti ^ 12 - 8607122* 11 * ti ^ 10 +1690685* 9 * ti ^ 8 -204931.5* 7 * ti^ 6 + 13882.46* 5 * ti^ 4 -435.1867* 3 * ti ^ 2 + 4.029506
  if (np.order>=26)
    L[26] <-  7389762 * 26 * ti ^ 25 -47091621 * 24 * ti ^ 23 + 132625380 * 22 * ti ^ 21  - 217279879 * 20 * ti ^ 19 + 229350983 * 18 * ti^ 17 - 163212560* 16 * ti^ 15 + 79615883* 14 * ti ^ 13 - 26538628* 12 * ti ^ 11 +5917397* 10 * ti ^ 9 -845342.4* 8 * ti^ 7 + 71726.02* 6 * ti^ 5 -3155.104* 4 * ti ^ 3 + 54.39834* 2 * ti
  
  return(L);
}


get_cl_meancurve <- function(data,clust){
  
  
  library(ggplot2)
  library(reshape2)
  library(scales)
  k=length( table(clust))
  fclmean <- c()
  yansee <- hue_pal()(k)
  for (i in 1:k) {
    
    fclm <- colMeans(data[[i]])
    
    
    fclmean <- rbind(fclmean,fclm)
  }
  
  rownames(fclmean) <- c(1:k)
  
  tt <- ggplot()
  for (i in 1:10) {
    
    brt <- data.frame(x=c(1:48),y=fclmean[i,])
    
    tt <- tt +geom_line(data = brt,mapping = aes(x,y),col=yansee[i])
    
    
    
  }
  
  
  clmeancurve_plot <- tt
  
  
  
  
  
  marker <- t(fclmean)
  n <- dim(marker)[1]
  marker_list <- list()
  name <- colnames(marker)
  
  
  library(glmnet)
  
  for (col in 1:length(marker[1,])) {
    #tim <- proc.time()
    
    
    m <- marker[,col]
    M <- marker[,-col]
    
    x_matrix <- M
    ridge1_cv <- cv.glmnet(x = x_matrix, y = m,type.measure = "mse",nfold = 10,alpha = 0)
    best_ridge_coef <- as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.min))[-1]
    fit_res <- cv.glmnet(x = x_matrix, y = m,type.measure = "mse",nfold = 10,alpha = 1,penalty.factor = 1/ abs(best_ridge_coef),keep = TRUE)
    best_alasso_coef1 <- coef(fit_res, s = fit_res$lambda.min)
    marker_list_one <- list()
    marker_list_one[[1]] <- name[col]#??һ???б???ֱ??qtl??????
    marker_list_one[[2]] <- best_alasso_coef1@Dimnames[[1]][best_alasso_coef1@i[-1]+1]#?ڶ????б??Ǽ???qtl??????
    marker_list_one[[3]] <- best_alasso_coef1@x[-1]#???????б??Ǳ?��ѡ??ϵ??
    
    marker_list[[col]] <- marker_list_one}
 
  for (i in 1:length(marker_list)) {    
    if(length(marker_list[[i]][[3]])>1){
      
      
      lsr <- quantile(abs(marker_list[[i]][[3]]))[2]
      marker_list[[i]][[2]] <-  marker_list[[i]][[2]][which(abs(marker_list[[i]][[3]])>lsr)]
      marker_list[[i]][[3]] <-   marker_list[[i]][[3]][which(abs(marker_list[[i]][[3]])>lsr)]}
    
  } 
  marker_list
  
  
  lasso_1 <- matrix(0,nrow=dim(marker)[2],ncol=dim(marker)[2])
  colnames(lasso_1) <- colnames(marker)
  rownames(lasso_1) <- colnames(marker)
  for(c in 1:k){
    
    for(d in 1:length(marker_list[[c]][[2]])){
      lasso_1[c,which(name== marker_list[[c]][[2]][d])]=0.5
      
    }
  }
  
  lasso=lasso_1
  for (i in 1:k) { lasso[i,i] <- 1
  
  }
  
  for(c in 1:k){
    
    for(d in 1:length(marker_list[[c]][[2]])){
      lasso_1[c,which(name== marker_list[[c]][[2]][d])]=1
      
    }
  }
  F2014_H.odee <- optim.parallel(connect=lasso_1,effect=t(fclmean),
                                 n.cores=1,proc=ode.optim,order=13,times=c(1:48),nstep=c(47))
  fc_res <- list(
    
  )
  fc_res[[1]] <-  clmeancurve_plot
  
  fc_res[[2]] <- F2014_H.odee
  names(fc_res)<-c("clmeancurve_plot","ode_res"
  )
  return(fc_res) 
  
  
  
  
  
  
  
  
}
fclall_res <- get_cl_meancurve(fclall,fitall)





