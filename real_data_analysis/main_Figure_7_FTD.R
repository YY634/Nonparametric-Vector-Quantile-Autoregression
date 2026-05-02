#################### This R script will produce the upper right panel (b) in Figure 7 of the main text. #####################

library(MASS)
library(transport)
library(doParallel)

load("./real_data_analysis/preprocessed.RData")
rm(CN_FTD, AD, CN_AD, dba_AD, dba_CNAD, dba_FTD, dba_CNFTD)
d <- 2

########## estimate the center-outward quantiles and medians of the stationary distribution (central panel)
samp <- matrix(0, 0, d)
for (i in 1:length(FTD)){  tmp <- FTD[[i]]; tmp <- tmp[1000:2000,];  samp <- rbind(samp, tmp) }
rm(tmp)

k_S <- 100; k_R <- 100; n <- k_S*k_R;  qgridn <- 120;  T <- dim(samp)[1]
sphere <- cbind( cos(2*pi*(0:(k_S-1))/k_S), sin(2*pi*(0:(k_S-1))/k_S) )
Ugrid <- NULL;  for (i in 1:k_R){ Ugrid <- rbind(Ugrid, (i/(k_R+1)) * sphere) }
rm(sphere)
a <- wpp(Ugrid, rep(1/n, n));  b <- wpp(samp, rep(1/T, T))
OTM <- transport(a, b, p = 2, method = "networkflow", fullreturn=FALSE, control = list(), threads=1)
OTplan <- matrix(0, nrow = dim(Ugrid)[1], ncol = dim(samp)[1])
for (i in 1:nrow(OTM)) {  OTplan[OTM$from[i], OTM$to[i]] <- OTM$mass[i]  } 
rm(a, b, OTM)
targets <- matrix(0, n, d)
for (i in 1:n){  targets[i,] <- n*OTplan[i,] %*% samp }

normas <- rep(0,n); for (i in 1:n){ normas[i]<-sqrt(sum(targets[i,]^2)) }; nsup <- max(normas) 
xxx <- Ugrid/nsup;  yyy <- targets/nsup
cij<-list(apply(xxx*yyy,1,sum))
cij<-do.call(cbind,rep(cij,n))
cij<-cij-xxx%*%t(yyy)
ind.diag<-cbind(1:n,1:n)
cij[ind.diag]<-rep(Inf,n)  # set the diagonal values of cij to be Inf
dkv<-matrix(Inf,nrow=n,ncol=n+1)
dkv[1,1]<-0
start.time <- Sys.time()
for (k in 2:(n+1)){
  aux.mat<-list(dkv[,(k-1)])
  aux.mat<-do.call(cbind,rep(aux.mat,n)) + cij
  dkv[,k]<-apply(aux.mat,2,min)
}
dndk<-list(dkv[,(n+1)])
dndk<-do.call(cbind,rep(dndk,n))-dkv[,1:n]
denom<-list(n:1)
denom<-do.call(rbind,rep(denom,n))
dndk<-dndk/denom
d.max<-apply(dndk,1,max)
mu.star<-min(d.max)
dkv.mu<-dkv-mu.star*matrix(1,ncol=1,nrow=n)%*%matrix(0:n,ncol=n+1,nrow=1)
shortest.distances<-apply(dkv.mu,1,min)
psi<-(-shortest.distances+shortest.distances[1])*nsup^2
e0<-abs(mu.star)*nsup^2
xxx <- Ugrid;  yyy <- targets
T.0 <- function(z){
  auxfun <- function(i){ return(sum(z*yyy[i,])-psi[i]) }
  scores<-sapply(1:n,auxfun)
  indice<-which.max(scores)
  return(yyy[indice,])
}
centr <- T.0(c(0,0))
disc.sphere <- cbind(cos(2*pi*(0:qgridn)/qgridn),sin(2*pi*(0:qgridn)/qgridn))
quant1 <- t(apply(0.2*disc.sphere,1,T.0))  
quant2 <- t(apply(0.4*disc.sphere,1,T.0))
quant3 <- t(apply(0.8*disc.sphere,1,T.0))

stationary <- list(centr, quant1, quant2, quant3)
#save(stationary, file = "./real_data_analysis/FTD_stationary.RData")


########## estimate the predictive quantiles and medians at different current values (surrounding panels)
Cond_Quantile_predict <- function(conds, h, qgridn, X, k_S, k_R, d=2){ 
  T <- dim(X)[1]-1;  n <- k_S*k_R;  tn <- dim(conds)[1]
  nc<-min(tn,detectCores());  registerDoParallel(nc)
  #generate uniform spherical grid denoted as Ugrid
  sphere <- cbind( cos(2*pi*(0:(k_S-1))/k_S), sin(2*pi*(0:(k_S-1))/k_S) )
  Ugrid <- NULL;  for (i in 1:k_R){ Ugrid <- rbind(Ugrid, (i/(k_R+1)) * sphere) }
  rm(sphere)
  quantiles <- foreach (t=1:tn) %dopar%{ 
    #compute the weights of the empirical conditional density
    cond <- conds[t,]  #is the conditional value
    wgts <- rep(0, T);  for (j in 1:T){  wgts[j] <- exp(-sum((X[j,]-cond)^2)/h^2)  }
    ord = sort(wgts, decreasing = TRUE, index.return=TRUE)$ix[1:n]
    Xtp <- X[ord+1,];  wgts <- wgts[ord]
    wgts <- wgts/sum(wgts)  
    # compute the empirical OT plan
    a <- wpp(Ugrid, rep(1/n, n));  b <- wpp(Xtp, wgts)
    OTM <- transport(a, b, p = 2, method = "networkflow", fullreturn=FALSE, control=list(), threads=1)
    OTplan <- matrix(0, nrow = dim(Ugrid)[1], ncol = dim(Xtp)[1])
    for (i in 1:nrow(OTM)) {  OTplan[OTM$from[i], OTM$to[i]] <- OTM$mass[i]  } 
    rm(wgts, a, b, OTM, i)
    # construct the empirical quantile map (Ugrid[i,] maps to target[i,]) 
    targets <- matrix(0, n, d)
    for (i in 1:n){  targets[i,] <- n*OTplan[i,] %*% Xtp }
    ##############################################
    ### computation of cyclically monotone interpolation
    ##############################################
    normas <- rep(0,n); for (i in 1:n){ normas[i]<-sqrt(sum(targets[i,]^2)) }; nsup <- max(normas) 
    xxx <- Ugrid/nsup;  yyy <- targets/nsup
    cij<-list(apply(xxx*yyy,1,sum))
    cij<-do.call(cbind,rep(cij,n))
    cij<-cij-xxx%*%t(yyy)
    ind.diag<-cbind(1:n,1:n)
    cij[ind.diag]<-rep(Inf,n)  # set the diagonal values of cij to be Inf
    dkv<-matrix(Inf,nrow=n,ncol=n+1)
    dkv[1,1]<-0
    start.time <- Sys.time()
    for (k in 2:(n+1)){
      aux.mat<-list(dkv[,(k-1)])
      aux.mat<-do.call(cbind,rep(aux.mat,n)) + cij
      dkv[,k]<-apply(aux.mat,2,min)
    }
    dndk<-list(dkv[,(n+1)])
    dndk<-do.call(cbind,rep(dndk,n))-dkv[,1:n]
    denom<-list(n:1)
    denom<-do.call(rbind,rep(denom,n))
    dndk<-dndk/denom
    d.max<-apply(dndk,1,max)
    mu.star<-min(d.max)
    dkv.mu<-dkv-mu.star*matrix(1,ncol=1,nrow=n)%*%matrix(0:n,ncol=n+1,nrow=1)
    shortest.distances<-apply(dkv.mu,1,min)
    psi<-(-shortest.distances+shortest.distances[1])*nsup^2
    e0<-abs(mu.star)*nsup^2
    xxx <- Ugrid;  yyy <- targets
    T.0 <- function(z){
      auxfun <- function(i){ return(sum(z*yyy[i,])-psi[i]) }
      scores<-sapply(1:n,auxfun)
      indice<-which.max(scores)
      return(c(t,yyy[indice,]))
    }
    #evaluate the quantile contours at three quantile levels
    centr <- T.0(c(0,0))
    disc.sphere<-cbind(cos(2*pi*(0:qgridn)/qgridn),sin(2*pi*(0:qgridn)/qgridn))
    quant1 <- t(apply(0.2*disc.sphere,1,T.0))  #quant1/2/3 is a qgridn*(1+d) matrix
    quant2 <- t(apply(0.4*disc.sphere,1,T.0))
    quant3 <- t(apply(0.8*disc.sphere,1,T.0))
    return( list(center=centr, quantile1=quant1, quantile2=quant2, quantile3=quant3) )
  }
  medians <- NULL; quant_list <- list()
  for (a in 1:tn){
    medians <- rbind(medians, quantiles[[a]]$center)
    quant_list[[a]] <- list(quantiles[[a]]$quantile1, quantiles[[a]]$quantile2, quantiles[[a]]$quantile3)
  } 
  return(list(medians=medians, quantiles=quant_list))  
}

k_S <- 100; k_R <- 100; qgridn <- 120 
X <- matrix(0, 0, d)
for (i in 1:length(FTD)){  tmp <- FTD[[i]]; L <- min(dim(tmp)[1], 100000); tmp <- tmp[1:L,];  X <- rbind(X, tmp) }
rm(tmp)
pairwisedist <- dist(X[1:20000,]);  h <- 0.2*sum(pairwisedist)/length(pairwisedist) 
rm(pairwisedist)
load("./real_data_analysis/FTD_conds.RData") 
prediction <- Cond_Quantile_predict(FTD_conds, h, qgridn, X, k_S, k_R, d=2)
#save(prediction, file = "./real_data_analysis/FTD_predictive.RData")



############################ Plot the above results (visualization) ###################################

medians <- prediction[[1]]; quantiles <- prediction[[2]]
library(plotrix)
layout(matrix(c(1,2,3, 4,5,6, 7,8,9), 3, 3, byrow = TRUE))
par(mai=c(0.3, 0.4, 0.1, 0.2))

for (i in 1:9){
  if (i < 5){
    medi <- medians[i,2:3]; Qtmp <- quantiles[[i]]
    quant1 <- Qtmp[[1]]; quant2 <- Qtmp[[2]]; quant3 <- Qtmp[[3]]
    plot(quant1[,2:3], type = "l", col = "brown", lwd = 2, xlim = c(-1, 1), ylim = c(-1, 1), asp=1, xaxt ="n", yaxt = "n", xlab = "", ylab = "") 
    xticks <- c(-1, -0.5, 0, 0.5, 1);  axis(1, at = xticks, labels = round(xticks, 1))
    yticks <- c(-1, -0.5, 0, 0.5, 1);  axis(2, at = yticks, labels = round(yticks, 1))
    lines(quant2[,2:3], col = 'tomato', lwd = 2)
    lines(quant3[,2:3], col = 'gold', lwd = 2)
    points(medi[1], medi[2], col = "darkred", pch = 16, cex=1) 
    rm(medi, Qtmp, quant1, quant2, quant3)
  }else if (i > 5){
    medi <- medians[i-1,2:3]; Qtmp <- quantiles[[i-1]]
    quant1 <- Qtmp[[1]]; quant2 <- Qtmp[[2]]; quant3 <- Qtmp[[3]]
    plot(quant1[,2:3], type = "l", col = 'brown', lwd = 2, xlim = c(-1, 1), ylim = c(-1, 1), asp=1, xaxt ="n", yaxt = "n", xlab = "", ylab = "")
    xticks <- c(-1, -0.5, 0, 0.5, 1);  axis(1, at = xticks, labels = round(xticks, 1))
    yticks <- c(-1, -0.5, 0, 0.5, 1);  axis(2, at = yticks, labels = round(yticks, 1))
    lines(quant2[,2:3], col = 'tomato', lwd = 2)
    lines(quant3[,2:3], col = 'gold', lwd = 2)
    points(medi[1], medi[2], col = "darkred", pch = 16, cex=1)
    rm(medi, Qtmp, quant1, quant2, quant3)
  }else{
    medi <- stationary[[1]]
    quant1 <- stationary[[2]]; quant2 <- stationary[[3]]; quant3 <- stationary[[4]]
    plot(quant1[,1:2], type = "l", col = 'brown', lwd = 2, xlim = c(-1, 1), ylim = c(-1, 1), asp=1, xaxt ="n", yaxt = "n", xlab = "", ylab = "") 
    xticks <- c(-1, -0.5, 0, 0.5, 1);  axis(1, at = xticks, labels = round(xticks, 1))
    yticks <- c(-1, -0.5, 0, 0.5, 1);  axis(2, at = yticks, labels = round(yticks, 1))
    lines(quant2[,1:2], col = 'tomato', lwd = 2)
    lines(quant3[,1:2], col = 'gold', lwd = 2)
    points(medi[1],medi[2], col = "darkred", pch = 16, cex=1)
    for (k in 1:8){ points(FTD_conds[k,1], FTD_conds[k,2], col = "seagreen3", pch = 17, cex=1.2) }
  }
}



