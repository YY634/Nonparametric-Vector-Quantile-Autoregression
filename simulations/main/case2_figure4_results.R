###### set "reproducibility_materilas" as the working directory

library(MASS)
library(mgcv)
library(Matrix)
library(transport)
library(doParallel)


#### simulate the data sample.
d <- 2; T <- 48000; burn <- 800000; samp <- matrix(0, T+1, d)  
set.seed(820) 
X1 <- runif(d, -1, 1)
for (t in 1:(burn+T)){
  v <- sqrt(sum(X1^2)) 
  c1 <- tanh((X1[1]+X1[2])/2)-0.5;  c2 <- cos((pi/10)*(X1[1]+X1[2])/(1+abs(X1[1]+X1[2])) ) 
  if (t > burn){ samp[(t-burn),] <- X1 } 
  inn <- runif(d, -1, 1)
  X2 <-  c(c1, c2) + v*inn/2
  X1 <- X2
} 
samp[T+1,] <- X2 
rm(X1, X2, c1, c2, v, inn, t, burn) 


################## estimating the quantiles and centers of the stationary distribution (the center panel in Figure 4) #################

#### estimating the stationary mean and covariance (quantiles) of a VAR model 
library(vars) 
var_model <- VAR(samp, p = 1, type = "both")
var_results <- var_model$varresult   
Y1_results <- var_results$y1;  Y2_results <- var_results$y2
A <- matrix(0, d, d)
test1 <- Y1_results$coefficients[1:2]; names(test1) <- NULL
test2 <- Y2_results$coefficients[1:2]; names(test2) <- NULL
A[1,] <- test1;  A[2,] <- test2
C <- c(Y1_results$coefficients[3], Y2_results$coefficients[3]);  names(C) <- NULL
mu_hat <- solve(diag(d) - A, C)
U_hat <- residuals(var_model)
Omega_hat <- crossprod(U_hat) / (T-d)
Kron <- kronecker(A, A)
Sigma_vec <- solve(diag(d*d) - Kron, as.vector(Omega_hat))
Sigma_hat <- matrix(Sigma_vec, nrow = d, ncol = d)
save(mu_hat, Sigma_hat, file = "./simulations/main/case2_VAR_stationary.RData")
rm(A, C, var_results, Y1_results, Y2_results, test1, test2, mu_hat, U_hat, Omega_hat, Kron, Sigma_hat, Sigma_vec, var_model)


#### compute the center-outward median and quantiles of the stationary distribution
select_ind <- seq(2, T+1, by =2); X <- samp[select_ind,]
k_S <- 300; k_R <- 80; n <- k_S*k_R;  qgridn <- 480  
sphere <- cbind( cos(2*pi*(0:(k_S-1))/k_S), sin(2*pi*(0:(k_S-1))/k_S) )
Ugrid <- NULL;  for (i in 1:k_R){ Ugrid <- rbind(Ugrid, (i/(k_R+1)) * sphere) }
rm(sphere)
# compute the empirical OT map
a <- pp(Ugrid);  b <- pp(X)
OTM <- transport(a, b, p=2, method='auctionbf', fullreturn=FALSE, control = list(), threads=1)
targets <- X[OTM$to,]
rm(a, b, OTM)
# computation of cyclically monotone interpolation
normas <- rep(0,n); for (i in 1:n){ normas[i]<-sqrt(sum(targets[i,]^2)) }; nsup <- max(normas) 
xxx <- Ugrid/nsup;  yyy <- targets/nsup
cij<-list(apply(xxx*yyy,1,sum))
cij<-do.call(cbind,rep(cij,n))
cij<-cij-xxx%*%t(yyy)
ind.diag<-cbind(1:n,1:n)
cij[ind.diag]<-rep(Inf,n) 
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
# obtain the center-outward median and quantiles of the stationary distribution. 
centr <- T.0(c(0,0))
disc.sphere <- cbind(cos(2*pi*(0:qgridn)/qgridn),sin(2*pi*(0:qgridn)/qgridn))
quant1 <- t(apply(0.2*disc.sphere,1,T.0))  
quant2 <- t(apply(0.4*disc.sphere,1,T.0))
quant3 <- t(apply(0.8*disc.sphere,1,T.0))
M <- colMeans(samp)
stationary <- list(M, centr, quant1, quant2, quant3)
save(stationary, file = "./simulations/main/case2_stationary.RData")




############### estimating the one-step ahead predictive quantiles and centers conditional on different current values (the surrounding panels in Figure 4)  #################
d <- 2; T <- 2000000; burn <- 10000; samp <- matrix(0, T+1, d)  
set.seed(820) 
X1 <- runif(d, -1, 1)
for (t in 1:(burn+T)){
  v <- sqrt(sum(X1^2)) 
  c1 <- tanh((X1[1]+X1[2])/2)-0.5;  c2 <- cos((pi/10)*(X1[1]+X1[2])/(1+abs(X1[1]+X1[2])) ) 
  if (t > burn){ samp[(t-burn),] <- X1 } 
  inn <- runif(d, -1, 1)
  X2 <-  c(c1, c2) + v*inn/2
  X1 <- X2
} 
samp[T+1,] <- X2
rm(X1, X2, c1, c2, v, inn, t, burn)

#### fitting a vector autoregressive (VAR) model which allows heterogeneous covariance
library(vars) 
var_model <- VAR(samp, p = 1, type = "both")
var_results <- var_model$varresult  
Y1_results <- var_results$y1;  Y2_results <- var_results$y2
A <- matrix(0, d, d)
test1 <- Y1_results$coefficients[1:2]; names(test1) <- NULL
test2 <- Y2_results$coefficients[1:2]; names(test2) <- NULL
A[1,] <- test1;  A[2,] <- test2
C <- c(Y1_results$coefficients[3], Y2_results$coefficients[3]);  names(C) <- NULL
load("./simulations/main/case2_conds.RData")
K <- dim(conds)[1];   VAR_centers <- matrix(0, K,2)
for (k in 1:K){ VAR_centers[k,] <- C + A %*% conds[k,] }
rm(A, C, var_results, Y1_results, Y2_results,  k, test1, test2)
resis <- residuals(var_model) 
response_list <- list();  resp_names <- c()
for(i in 1:d) {
  for(j in 1:i) {
    response_list[[length(response_list) + 1]] <- (resis[, i] * resis[, j])
    resp_names <- c(resp_names, paste0("Sigma_", i, "_", j))
  }
}
response_mat <- do.call(cbind, response_list)   
colnames(response_mat) <- resp_names;  rm(response_list, resp_names)
# fit each element of the covariance matrix of the residual as a nonparametric function of the previous observation
M <- dim(response_mat)[2];  hatSigs <- matrix(0, K, M)   
data_train <- as.data.frame(samp[1:T,])
for (m in 1:M){
  y <- response_mat[,m];  data_train$y <- y  
  fit <- gam(y ~ s(V1, V2), data = data_train)
  hatSigs[, m] <- predict(fit, newdata = as.data.frame(conds))
}
VAR_SIGMs <- list()
for (k in 1:K){
  sigM <- matrix(0, 2, 2); sigM[1,1] <- hatSigs[k,1]
  sigM[2,1] <- sigM[1,2] <- hatSigs[k,2];  sigM[2,2] <- hatSigs[k,3]
  mtemp <- nearPD(sigM, corr = FALSE)$mat
  VAR_SIGMs[[k]] <- as.matrix(mtemp)
}  
save(VAR_centers, VAR_SIGMs, file = "./simulations/main/case2_VAR_predictive.RData")



#### compute the center-outward predictive medians and quantiles
Cond_Quantile_predict <- function(conds, h, qgridn, X, k_S, k_R, d=2){ 
  T <- dim(X)[1]-1;  n <- k_S*k_R;  tn <- dim(conds)[1]
  nc<-min(tn,detectCores());  registerDoParallel(nc)
  # generate uniform spherical grid denoted as Ugrid
  sphere <- cbind( cos(2*pi*(0:(k_S-1))/k_S), sin(2*pi*(0:(k_S-1))/k_S) )
  Ugrid <- NULL;  for (i in 1:k_R){ Ugrid <- rbind(Ugrid, (i/(k_R+1)) * sphere) }
  rm(sphere)
  
  quantiles <- foreach (t=1:tn) %dopar%{ 
    # compute the weights of the empirical conditional density
    cond <- conds[t,] 
    wgts <- rep(0, T);  for (j in 1:T){  wgts[j] <- exp(-sum((X[j,]-cond)^2)/h^2)  }
    ord = sort(wgts, decreasing = TRUE, index.return=TRUE)$ix[1:n]
    Xtp <- X[ord+1,];  wgts <- wgts[ord]
    wgts <- wgts/sum(wgts)  
    # compute the empirical OT plan
    a <- wpp(Ugrid, rep(1/n, n));  b <- wpp(Xtp, wgts)
    OTM <- transport(a, b, p = 2, method = "networkflow", fullreturn=FALSE, control=list(), threads=1)
    OTplan <- matrix(0, nrow = dim(Ugrid)[1], ncol = dim(Xtp)[1])
    for (i in 1:nrow(OTM)) {
      OTplan[OTM$from[i], OTM$to[i]] <- OTM$mass[i]
    } 
    rm(wgts, a, b, OTM, i)
    # construct the empirical quantile map (Ugrid[i,] maps to target[i,]) 
    targets <- matrix(0, n, d)
    for (i in 1:n){  targets[i,] <- n*OTplan[i,] %*% Xtp }
    # computation of cyclically monotone interpolation
    normas <- rep(0,n); for (i in 1:n){ normas[i]<-sqrt(sum(targets[i,]^2)) }; nsup <- max(normas) 
    xxx <- Ugrid/nsup;  yyy <- targets/nsup
    cij<-list(apply(xxx*yyy,1,sum))
    cij<-do.call(cbind,rep(cij,n))
    cij<-cij-xxx%*%t(yyy)
    ind.diag<-cbind(1:n,1:n)
    cij[ind.diag]<-rep(Inf,n)  
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
    # evaluate the quantile contours at three quantile levels
    centr <- T.0(c(0,0))
    disc.sphere<-cbind(cos(2*pi*(0:qgridn)/qgridn),sin(2*pi*(0:qgridn)/qgridn))
    quant1 <- t(apply(0.2*disc.sphere,1,T.0)) 
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

pairwisedist <- dist(samp[1:80000,]);  h <- 0.4*sum(pairwisedist)/length(pairwisedist)
k_S <- 200; k_R <- 80; qgridn <- 240 
#load("./simulations/main/case2_conds.RData")
prediction <- Cond_Quantile_predict(conds, h, qgridn, samp, k_S, k_R, d=2)
save(prediction, file = "./simulations/main/case2_predictive.RData")






