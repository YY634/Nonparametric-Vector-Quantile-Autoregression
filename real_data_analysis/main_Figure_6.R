
library(MASS)
library(transport)
library(doParallel)

# cond is the value that conditional on; h is the kernel width; 
# Ugrid is the reference measure, empirical grids which is a matrix of n*d, n=T-1
# tims is a vector containing time points at which to compute the quantiles
Cond_Quantile <- function(tims, conds, X, h, qgridn, k_S, k_R, d=2){ 
  T <- dim(X)[1]-1;  n <- k_S*k_R;  tn <- length(tims)
  nc<-min(tn,detectCores());  registerDoParallel(nc)
  #generate uniform spherical grid denoted as Ugrid
  sphere <- cbind( cos(2*pi*(0:(k_S-1))/k_S), sin(2*pi*(0:(k_S-1))/k_S) )
  Ugrid <- NULL;  for (i in 1:k_R){ Ugrid <- rbind(Ugrid, (i/(k_R+1)) * sphere) }
  rm(sphere)
  quantiles <- foreach (t=1:tn) %dopar%{ 
    #compute the weights of the empirical conditional density
    cond <- conds[t,] 
    wgts <- rep(0, T);  for (j in 1:T){  wgts[j] <- exp(-sum((X[j,]-cond)^2)/h^2)  }
    ord = sort(wgts, decreasing = TRUE, index.return=TRUE)$ix[1:n]
    Xtp <- X[ord+1,];  wgts <- wgts[ord]
    wgts <- wgts/sum(wgts)  
    # compute the empirical OT plan
    a <- wpp(Ugrid, rep(1/n, n));  b <- wpp(Xtp, wgts)
    OTM <- transport(a, b, p = 2, method = "networkflow", fullreturn=FALSE, control = list(), threads=1)
    OTplan <- matrix(0, nrow = dim(Ugrid)[1], ncol = dim(Xtp)[1])
    for (i in 1:nrow(OTM)) { OTplan[OTM$from[i], OTM$to[i]] <- OTM$mass[i] } 
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
      scores <- sapply(1:n,auxfun)
      indice <- which.max(scores)
      return(c(tims[t],yyy[indice,]))
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


load("./real_data_analysis/preprocessed_data/consensus_series.RData")
d <- 2;  k_S <- 100;  k_R <- 100;  qgridn <- 120;  tims <- seq(1000, 12000, by=1000) 

load("./real_data_analysis/preprocessed_data/dataset_CNFTD.RData") 
conds <- dba_CNFTD[tims,];  X <- matrix(0, 0, d)
for (i in 1:length(CN_FTD)){  tmp <- CN_FTD[[i]]; L <- min(dim(tmp)[1], 100000); tmp <- tmp[1:L,];  X <- rbind(X, tmp) }
rm(CN_FTD, tmp, dba_CNFTD)
pairwisedist <- dist(X[1:20000,]);  h <- 0.2*sum(pairwisedist)/length(pairwisedist) 
rm(pairwisedist)
upperleft <- Cond_Quantile(tims, conds, X, h, qgridn, k_S, k_R, d=2)
#save(upperleft, file = "./real_data_analysis/Figure6_CNFTD.RData")


load("./real_data_analysis/preprocessed_data/dataset_FTD.RData") 
conds <- dba_FTD[tims,];  X <- matrix(0, 0, d)
for (i in 1:length(FTD)){  tmp <- FTD[[i]]; L <- min(dim(tmp)[1], 100000); tmp <- tmp[1:L,];  X <- rbind(X, tmp) }
rm(FTD, tmp, dba_FTD)
pairwisedist <- dist(X[1:20000,]);  h <- 0.2*sum(pairwisedist)/length(pairwisedist) 
rm(pairwisedist)
upperright <- Cond_Quantile(tims, conds, X, h, qgridn, k_S, k_R, d=2)
#save(upperright, file = "./real_data_analysis/Figure6_FTD.RData")


load("./real_data_analysis/preprocessed_data/dataset_CNAD.RData") 
conds <- dba_CNAD[tims,];  X <- matrix(0, 0, d)
for (i in 1:length(CN_AD)){  tmp <- CN_AD[[i]]; L <- min(dim(tmp)[1], 100000); tmp <- tmp[1:L,];  X <- rbind(X, tmp) }
rm(CN_AD, tmp, dba_CNAD)
pairwisedist <- dist(X[1:20000,]);  h <- 0.2*sum(pairwisedist)/length(pairwisedist) 
rm(pairwisedist)
lowerleft <- Cond_Quantile(tims, conds, X, h, qgridn, k_S, k_R, d=2)
#save(lowerleft, file = "./real_data_analysis/Figure6_CNAD.RData")


load("./real_data_analysis/preprocessed_data/dataset_AD.RData") 
conds <- dba_AD[tims,];  X <- matrix(0, 0, d)
for (i in 1:length(AD)){  tmp <- AD[[i]]; L <- min(dim(tmp)[1], 100000); tmp <- tmp[1:L,];  X <- rbind(X, tmp) }
rm(AD, tmp, dba_AD)
pairwisedist <- dist(X[1:20000,]);  h <- 0.2*sum(pairwisedist)/length(pairwisedist) 
rm(pairwisedist)
lowerright <- Cond_Quantile(tims, conds, X, h, qgridn, k_S, k_R, d=2)
#save(lowerright, file = "./real_data_analysis/Figure6_AD.RData")






##################################### plot the above results ############################################
grid3d_faces <- function(nx, ny, nz, x_plane = "min", y_plane = "max", z_plane = "min", col = "gray", lwd = 3) {
  lim <- par3d("bbox")
  if (is.null(lim) || length(lim) != 6) stop("No rgl scene found. Run open3d(); plot3d(...) first.")
  xlim <- lim[1:2]; ylim <- c(-1, 1); zlim <- c(-1,1) #ylim <- lim[3:4]; zlim <- lim[5:6]
  
  plane_pos <- function(pos, range) {
    if (is.null(pos)) return(NULL)
    if (is.character(pos)) {
      pos <- match.arg(pos, choices = c("min","max"))
      return(if (pos == "min") range[1] else range[2])
    }
    if (is.numeric(pos) && length(pos) == 1) {
      return(max(min(pos, range[2]), range[1])) 
    }
    stop("Plane must be 'min','max', numeric, or NULL")
  }
  x0 <- plane_pos(x_plane, xlim)  # YZ plane at x = x0
  y0 <- plane_pos(y_plane, ylim)  # XZ plane at y = ymax
  z0 <- plane_pos(z_plane, zlim)  # XY plane at z = z0
  x.at <- seq(xlim[1], xlim[2], length.out = nx)
  y.at <- seq(ylim[1], ylim[2], length.out = ny)
  z.at <- seq(zlim[1], zlim[2], length.out = nz)
  # XY plane at z = z0
  if (!is.null(z0)) {
    for (xv in x.at) lines3d(c(xv, xv), c(ylim[1], ylim[2]), c(z0, z0), col = col, lwd = lwd) # vertical in y
    for (yv in y.at) lines3d(c(xlim[1], xlim[2]), c(yv, yv), c(z0, z0), col = col, lwd = lwd) # horizontal in x
  }
  # XZ plane at y = ymax
  if (!is.null(y0)) {
    for (xv in x.at) lines3d(c(xv, xv), c(y0, y0), c(zlim[1], zlim[2]), col = col, lwd = lwd) # vertical in z
    for (zv in z.at) lines3d(c(xlim[1], xlim[2]), c(y0, y0), c(zv, zv), col = col, lwd = lwd) # horizontal in x
  }
  # YZ plane at x = x0
  if (!is.null(x0)) {
    for (yv in y.at) lines3d(c(x0, x0), c(yv, yv), c(zlim[1], zlim[2]), col = col, lwd = lwd) # vertical in z
    for (zv in z.at) lines3d(c(x0, x0), c(ylim[1], ylim[2]), c(zv, zv), col = col, lwd = lwd) # horizontal in y
  }
}

options(rgl.useNULL = TRUE)
options(rgl.printRglwidget = TRUE)
library(rgl)
library(MASS)
library(htmlwidgets)
library(webshot2)
load("./real_data_analysis/3D_view.RData")

d <- 2; tims <- seq(1000, 12000, by=1000)/1000; tn <- length(tims)

medians <- upperleft[[1]];  quant_list <- upperleft[[2]]
clear3d()
par3d(userMatrix = rotation_matrix);aspect3d(1,2.5,2.5);  highlevel()
t <- seq(0, 1, length.out = 10)
for (i in 1:(tn-1)){
  x_line <- (1 - t) * tims[i] + t * tims[i+1]
  y_line <- (1 - t) * medians[i,2] + t * medians[i+1,2]
  z_line <- (1 - t) * medians[i,3] + t * medians[i+1,3]
  for (h in seq(1, length(t)-1, by = 2)) { lines3d(x_line[h:(h+1)], y_line[h:(h+1)], z_line[h:(h+1)], col = "seagreen3", lwd = 8) }
}
for (i in 1:tn){
  temp <- quant_list[[i]];  Q1 <- temp[[1]];  Q2 <- temp[[2]];  Q3 <- temp[[3]]
  points3d(x=tims[i], y=medians[i,2], z=medians[i,3], zlim=c(-1, 1), col = "seagreen3", size=16)
  lines3d(x=tims[i], y=Q3[,2], z=Q3[,3], col='gold', add=TRUE, lwd=8)
  lines3d(x=tims[i], y=Q2[,2], z=Q2[,3], col='tomato', add=TRUE, lwd=8)
  lines3d(x=tims[i], y=Q1[,2], z=Q1[,3], col='brown', add=TRUE, lwd=8)
}
par3d(windowRect=c(100,100,700,700))
grid3d_faces(nx = 7, ny = 5, nz = 5, x_plane = "min", y_plane = "max", z_plane = "min", col = "gray", lwd = 3)
x_ticks <- seq(2, 10, length.out = 5); axis3d("x-", at = x_ticks, labels = round(x_ticks, 1), cex = 3, lwd = 3) 
y_ticks <- seq(-1, 1, length.out = 3); axis3d("y+", at = y_ticks, labels = round(y_ticks, 1), cex = 3, lwd = 3)
z_ticks <- seq(-1, 1, length.out = 3); axis3d("z-", at = z_ticks, labels = round(z_ticks, 1), cex = 3, lwd = 3)
rglwidget() 
snapshot3d("./real_data_analysis/Figure6_CNFTD.png", width = 3000, height = 3000)  


medians <- upperright[[1]];  quant_list <- upperright[[2]]
clear3d()
par3d(userMatrix = rotation_matrix);aspect3d(1,2.5,2.5);  highlevel()
t <- seq(0, 1, length.out = 10)
for (i in 1:(tn-1)){
  x_line <- (1 - t) * tims[i] + t * tims[i+1]
  y_line <- (1 - t) * medians[i,2] + t * medians[i+1,2]
  z_line <- (1 - t) * medians[i,3] + t * medians[i+1,3]
  for (h in seq(1, length(t)-1, by = 2)) { lines3d(x_line[h:(h+1)], y_line[h:(h+1)], z_line[h:(h+1)], col = "seagreen3", lwd = 8) }
}
for (i in 1:tn){
  temp <- quant_list[[i]];  Q1 <- temp[[1]];  Q2 <- temp[[2]];  Q3 <- temp[[3]]
  points3d(x=tims[i], y=medians[i,2], z=medians[i,3], zlim=c(-1, 1), col = "seagreen3", size=16)
  lines3d(x=tims[i], y=Q3[,2], z=Q3[,3], col='gold', add=TRUE, lwd=8)
  lines3d(x=tims[i], y=Q2[,2], z=Q2[,3], col='tomato', add=TRUE, lwd=8)
  lines3d(x=tims[i], y=Q1[,2], z=Q1[,3], col='brown', add=TRUE, lwd=8)
}
par3d(windowRect=c(100,100,700,700))
grid3d_faces(nx = 7, ny = 5, nz = 5, x_plane = "min", y_plane = "max", z_plane = "min", col = "gray", lwd = 3)
x_ticks <- seq(2, 10, length.out = 5); axis3d("x-", at = x_ticks, labels = round(x_ticks, 1), cex = 3, lwd = 3) 
y_ticks <- seq(-1, 1, length.out = 3); axis3d("y+", at = y_ticks, labels = round(y_ticks, 1), cex = 3, lwd = 3)
z_ticks <- seq(-1, 1, length.out = 3); axis3d("z-", at = z_ticks, labels = round(z_ticks, 1), cex = 3, lwd = 3)
rglwidget() 
snapshot3d("./real_data_analysis/Figure6_FTD.png", width = 3000, height = 3000)  


medians <- lowerleft[[1]];  quant_list <- lowerleft[[2]]
clear3d()
par3d(userMatrix = rotation_matrix);aspect3d(1,2.5,2.5);  highlevel()
t <- seq(0, 1, length.out = 10)
for (i in 1:(tn-1)){
  x_line <- (1 - t) * tims[i] + t * tims[i+1]
  y_line <- (1 - t) * medians[i,2] + t * medians[i+1,2]
  z_line <- (1 - t) * medians[i,3] + t * medians[i+1,3]
  for (h in seq(1, length(t)-1, by = 2)) { lines3d(x_line[h:(h+1)], y_line[h:(h+1)], z_line[h:(h+1)], col = "seagreen3", lwd = 8) }
}
for (i in 1:tn){
  temp <- quant_list[[i]];  Q1 <- temp[[1]];  Q2 <- temp[[2]];  Q3 <- temp[[3]]
  points3d(x=tims[i], y=medians[i,2], z=medians[i,3], zlim=c(-1, 1), col = "seagreen3", size=16)
  lines3d(x=tims[i], y=Q3[,2], z=Q3[,3], col='gold', add=TRUE, lwd=8)
  lines3d(x=tims[i], y=Q2[,2], z=Q2[,3], col='tomato', add=TRUE, lwd=8)
  lines3d(x=tims[i], y=Q1[,2], z=Q1[,3], col='brown', add=TRUE, lwd=8)
}
par3d(windowRect=c(100,100,700,700))
grid3d_faces(nx = 7, ny = 5, nz = 5, x_plane = "min", y_plane = "max", z_plane = "min", col = "gray", lwd = 3)
x_ticks <- seq(2, 10, length.out = 5); axis3d("x-", at = x_ticks, labels = round(x_ticks, 1), cex = 3, lwd = 3) 
y_ticks <- seq(-1, 1, length.out = 3); axis3d("y+", at = y_ticks, labels = round(y_ticks, 1), cex = 3, lwd = 3)
z_ticks <- seq(-1, 1, length.out = 3); axis3d("z-", at = z_ticks, labels = round(z_ticks, 1), cex = 3, lwd = 3)
rglwidget() 
snapshot3d("./real_data_analysis/Figure6_CNAD.png", width = 3000, height = 3000)  


medians <- lowerright[[1]];  quant_list <- lowerright[[2]]
clear3d()
par3d(userMatrix = rotation_matrix);aspect3d(1,2.5,2.5);  highlevel()
t <- seq(0, 1, length.out = 10)
for (i in 1:(tn-1)){
  x_line <- (1 - t) * tims[i] + t * tims[i+1]
  y_line <- (1 - t) * medians[i,2] + t * medians[i+1,2]
  z_line <- (1 - t) * medians[i,3] + t * medians[i+1,3]
  for (h in seq(1, length(t)-1, by = 2)) { lines3d(x_line[h:(h+1)], y_line[h:(h+1)], z_line[h:(h+1)], col = "seagreen3", lwd = 8) }
}
for (i in 1:tn){
  temp <- quant_list[[i]];  Q1 <- temp[[1]];  Q2 <- temp[[2]];  Q3 <- temp[[3]]
  points3d(x=tims[i], y=medians[i,2], z=medians[i,3], zlim=c(-1, 1), col = "seagreen3", size=16)
  lines3d(x=tims[i], y=Q3[,2], z=Q3[,3], col='gold', add=TRUE, lwd=8)
  lines3d(x=tims[i], y=Q2[,2], z=Q2[,3], col='tomato', add=TRUE, lwd=8)
  lines3d(x=tims[i], y=Q1[,2], z=Q1[,3], col='brown', add=TRUE, lwd=8)
}
par3d(windowRect=c(100,100,700,700))
grid3d_faces(nx = 7, ny = 5, nz = 5, x_plane = "min", y_plane = "max", z_plane = "min", col = "gray", lwd = 3)
x_ticks <- seq(2, 10, length.out = 5); axis3d("x-", at = x_ticks, labels = round(x_ticks, 1), cex = 3, lwd = 3) 
y_ticks <- seq(-1, 1, length.out = 3); axis3d("y+", at = y_ticks, labels = round(y_ticks, 1), cex = 3, lwd = 3)
z_ticks <- seq(-1, 1, length.out = 3); axis3d("z-", at = z_ticks, labels = round(z_ticks, 1), cex = 3, lwd = 3)
rglwidget() 
snapshot3d("./real_data_analysis/Figure6_AD.png", width = 3000, height = 3000)  



