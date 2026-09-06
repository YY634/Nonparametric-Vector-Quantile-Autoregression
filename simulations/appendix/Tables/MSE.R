
################################## set "reproducibility_materilas" as the working directory #########################################
setwd("~/Desktop/reproducibility_materials/")
#The MSE is computed as the mean square discrepancy between the estimated quantile and theoretical quantile on a common set of grid of the uniform spherical distribution. 



#Appendix Table 1 (case 1):
load("./simulations/appendix/Tables/MSE_case1.RData")
TM <- theoretic[[1]];  TQ <- theoretic[[2]]

# compute the MSEs of the estimated quantiles with T=40000 and h=0.5
M <- T40000h05[[1]];  Q <- T40000h05[[2]]
M_sqe <- rowSums((M[,2:3] - TM[,2:3])^2);  M_mse <- mean(M_sqe);  rm(M_sqe)
Q1_sqe <- 0;  Q2_sqe <- 0;  Q3_sqe <- 0
for (i in 1:length(TQ)){
  T_tmp <- TQ[[i]];  Q_tmp <- Q[[i]]
  TQ1 <- T_tmp[[1]]; Q1 <- Q_tmp[[1]]; Q1_sqe <- Q1_sqe + mean( rowSums( (TQ1[,2:3] - Q1[,2:3])^2 ) )
  TQ2 <- T_tmp[[2]]; Q2 <- Q_tmp[[2]]; Q2_sqe <- Q2_sqe + mean( rowSums( (TQ2[,2:3] - Q2[,2:3])^2 ) )
  TQ3 <- T_tmp[[3]]; Q3 <- Q_tmp[[3]]; Q3_sqe <- Q3_sqe + mean( rowSums( (TQ3[,2:3] - Q3[,2:3])^2 ) )
}
Q1_mse <- Q1_sqe/10;  Q2_mse <- Q2_sqe/10;  Q3_mse <- Q3_sqe/10;  rm(Q1_sqe, Q2_sqe, Q3_sqe)
MSE_T40000h05 <- c(M_mse, Q1_mse, Q2_mse, Q3_mse)
rm(M, Q, M_mse, i, T_tmp, Q_tmp, TQ1, TQ2, TQ3, Q1, Q2, Q3, Q1_mse, Q2_mse, Q3_mse)

# compute the MSEs of the estimated quantiles with T=80000 and h=0.5
M <- T80000h05[[1]];  Q <- T80000h05[[2]]
M_sqe <- rowSums((M[,2:3] - TM[,2:3])^2);  M_mse <- mean(M_sqe);  rm(M_sqe)
Q1_sqe <- 0;  Q2_sqe <- 0;  Q3_sqe <- 0
for (i in 1:length(TQ)){
  T_tmp <- TQ[[i]];  Q_tmp <- Q[[i]]
  TQ1 <- T_tmp[[1]]; Q1 <- Q_tmp[[1]]; Q1_sqe <- Q1_sqe + mean( rowSums( (TQ1[,2:3] - Q1[,2:3])^2 ) )
  TQ2 <- T_tmp[[2]]; Q2 <- Q_tmp[[2]]; Q2_sqe <- Q2_sqe + mean( rowSums( (TQ2[,2:3] - Q2[,2:3])^2 ) )
  TQ3 <- T_tmp[[3]]; Q3 <- Q_tmp[[3]]; Q3_sqe <- Q3_sqe + mean( rowSums( (TQ3[,2:3] - Q3[,2:3])^2 ) )
}
Q1_mse <- Q1_sqe/10;  Q2_mse <- Q2_sqe/10;  Q3_mse <- Q3_sqe/10;  rm(Q1_sqe, Q2_sqe, Q3_sqe)
MSE_T80000h05 <- c(M_mse, Q1_mse, Q2_mse, Q3_mse)
rm(M, Q, M_mse, i, T_tmp, Q_tmp, TQ1, TQ2, TQ3, Q1, Q2, Q3, Q1_mse, Q2_mse, Q3_mse)

# compute the MSEs of the estimated quantiles with T=400000 and h=0.5
M <- T400000h05[[1]];  Q <- T400000h05[[2]]
M_sqe <- rowSums((M[,2:3] - TM[,2:3])^2);  M_mse <- mean(M_sqe);  rm(M_sqe)
Q1_sqe <- 0;  Q2_sqe <- 0;  Q3_sqe <- 0
for (i in 1:length(TQ)){
  T_tmp <- TQ[[i]];  Q_tmp <- Q[[i]]
  TQ1 <- T_tmp[[1]]; Q1 <- Q_tmp[[1]]; Q1_sqe <- Q1_sqe + mean( rowSums( (TQ1[,2:3] - Q1[,2:3])^2 ) )
  TQ2 <- T_tmp[[2]]; Q2 <- Q_tmp[[2]]; Q2_sqe <- Q2_sqe + mean( rowSums( (TQ2[,2:3] - Q2[,2:3])^2 ) )
  TQ3 <- T_tmp[[3]]; Q3 <- Q_tmp[[3]]; Q3_sqe <- Q3_sqe + mean( rowSums( (TQ3[,2:3] - Q3[,2:3])^2 ) )
}
Q1_mse <- Q1_sqe/10;  Q2_mse <- Q2_sqe/10;  Q3_mse <- Q3_sqe/10;  rm(Q1_sqe, Q2_sqe, Q3_sqe)
MSE_T400000h05 <- c(M_mse, Q1_mse, Q2_mse, Q3_mse)
rm(M, Q, M_mse, i, T_tmp, Q_tmp, TQ1, TQ2, TQ3, Q1, Q2, Q3, Q1_mse, Q2_mse, Q3_mse)

# compute the MSEs of the estimated quantiles with T=400000 and h=0.1
M <- T400000h01[[1]];  Q <- T400000h01[[2]]
M_sqe <- rowSums((M[,2:3] - TM[,2:3])^2);  M_mse <- mean(M_sqe);  rm(M_sqe)
Q1_sqe <- 0;  Q2_sqe <- 0;  Q3_sqe <- 0
for (i in 1:length(TQ)){
  T_tmp <- TQ[[i]];  Q_tmp <- Q[[i]]
  TQ1 <- T_tmp[[1]]; Q1 <- Q_tmp[[1]]; Q1_sqe <- Q1_sqe + mean( rowSums( (TQ1[,2:3] - Q1[,2:3])^2 ) )
  TQ2 <- T_tmp[[2]]; Q2 <- Q_tmp[[2]]; Q2_sqe <- Q2_sqe + mean( rowSums( (TQ2[,2:3] - Q2[,2:3])^2 ) )
  TQ3 <- T_tmp[[3]]; Q3 <- Q_tmp[[3]]; Q3_sqe <- Q3_sqe + mean( rowSums( (TQ3[,2:3] - Q3[,2:3])^2 ) )
}
Q1_mse <- Q1_sqe/10;  Q2_mse <- Q2_sqe/10;  Q3_mse <- Q3_sqe/10;  rm(Q1_sqe, Q2_sqe, Q3_sqe)
MSE_T400000h01 <- c(M_mse, Q1_mse, Q2_mse, Q3_mse)
rm(M, Q, M_mse, i, T_tmp, Q_tmp, TQ1, TQ2, TQ3, Q1, Q2, Q3, Q1_mse, Q2_mse, Q3_mse)

# compute the MSEs of the estimated quantiles with T=400000 and h=1.5
M <- T400000h15[[1]];  Q <- T400000h15[[2]]
M_sqe <- rowSums((M[,2:3] - TM[,2:3])^2);  M_mse <- mean(M_sqe);  rm(M_sqe)
Q1_sqe <- 0;  Q2_sqe <- 0;  Q3_sqe <- 0
for (i in 1:length(TQ)){
  T_tmp <- TQ[[i]];  Q_tmp <- Q[[i]]
  TQ1 <- T_tmp[[1]]; Q1 <- Q_tmp[[1]]; Q1_sqe <- Q1_sqe + mean( rowSums( (TQ1[,2:3] - Q1[,2:3])^2 ) )
  TQ2 <- T_tmp[[2]]; Q2 <- Q_tmp[[2]]; Q2_sqe <- Q2_sqe + mean( rowSums( (TQ2[,2:3] - Q2[,2:3])^2 ) )
  TQ3 <- T_tmp[[3]]; Q3 <- Q_tmp[[3]]; Q3_sqe <- Q3_sqe + mean( rowSums( (TQ3[,2:3] - Q3[,2:3])^2 ) )
}
Q1_mse <- Q1_sqe/10;  Q2_mse <- Q2_sqe/10;  Q3_mse <- Q3_sqe/10;  rm(Q1_sqe, Q2_sqe, Q3_sqe)
MSE_T400000h15 <- c(M_mse, Q1_mse, Q2_mse, Q3_mse)
rm(M, Q, M_mse, i, T_tmp, Q_tmp, TQ1, TQ2, TQ3, Q1, Q2, Q3, Q1_mse, Q2_mse, Q3_mse)

# compute the MSEs of the estimated quantiles with T=400000 and h=3.0
M <- T400000h30[[1]];  Q <- T400000h30[[2]]
M_sqe <- rowSums((M[,2:3] - TM[,2:3])^2);  M_mse <- mean(M_sqe);  rm(M_sqe)
Q1_sqe <- 0;  Q2_sqe <- 0;  Q3_sqe <- 0
for (i in 1:length(TQ)){
  T_tmp <- TQ[[i]];  Q_tmp <- Q[[i]]
  TQ1 <- T_tmp[[1]]; Q1 <- Q_tmp[[1]]; Q1_sqe <- Q1_sqe + mean( rowSums( (TQ1[,2:3] - Q1[,2:3])^2 ) )
  TQ2 <- T_tmp[[2]]; Q2 <- Q_tmp[[2]]; Q2_sqe <- Q2_sqe + mean( rowSums( (TQ2[,2:3] - Q2[,2:3])^2 ) )
  TQ3 <- T_tmp[[3]]; Q3 <- Q_tmp[[3]]; Q3_sqe <- Q3_sqe + mean( rowSums( (TQ3[,2:3] - Q3[,2:3])^2 ) )
}
Q1_mse <- Q1_sqe/10;  Q2_mse <- Q2_sqe/10;  Q3_mse <- Q3_sqe/10;  rm(Q1_sqe, Q2_sqe, Q3_sqe)
MSE_T400000h30 <- c(M_mse, Q1_mse, Q2_mse, Q3_mse)
rm(M, Q, M_mse, i, T_tmp, Q_tmp, TQ1, TQ2, TQ3, Q1, Q2, Q3, Q1_mse, Q2_mse, Q3_mse)

MSE_table <- rbind(MSE_T40000h05, MSE_T80000h05, MSE_T400000h05, MSE_T400000h01, MSE_T400000h15, MSE_T400000h30)
colnames(MSE_table) <- c("median", "quantile1", "quantile2", "quantile3")
#library(knitr); kable(MSE_table, format = "latex")
rm(list=ls())



#Appendix Table 2 (case 2):
load("./simulations/appendix/Tables/MSE_case2.RData")
TM <- theoretic[[1]];  TQ <- theoretic[[2]];  rm(theoretic)

# compute the MSEs of the estimated quantiles with T=20000 and h=0.4
M <- T20000h04[[1]];  Q <- T20000h04[[2]]
M_sqe <- rowSums((M[,2:3] - TM[,2:3])^2);  M_mse <- mean(M_sqe);  rm(M_sqe)
Q1_sqe <- 0;  Q2_sqe <- 0;  Q3_sqe <- 0
for (i in 1:length(TQ)){
  T_tmp <- TQ[[i]];  Q_tmp <- Q[[i]]
  TQ1 <- T_tmp[[1]]; Q1 <- Q_tmp[[1]]; Q1_sqe <- Q1_sqe + mean( rowSums( (TQ1[,2:3] - Q1[,2:3])^2 ) )
  TQ2 <- T_tmp[[2]]; Q2 <- Q_tmp[[2]]; Q2_sqe <- Q2_sqe + mean( rowSums( (TQ2[,2:3] - Q2[,2:3])^2 ) )
  TQ3 <- T_tmp[[3]]; Q3 <- Q_tmp[[3]]; Q3_sqe <- Q3_sqe + mean( rowSums( (TQ3[,2:3] - Q3[,2:3])^2 ) )
}
Q1_mse <- Q1_sqe/10;  Q2_mse <- Q2_sqe/10;  Q3_mse <- Q3_sqe/10;  rm(Q1_sqe, Q2_sqe, Q3_sqe)
MSE_T20000h04 <- c(M_mse, Q1_mse, Q2_mse, Q3_mse)
rm(M, Q, M_mse, i, T_tmp, Q_tmp, TQ1, TQ2, TQ3, Q1, Q2, Q3, Q1_mse, Q2_mse, Q3_mse)

# compute the MSEs of the estimated quantiles with T=80000 and h=0.4
M <- T80000h04[[1]];  Q <- T80000h04[[2]]
M_sqe <- rowSums((M[,2:3] - TM[,2:3])^2);  M_mse <- mean(M_sqe);  rm(M_sqe)
Q1_sqe <- 0;  Q2_sqe <- 0;  Q3_sqe <- 0
for (i in 1:length(TQ)){
  T_tmp <- TQ[[i]];  Q_tmp <- Q[[i]]
  TQ1 <- T_tmp[[1]]; Q1 <- Q_tmp[[1]]; Q1_sqe <- Q1_sqe + mean( rowSums( (TQ1[,2:3] - Q1[,2:3])^2 ) )
  TQ2 <- T_tmp[[2]]; Q2 <- Q_tmp[[2]]; Q2_sqe <- Q2_sqe + mean( rowSums( (TQ2[,2:3] - Q2[,2:3])^2 ) )
  TQ3 <- T_tmp[[3]]; Q3 <- Q_tmp[[3]]; Q3_sqe <- Q3_sqe + mean( rowSums( (TQ3[,2:3] - Q3[,2:3])^2 ) )
}
Q1_mse <- Q1_sqe/10;  Q2_mse <- Q2_sqe/10;  Q3_mse <- Q3_sqe/10;  rm(Q1_sqe, Q2_sqe, Q3_sqe)
MSE_T80000h04 <- c(M_mse, Q1_mse, Q2_mse, Q3_mse)
rm(M, Q, M_mse, i, T_tmp, Q_tmp, TQ1, TQ2, TQ3, Q1, Q2, Q3, Q1_mse, Q2_mse, Q3_mse)

# compute the MSEs of the estimated quantiles with T=400000 and h=0.4
M <- T400000h04[[1]];  Q <- T400000h04[[2]]
M_sqe <- rowSums((M[,2:3] - TM[,2:3])^2);  M_mse <- mean(M_sqe);  rm(M_sqe)
Q1_sqe <- 0;  Q2_sqe <- 0;  Q3_sqe <- 0
for (i in 1:length(TQ)){
  T_tmp <- TQ[[i]];  Q_tmp <- Q[[i]]
  TQ1 <- T_tmp[[1]]; Q1 <- Q_tmp[[1]]; Q1_sqe <- Q1_sqe + mean( rowSums( (TQ1[,2:3] - Q1[,2:3])^2 ) )
  TQ2 <- T_tmp[[2]]; Q2 <- Q_tmp[[2]]; Q2_sqe <- Q2_sqe + mean( rowSums( (TQ2[,2:3] - Q2[,2:3])^2 ) )
  TQ3 <- T_tmp[[3]]; Q3 <- Q_tmp[[3]]; Q3_sqe <- Q3_sqe + mean( rowSums( (TQ3[,2:3] - Q3[,2:3])^2 ) )
}
Q1_mse <- Q1_sqe/10;  Q2_mse <- Q2_sqe/10;  Q3_mse <- Q3_sqe/10;  rm(Q1_sqe, Q2_sqe, Q3_sqe)
MSE_T400000h04 <- c(M_mse, Q1_mse, Q2_mse, Q3_mse)
rm(M, Q, M_mse, i, T_tmp, Q_tmp, TQ1, TQ2, TQ3, Q1, Q2, Q3, Q1_mse, Q2_mse, Q3_mse)

# compute the MSEs of the estimated quantiles with T=400000 and h=0.2
M <- T400000h02[[1]];  Q <- T400000h02[[2]]
M_sqe <- rowSums((M[,2:3] - TM[,2:3])^2);  M_mse <- mean(M_sqe);  rm(M_sqe)
Q1_sqe <- 0;  Q2_sqe <- 0;  Q3_sqe <- 0
for (i in 1:length(TQ)){
  T_tmp <- TQ[[i]];  Q_tmp <- Q[[i]]
  TQ1 <- T_tmp[[1]]; Q1 <- Q_tmp[[1]]; Q1_sqe <- Q1_sqe + mean( rowSums( (TQ1[,2:3] - Q1[,2:3])^2 ) )
  TQ2 <- T_tmp[[2]]; Q2 <- Q_tmp[[2]]; Q2_sqe <- Q2_sqe + mean( rowSums( (TQ2[,2:3] - Q2[,2:3])^2 ) )
  TQ3 <- T_tmp[[3]]; Q3 <- Q_tmp[[3]]; Q3_sqe <- Q3_sqe + mean( rowSums( (TQ3[,2:3] - Q3[,2:3])^2 ) )
}
Q1_mse <- Q1_sqe/10;  Q2_mse <- Q2_sqe/10;  Q3_mse <- Q3_sqe/10;  rm(Q1_sqe, Q2_sqe, Q3_sqe)
MSE_T400000h02 <- c(M_mse, Q1_mse, Q2_mse, Q3_mse)
rm(M, Q, M_mse, i, T_tmp, Q_tmp, TQ1, TQ2, TQ3, Q1, Q2, Q3, Q1_mse, Q2_mse, Q3_mse)

# compute the MSEs of the estimated quantiles with T=400000 and h=1.2
M <- T400000h12[[1]];  Q <- T400000h12[[2]]
M_sqe <- rowSums((M[,2:3] - TM[,2:3])^2);  M_mse <- mean(M_sqe);  rm(M_sqe)
Q1_sqe <- 0;  Q2_sqe <- 0;  Q3_sqe <- 0
for (i in 1:length(TQ)){
  T_tmp <- TQ[[i]];  Q_tmp <- Q[[i]]
  TQ1 <- T_tmp[[1]]; Q1 <- Q_tmp[[1]]; Q1_sqe <- Q1_sqe + mean( rowSums( (TQ1[,2:3] - Q1[,2:3])^2 ) )
  TQ2 <- T_tmp[[2]]; Q2 <- Q_tmp[[2]]; Q2_sqe <- Q2_sqe + mean( rowSums( (TQ2[,2:3] - Q2[,2:3])^2 ) )
  TQ3 <- T_tmp[[3]]; Q3 <- Q_tmp[[3]]; Q3_sqe <- Q3_sqe + mean( rowSums( (TQ3[,2:3] - Q3[,2:3])^2 ) )
}
Q1_mse <- Q1_sqe/10;  Q2_mse <- Q2_sqe/10;  Q3_mse <- Q3_sqe/10;  rm(Q1_sqe, Q2_sqe, Q3_sqe)
MSE_T400000h12 <- c(M_mse, Q1_mse, Q2_mse, Q3_mse)
rm(M, Q, M_mse, i, T_tmp, Q_tmp, TQ1, TQ2, TQ3, Q1, Q2, Q3, Q1_mse, Q2_mse, Q3_mse)

# compute the MSEs of the estimated quantiles with T=400000 and h=3.0
M <- T400000h30[[1]];  Q <- T400000h30[[2]]
M_sqe <- rowSums((M[,2:3] - TM[,2:3])^2);  M_mse <- mean(M_sqe);  rm(M_sqe)
Q1_sqe <- 0;  Q2_sqe <- 0;  Q3_sqe <- 0
for (i in 1:length(TQ)){
  T_tmp <- TQ[[i]];  Q_tmp <- Q[[i]]
  TQ1 <- T_tmp[[1]]; Q1 <- Q_tmp[[1]]; Q1_sqe <- Q1_sqe + mean( rowSums( (TQ1[,2:3] - Q1[,2:3])^2 ) )
  TQ2 <- T_tmp[[2]]; Q2 <- Q_tmp[[2]]; Q2_sqe <- Q2_sqe + mean( rowSums( (TQ2[,2:3] - Q2[,2:3])^2 ) )
  TQ3 <- T_tmp[[3]]; Q3 <- Q_tmp[[3]]; Q3_sqe <- Q3_sqe + mean( rowSums( (TQ3[,2:3] - Q3[,2:3])^2 ) )
}
Q1_mse <- Q1_sqe/10;  Q2_mse <- Q2_sqe/10;  Q3_mse <- Q3_sqe/10;  rm(Q1_sqe, Q2_sqe, Q3_sqe)
MSE_T400000h30 <- c(M_mse, Q1_mse, Q2_mse, Q3_mse)
rm(M, Q, M_mse, i, T_tmp, Q_tmp, TQ1, TQ2, TQ3, Q1, Q2, Q3, Q1_mse, Q2_mse, Q3_mse)

MSE_table <- rbind(MSE_T20000h04, MSE_T80000h04, MSE_T400000h04, MSE_T400000h02, MSE_T400000h12, MSE_T400000h30)
colnames(MSE_table) <- c("median", "quantile1", "quantile2", "quantile3")
#library(knitr); kable(MSE_table, format = "latex")
rm(list=ls())




#Appendix Table 3 (case 3):
load("./simulations/appendix/Tables/MSE_case3.RData")
TM <- theoretic[[1]];  TQ <- theoretic[[2]];  rm(theoretic)

# compute the MSEs of the estimated quantiles with T=20000 and h=0.1
M <- T20000h01[[1]];  Q <- T20000h01[[2]]
M_sqe <- rowSums((M[,2:3] - TM[,2:3])^2);  M_mse <- mean(M_sqe);  rm(M_sqe)
Q1_sqe <- 0;  Q2_sqe <- 0;  Q3_sqe <- 0
for (i in 1:length(TQ)){
  T_tmp <- TQ[[i]];  Q_tmp <- Q[[i]]
  TQ1 <- T_tmp[[1]]; Q1 <- Q_tmp[[1]]; Q1_sqe <- Q1_sqe + mean( rowSums( (TQ1[,2:3] - Q1[,2:3])^2 ) )
  TQ2 <- T_tmp[[2]]; Q2 <- Q_tmp[[2]]; Q2_sqe <- Q2_sqe + mean( rowSums( (TQ2[,2:3] - Q2[,2:3])^2 ) )
  TQ3 <- T_tmp[[3]]; Q3 <- Q_tmp[[3]]; Q3_sqe <- Q3_sqe + mean( rowSums( (TQ3[,2:3] - Q3[,2:3])^2 ) )
}
Q1_mse <- Q1_sqe/10;  Q2_mse <- Q2_sqe/10;  Q3_mse <- Q3_sqe/10;  rm(Q1_sqe, Q2_sqe, Q3_sqe)
MSE_T20000h01 <- c(M_mse, Q1_mse, Q2_mse, Q3_mse)
rm(M, Q, M_mse, i, T_tmp, Q_tmp, TQ1, TQ2, TQ3, Q1, Q2, Q3, Q1_mse, Q2_mse, Q3_mse)

# compute the MSEs of the estimated quantiles with T=40000 and h=0.1
M <- T40000h01[[1]];  Q <- T40000h01[[2]]
M_sqe <- rowSums((M[,2:3] - TM[,2:3])^2);  M_mse <- mean(M_sqe);  rm(M_sqe)
Q1_sqe <- 0;  Q2_sqe <- 0;  Q3_sqe <- 0
for (i in 1:length(TQ)){
  T_tmp <- TQ[[i]];  Q_tmp <- Q[[i]]
  TQ1 <- T_tmp[[1]]; Q1 <- Q_tmp[[1]]; Q1_sqe <- Q1_sqe + mean( rowSums( (TQ1[,2:3] - Q1[,2:3])^2 ) )
  TQ2 <- T_tmp[[2]]; Q2 <- Q_tmp[[2]]; Q2_sqe <- Q2_sqe + mean( rowSums( (TQ2[,2:3] - Q2[,2:3])^2 ) )
  TQ3 <- T_tmp[[3]]; Q3 <- Q_tmp[[3]]; Q3_sqe <- Q3_sqe + mean( rowSums( (TQ3[,2:3] - Q3[,2:3])^2 ) )
}
Q1_mse <- Q1_sqe/10;  Q2_mse <- Q2_sqe/10;  Q3_mse <- Q3_sqe/10;  rm(Q1_sqe, Q2_sqe, Q3_sqe)
MSE_T40000h01 <- c(M_mse, Q1_mse, Q2_mse, Q3_mse)
rm(M, Q, M_mse, i, T_tmp, Q_tmp, TQ1, TQ2, TQ3, Q1, Q2, Q3, Q1_mse, Q2_mse, Q3_mse)

# compute the MSEs of the estimated quantiles with T=400000 and h=0.1
M <- T400000h01[[1]];  Q <- T400000h01[[2]]
M_sqe <- rowSums((M[,2:3] - TM[,2:3])^2);  M_mse <- mean(M_sqe);  rm(M_sqe)
Q1_sqe <- 0;  Q2_sqe <- 0;  Q3_sqe <- 0
for (i in 1:length(TQ)){
  T_tmp <- TQ[[i]];  Q_tmp <- Q[[i]]
  TQ1 <- T_tmp[[1]]; Q1 <- Q_tmp[[1]]; Q1_sqe <- Q1_sqe + mean( rowSums( (TQ1[,2:3] - Q1[,2:3])^2 ) )
  TQ2 <- T_tmp[[2]]; Q2 <- Q_tmp[[2]]; Q2_sqe <- Q2_sqe + mean( rowSums( (TQ2[,2:3] - Q2[,2:3])^2 ) )
  TQ3 <- T_tmp[[3]]; Q3 <- Q_tmp[[3]]; Q3_sqe <- Q3_sqe + mean( rowSums( (TQ3[,2:3] - Q3[,2:3])^2 ) )
}
Q1_mse <- Q1_sqe/10;  Q2_mse <- Q2_sqe/10;  Q3_mse <- Q3_sqe/10;  rm(Q1_sqe, Q2_sqe, Q3_sqe)
MSE_T400000h01 <- c(M_mse, Q1_mse, Q2_mse, Q3_mse)
rm(M, Q, M_mse, i, T_tmp, Q_tmp, TQ1, TQ2, TQ3, Q1, Q2, Q3, Q1_mse, Q2_mse, Q3_mse)

# compute the MSEs of the estimated quantiles with T=400000 and h=0.03
M <- T400000h003[[1]];  Q <- T400000h003[[2]]
M_sqe <- rowSums((M[,2:3] - TM[,2:3])^2);  M_mse <- mean(M_sqe);  rm(M_sqe)
Q1_sqe <- 0;  Q2_sqe <- 0;  Q3_sqe <- 0
for (i in 1:length(TQ)){
  T_tmp <- TQ[[i]];  Q_tmp <- Q[[i]]
  TQ1 <- T_tmp[[1]]; Q1 <- Q_tmp[[1]]; Q1_sqe <- Q1_sqe + mean( rowSums( (TQ1[,2:3] - Q1[,2:3])^2 ) )
  TQ2 <- T_tmp[[2]]; Q2 <- Q_tmp[[2]]; Q2_sqe <- Q2_sqe + mean( rowSums( (TQ2[,2:3] - Q2[,2:3])^2 ) )
  TQ3 <- T_tmp[[3]]; Q3 <- Q_tmp[[3]]; Q3_sqe <- Q3_sqe + mean( rowSums( (TQ3[,2:3] - Q3[,2:3])^2 ) )
}
Q1_mse <- Q1_sqe/10;  Q2_mse <- Q2_sqe/10;  Q3_mse <- Q3_sqe/10;  rm(Q1_sqe, Q2_sqe, Q3_sqe)
MSE_T400000h003 <- c(M_mse, Q1_mse, Q2_mse, Q3_mse)
rm(M, Q, M_mse, i, T_tmp, Q_tmp, TQ1, TQ2, TQ3, Q1, Q2, Q3, Q1_mse, Q2_mse, Q3_mse)

# compute the MSEs of the estimated quantiles with T=400000 and h=1.0
M <- T400000h10[[1]];  Q <- T400000h10[[2]]
M_sqe <- rowSums((M[,2:3] - TM[,2:3])^2);  M_mse <- mean(M_sqe);  rm(M_sqe)
Q1_sqe <- 0;  Q2_sqe <- 0;  Q3_sqe <- 0
for (i in 1:length(TQ)){
  T_tmp <- TQ[[i]];  Q_tmp <- Q[[i]]
  TQ1 <- T_tmp[[1]]; Q1 <- Q_tmp[[1]]; Q1_sqe <- Q1_sqe + mean( rowSums( (TQ1[,2:3] - Q1[,2:3])^2 ) )
  TQ2 <- T_tmp[[2]]; Q2 <- Q_tmp[[2]]; Q2_sqe <- Q2_sqe + mean( rowSums( (TQ2[,2:3] - Q2[,2:3])^2 ) )
  TQ3 <- T_tmp[[3]]; Q3 <- Q_tmp[[3]]; Q3_sqe <- Q3_sqe + mean( rowSums( (TQ3[,2:3] - Q3[,2:3])^2 ) )
}
Q1_mse <- Q1_sqe/10;  Q2_mse <- Q2_sqe/10;  Q3_mse <- Q3_sqe/10;  rm(Q1_sqe, Q2_sqe, Q3_sqe)
MSE_T400000h10 <- c(M_mse, Q1_mse, Q2_mse, Q3_mse)
rm(M, Q, M_mse, i, T_tmp, Q_tmp, TQ1, TQ2, TQ3, Q1, Q2, Q3, Q1_mse, Q2_mse, Q3_mse)

# compute the MSEs of the estimated quantiles with T=400000 and h=2.0
M <- T400000h20[[1]];  Q <- T400000h20[[2]]
M_sqe <- rowSums((M[,2:3] - TM[,2:3])^2);  M_mse <- mean(M_sqe);  rm(M_sqe)
Q1_sqe <- 0;  Q2_sqe <- 0;  Q3_sqe <- 0
for (i in 1:length(TQ)){
  T_tmp <- TQ[[i]];  Q_tmp <- Q[[i]]
  TQ1 <- T_tmp[[1]]; Q1 <- Q_tmp[[1]]; Q1_sqe <- Q1_sqe + mean( rowSums( (TQ1[,2:3] - Q1[,2:3])^2 ) )
  TQ2 <- T_tmp[[2]]; Q2 <- Q_tmp[[2]]; Q2_sqe <- Q2_sqe + mean( rowSums( (TQ2[,2:3] - Q2[,2:3])^2 ) )
  TQ3 <- T_tmp[[3]]; Q3 <- Q_tmp[[3]]; Q3_sqe <- Q3_sqe + mean( rowSums( (TQ3[,2:3] - Q3[,2:3])^2 ) )
}
Q1_mse <- Q1_sqe/10;  Q2_mse <- Q2_sqe/10;  Q3_mse <- Q3_sqe/10;  rm(Q1_sqe, Q2_sqe, Q3_sqe)
MSE_T400000h20 <- c(M_mse, Q1_mse, Q2_mse, Q3_mse)
rm(M, Q, M_mse, i, T_tmp, Q_tmp, TQ1, TQ2, TQ3, Q1, Q2, Q3, Q1_mse, Q2_mse, Q3_mse)

MSE_table <- rbind(MSE_T20000h01, MSE_T40000h01, MSE_T400000h01, MSE_T400000h003, MSE_T400000h10, MSE_T400000h20)
colnames(MSE_table) <- c("median", "quantile1", "quantile2", "quantile3")
#library(knitr); kable(MSE_table, format = "latex")
rm(list=ls())




