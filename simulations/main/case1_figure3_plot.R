
library(ellipse)
load("./simulations/main/case1_conds.RData")
load("./simulations/main/case1_predictive.RData")
load("./simulations/main/case1_VAR_predictive.RData")
load("./simulations/main/case1_stationary.RData")
load("./simulations/main/case1_VAR_stationary.RData")
load("./simulations/main/plot_case1_figure3.RData")

medians <- prediction[[1]]; quantiles <- prediction[[2]]

library(plotrix)
layout(matrix(c(1,2,3,1,2,3,4,5,6,4,5,6,7,8,9,7,8,9,10,10,10), 7, 3, byrow = TRUE))
par(mai=c(0.2, 0.3, 0.2, 0.2))

for (i in 1:9){
  VAR_tmp <- matrix(0, nrow = dim(VAR_centers)[1], ncol = dim(VAR_centers)[2])
  if (i < 5){
    medi <- medians[i,2:3]; Qtmp <- quantiles[[i]]
    quant1 <- Qtmp[[1]]; quant2 <- Qtmp[[2]]; quant3 <- Qtmp[[3]]
    plot(quant2[,2:3], type = "l", col = "orange", lwd = 2, xlim = c(-2.5, 1.5), ylim = c(0, 4), asp=1,  xaxt ="n", yaxt = "n", xlab = "", ylab = "") 
    xticks <- c(-2, -1, 0, 1, 2);  axis(1, at = xticks, labels = round(xticks, 0))
    yticks <- c(0, 1, 2, 3, 4);  axis(2, at = yticks, labels = round(yticks, 0))
    lines(quant3[,2:3], col = 'orange', lwd = 2)
    points(medi[1], medi[2], col = "darkred", pch = 16, cex=1)   
    
    VAR_tmp[i, ] <- VAR_centers [i, ];  VAR_sig <- VAR_SIGMs[[i]]
    lines(ellipse(VAR_sig, centre = VAR_tmp[i,], level = 0.4), type="l", col="lightslateblue")
    lines(ellipse(VAR_sig, centre = VAR_tmp[i,], level = 0.8), type="l", col="lightslateblue")
    points(VAR_tmp[i,1], VAR_tmp[i,2], col = "blue", pch = 12, cex=1)
    rm(medi, Qtmp, quant1, quant2, quant3, VAR_tmp, VAR_sig)
  }else if (i > 5){
    medi <- medians[i-1,2:3]; Qtmp <- quantiles[[i-1]]
    quant1 <- Qtmp[[1]]; quant2 <- Qtmp[[2]]; quant3 <- Qtmp[[3]]
    plot(quant2[,2:3], type = "l", col = "orange", lwd = 2, xlim = c(-2.5, 1.5), ylim = c(0, 4),  asp=1,  xaxt ="n", yaxt = "n", xlab = "", ylab = "") 
    xticks <- c(-2, -1, 0, 1, 2);  axis(1, at = xticks, labels = round(xticks, 0))
    yticks <- c(0, 1, 2, 3, 4);  axis(2, at = yticks, labels = round(yticks, 0))
    lines(quant3[,2:3], col = 'orange', lwd = 2)
    points(medi[1], medi[2], col = "darkred", pch = 16, cex=1)
    
    VAR_tmp[i-1, ] <- VAR_centers[i-1, ];   VAR_sig <- VAR_SIGMs[[i-1]]
    lines(ellipse(VAR_sig, centre = VAR_tmp[i-1,], level = 0.4), type="l", col="lightslateblue")
    lines(ellipse(VAR_sig, centre = VAR_tmp[i-1,], level = 0.8), type="l", col="lightslateblue")
    points(VAR_tmp[i-1,1], VAR_tmp[i-1,2], col = "blue", pch = 12, cex=1)
    rm(medi, Qtmp, quant1, quant2, quant3, VAR_tmp, VAR_sig)
  }else {
    #M <- stationary[[1]] # M is the sample mean of the time series
    medi <- stationary[[2]]  
    quant1 <- stationary[[3]]; quant2 <- stationary[[4]]; quant3 <- stationary[[5]]
    plot(quant2[,1:2], type = "l", col = 'orange', lwd = 2, xlim = c(-2.5, 1.5), ylim = c(0, 4), xaxt ="n", yaxt = "n",  asp=1, xlab = "", ylab = "") 
    xticks <- c(-2, -1, 0, 1, 2);  axis(1, at = xticks, labels = round(xticks, 0))
    yticks <- c(0, 1, 2, 3, 4);  axis(2, at = yticks, labels = round(yticks, 0))
    lines(quant3[,1:2], col = 'orange', lwd = 2)
    points(medi[1], medi[2], col = "darkred", pch = 16, cex=1)
    
    mu_hat <- mu_hat; Sigma_hat <- Sigma_hat
    lines(ellipse(Sigma_hat, centre = mu_hat, level = 0.4), type="l", col="lightslateblue")
    lines(ellipse(Sigma_hat, centre = mu_hat, level = 0.8), type="l", col="lightslateblue")
    points(mu_hat[1], mu_hat[2], col = "blue", pch = 12, cex=1)
    
    for (k in 1:8){ points(conds[k,1], conds[k,2], col = "seagreen3", pch = 17, cex=1.2) }
  }
}


par(mai=c(0,0,0.2,0))
plot.new()
legend(x="top", bty="n", c("predictive center-outward median", "VAR predicted mean", "center-outward predictive quantiles", 
                                "Gaussian elliptic quantiles"),
       pch = c(16, 12, NA, NA), lty=c(0,0,1,1),
       pt.bg=c("darkred", "blue", "orange", "lightslateblue"),
       col = c("darkred", "blue", "orange", "lightslateblue"), text.font=2, ncol = 2, x.intersp = 0.4, inset = c(0, 0)) 

#dev.off()





