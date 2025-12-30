
library(genlasso)
library(dplyr)
library(changepoint)
library(wbs)

# Define predict function given beta and x
pred_fcn <- function(beta, X) {
  if (is.matrix(beta) && ncol(beta) == 1L) {
    return(as.numeric(X %*% beta[, 1L, drop = FALSE]))
  }
  if (is.matrix(beta) && ncol(beta) > 1L) {
    stop("pred_fcn expects a coefficient vector or a single-column matrix.")
  }
  as.numeric(X %*% beta)
}

# Define BIC function
bic_fcn <- function(beta, y, X, tol = 1e-5) {
  n  <- length(y)
  L  <- ncol(beta)
  bic <- numeric(L)
  numberknots <- integer(L)
  
  for (k in seq_len(L)) {
    # changepoints from mean block only (exclude slope in row 1)
    knots <- sum(abs(diff(beta[-1L, k], differences = 1L)) > tol)
    
    # predictions via design matrix (prevents length/recycling issues)
    pred <- as.numeric(X %*% beta[, k, drop = FALSE])
    
    ss <- sum((y - pred)^2)
    # df = slope(1) + segments(K+1) = K + 2
    bic[k] <- n * log(ss / n) + (knots + 2L) * log(n)
    numberknots[k] <- knots
  }
  list(bic = bic, knots = numberknots)
}

# Define function to predict using PELT
yhat_pelt = function(x){
  out = rep(x@param.est$mean[1],x@cpts[1])
  if(length(x@cpts)==1){
    return(out)
  }
  else{
    for(k in 2:length(x@param.est$mean)){
      out = c(out,
              rep(x@param.est$mean[k],x@cpts[k]-x@cpts[k-1]))
    } 
  }
  return(out)
}

# Define function to predict using WBS
yhat_wbs = function(x){
  cpts = c(as.numeric(sort(changepoints(x)$cpt.ic$bic.penalty)),length(y)) 
  out = rep(mean(y[1:cpts[1]]),cpts[1])
  if(length(cpts)==1){
    return(out)
  }
  else{
    for(k in 2:(length(cpts))){
      out = c(out,
              rep(mean(y[(cpts[k-1]+1):cpts[k]]),cpts[k]-cpts[k-1]))
    }
  }
  return(out)
}

# Define BIC penalty for PELT, BS
bic_pelt = function(x,y,n){
  n*log(sum((yhat_pelt(x)-y)^2)/n) + (2*length(x@cpts)-1)*log(n)
}

# Define BIC penalty for WBS
bic_wbs = function(x,y,n){
  n*log(sum((yhat_wbs(x)-y)^2)/n) + 
    (2*(length(changepoints(x)$cpt.ic$bic.penalty)+1)-1)*log(n)
}

#############
## Adaptive / Iteratively Reweighted Fused Lasso
#############

# Repeat 20 iterations on 50 pc's
n = 1000
x = 1:n
true_model = rep(c(0,2,0,2),each=250) + 0.005*x
# Include number of changepoints
out_df = as.data.frame(matrix(NA,20,57))
cnames = c("true_model","true_slope","true_nknots","bic_true_model",
           "flasso","slope_flasso","nknots_flasso","bic_flasso",
           "adaflasso","slope_adaflasso","nknots_adaflasso","bic_adaflasso",
           paste(c("iter","slope_iter","nknots_iter","bic_iter"),
                 rep(2:10,each=4),"flasso",sep=""),
           "PELT","nknots_PELT","bic_PELT","knots_BS","nknots_BS","bic_BS",
           "knots_WBS","nknots_WBS","bic_WBS")
colnames(out_df)=cnames

# Set design matrix
Xstar = cbind(x,diag(n))
X = Xstar[,-2]
X[1,2] = 1

if(interactive()){
  X[1:5,1:5]
}

D0_mean <- diff(diag(n-1), differences = 1)      # (n-1) x n
D0      <- cbind(0, D0_mean)                   # (n-1) x (n+1); zero column for slope

for(i in 1:20){
  y = rnorm(rep(0,n)) + 0.005*x
  true_bic = n*log(sum((true_model-y)^2)/n)+8*log(n)
  out = c("251,501,751",.005,3,true_bic)    # Repeat 11 iteratively reweighted estimations
  for(k in 1:11){
    
    # Set initial D
    if (k == 1) {
      D <- D0
    } else {
      wt <- abs(diff(beta_hat$beta[-1, loc], differences = 1)) + 1e-16   # length n-1 (on means only)
      D <- cbind(0, diag(1 / wt) %*% D0_mean)    # row-scale the mean part only
    }
    
    # Run Genlasso
    fit = genlasso(y, X=X, D=D,maxsteps=25)
    beta_hat = coef(fit)
    
    # Find BIC for all models
    bic_find = bic_fcn(beta_hat$beta,y,X)
    
    # Best model by BIC
    bic = bic_find$bic
    loc = which.min(bic)
    minbic = bic[loc]
    
    # Location of knots found by best model
    knots = paste(which(abs(diff(beta_hat$beta[-1,loc],
                                 differences=1))>1e-9)+2,
                  collapse=",")
    numberknots = length(which(abs(diff(beta_hat$beta[-1,loc],
                                        differences=1))>1e-9))
    slope = as.numeric(beta_hat$beta[1,loc])
    
    out = c(out,knots,slope,numberknots,minbic)
    
    cat(sprintf("Simulation:%d; Iteration:%d; Changepoints:%d; Locations:%s; Slope:%g\n",
                i, k, numberknots,
                if (numberknots) knots else "none",
                round(slope, 4)))
    
    if (interactive()) {
      # changepoints from mean block only (exclude slope)
      knot_locs <- which(abs(diff(beta_hat$beta[-1, loc], differences = 1L)) > 1e-9) + 2L
      
      # fitted values via matrix multiplication
      yhat_irfl <- pred_fcn(beta_hat$beta[, loc, drop = FALSE], X)
      
      plot(y, type = "l", main = sprintf("IRFL Iteration %d", k),
           xlab = "Index", ylab = "Value")
      lines(yhat_irfl, lwd = 2, col = 2)
      if (length(knot_locs)) abline(v = knot_locs, lty = 3)
    }
  }
  
  # Now add PELT 
  cpt = cpt.mean(y, method="PELT",penalty="BIC")
  pelt_bic = bic_pelt(cpt,y,1000)
  pelt_nknots = length(cpt@cpts)-1
  pelt_knots = paste0(cpt@cpts[-length(cpt@cpts)]+1,collapse=",")
  
  # Now add BS
  bs = cpt.mean(y, method="BinSeg",penalty="BIC",Q=25)
  bs_bic = bic_pelt(bs,y,1000)
  bs_nknots = length(bs@cpts)-1
  bs_knots = paste0(bs@cpts[-length(bs@cpts)]+1,collapse=",")
  
  # Now add WBS
  wbs = wbs(y)
  wbs_cp = sort(changepoints(wbs)$cpt.ic$bic.penalty)
  wbs_bic = bic_wbs(wbs,y,1000)
  wbs_nknots = length(wbs_cp)
  wbs_knots = paste0(sort(changepoints(wbs)$cpt.ic$bic.penalty+1),collapse=",")
  
  # Update the dataframe
  out = c(out,
          pelt_knots,pelt_nknots,pelt_bic,
          bs_knots,bs_nknots,bs_bic,
          wbs_knots,wbs_nknots,wbs_bic)
  
  # Write to out_df
  out_df[i,] = out
}
write.csv(out_df,"out.csv")