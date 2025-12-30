#!/usr/bin/env Rscript
.libPaths(c("/work/xueheng/mgrantham/r_lib", .libPaths()))

# ----- knobs to tweak -----
PHI   <- 0.0     # set phi (e.g., 0, 0.3, -0.3, 0.7, -0.7)
N     <- 1000
SIMS  <- 20
ITERS <- 11

###############
##   Setup
###############

suppressPackageStartupMessages({
  library(genlasso)
  library(dplyr)
  library(Matrix)      # Diagonal()
  library(breakfast)   # WCM
  library(AR1seg)      # AR1seg
  library(changepoint) # PELT, BS
  library(wbs)
})

# Define function to create mean vector for BIC calculation
find_means = function(y,changepoints){
  changepoints = c(0, changepoints, length(y))
  
  # Initialize the mean vector
  mean_vector = numeric(length(y))
  
  # Loop through each segment defined by the changepoints
  for (i in 1:(length(changepoints) - 1)) {
    start = changepoints[i] + 1  # Start index of the segment
    end = changepoints[i + 1]    # End index of the segment
    segment_mean = mean(y[start:end])  # Compute the mean of the segment
    mean_vector[start:end] = segment_mean  # Assign mean to the segment
  }
  
  return(mean_vector)
}

# Define predict function given beta and x
pred_fcn = function(beta,y){
  phi = beta[1]
  pred = c()
  pred[1] = beta[2]
  for(j in 2:length(beta)){
    pred[j] = phi*y[j-1] + beta[j]
  }
  return(pred)
}

# Define BIC function
bic_fcn = function(beta,y){
  n = length(y)
  bic = c()
  numberknots = c()
  for(k in 1:ncol(beta)){
    nknots = sum(abs(diff(beta[-1,k],differences=1))>1e-4)
    knots = which(abs(diff(beta[-1,k],differences=1))>1e-4)
    means = find_means(y,knots)
    phi = beta[1,k]
    pred = c()
    pred[1] = beta[2,k]
    for(j in 2:length(beta[,k])){
      pred[j] = phi*y[j-1] + means[j]
    }
    ss = sum((pred-y)^2)
    bic[k] = n/2*log(ss/n) + (nknots+1)*log(n)
    numberknots[k] = length(knots)
  }
  return(list(bic=bic,knots=numberknots))
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
  n/2*log(sum((yhat_pelt(x)-y)^2)/n) + (length(x@cpts)-1)*log(n)
}

# Define BIC penalty for WBS
bic_wbs = function(x,y,n){
  n/2*log(sum((yhat_wbs(x)-y)^2)/n) + 
    (length(changepoints(x)$cpt.ic$bic.penalty)+1)*log(n)
}

# Define function to predict using PELT
generate_meanvec = function(tau,means){
  out = rep(means[1],tau[2]-tau[1])
  if(length(means)==1){
    return(out)
  }
  else{
    for(k in 2:length(means)){
      out = c(out,
              rep(means[k],tau[k+1]-tau[k]))
    } 
  }
  return(out)
}

# Define function to predict using WCM
yhat_wcm = function(y){
  sol = sol.wcm(y,M=1000, Q=50)
  mod = model.gsa(sol,p.max=1,pen = log(length(sol$x)))
  meanvec = mod$est
  data = y-meanvec
  
  # Set up YW
  y_lag = data[-n] 
  y_curr = data[-1]
  
  # Estimate phi using MLE formula for AR(1)
  phi = sum(y_lag * y_curr) / sum(y_lag^2)
  
  # Estimate y using this value of phi
  pred = c()
  pred[1] = meanvec[1]
  for(j in 2:length(meanvec)){
    pred[j] = phi*y[j-1] + meanvec[j]
  }
  return(list(phi=phi,yhat=pred))
}

# Define BIC penalty for WCM
bic_wcm = function(y){
  
  sol = sol.wcm(y,M=1000, Q=50)
  mod = model.gsa(sol,p.max=1,pen = log(length(sol$x)))
  n = length(mod$est)
  
  
  n/2*log(sum((yhat_wcm(y)$yhat-y)^2)/n) + 
    (length(mod$cpts)+1)*log(n)
}

# Define predict function for AR1seg given y
yhat_ar1seg = function(y){
  res = AR1seg_func(y)
  a = c(1,res$PPSelectedBreaks[1:(res$PPselected-1)]+1)
  b = res$PPSelectedBreaks[1:(res$PPselected)]
  Bounds = cbind(a,b)
  if (res$PPselected != 1) {
    mean = rep(res$PPmean, Bounds[,2] - Bounds[,1] + 1)
  } else {
    mean = rep(mean(y),length(y))
  }
  phi = res$rho
  pred = c()
  pred[1] = mean[1]
  for(j in 2:length(mean)){
    pred[j] = phi*y[j-1] + mean[j]
  }
  return(pred)
}

# Define BIC penalty for AR1seg
bic_ar1seg = function(y){
  res=AR1seg_func(y)
  n = length(y)
  n/2*log(sum((yhat_ar1seg(y)-y)^2)/n) + 
    (res$PPselected)*log(n)
}

# --- main ---

n = N
true_model = rep(c(0,2,0,2),each=250) 

# Include number of changepoints
cnames = c("true_model","true_phi","true_nknots","bic_true_model",
           "flasso","phi_flasso","phi_yw_flasso","nknots_flasso","bic_flasso",
           "adaflasso","phi_adaflasso","phi_yw_adaflasso",
           "nknots_adaflasso","bic_adaflasso",
           paste(c("iter","phi_iter","phi_yw_iter","nknots_iter","bic_iter"),
                 rep(2:(ITERS-1),each=5),"flasso",sep=""),
           "knots_PELT","nknots_PELT","bic_PELT",
           "knots_BS","nknots_BS","bic_BS",
           "knots_WBS","nknots_WBS","bic_WBS",
           "knots_WCM","nknots_WCM","bic_WCM","phi_WCM",
           "knots_AR1seg","nknots_AR1seg","bic_AR1seg","phi_AR1seg")
out_df = as.data.frame(matrix(NA,SIMS,length(cnames)))
colnames(out_df) <- cnames

for(i in 1:SIMS){
  y = rnorm(1)
  phi = PHI
  for(j in 2:n){
    y[j] = phi*y[j-1]+rnorm(1)
  }
  y = y + true_model
  yhat = c()
  yhat[1] = true_model[1]
  for(k in 2:length(y)){
    yhat[k] = phi*y[k-1] + true_model[k]
  }
  true_bic = n*log(sum((yhat-y)^2)/n)+4*log(n)
  out = c("251,501,751",phi,3,true_bic)    
  
  # Set design matrix
  lessalpha = rbind( c(1,rep(0,n-2)),diag(n-1))
  X = cbind(c(0,y[1:(n-1)]),lessalpha)
  
  # build once, outside the loop
  base_Dprime <- diff(diag(n - 1), differences = 1)  

  for (k in seq_len(ITERS)) {
    if (k == 1) {
      Dprime <- base_Dprime
    } else {
      delta <- diff(beta_hat$beta[-1, loc], differences = 1)  # length n-2
      wt    <- 1 / (abs(delta)+ 1e-16)                        # adaptive weights
      Dprime <- diag(wt) %*% base_Dprime                      # row-scale via diag
    }
    
    D <- cbind(0, Dprime)  # keep phi unpenalized
    
    # Run Genlasso
    fit = genlasso(y, X=X, D=D,maxsteps=500)
    beta_hat = coef(fit)
    
    # Find BIC for all models
    bic_find = bic_fcn(beta_hat$beta,y)
    
    # Best model by BIC
    bic = bic_find$bic
    loc = which.min(bic)
    minbic = bic[loc]
    
    # Location of knots found by best model
    knots = paste(which(abs(diff(beta_hat$beta[-1,loc],
                                 differences=1))>1e-4)+2,
                  collapse=",")
    
    knot_locs = which(abs(diff(beta_hat$beta[-1,loc],
                               differences=1))>1e-4)+2
    
    numberknots = length(which(abs(diff(beta_hat$beta[-1,loc],
                                        differences=1))>1e-4))
    phi = as.numeric(beta_hat$beta[1,loc])
    
    # Re-estimate phi using YW after de-meaning
    tau = as.numeric(c(0,as.numeric(strsplit(knots,split=",")[[1]])-1,1000))
    
    # Estimate means using OLS design matrix X
    dm = matrix(NA,nrow=1000,ncol=length(tau)-1)
    for(id in 2:length(tau)){
      dm[,id-1] = c(rep(0,tau[id-1]),
                    rep(1,tau[id]-tau[id-1]),
                    rep(0,1000-tau[id]))
    }
    
    # Initial estimation of means
    mu_hat = solve(t(dm)%*%dm)%*%t(dm)%*%y
    
    # Generate vector of means by length in order to de-mean
    meanvec = generate_meanvec(tau,mu_hat)
    
    # De-mean the data
    data = y-meanvec
    
    # Set up YW
    y_lag = data[-n] 
    y_curr = data[-1]
    
    # Estimate phi using MLE formula for AR(1)
    phi_yw = sum(y_lag * y_curr) / sum(y_lag^2)
    
    out = c(out,knots,phi,phi_yw,numberknots,minbic)
    cat(sprintf("Simulation:%d; Iteration:%d; Changepoints:%d; Locations:%s; Phi:%g; Phi_YW:%g\n",
                i, k, numberknots,
                if (numberknots) knots else "none",
                round(phi, 4), round(phi_yw, 4)))
    
    if (interactive()){
      yhat_irfl <- pred_fcn(beta_hat$beta[, loc], y)
      plot(y, type = "l", main = paste("IRFL Iteration", k), 
           xlab = "Index", 
           ylab = "Value")
      lines(yhat_irfl, col = 2, lwd = 2)
      if (numberknots) abline(v = knot_locs, lty = 3)
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
  wbs_bic = bic_wbs(wbs,y,1000)
  wbs_cp = sort(changepoints(wbs)$cpt.ic$bic.penalty)
  wbs_nknots = length(wbs_cp)
  wbs_knots = paste0(sort(changepoints(wbs)$cpt.ic$bic.penalty+1),collapse=",")
  
  # Now add WCM
  sol = sol.wcm(y,M=1000, Q=50)
  wcm = model.gsa(sol,p.max=1,pen = log(length(sol$x)))
  wcm_cp = wcm$cpts+1
  wcm_nknots = length(wcm$cpts)
  wcm_bic = bic_wcm(y)
  wcm_knots = paste0(wcm_cp,collapse=",")
  wcm_phi = yhat_wcm(y)$phi
  
  # Now add AR1seg
  res = AR1seg_func(y)
  ar1seg_cp = as.vector(res$PPSelectedBreaks[-res$PPselected])+1
  ar1seg_nknots = res$PPselected-1
  ar1seg_bic = bic_ar1seg(y)
  ar1seg_knots = paste0(ar1seg_cp,collapse=",")
  ar1seg_phi = res$rho
  
  # Update the dataframe
  out = c(out,
          pelt_knots,pelt_nknots,pelt_bic,
          bs_knots,bs_nknots,bs_bic,
          wbs_knots,wbs_nknots,wbs_bic,
          wcm_knots,wcm_nknots,wcm_bic,wcm_phi,
          ar1seg_knots,ar1seg_nknots,ar1seg_bic,ar1seg_phi)
  
  # Write to out_df
  out_df[i,] = out
}

write.csv(out_df,"out.csv", row.names=FALSE)

