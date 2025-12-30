library(genlasso)
library(dplyr)
library(changepoint)
library(wbs)

# Define BIC function
bic_fcn <- function(beta, y){
  n <- length(y)
  p <- ncol(beta)
  bic <- numeric(p)
  numberknots <- integer(p)
  for(k in seq_len(p)){
    knots <- sum(abs(diff(beta[, k], 1)) > 1e-5)
    ss <- sum((beta[, k] - y)^2)
    bic[k] <- n * log( max(ss/n, .Machine$double.eps) ) + (2*knots + 1) * log(n)
    numberknots[k] <- knots
  }
  list(bic = bic, knots = numberknots)
}

# Define function to predict using PELT
yhat_pelt = function(x){
  # Robust to zero changepoints and accessor differences
  cp <- tryCatch(cpts(x), error = function(e) integer(0))
  mu <- tryCatch(param.est(x)$mean, error = function(e) numeric(0))
  
  # If no changepoints are reported, fall back to one constant segment
  if (length(cp) == 0L) {
    n <- tryCatch(length(x@data.set), error = function(e) 0L)
    m <- if (length(mu) >= 1L) mu[1L] else 0
    return(rep(m, n))
  }
  
  # Ensure terminal index present; if not, assume last cp is terminal
  cp <- sort(as.integer(cp))
  seg_len <- diff(c(0L, cp))
  
  # Guard against any nonpositive lengths (rare but safer)
  if (any(seg_len <= 0L)) {
    n <- tryCatch(length(x@data.set), error = function(e) sum(pmax(seg_len, 0L)))
    m <- if (length(mu) >= 1L) mu[1L] else 0
    return(rep(m, n))
  }
  
  # Align means with segments (some versions already match; otherwise recycle last)
  if (length(mu) < length(seg_len)) {
    mu <- c(mu, rep(tail(mu, 1L), length(seg_len) - length(mu)))
  }
  rep(mu[seq_along(seg_len)], times = seg_len)
}

# Define function to predict using WBS
yhat_wbs <- function(x, y){
  n <- length(y)
  cp <- tryCatch(sort(as.numeric(changepoints(x)$cpt.ic$bic.penalty)),
                 error = function(e) numeric(0))
  
  # If no changepoints, constant fit
  if (length(cp) == 0L) return(rep(mean(y), n))
  
  # Append terminal index and build segmentwise means safely
  cp <- c(cp, n)
  out <- rep(mean(y[1:cp[1]]), cp[1])
  if (length(cp) > 1L) {
    for (k in 2:length(cp)) {
      out <- c(out, rep(mean(y[(cp[k-1] + 1):cp[k]]), cp[k] - cp[k-1]))
    }
  }
  out
}

# Define BIC penalty for PELT, BS
bic_pelt = function(x,y,n){
  n*log(sum((yhat_pelt(x)-y)^2)/n) + (2*length(x@cpts)-1)*log(n)
}

# Define BIC penalty for WBS
# --- replace ONLY this function ---
bic_wbs <- function(w, y){
  n <- length(y)
  # Extract CPs robustly
  cp <- tryCatch(sort(as.numeric(changepoints(w)$cpt.ic$bic.penalty)),
                 error = function(e) numeric(0))
  ncp <- length(cp)
  
  # If no CPs, force identical handling to PELT/BS zero-CP case
  if (ncp == 0L) {
    fitvec <- rep(mean(y), n)
    return(n * log( max(sum((fitvec - y)^2)/n, .Machine$double.eps) ) +
             (2 * 0L + 1L) * log(n))
  }
  
  # Otherwise, use segmentwise means as before
  fitvec <- yhat_wbs(w, y)
  # Guard against any accidental length mismatch
  if (length(fitvec) != n) {
    fitvec <- rep(mean(y), n)
    ncp <- 0L
  }
  
  n * log( max(sum((fitvec - y)^2)/n, .Machine$double.eps) ) +
    (2 * ncp + 1L) * log(n)
}



#############
## Adaptive / Iteratively Reweighted Fused Lasso
#############

# Repeat 20 iterations on 50 pc's
n = 1000
true_model = rep(0,n)
# Include number of changepoints
out_df = as.data.frame(matrix(NA,20,45))
cnames = c("true_model","true_nknots","bic_true_model",
           "flasso","nknots_flasso","bic_flasso",
           "adaflasso","nknots_adaflasso","bic_adaflasso",
           paste(c("iter","nknots_iter","bic_iter"),rep(2:10,each=3),"flasso",sep=""),
           "PELT","nknots_PELT","bic_PELT","knots_BS","nknots_BS","bic_BS",
           "knots_WBS","nknots_WBS","bic_WBS")
colnames(out_df)=cnames

for(i in 1:20){
  y = rnorm(n,0)
  true_bic = n*log(sum((true_model-y)^2)/(n))+7*log(n)
  out = c("251,501,751",3,true_bic)    # Repeat 11 iteratively reweighted estimations
  D0 <- diff(diag(n), differences = 1)
  for(k in 1:11){
    if(k == 1){
      D <- D0
    } else {
      wt <- abs(diff(beta_hat$beta[, which.min(bic)]))+1e-16
      D <- diag(1/wt) %*% D0
    }
    fit <- genlasso(y, D = D, maxsteps=500)
    beta_hat <- coef(fit)
    bic_find <- bic_fcn(beta_hat$beta, y)
    bic <- bic_find$bic
    loc = which.min(bic)
    minbic = bic[loc]
    
    # Location of knots found by best model
    knots = paste(which(abs(diff(beta_hat$beta[,loc], 
                                   differences=1))>1e-5)+1,
                   collapse=",")
    knot_locs = which(abs(diff(beta_hat$beta[, loc])) > 1e-5) + 1
    numberknots = length(which(abs(diff(beta_hat$beta[,loc], 
                                    differences=1))>1e-5))
    

    # Console line
    cat(sprintf("Simulation:%d; Iteration:%d; Changepoints:%d; Locations:%s\n",
                i, k, numberknots,
                if (numberknots) paste(knot_locs, collapse = ",") else "none"))
    
    # Plot 
    if (interactive()) {
      yhat <- beta_hat$beta[, loc]
      plot(y, type = "l",
           main = paste("Iteration", k),
           xlab = "Index", ylab = "Value")
      lines(yhat, lwd = 2, col = 2)
      if (numberknots) abline(v = knot_locs, lty = 3)
    }
  
    out = c(out,knots,numberknots,minbic)
    
    }
  # Now add PELT 
  cpt = cpt.mean(y, method="PELT",penalty="BIC")
  pelt_bic = bic_pelt(cpt,y,1000)
  pelt_nknots = length(cpt@cpts)-1
  pelt_knots = paste0(cpt@cpts[-length(cpt@cpts)]+1,collapse=",")

  # Now add BS
  bs = cpt.mean(y, method="BinSeg",penalty="BIC")
  bs_bic = bic_pelt(bs,y,1000)
  bs_nknots = length(bs@cpts)-1
  bs_knots = paste0(bs@cpts[-length(bs@cpts)]+1,collapse=",")

  # Now add WBS
  wbs = wbs(y)
  wbs_cp = sort(changepoints(wbs)$cpt.ic$bic.penalty)
  wbs_bic = bic_wbs(wbs,y)
  wbs_knots = paste0(sort(changepoints(wbs)$cpt.ic$bic.penalty+1),collapse=",")
  wbs_nknots = length(wbs_cp)
  
  # Update the dataframe
  out = c(out,
          pelt_knots,pelt_nknots,pelt_bic,
          bs_knots,bs_nknots,bs_bic,
          wbs_knots,wbs_nknots,wbs_bic)
  
  # Write to out_df
  out_df[i,] = out
}
write.csv(out_df,"out.csv")

