library(Matrix)
library(genlasso)

# Need 1-based modular function
mod12 <- function(x) {
  (x - 1) %% 12 + 1
}

# Vectorized predict function
pred_fcn = function(beta_mat,s) {
  n = nrow(beta_mat)
  p = ncol(beta_mat)
  
  # Preallocate prediction matrix
  pred = beta_mat
  s_mat = matrix(NA,s,ncol(beta_mat))
  s_mat[1:(s-1),]=beta_mat[1:(s-1),]
  s_mat[s,]=-colSums(beta_mat[1:(s-1),])
  
  pred[1:s, ] = matrix(rep(beta_mat[s, ], s), 
                       nrow = s, byrow = TRUE) + s_mat[1:s,]
  
  pred[(s+1):n, ] = beta_mat[(s+1):n, ] + s_mat[mod12((s+1):n),]
  
  return(pred)
}

# Vectorized deseasoned predict function
pred_fcn_deseasoned = function(beta_mat, s){
  n = nrow(beta_mat)
  p = ncol(beta_mat)
  pred = matrix(NA, n, p)
  pred[1:s, ] = matrix(rep(beta_mat[s, ], s), nrow = s, byrow = TRUE)
  pred[(s+1):n, ] = beta_mat[(s+1):n, ]
  return(pred)
}

# Vectorized BIC computation
bic_fcn = function(beta, y, s) {
  n = length(y)
  p = ncol(beta)
  
  # Compute number of knots per column (only from position `s` onward)
  dbeta = diff(beta[s:n, , drop=FALSE], differences=1)
  knots = colSums(abs(dbeta) > 1e-5)
  
  # Vectorized predictions
  preds = pred_fcn(beta,s)
  
  # Compute residuals and sum of squares
  resids = preds - matrix(y, nrow=length(y), ncol=p)
  ss = colSums(resids^2)
  
  # (3) BIC computation with numerical safety
  bic_vals = n * log(pmax(ss / n, .Machine$double.eps)) + (2 * knots + 12) * log(n)
  
  return(list(bic = bic_vals, knots = knots))
}

# Set parameters
s = 12       # Season length (e.g., months)
T = 100      # Number of seasons (e.g., years)
n = s * T    # Total number of time points

# Create season index vector
S_block = rbind(diag(s-1),rep(-1,s-1))
season_dummies <- do.call(rbind, replicate(T, S_block, simplify = FALSE))

# Construct full design matrix: [seasonal dummies | identity]
X = cbind(season_dummies, Diagonal(n)[,s:n])  # (n x n)

# (1) Link the entire first seasonal block (rows 1..s) to the first mean column
X[1:s, s] = 1

# Visual Inspection if in manual mode
if(interactive()){
  X[1:(s+3),1:(s+3)]
}

m = n - s + 1  # number of mean parameters

# Fused lasso penalty on mu_s to mu_n
D_mu = sparseMatrix(
  i = rep(1:(n - s), each = 2),
  j = as.vector(t(cbind(s:(n - 1), (s + 1):n))),
  x = rep(c(-1, 1), times = n - s),
  dims = c(n - s, n)
)

# Visual Inspection if in manual mode
if(interactive()){
  D_mu[1:(s+3),1:(s+3)]
}
if(interactive()){
  dim(D_mu) # should be n-s x n
}

# Simulate mu
mu = rep(0, n)

# Set up output out_df
out_df = as.data.frame(matrix(NA,20,32))
cnames = c("true_model","true_nknots","bic_true_model","true season",
           "flasso","nknots_flasso","bic_flasso","season_flasso",
           "adaflasso","nknots_adaflasso","bic_adaflasso","season_adaflasso",
           paste(c("iter","nknots_iter","bic_iter","season_iter"),
                 rep(2:6,each=4),"flasso",sep=""))
colnames(out_df)=cnames

for(i in 1:20){
  
  # Simulate seasonal component with scaled brownian bridge
  Z = rnorm(s)
  W = cumsum(Z) 
  t_grid = (1:s) / s
  B = W - t_grid * W[s]  
  B = B - mean(B)
  seasonal_effect = rep(B, T)
  y = mu + seasonal_effect + rnorm(n)
  
  # For the no-trend model with 11 seasonal params and 0 changepoints : 11 
  true_param_count <- 11
  true_bic = n * log(sum((mu + seasonal_effect - y)^2) / n) + true_param_count * log(n)
  
  out = c("",3,true_bic,paste(B,collapse=","))    # Repeat 7 iteratively reweighted estimations
  
  for(k in 1:7){
    
    # Set initial D
    if(k==1){
      D = D_mu  # size: (n-s) × n
    }
    
    # (2) Row-wise reweighting using a diagonal matrix (numerically safe with floor)
    if(k > 1){
      wt = abs(diff(beta_hat$beta[s:n, loc]))
      wt = pmax(wt, 1e-16)
      D  = Diagonal(x = 1/wt) %*% D_mu
    }
    
    # Run Genlasso
    suppressWarnings({
      fit = genlasso(y, X=as.matrix(X), D=D, maxsteps = 500)
    })
    beta_hat = coef(fit)
    
    # Find BIC for all models
    bic_find = bic_fcn(beta_hat$beta,y,s)
    
    # Best model by BIC
    bic = bic_find$bic
    loc = which.min(bic)
    minbic = bic[loc]
    
    # Location of knots found by best model
    knots = paste(which(abs(diff(beta_hat$beta[s:n,loc], 
                                 differences=1))>1e-5)+12,
                  collapse=",")
    
    knot_locs = which(abs(diff(beta_hat$beta[s:n,loc], 
                               differences=1))>1e-5) + s
    
    numberknots = length(knot_locs)
    season = paste(beta_hat$beta[1:(s-1),loc],collapse=",")
    
    out = c(out,knots,numberknots,minbic,season)
    
    loc_str <- if (numberknots > 0) paste(knot_locs, collapse = ",") else "none"
    cat("Simulation: ", i, "; Iteration: ", k,
        "; Changepoints: ", numberknots, "; Locations: ", loc_str, "\n", sep = "")
    
    if(interactive()){
      yhat = pred_fcn_deseasoned(beta_hat$beta, s)[, loc]
      plot(y, type = "l", main = paste("Iteration ", k),
           xaxt = "n", xlab = "Index")
      lines(yhat,type="l",col="red")
      abline(v = knot_locs, lty = "dotted")
    }
  }
  out_df[i,] = out
}
write.csv(out_df,"out.csv")
