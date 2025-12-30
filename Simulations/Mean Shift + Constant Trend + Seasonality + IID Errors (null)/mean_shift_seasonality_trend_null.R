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
  s_mat = matrix(NA,s,ncol(beta_mat))
  time_index = matrix(1:n,nrow=n)
  pred = matrix(1:n, nrow = n, ncol = 1) %*% matrix(beta_mat[1,], nrow = 1)
  
  s_mat[1:(s-1),]=beta_mat[2:(s),]
  s_mat[s,]=-colSums(beta_mat[2:(s),])
  
  pred[1:(s-1),] = pred[1:(s-1),] + s_mat[1:(s-1),] +
    matrix(rep(beta_mat[s+1, ], s-1), 
           nrow = s-1, byrow = TRUE)
  
  pred[s:n, ] = pred[s:n, ] + beta_mat[s:n, ] + s_mat[mod12(s:n), ]
  
  
  return(pred)
}

# Vectorized deseasoned predict function
pred_fcn_deseasoned = function(beta_mat, s){
  n = nrow(beta_mat)
  p = ncol(beta_mat)
  
  # Preallocate prediction matrix
  time_index = matrix(1:n,nrow=n)
  pred = matrix(1:n, nrow = n, ncol = 1) %*% matrix(beta_mat[1,], nrow = 1)
  s_mat = matrix(NA,s,ncol(beta_mat))
  s_mat[1:(s-1),]=beta_mat[2:(s),]
  s_mat[s,]=-colSums(beta_mat[2:(s),])
  
  pred[1:(s-1),] = pred[1:(s-1),] +
    matrix(rep(beta_mat[s+1, ], s-1), 
           nrow = s-1, byrow = TRUE)
  
  pred[s:n,]= pred[s:n,] + beta_mat[s:n,]
  return(pred)
}

# Vectorized BIC computation
bic_fcn = function(beta, y, s) {
  n = length(y)
  p = ncol(beta)
  
  # Compute number of knots per column (only from position `s` onward)
  dbeta = diff(beta[s:n, , drop=FALSE], differences=1)
  knots = colSums(abs(dbeta) > 1e-3)
  
  # Vectorized predictions
  preds = pred_fcn(beta,s)
  
  # Compute residuals and sum of squares
  resids = preds - matrix(y, nrow=length(y), ncol=p)
  ss = colSums(resids^2)
  
  # BIC computation
  bic_vals = n * log(pmax(ss / n, .Machine$double.eps)) + (2 * knots + 13) * log(n)
  
  return(list(bic = bic_vals, knots = knots))
}

# Set parameters
s = 12       # Season length (e.g., months)
T = 100      # Number of seasons (e.g., years)
n = s * T    # Total number of time points
alpha = .005 # Trend

# Trend column
trend_column <- matrix(1:n, nrow = n, ncol = 1)
trend_affect = as.vector(alpha * trend_column)

# Create season index vector
S_block = rbind(diag(s - 1), rep(-1, s - 1))
season_dummies <- do.call(rbind, replicate(T, S_block, simplify = FALSE))

# Construct full design matrix: [trend | seasonal dummies | identity]
X = cbind(trend_column, season_dummies, as.matrix(Diagonal(n)[, (s + 1):n]))  # (n x n)
X = as.matrix(X)   # ensure base matrix to avoid non-conformable errors
storage.mode(X) <- "double"

# Make first seasonal block depend on first mean term
X[1:s, s + 1] = 1

# Visual Inspection if in manual mode
if(interactive()){
  X[1:(s+3),1:(s+3)]
}

m = n - s   # number of mean parameters

# Fused lasso penalty on mu_s to mu_n
D_mu <- sparseMatrix(
  i = rep(1:(m - 1), each = 2),
  j = as.vector(t(cbind((s + 1):(s + m - 1), (s + 2):(s + m)))),
  x = rep(c(-1, 1), times = m - 1),
  dims = c(m - 1, n)
)

# Visual Inspection if in manual mode
if(interactive()){
  D_mu[1:(s+3),1:(s+3)]
}
if(interactive()){
  dim(D_mu) # should be n-s-1 x n
}

# Simulate mu
mu = rep(0, n)

# Set up output out_df
out_df = as.data.frame(matrix(NA,20,40))
cnames = c("true_model","true_nknots","bic_true_model","true season", "true alpha",
           "flasso","nknots_flasso","bic_flasso","season_flasso", "alpha_flasso",
           "adaflasso","nknots_adaflasso","bic_adaflasso","season_adaflasso", "alpha_adaflasso",
           paste(c("iter","nknots_iter","bic_iter","season_iter", "alpha_iter"),
                 rep(2:6,each=5),"flasso",sep=""))
colnames(out_df)=cnames

for(i in 1:20){
  
  # Simulate seasonal component with scaled brownian bridge
  Z = rnorm(s)
  W = cumsum(Z) 
  t_grid = (1:s) / s
  B = W - t_grid * W[s]  
  B = B - mean(B)
  seasonal_effect = rep(B, T)
  y = mu + seasonal_effect + trend_affect + rnorm(n)
  true_bic = n*log(sum((mu + seasonal_effect+trend_affect -y)^2)/(n))+19*log(n)
  out = c("301,601,901",3,true_bic,paste(B,collapse=","),alpha)    # Repeat 7 iteratively reweighted estimations
  
  for(k in 1:7){
    
    # Set initial D
    if(k==1){
      D = D_mu  # size: (n-s-1) × n
    }
    
    # Weights
    if(k >1){
      wt =abs(diff( beta_hat$beta[(s+1):n,loc] ))+1e-16
      D = Diagonal(x = 1/wt) %*% D_mu   # (m-1 x m-1) %*% (m-1 x n) -> (m-1 x n)
    }
    
    # Run Genlasso
    suppressWarnings({
      fit = genlasso(y, X=X,D=D,maxsteps = 900)
    })
    beta_hat = coef(fit)
    
    # Find BIC for all models
    bic_find = bic_fcn(beta_hat$beta,y,s)
    
    # Best model by BIC
    bic = bic_find$bic
    loc = which.min(bic)
    minbic = bic[loc]
    
    # Location of knots found by best model
    knots = paste(which(abs(diff(beta_hat$beta[(s+1):n,loc], 
                                 differences=1))>1e-3)+13,
                  collapse=",")
    
    knot_locs = which(abs(diff(beta_hat$beta[(s+1):n,loc], 
                               differences=1))>1e-3) + s+1
    
    numberknots = length(knot_locs)
    
    season = paste(beta_hat$beta[2:s,loc],collapse=",")
    
    slope = beta_hat$beta[1,loc]
    
    out = c(out,knots,numberknots,minbic,season, slope)
    
    cat("Simulation: ",i,";  Iteration: ",k,";  Changepoints: ",
        numberknots,"; Slope: ",round(slope,6),"\n",sep="")
    
    if(interactive()){
      yhat = pred_fcn_deseasoned(beta_hat$beta, s)[, loc]
      plot(y, type = "l", main = paste("Iteration ", k))
      lines(yhat,type="l",col="red")
      abline(v = knot_locs, lty = "dotted")
    }
    
  }
  out_df[i,] = out
}
write.csv(out_df,"out.csv")