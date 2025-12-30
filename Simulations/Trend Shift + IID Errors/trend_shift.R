library(genlasso)
library(dplyr)
library(cpop)
library(Matrix)

# Threshold for detecting changepoints
THR <- 1e-4

# ---------------------------------------------------------------
# BIC function for generalized lasso trend filtering (TF(1))
# df = K + 2, where K is the number of interior knots
# ---------------------------------------------------------------
bic_fcn <- function(beta, y) {
  n  <- length(y)
  p  <- ncol(beta)
  bic <- numeric(p)
  numberknots <- integer(p)
  
  for (k in 1:p) {
    K   <- sum(abs(diff(beta[, k], differences = 2)) > THR)
    ss  <- sum((beta[, k] - y)^2)
    df  <- K + 2
    bic[k] <- n * log(pmax(ss / n, .Machine$double.eps)) + df * log(n)
    numberknots[k] <- K
  }
  
  list(bic = bic, knots = numberknots)
}

# ---------------------------------------------------------------
# Compute fitted values from a CPOP model
# ---------------------------------------------------------------
yhat_cpop <- function(x) {
  fit <- fitted(x)
  p   <- nrow(fit)
  n   <- fit[p, 3]
  yhat <- numeric(n)
  
  for (i in 1:p) {
    idx <- fit[i, 1]:(fit[i, 3] - 1)
    m   <- fit[i, 5]
    b   <- fit[i, 6]
    yhat[idx] <- m * idx + b
  }
  
  yhat[n] <- m * n + b
  yhat
}

# ---------------------------------------------------------------
# BIC for CPOP model using the same df convention as TF(1)
# df = K + 2, where K is the number of interior knots
# ---------------------------------------------------------------
bic_cpop <- function(x, y, n) {
  r  <- yhat_cpop(x) - y
  K  <- length(x@changepoints) - 2
  df <- K + 2
  n * log(pmax(sum(r^2) / n, .Machine$double.eps)) + df * log(n)
}

# ---------------------------------------------------------------
# Simulation setup
# ---------------------------------------------------------------
n <- 1000
true_model <- c(1:250, 251:2, 1:250, 251:2) / 100

out_df <- as.data.frame(matrix(NA, 20, 39))
cnames <- c(
  "true_model", "true_nknots", "bic_true_model",
  "flasso", "nknots_flasso", "bic_flasso",
  "adaflasso", "nknots_adaflasso", "bic_adaflasso",
  paste(c("iter", "nknots_iter", "bic_iter"), rep(2:10, each = 3), "flasso", sep = ""),
  "cpop", "nknots_cpop", "bic_cpop"
)
colnames(out_df) <- cnames

# Base second-difference operator for trend filtering (sparse)
D0 <- diff(Diagonal(n), differences = 2)

# ---------------------------------------------------------------
# Simulation loop
# ---------------------------------------------------------------
for (i in 1:20) {
  y <- rnorm(n, true_model)
  
  # True-model BIC with consistent df
  K_true  <- 3
  df_true <- K_true + 2
  true_ss  <- sum((true_model - y)^2)
  true_bic <- n * log(pmax(true_ss / n, .Machine$double.eps)) + df_true * log(n)
  
  out <- c("252,502,752", K_true, true_bic)
  
  for (k in 1:11) {
    if (k == 1) {
      D <- D0
    } else {
      wt <- abs(diff(beta_hat$beta[, loc], differences = 2)) + 1e-16
      D  <- Diagonal(x = 1 / wt) %*% D0
    }
    
    fit <- genlasso(y, D = D, maxsteps = 800)
    beta_hat <- coef(fit)
    
    bic_find <- bic_fcn(beta_hat$beta, y)
    bic <- bic_find$bic
    loc <- which.min(bic)
    minbic <- bic[loc]
    
    cp_idx <- which(abs(diff(beta_hat$beta[, loc], differences = 2)) > THR)
    numberknots <- length(cp_idx)
    knot_locs   <- cp_idx + 1
    knots_str   <- if (numberknots > 0) paste(knot_locs, collapse = ",") else "none"
    
    out <- c(out, knots_str, numberknots, minbic)
    
    cat("Simulation: ", i, "; Iteration: ", k,
        "; Changepoints: ", numberknots, "; Locations: ", knots_str, "\n", sep = "")
    
    if (interactive()) {
      plot(y, type = "l", main = paste("Iteration", k),
           xlab = "Index", ylab = "Value")
      lines(beta_hat$beta[, loc], col = 2, lwd = 2)
      if (numberknots > 0) abline(v = knot_locs, lty = 3)
    }
  }
  
  suppressMessages(
    suppressWarnings(
      cpop_fit <- cpop(y, 1:n)
    )
  )
  cpop_bic     <- bic_cpop(cpop_fit, y, n)
  cpop_nknots  <- length(cpop_fit@changepoints) - 2
  cpop_knots   <- if (cpop_nknots > 0) {
    paste0(cpop_fit@changepoints[-c(1, length(cpop_fit@changepoints))] + 1, collapse = ",")
  } else {
    "none"
  }
  
  out <- c(out, cpop_knots, cpop_nknots, cpop_bic)
  out_df[i, ] <- out
}

write.csv(out_df, "out.csv", row.names = FALSE)
