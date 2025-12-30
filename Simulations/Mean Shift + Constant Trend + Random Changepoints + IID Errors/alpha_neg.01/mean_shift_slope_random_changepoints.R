
library(genlasso)
library(dplyr)
library(changepoint)
library(wbs)

# ---- choose one of: -0.01, -0.005, 0.005, 0.01 ----
SLOPE <- -0.01

# ---------- helpers ----------
pred_fcn <- function(beta, X) as.numeric(X %*% matrix(beta, ncol = 1L))

bic_fcn <- function(beta, y, X, tol = 1e-5) {
  n  <- length(y)
  L  <- ncol(beta)
  bic <- numeric(L)
  numberknots <- integer(L)
  for (k in seq_len(L)) {
    K    <- sum(abs(diff(beta[-1L, k])) > tol)
    pred <- as.numeric(X %*% beta[, k, drop = FALSE])
    sse  <- sum((y - pred)^2)
    bic[k] <- n * log(sse / n) + (K + 2L) * log(n)
    numberknots[k] <- K
  }
  list(bic = bic, knots = numberknots)
}

generate_cp_vector <- function(N, m, d) {
  if (N - (m + 1) * d < m) stop("Not possible: N - (m+1)d < m. Decrease m or d.")
  R <- N - (m + 1) * d
  random_parts <- function(total, n) {
    for (tries in 1:1000) {
      parts <- runif(n); parts <- parts / sum(parts)
      parts <- round(parts * total)
      diff  <- total - sum(parts)
      if (diff != 0) {
        idx <- sample(seq_len(n), size = abs(diff), replace = TRUE)
        parts[idx] <- parts[idx] + sign(diff)
      }
      if (all(parts >= 0)) return(parts)
    }
    stop("Failed to allocate random parts after 1000 tries.")
  }
  r_values <- random_parts(R, m + 1)
  semi <- r_values + d
  edges <- c(0, cumsum(semi))
  c(edges[1:(length(edges) - 1)], N)
}

generate_means <- function(cp.vec) {
  v <- numeric(length(cp.vec) - 1)
  for (i in 2:(length(cp.vec) - 1)) {
    step <- runif(1, min = 1.5, max = 2.5)
    direction <- sample(c(-1, 1), 1)
    v[i] <- v[i - 1] + direction * step
  }
  means <- numeric(0)
  for (k in 2:length(cp.vec)) {
    means <- c(means, rep(v[k - 1], cp.vec[k] - cp.vec[k - 1]))
  }
  means
}

# ---------- main ----------
n <- 1000L
x <- seq_len(n)
# --- Design matrix: (n x n) ---
Xstar <- cbind(x, diag(n))          # n x (n+1)
X     <- Xstar[, -2, drop = FALSE]  # remove col 2 -> n x n
X[1, 2] <- 1
if(interactive()){
  X[1:5,1:5]
}

cnames <- c(
  "true_model","true_slope","true_nknots","true_means","bic_true_model",
  "flasso","slope_flasso","slope_OLS_flasso","nknots_flasso","bic_flasso",
  "adaflasso","slope_adaflasso","slope_OLS_adaflasso","nknots_adaflasso","bic_adaflasso",
  paste(c("iter","slope_iter","slope_OLS_iter","nknots_iter","bic_iter"),
        rep(2:10, each = 5), "flasso", sep = "")
)
out_df <- as.data.frame(matrix(NA, 20, length(cnames)))
colnames(out_df) <- cnames

for (i in 1:20) {
  m <- sample(1:7, 1)
  g <- generate_cp_vector(n, m, d = 100)
  cpts <- g[-c(1, length(g))] + 1
  true_model <- generate_means(g)
  true_means <- round(unique(true_model), 4)
  y <- rnorm(n, true_model) + SLOPE * x
  
  true_sse <- sum((true_model + SLOPE * x - y)^2)
  true_bic <- n * log(true_sse / n) + (m + 2L) * log(n)
  cpts_string <- paste0(cpts, collapse = ",")
  
  out <- c(cpts_string, SLOPE, m, paste0(true_means, collapse = ","), true_bic)
  
  D_mean0 <- diff(diag(n), differences = 1L)
  
  for (k in 1:11) {
    if (k == 1) {
      Dprime <- diff(diag(n - 1L), differences = 1L)      # (n-2) x (n-1)
      D      <- cbind(0, Dprime)                           # (n-2) x n
    } else {
      wt_raw <- abs(diff(beta_hat$beta[-1L, loc]))         # length n-2 (mean block diffs)
      wts    <- 1 / (wt_raw + 1e-16)
      Dprime <- diag(wts, nrow = length(wts)) %*% diff(diag(n - 1L), 1L)
      D      <- cbind(0, Dprime)
    }
    
    fit <- genlasso(y, X = X, D = D, maxsteps = 500)
    beta_hat <- coef(fit)
    
    bic_find <- bic_fcn(beta_hat$beta, y, X)
    bic <- bic_find$bic
    loc <- which.min(bic)
    minbic <- bic[loc]
    
    knots <- paste(which(abs(diff(beta_hat$beta[-1L, loc])) > 1e-5) + 2L, collapse = ",")
    numberknots <- if (nzchar(knots)) length(strsplit(knots, ",")[[1]]) else 0L
    
    slope_from_FL <- as.numeric(beta_hat$beta[1L, loc])
    
    tau <- as.numeric(c(0, if (numberknots) as.integer(strsplit(knots, ",")[[1]]) - 1 else NULL, n))
    dm <- matrix(0, nrow = n, ncol = length(tau) - 1L)
    for (j in 2:length(tau)) dm[(tau[j - 1] + 1):tau[j], j - 1] <- 1
    X_OLS <- cbind(x, dm)
    beta_hat_OLS <- solve(t(X_OLS) %*% X_OLS, t(X_OLS) %*% y)
    slope_OLS <- as.numeric(beta_hat_OLS[1L])
    
    vals <- c(knots, slope_from_FL, slope_OLS, numberknots, minbic)
    out <- c(out, vals)
    
    cat(sprintf("Simulation:%d; Iteration:%d; Changepoints:%d; Locations:%s; Slope_FL:%g; Slope_OLS:%g\n",
                i, k, numberknots, if (numberknots) knots else "none",
                round(slope_from_FL, 4), round(slope_OLS, 4)))
    
    if (interactive()) {
      knot_locs <- if (numberknots) as.integer(strsplit(knots, ",")[[1]]) else integer(0)
      yhat_fl <- pred_fcn(beta_hat$beta[, loc, drop = FALSE], X)
      plot(y, type = "l", main = sprintf("Slope %.3f | Sim %d | Iter %d", SLOPE, i, k),
           xlab = "Index", ylab = "Value")
      lines(yhat_fl, col = 2, lwd = 2)
      if (length(knot_locs)) abline(v = knot_locs, lty = 3)
    }
  }
  
  out_df[i, ] <- out
}

write.csv(out_df, "out.csv", row.names = FALSE)
