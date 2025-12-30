library(genlasso)
library(dplyr)
library(cpop)
library(not)
library(Matrix)

THR <- 1e-4

# --------------------------
# Utilities
# --------------------------

# Sample m interior changepoints in {2,...,n-1} with minimum spacing d
sample_cp_minsep <- function(n, m, d) {
  if (m == 0) return(integer(0))
  if (n - 2 < m * d + 1) stop("Infeasible: increase n or reduce m/d.")
  repeat {
    tau <- sort(sample(2:(n - 1), m, replace = FALSE))
    ok <- (tau[1] - 1) >= d && all(diff(tau) >= d) && (n - tau[length(tau)]) >= d
    if (ok) return(tau)
  }
}

# Continuous piecewise-linear signal with user-specified slopes per segment
# changepoints are interior indices in {2,...,n-1}
generate_signal_pl <- function(n, changepoints, slopes) {
  if (length(slopes) != length(changepoints) + 1) {
    stop("length(slopes) must be one more than the number of changepoints.")
  }
  tau <- sort(unique(changepoints))
  seg_starts <- c(1, tau + 1)
  seg_ends   <- c(tau, n)
  
  # Build slope per step and integrate
  slope_t <- numeric(n - 1)
  for (s in seq_along(slopes)) {
    a <- seg_starts[s]; b <- seg_ends[s]
    if (a < b) slope_t[a:(b - 1)] <- slopes[s]
  }
  y <- numeric(n)
  for (t in 1:(n - 1)) y[t + 1] <- y[t] + slope_t[t]
  y
}

# --------------------------
# BIC functions
# --------------------------

# BIC along a genlasso path for TF(1): df = K + 2
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

# Fitted values from a CPOP model
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

# BIC for CPOP with df = K + 2
bic_cpop <- function(x, y, n) {
  r  <- yhat_cpop(x) - y
  K  <- length(x@changepoints) - 2
  df <- K + 2
  n * log(pmax(sum(r^2) / n, .Machine$double.eps)) + df * log(n)
}

# Fitted values from NOT (piecewise-linear with continuity)
yhat_not <- function(obj) {
  x   <- 1:obj$n
  y   <- obj$x
  cps <- features(obj)$cpt
  if (length(cps) == 0) {
    dm <- cbind(1, x)
  } else {
    dm <- cbind(1, x, sapply(cps, function(tau) pmax(0, x - tau)))
  }
  coef <- qr.solve(dm, y)
  drop(dm %*% coef)
}

# BIC for NOT with df = K + 2
bic_not <- function(obj, y) {
  n   <- obj$n
  K   <- length(features(obj)$cpt)
  df  <- K + 2
  r   <- yhat_not(obj) - y
  n * log(pmax(sum(r^2) / n, .Machine$double.eps)) + df * log(n)
}

# --------------------------
# Simulation setup
# --------------------------

n <- 1000
d <- 100

out_df <- as.data.frame(matrix(NA, 20, 42))
cnames <- c(
  "true_model","true_nknots","bic_true_model",
  "flasso","nknots_flasso","bic_flasso",
  "adaflasso","nknots_adaflasso","bic_adaflasso",
  paste(c("iter","nknots_iter","bic_iter"), rep(2:10, each = 3), "flasso", sep = ""),
  "cpop","nknots_cpop","bic_cpop",
  "not","nknots_not","bic_not"
)
colnames(out_df) <- cnames

# Sparse base second-difference operator (TF order 1)
D0 <- diff(Diagonal(n), differences = 2)

# --------------------------
# Simulation loop
# --------------------------

for (i in 1:20) {
  # Random number of interior changepoints and locations
  m   <- sample(1:7, 1)
  tau <- sample_cp_minsep(n, m, d)
  
  # Accelerating slopes across segments (monotone increasing)
  slopes <- seq(0.005, 0.02, length.out = m + 1)
  
  # True signal and data
  true_model <- generate_signal_pl(n, tau, slopes)
  y <- rnorm(n, true_model)
  
  # True-model BIC with df = m + 2
  df_true  <- m + 2
  true_ss  <- sum((true_model - y)^2)
  true_bic <- n * log(pmax(true_ss / n, .Machine$double.eps)) + df_true * log(n)
  
  out <- c(paste0(tau, collapse = ","), m, true_bic)
  
  for (k in 1:11) {
    if (k == 1) {
      D <- D0
    } else {
      wt <- abs(diff(beta_hat$beta[, loc], differences = 2)) + 1e-16
      D  <- Diagonal(x = 1 / wt) %*% D0
    }
    
    fit <- genlasso(y, D = D, maxsteps = n)
    beta_hat <- coef(fit)
    
    # BIC over the path and model choice
    bic_find <- bic_fcn(beta_hat$beta, y)
    bic <- bic_find$bic
    loc <- which.min(bic)
    minbic <- bic[loc]
    
    # Knot extraction with unified threshold
    cp_idx      <- which(abs(diff(beta_hat$beta[, loc], differences = 2)) > THR)
    numberknots <- length(cp_idx)
    knot_locs   <- cp_idx + 1
    knots_str   <- if (numberknots > 0) paste(knot_locs, collapse = ",") else "none"
    
    out <- c(out, knots_str, numberknots, minbic)
    
    cat("Simulation: ", i, ", Iteration: ", k,
        "; Changepoints: ", numberknots, "; Locations: ", knots_str, "\n", sep = "")
    
    if (interactive()) {
      # --- BIC plot ---
      plot(bic, type = "l", lwd = 2,
           main = paste0("BIC Plot\nSimulation: ", i, "\nIteration: ", k))
      abline(v = loc, lty = "dotted", col = "red")
      
      # --- Model plot with clipped vertical lines and outside legend ---
      op <- par(no.readonly = TRUE)
      par(mar = c(5, 4, 4, 12))
      par(xpd = FALSE)
      
      plot(y, type = "l",
           main = paste0("Model Plot\nSimulation: ", i, "\nIteration: ", k),
           xlab = "Index", ylab = expression(y[t]), lwd = 1.5)
      lines(beta_hat$beta[, loc], col = "red",  lwd = 2)
      lines(1:n, true_model,     col = "blue", lwd = 2)
      
      # Vertical changepoint lines (IRFL estimate)
      if (numberknots > 0) {
        usr <- par("usr")
        segments(x0 = knot_locs, y0 = usr[3],
                 x1 = knot_locs, y1 = usr[4],
                 lty = 3, col = "gray40")
      }
      
      # True changepoints (vertical blue dashed lines)
      if (length(tau) > 0) {
        usr <- par("usr")
        segments(x0 = tau, y0 = usr[3],
                 x1 = tau, y1 = usr[4],
                 lty = 2, col = "blue")
      }
      
      # Legend outside the panel on the right, boxed
      usr <- par("usr")
      par(xpd = NA)
      legend(x = usr[2] + 0.06 * diff(usr[1:2]),
             y = usr[4],
             legend = c("True Model", "IRFL Estimate", "Estimated CPs", "True CPs"),
             col = c("blue", "red", "gray40", "blue"),
             lty = c(1, 1, 3, 2), lwd = 2,
             bty = "o", box.lwd = 0.8,
             xjust = 0, yjust = 1)
      par(op)
    }
    
  }
  
  # CPOP
  cpop_fit     <- suppressMessages(suppressWarnings(cpop(y, 1:n)))
  cpop_bic     <- bic_cpop(cpop_fit, y, n)
  cpop_nknots  <- length(cpop_fit@changepoints) - 2
  cpop_knots   <- if (cpop_nknots > 0) {
    paste0(cpop_fit@changepoints[-c(1, length(cpop_fit@changepoints))] + 1, collapse = ",")
  } else {
    "none"
  }
  
  # NOT
  not_fit     <- not(y, contrast = "pcwsLinContMean")
  not_bic     <- bic_not(not_fit, y)
  not_nknots  <- length(features(not_fit)$cpt)
  not_knots   <- if (not_nknots > 0) paste0(features(not_fit)$cpt + 1, collapse = ",") else "none"
  
  out <- c(out,
           cpop_knots, cpop_nknots, cpop_bic,
           not_knots,  not_nknots,  not_bic)
  out_df[i, ] <- out
}

write.csv(out_df, "out.csv", row.names = FALSE)
