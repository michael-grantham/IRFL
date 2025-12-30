library(genlasso)
library(dplyr)
library(cpop)
library(Matrix)

THR <- 1e-4

# --------------------------
# Utilities
# --------------------------

sample_cp_minsep <- function(n, m, d) {
  if (m == 0) return(integer(0))
  if (n - 2 < m * d + 1) stop("Infeasible: increase n or reduce m/d.")
  repeat {
    tau <- sort(sample(2:(n - 1), m, replace = FALSE))
    ok <- TRUE
    if ((tau[1] - 1) < d) ok <- FALSE
    if (ok && any(diff(tau) < d)) ok <- FALSE
    if (ok && (n - tau[length(tau)]) < d) ok <- FALSE
    if (ok) return(tau)
  }
}

generate_signal_alt <- function(n, tau, slope_mag = 0.01, start_sign = +1) {
  tau <- sort(unique(tau))
  m <- length(tau)
  seg_starts <- c(1, tau + 1)
  seg_ends   <- c(tau, n)
  seg_signs <- start_sign * ((-1)^(0:m))
  seg_slopes <- slope_mag * seg_signs
  
  slope_t <- numeric(n - 1)
  for (s in seq_along(seg_slopes)) {
    a <- seg_starts[s]; b <- seg_ends[s]
    if (a < b) slope_t[a:(b - 1)] <- seg_slopes[s]
  }
  y <- numeric(n)
  for (t in 1:(n - 1)) y[t + 1] <- y[t] + slope_t[t]
  y
}

# --------------------------
# BIC functions
# --------------------------

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

bic_cpop <- function(x, y, n) {
  r  <- yhat_cpop(x) - y
  K  <- length(x@changepoints) - 2
  df <- K + 2
  n * log(pmax(sum(r^2) / n, .Machine$double.eps)) + df * log(n)
}

# --------------------------
# Simulation setup
# --------------------------

n <- 1000
d <- 100

out_df <- as.data.frame(matrix(NA, 20, 39))
cnames <- c(
  "true_model","true_nknots","bic_true_model",
  "flasso","nknots_flasso","bic_flasso",
  "adaflasso","nknots_adaflasso","bic_adaflasso",
  paste(c("iter","nknots_iter","bic_iter"), rep(2:10, each = 3), "flasso", sep = ""),
  "cpop","nknots_cpop","bic_cpop"
)
colnames(out_df) <- cnames

D0 <- diff(Diagonal(n), differences = 2)

# --------------------------
# Simulation loop
# --------------------------

for (i in 1:20) {
  m <- sample(1:7, 1)
  tau <- sample_cp_minsep(n, m, d)
  true_model <- generate_signal_alt(n, tau, slope_mag = 0.01, start_sign = +1)
  
  y <- rnorm(n, true_model)
  df_true <- m + 2
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
    bic_find <- bic_fcn(beta_hat$beta, y)
    bic <- bic_find$bic
    loc <- which.min(bic)
    minbic <- bic[loc]
    
    cp_idx <- which(abs(diff(beta_hat$beta[, loc], differences = 2)) > THR)
    numberknots <- length(cp_idx)
    knot_locs <- cp_idx + 1
    knots_str <- if (numberknots > 0) paste(knot_locs, collapse = ",") else "none"
    
    out <- c(out, knots_str, numberknots, minbic)
    
    cat("Simulation: ", i, ", Iteration: ", k,
        "; Changepoints: ", numberknots, "; Locations: ", knots_str, "\n", sep = "")

    if (interactive()) {
      ## --- BIC plot (unchanged) ---
      plot(bic, type = "l", lwd = 2,
           main = paste0("BIC Plot\nSimulation: ", i, "\nIteration: ", k))
      abline(v = loc, lty = "dotted", col = "red")
      
      ## --- Model plot with clipped vertical lines and outside legend ---
      op <- par(no.readonly = TRUE)
      par(mar = c(5, 4, 4, 12))         # extra right margin for legend
      par(xpd = FALSE)                  # ensure panel clipping for segments
      
      plot(y, type = "l",
           main = paste0("Model Plot\nSimulation: ", i, "\nIteration: ", k),
           xlab = "Index", ylab = expression(y[t]), lwd = 1.5)
      lines(beta_hat$beta[, loc], col = "red",  lwd = 2)
      lines(1:n, true_model,          col = "blue", lwd = 2)
      
      # vertical changepoint lines, clipped to panel
      if (numberknots > 0) {
        usr <- par("usr")  # c(xmin,xmax,ymin,ymax)
        segments(x0 = knot_locs, y0 = usr[3],
                 x1 = knot_locs, y1 = usr[4],
                 lty = 3, col = "gray40")
      }
      
      # legend outside the panel on the right, in its own box
      usr <- par("usr")
      par(xpd = NA)  # allow drawing in the right margin
      legend(x = usr[2] + 0.06 * diff(usr[1:2]),  # a bit to the right of panel
             y = usr[4],                           # top-aligned
             legend = c("True Model", "IRFL Estimate"),
             col = c("blue", "red"),
             lty = 1, lwd = 2,
             bty = "o", box.lwd = 0.8,
             xjust = 0, yjust = 1)
      par(op)  # restore all settings
    }
    
    
  
  }
  
  cpop_fit     <- suppressMessages(suppressWarnings(cpop(y, 1:n)))
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
