###############################
# 0) Setup
###############################

library(Matrix)
library(genlasso)


###############################
# 1) Small Utilities
###############################

# 1-based modulo-12 helper
mod12 = function(x) {
  (x - 1) %% 12 + 1
}

# AR(1) correlation matrix: R_{ij} = phi^{|i-j|}
build_R = function(n, phi) {
  i = row(matrix(0, n, n))
  j = col(matrix(0, n, n))
  phi^abs(i - j)
}

# Whitening factor via Cholesky of R^{-1}
make_whiter = function(n, phi) {
  R = build_R(n, phi)
  chol(solve(R))
}

# Yule–Walker estimate of phi from residuals
phi_yw = function(preds, y){
  eps = y - preds
  n = length(y)
  sum(eps[2:n] * eps[1:(n-1)]) / sum(eps[1:(n-1)]^2)
}


###############################
# 2) Design Matrices
###############################

# Construct design matrix: [linear | quadratic | seasonal | segment means]
build_DM = function(n, s, cpts, season_dummies = NULL) {
  t = 1:n
  lin_col  = matrix(t,  ncol = 1)
  quad_col = matrix(t^2,ncol = 1)
  
  # Seasonal block (sum-to-zero per cycle)
  if (is.null(season_dummies)) {
    S_block = rbind(diag(s - 1), rep(-1, s - 1))
    season_dummies = do.call(rbind, replicate(ceiling(n / s), list(S_block)))
    season_dummies = season_dummies[1:n, ]
  }
  
  # Segment indicators
  cpts = sort(unique(cpts))
  breaks = c(0, cpts, n)
  m = length(breaks) - 1
  segment_mat = matrix(0, nrow = n, ncol = m)
  for (i in 1:m) segment_mat[(breaks[i] + 1):breaks[i + 1], i] = 1
  
  cbind(lin_col, quad_col, season_dummies, segment_mat)
}


###############################
# 3) Prediction & Criteria
###############################

# Predicted values for each column of beta_mat
pred_fcn = function(beta_mat, X) {
  as.matrix(X) %*% beta_mat
}

# Predicted values with seasonal coefficients set to zero
pred_fcn_deseasoned = function(beta_mat, X, s) {
  beta_mat[(3:(s+1)),] = 0
  as.matrix(X) %*% beta_mat
}

# BIC for each column of beta_mat under AR(1) errors
bic_fcn = function(beta_mat, X, y, phi) {
  n = nrow(beta_mat)
  p = ncol(beta_mat)
  
  # Changepoints inferred from jumps in fused means
  dbeta = diff(beta_mat[(s + 2):n, , drop = FALSE], differences = 1)
  knots = colSums(abs(dbeta) > 1e-5)
  
  preds  = pred_fcn(beta_mat, X)
  resids = preds - matrix(y, nrow = n, ncol = p)
  ss     = colSums(resids^2)
  
  # Closed-form log|R| for AR(1) (up to constants)
  log_det_R = (n - 1) * log(1 - phi^2)
  
  # df: 2 per changepoint + 14 fixed (trend + seasonal)
  df = 2 * knots + 14
  
  n * log(ss / n) + log_det_R + df * log(n) |>
    (\(bic_vals) list(bic = bic_vals, knots = knots))()
}


###############################
# 4) GLS with Known phi
###############################

gls_with_known_phi = function(y, X, phi, df_correction) {
  n = length(y)
  R = build_R(n, phi)
  R_inv = solve(R)
  
  Xt_Rinv = t(X) %*% R_inv
  beta_hat = solve(Xt_Rinv %*% X, Xt_Rinv %*% y)
  
  y_hat = as.numeric(X %*% beta_hat)
  resid = y - y_hat
  quad_form = as.numeric(t(resid) %*% R_inv %*% resid)
  
  log_det_R = (n - 1) * log(1 - phi^2)
  
  BIC = n * log(quad_form / n) + log_det_R + df_correction * log(n)
  list(beta = beta_hat, yhat = y_hat, BIC = BIC, resid = resid)
}


###############################
# 5) IRFL Solver (genlasso) + Reweighting
###############################

run_irfl = function(X, y, s, phi, niter = 10, maxsteps = 200, verbose = TRUE){
  BIC_list     = list()
  yhat_list    = list()
  yhat_d_list  = list()
  cpts_list    = list()
  beta_list    = list()
  
  n = length(y)
  p = ncol(X)
  m = p - (s + 1)  # number of fused mean parameters
  
  # First-order difference on fused means
  delta = sparseMatrix(
    i = rep(1:(m - 1), 2),
    j = c(1:(m - 1), 2:m),
    x = rep(c(-1, 1), each = m - 1),
    dims = c(m - 1, m)
  )
  
  for (k in 1:niter) {
    if (k == 1) {
      # Initial uniform penalty on fused means
      D = cbind(Matrix(0, nrow = m - 1, ncol = s + 1), delta)
    } else {
      # Iterative reweighting based on previous selection
      wt = abs(diff(beta_mat[(s + 2):n, loc]))
      wt[wt == 0] = 1e-16
      D = cbind(Matrix(0, nrow = m - 1, ncol = s + 1), Diagonal(x = 1 / wt) %*% delta)
    }
    
    fit = suppressWarnings(
      genlasso(y, X = as.matrix(X), D = D, maxsteps = maxsteps)
    )
    
    beta_mat = coef(fit)$beta
    bic_res  = bic_fcn(beta_mat, X, y, phi)
    bic      = bic_res$bic
    loc      = which.min(bic)
    minbic   = bic[loc]
    
    cpts   = which(abs(diff(beta_mat[(s + 2):n, loc])) > 1e-5) + (s + 2)
    n_cpts = length(cpts)
    yhat_d = pred_fcn_deseasoned(beta_mat, X, s)[, loc]
    yhat   = pred_fcn(beta_mat, X)[, loc]
    
    beta_list[[k]]   = beta_mat[, loc]
    BIC_list[[k]]    = minbic
    yhat_list[[k]]   = yhat
    yhat_d_list[[k]] = yhat_d
    cpts_list[[k]]   = cpts
    
    if (isTRUE(verbose)) {
      cat("Iteration: ", k, "; Changepoints: ", n_cpts,
          "; BIC: ", round(minbic, 6), "\n", sep = "")
      
      plot(y, type = "l", main = paste("Iteration ", k),
           xaxt = "n", xlab = "Year")
      lines(yhat_d, type = "l", col = "red")
      abline(v = cpts, lty = "dotted")
      
      # X-axis ticks every 60 observations (requires `year` in scope)
      tick_indices = seq(1, length(year), by = 60)
      tick_labels  = year[tick_indices]
      axis(1, at = tick_indices, labels = tick_labels)
    }
  }
  
  D_mu = cbind(Matrix(0, nrow = m - 1, ncol = s + 1), delta)
  
  list(
    y = y,
    X = X,
    Dmu = D_mu,
    beta = beta_list,
    BIC = BIC_list,
    yhat = yhat_list,
    yhat_d = yhat_d_list,
    cpts = cpts_list
  )
}


###############################
# 6) Data: Mauna Loa
###############################

dat  = read.table("MaunaLoa.txt", header = TRUE, sep = "", strip.white = TRUE)
y    = dat$avg
year = dat$year
n    = length(y)
s    = 12

# Base trend terms
quad_col = Matrix((1:n)^2, nrow = n, ncol = 1)
lin_col  = Matrix(1:n,     nrow = n, ncol = 1)

# Seasonal block (sum-to-zero per cycle)
S_block = rbind(diag(s - 1), rep(-1, s - 1))
season_dummies = do.call(rbind, replicate(ceiling(n / s), S_block, simplify = FALSE))
season_dummies = season_dummies[1:n, ]

# Full design: [linear | quadratic | seasonal | fused means]
X = cbind(lin_col, quad_col, season_dummies, Diagonal(n)[, (s + 2):n])
X[1:(s + 1), s + 2] = 1  # shared mean for t = 1..(s+1)

# Quick inspection of a small block
X[1:(s + 3), 1:(s + 3)]


###############################
# 7) phi Grid Search + IRFL (whitened)
###############################

result   = list()
phi_grid = seq(0.5, 0.90, length.out = 17)

for (phi in phi_grid) {
  phi_str = formatC(phi, format = "f", digits = 4)
  
  # Whitening transform
  L = make_whiter(n, phi)
  
  # IRFL on whitened data
  ds   = run_irfl(L %*% X, L %*% y, s, phi, niter = 8, verbose = FALSE)
  cpts = ds$cpts[[length(ds$cpts)]] - 1
  
  # GLS evaluation at this phi
  X_final   = build_DM(n, s, cpts)
  gls_result = gls_with_known_phi(y, X_final, phi, df_correction = 2 * length(cpts) + 14)
  yhat      = gls_result$yhat
  BIC       = gls_result$BIC
  cat("Final BIC with \u03C6 =", round(phi, 4), ":", BIC, "\n")
  
  # Store results
  result[[phi_str]] = list(BIC = BIC, cpts = cpts)
  
  # Visualization for this phi
  plot(y, type = "l", xaxt = "n", xlab = "Year")
  title(main = bquote(phi == .(round(phi, 5))), line = 2)
  title(main = bquote(BIC == .(round(BIC, 5))), line = 1)
  lines(yhat, type = "l", col = "red")
  abline(v = cpts, lty = "dotted")
  tick_indices = seq(1, length(year), by = 60)
  tick_labels  = year[tick_indices]
  axis(1, at = tick_indices, labels = tick_labels)
}

# Optional persistence
# saveRDS(result, "result.RDS")
if (file.exists("result.RDS")) result <- readRDS("result.RDS")

BIC_values = sapply(result, function(x) x$BIC)
best_phi   = names(BIC_values)[which.min(BIC_values)]
result[[best_phi]]$cpts
# Example reported: phi = 0.7000; cpts = 218 254 291 402 687


###############################
# 8) BIC for Fixed cpts (EM for phi within GLS)
###############################

estimate_bic_given_cpts <- function(y, cpts, s, phi_init = 0.3, tol = 1e-6, max_iter = 20, df_base = 14) {
  n        <- length(y)
  cpts     <- sort(unique(cpts))
  phi      <- phi_init
  phi_prev <- NA
  iter     <- 0
  
  X <- build_DM(n, s, cpts)
  
  # Alternates GLS for beta | phi and Yule–Walker update for phi | beta
  repeat {
    iter <- iter + 1
    R      <- phi^abs(outer(1:n, 1:n, "-"))
    R_inv  <- solve(R)
    Xt_Rinv <- crossprod(X, R_inv)
    beta_hat <- solve(Xt_Rinv %*% X, Xt_Rinv %*% y)
    y_hat    <- as.numeric(X %*% beta_hat)
    
    eps     <- y - y_hat
    phi_new <- sum(eps[2:n] * eps[1:(n-1)]) / sum(eps[1:(n-1)]^2)
    phi_new <- max(min(phi_new, 0.98), 0.01)
    
    if (!is.na(phi_prev) && abs(phi_new - phi_prev) < tol) break
    if (iter >= max_iter) break
    phi_prev <- phi
    phi      <- phi_new
  }
  
  resid     <- y - y_hat
  quad_form <- as.numeric(crossprod(resid, R_inv %*% resid))
  chol_R    <- chol(R)
  log_det_R <- 2 * sum(log(diag(chol_R)))
  
  df  <- 2 * length(cpts) + df_base
  BIC <- log_det_R + n * log(quad_form / n) + df * log(n)
  
  # Deseasoned prediction
  beta_ds         <- beta_hat
  beta_ds[(3:(s+1))] <- 0
  yhat_ds         <- as.numeric(X %*% beta_ds)
  
  list(BIC = BIC, phi = phi, yhat = y_hat, yhat_ds = yhat_ds,
       beta = beta_hat, resid = resid, X = X, cpts = cpts, iter = iter)
}


###############################
# 9) Reference cpts and Comparison
###############################

IRFL_cpts  <- result[[best_phi]]$cpts

# Print changepoint dates (monthly index → calendar)
print_cpt_dates <- function(cpts, start_year = 1959, start_month = 1) {
  for (idx in cpts) {
    total_months <- start_month - 1 + idx
    year <- start_year + (total_months - 1) %/% 12
    month <- (total_months - 1) %% 12 + 1
    cat(sprintf("Index %3d => %s %d\n", idx, month.name[month], year))
  }
}

print_cpt_dates(IRFL_cpts)

# BIC/phi under GLS with fixed cpts from IRFL
res_IRFL <- estimate_bic_given_cpts(y, IRFL_cpts, s = 12)
cat("IRFL BIC:",  res_IRFL$BIC,  "   φ:", round(res_IRFL$phi,  4), "\n")


###############################
# 10) GLS Cross-Check via nlme
###############################

library(nlme)

gls_analysis = function(n, s, cpts){
  X_final = build_DM(n, s, cpts)
  X_df <- as.data.frame(as.matrix(X_final))
  names(X_df) <- c("Slope",
                   "Quad",
                   paste0("s", 1:11),
                   paste0("Mean", 1:(ncol(X_final) - 13)))
  df <- cbind(y = y, X_df)
  
  fit_gls <- gls(
    y ~ . - 1,
    data = df,
    correlation = corAR1(form = ~ 1),
    method = "REML"
  )
  sum = summary(fit_gls)
  BIC = sum$BIC
  phi = coef(fit_gls$modelStruct$corStruct, unconstrained = FALSE)
  list(BIC = BIC, phi = phi, fit = fit_gls)
}

gls_IRFL = gls_analysis(n, s, IRFL_cpts)

gls_IRFL$BIC
summary(gls_IRFL$fit)

# Extract intervals and standard error for AR(1) parameter
phi_ci <- as.data.frame(intervals(gls_IRFL$fit)$corStruct)
phi_est <- phi_ci$est.
phi_lo  <- phi_ci$lower
phi_hi  <- phi_ci$upper
phi_moe <- (phi_hi - phi_lo) / 2

cat(sprintf("φ̂ = %.6f  MOE(φ̂) = %.6f   95%% CI: [%.6f, %.6f]\n",
            phi_est, phi_moe, phi_lo, phi_hi))


###############################
# 11) Tangent Slopes and Margins of Error (Yearly Scale)
###############################

tt   <- summary(gls_IRFL$fit)$tTable
b1_m <- tt[1, "Value"]      # monthly slope
b2_m <- tt[2, "Value"]      # monthly quadratic
p_b1 <- tt[1, "p-value"]
p_b2 <- tt[2, "p-value"]

# Convert to yearly coefficients
b1_y <- 12  * b1_m
b2_y <- 144 * b2_m

cat(sprintf("Yearly slope   β1(year): % .6f   (p = %.3g)\n", b1_y, p_b1))
cat(sprintf("Yearly quad    β2(year²): % .6f   (p = %.3g)\n", b2_y, p_b2))

# Covariance of (β1, β2) in monthly units → convert to yearly scale
V_month <- vcov(gls_IRFL$fit)[1:2, 1:2]
D <- diag(c(12, 144))
V_year <- D %*% V_month %*% D
b_year <- c(b1_y, b2_y)

# Tangent (yearly rate) at start (1959) and end (2025)
year_start <- 0
year_end   <- 66  # 1959 + 66 = 2025

rate_at <- function(t, b_year, V_year) {
  g <- c(1, 2 * t)
  rate <- sum(g * b_year)
  se   <- sqrt(as.numeric(t(g) %*% V_year %*% g))
  moe  <- qnorm(0.975) * se
  list(rate = rate, se = se, moe = moe)
}

r0  <- rate_at(year_start, b_year, V_year)
r66 <- rate_at(year_end,   b_year, V_year)

cat(sprintf("\nYearly tangent at start (1959): %.4f ± %.4f ppm/year\n", 
            r0$rate, r0$moe))
cat(sprintf("Yearly tangent at end   (2025): %.4f ± %.4f ppm/year\n", 
            r66$rate, r66$moe))


###############################
# 12) Plots (IRFL deseasoned fit)
###############################

library(ggplot2)
yhat_ds = res_IRFL$yhat_ds

df_plot <- data.frame(
  index   = seq_along(y),
  y       = y,
  yhat_ds = yhat_ds,
  year    = year
)

tick_indices = seq(1, length(year), by = 60)
tick_labels  = year[tick_indices]

p <- ggplot(df_plot, aes(x = index)) +
  geom_line(aes(y = y)) +
  geom_line(aes(y = yhat_ds), color = "red") +
  geom_vline(xintercept = IRFL_cpts, linetype = "dotted") +
  scale_x_continuous(name = NULL, breaks = tick_indices, labels = tick_labels) +
  ylab(expression(CO[2]~"PPM")) +
  theme_minimal()

print(p)
ggsave(filename = "MaunaLoaDS.pdf", plot = p, device = "pdf", width = 14 * .5, height = 6 * .5)


###############################
# 13) Plots (raw data)
###############################

month_index   <- 1:n
tick_positions <- seq(1, n, by = 60)
tick_labels    <- year[tick_positions]

df <- data.frame(Month = month_index, CO2 = y)

dataplot = ggplot(df, aes(x = Month, y = CO2)) +
  geom_line() +
  scale_x_continuous(breaks = tick_positions, labels = tick_labels) +
  scale_y_continuous(breaks = seq(floor(min(y) / 20) * 20, ceiling(max(y) / 20) * 20, by = 20)) +
  labs(x = NULL, y = expression("CO"[2]*" (ppm)")) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_blank())

print(dataplot)
ggsave(filename = "MaunaLoa.pdf", plot = dataplot, device = "pdf", width = 14 * .5, height = 6 * .5)


###############################
# 14) Local REML Grid Search Around Init cpts
###############################

search_best_cpts_REML <- function(y, s, init_cpts, width = 10, spacing = 5) {
  library(nlme)
  library(progress)
  
  n <- length(y)
  m <- length(init_cpts)
  
  # Grid per changepoint (centered at init, step = spacing, window = width)
  grid_list <- lapply(init_cpts, function(cpt) seq(cpt - width, cpt + width, by = spacing))
  cpt_grid <- expand.grid(grid_list)
  
  # Ensure the exact init configuration is present
  if (!any(apply(cpt_grid, 1, function(row) all(sort(as.integer(row)) == sort(init_cpts))))) {
    cpt_grid <- rbind(cpt_grid, init_cpts)
  }
  
  total_models <- nrow(cpt_grid)
  pb <- progress_bar$new(format = "Evaluating [:bar] :percent ETA: :eta",
                         total = total_models, clear = FALSE, width = 60)
  
  best_bic <- Inf
  best_fit <- NULL
  best_cpts <- NULL
  
  for (i in 1:total_models) {
    cpts <- sort(unique(as.integer(cpt_grid[i, ])))
    
    # Design and names compatible with gls_analysis
    X <- build_DM(n = n, s = s, cpts = cpts)
    X_df <- as.data.frame(as.matrix(X))
    colnames(X_df) <- c("Slope", "Quad", paste0("s", 1:(s - 1)), paste0("Mean", 1:(ncol(X) - 1 - 1 - (s - 1))))
    df <- cbind(y = y, X_df)
    
    # GLS REML fit
    try({
      fit <- gls(y ~ . - 1, data = df, correlation = corAR1(form = ~ 1), method = "REML", control = glsControl(singular.ok = TRUE))
      curr_bic <- BIC(fit)
      if (curr_bic < best_bic) {
        best_bic  <- curr_bic
        best_fit  <- fit
        best_cpts <- cpts
      }
    }, silent = TRUE)
    
    pb$tick()
  }
  
  list(best_fit = best_fit, best_cpts = best_cpts, best_bic = best_bic,
       best_phi = coef(best_fit$modelStruct$corStruct, unconstrained = FALSE))
}

# Local search around nominal IRFL cpts 
init_cpts   <- c(218, 254, 291, 402, 687)
result_local <- search_best_cpts_REML(y, s = 12, init_cpts = init_cpts, width = 12)

summary(result_local$best_fit)
result_local$best_bic
gls_IRFL$BIC
print_cpt_dates(result_local$best_cpts)
print_cpt_dates(IRFL_cpts)
