###############################################################################
# IRFL image denoising (weighted TV) + edge overlays
###############################################################################

########## ============== LIBRARIES ============== ##########
library(Matrix)
library(png)
library(jpeg)
library(grid)
library(gridExtra)

###############################################################################
# SECTION 1: IMAGE I/O (CONSISTENT ORIENTATION)
###############################################################################

# Read a JPG as a grayscale matrix in [0,1], with [row, col] = [y, x] (top-left origin).
read_jpg_gray <- function(path) {
  img <- readJPEG(path) # array [row, col, channel] in [0,1]
  if (length(dim(img)) == 2L) return(img)
  
  # Luma transform
  0.299 * img[,,1] + 0.587 * img[,,2] + 0.114 * img[,,3]
}

# Write a grayscale matrix (values arbitrary) to PNG, normalizing to [0,1].
write_png_gray <- function(mat, path) {
  min_val <- min(mat)
  max_val <- max(mat)
  
  if (max_val - min_val < 1e-12) {
    norm <- matrix(0, nrow(mat), ncol(mat))
  } else {
    norm <- (mat - min_val) / (max_val - min_val)
  }
  
  writePNG(norm, target = path)
}

# Read a grayscale PNG (or RGB PNG) as a grayscale matrix in [0,1].
read_png_gray <- function(path) {
  img <- readPNG(path)
  if (length(dim(img)) == 2L) return(img)
  0.299 * img[,,1] + 0.587 * img[,,2] + 0.114 * img[,,3]
}

###############################################################################
# SECTION 2: SPARSE DIFFERENCE OPERATORS
###############################################################################

# Build forward differences on a nrow x ncol grid.
# Dx: horizontal differences (right - current) for all valid edges.
# Dy: vertical differences   (down  - current) for all valid edges.
#
# Performance note:
# This mirrors the "fast" script’s approach (lists + unlist), which is typically
# very efficient in base R for constructing sparse operator triplets.
get_diff_operators <- function(nrow, ncol) {
  pixel_index <- function(i, j) (j - 1L) * nrow + i
  
  rows_i <- list()
  cols_j <- list()
  vals_x <- list()
  horiz_flag <- logical(0)
  
  counter <- 1L
  
  for (j in 1L:ncol) {
    for (i in 1L:nrow) {
      curr <- pixel_index(i, j)
      
      if (i < nrow) {
        down <- pixel_index(i + 1L, j)
        rows_i[[counter]] <- c(counter, counter)
        cols_j[[counter]] <- c(curr, down)
        vals_x[[counter]] <- c(-1, 1)
        horiz_flag[counter] <- FALSE
        counter <- counter + 1L
      }
      
      if (j < ncol) {
        right <- pixel_index(i, j + 1L)
        rows_i[[counter]] <- c(counter, counter)
        cols_j[[counter]] <- c(curr, right)
        vals_x[[counter]] <- c(-1, 1)
        horiz_flag[counter] <- TRUE
        counter <- counter + 1L
      }
    }
  }
  
  D_all <- sparseMatrix(
    i = unlist(rows_i),
    j = unlist(cols_j),
    x = unlist(vals_x),
    dims = c(counter - 1L, nrow * ncol)
  )
  
  list(
    Dx = D_all[which(horiz_flag), , drop = FALSE],
    Dy = D_all[which(!horiz_flag), , drop = FALSE]
  )
}

###############################################################################
# SECTION 3: WEIGHTED TV VIA ADMM (ROBUST STOPPING + ADAPTIVE RHO)
###############################################################################

# Solve:  minimize_x  0.5 ||x - y||^2 + lambda * ( ||W_x Dx x||_1 + ||W_y Dy x||_1 )
# where W_x = diag(wx), W_y = diag(wy) are iteration-dependent weights.
#
# This imitates the "fast" script:
# - explicitly forms D = rbind(Wx Dx, Wy Dy)
# - uses sparse Cholesky factorization of (I + rho D'D)
# - uses Boyd-style residuals and stopping
tv_admm <- function(
    y, Dx, Dy, wx, wy, lambda,
    rho = 1,
    maxit = 200,
    abs_tol = 1e-7,
    rel_tol = 1e-4,
    alpha = 1.5,     # over-relaxation
    mu = 10,         # residual ratio threshold
    tau_incr = 2,
    tau_decr = 2,
    x_init = NULL
) {
  n <- length(y)
  
  # Weighted difference operator (as in the fast script)
  Dx_w <- Diagonal(x = wx) %*% Dx
  Dy_w <- Diagonal(x = wy) %*% Dy
  D <- rbind(Dx_w, Dy_w)
  m <- nrow(D)
  Dt <- t(D)
  
  # Factor for x-update: (I + rho D'D)
  AtA <- Diagonal(n) + rho * crossprod(D)
  AtA_chol <- Cholesky(AtA)
  
  # Warm start
  x <- if (!is.null(x_init)) as.numeric(x_init) else as.numeric(y)
  z <- as.numeric(D %*% x)
  u <- rep(0, m)  # scaled dual variable
  
  for (iter in 1:maxit) {
    # x-update: solve (I + rho D'D) x = y + rho D'(z - u)
    rhs <- y + rho * as.numeric(Dt %*% (z - u))
    x <- as.numeric(solve(AtA_chol, rhs))
    
    # z-update (soft-threshold) with over-relaxation
    Dx_x <- as.numeric(D %*% x)
    z_old <- z
    v <- alpha * Dx_x + (1 - alpha) * z_old + u
    z <- pmax(0, v - lambda / rho) - pmax(0, -v - lambda / rho)
    
    # u-update (scaled dual)
    r <- Dx_x - z
    u <- u + r
    
    # Residual norms (Boyd et al.)
    r_norm <- sqrt(sum(r * r))
    dual_vec <- rho * as.numeric(Dt %*% (z - z_old))
    s_norm <- sqrt(sum(dual_vec * dual_vec))
    
    # Stopping thresholds (Boyd et al.)
    Dx_norm <- sqrt(sum(Dx_x * Dx_x))
    z_norm  <- sqrt(sum(z * z))
    rho_Dt_u_norm <- sqrt(sum((as.numeric(rho * (Dt %*% u)))^2))
    
    eps_pri  <- sqrt(m) * abs_tol + rel_tol * max(Dx_norm, z_norm)
    eps_dual <- sqrt(n) * abs_tol + rel_tol * rho_Dt_u_norm
    
    if (r_norm < eps_pri && s_norm < eps_dual) break
    
    # Adaptive rho with invariant scaling of the (unscaled) dual y = rho*u
    rho_new <- rho
    if (r_norm > mu * s_norm) {
      rho_new <- rho * tau_incr
    } else if (s_norm > mu * r_norm) {
      rho_new <- rho / tau_decr
    }
    
    if (rho_new != rho) {
      u <- (rho / rho_new) * u
      rho <- rho_new
      
      # Refactor since AtA depends on rho
      AtA <- Diagonal(n) + rho * crossprod(D)
      AtA_chol <- Cholesky(AtA)
    }
  }
  
  x
}

###############################################################################
# SECTION 4: IRFL LOOP (REWEIGHTED TV)
###############################################################################

# Lambda generator:
# - Provide a numeric vector directly, OR
# - Provide a function(img_matrix) -> numeric vector.
make_lambda <- function(lambda_spec, img_matrix) {
  if (is.function(lambda_spec)) return(lambda_spec(img_matrix))
  as.numeric(lambda_spec)
}

# IRFL denoise + save every iteration
#
# Performance alignment with the fast script:
# - beta initialized once per lambda
# - Dx, Dy computed once per image size (this is strictly better than the fast script)
irfl_denoise_and_save <- function(
    img_matrix,
    label,
    lambda_spec,
    num_iters = 5,
    eps_w = 1e-5,
    admm = list(rho = 1, maxit = 200, abs_tol = 1e-7, rel_tol = 1e-4, alpha = 1.5, mu = 10)
) {
  nr <- nrow(img_matrix)
  nc <- ncol(img_matrix)
  
  diff_ops <- get_diff_operators(nr, nc)
  Dx <- diff_ops$Dx
  Dy <- diff_ops$Dy
  
  y <- as.vector(img_matrix)
  lambda_vals <- make_lambda(lambda_spec, img_matrix)
  
  for (lambda in lambda_vals) {
    beta <- y
    
    for (iter in 1:num_iters) {
      dx <- as.numeric(Dx %*% beta)
      dy <- as.numeric(Dy %*% beta)
      
      wx <- 1 / (abs(dx) + eps_w)
      wy <- 1 / (abs(dy) + eps_w)
      
      beta <- tv_admm(
        y = y, Dx = Dx, Dy = Dy,
        wx = wx, wy = wy,
        lambda = lambda,
        rho = admm$rho,
        maxit = admm$maxit,
        abs_tol = admm$abs_tol,
        rel_tol = admm$rel_tol,
        alpha = admm$alpha,
        mu = admm$mu,
        x_init = beta
      )
      
      beta_mat <- matrix(beta, nrow = nr, ncol = nc)
      out_png <- sprintf("%s_lambda_%g_iter_%d.png", label, lambda, iter)
      write_png_gray(beta_mat, out_png)
    }
  }
  
  invisible(TRUE)
}

###############################################################################
# SECTION 5: EDGE THRESHOLDING + OVERLAY EXPORT
###############################################################################

# Compute absolute neighbor differences (horizontal + vertical) from a matrix.
neighbor_abs_diffs <- function(mat) {
  dv <- abs(mat[-1, , drop = FALSE] - mat[-nrow(mat), , drop = FALSE])
  dh <- abs(mat[, -1, drop = FALSE] - mat[, -ncol(mat), drop = FALSE])
  c(dv, dh)
}

# Gap-based threshold selection in log-space.
auto_threshold_from_gap <- function(gray_mat, bins = 200, plot = FALSE) {
  diffs <- neighbor_abs_diffs(gray_mat)
  diffs <- diffs[diffs > 1e-12]
  if (length(diffs) < 2L) return(0)
  
  sorted <- sort(diffs)
  log_sorted <- log(sorted)
  gaps <- diff(log_sorted)
  k <- which.max(gaps)
  threshold <- mean(sorted[c(k, k + 1L)])
  
  if (plot) {
    hist(diffs, breaks = bins, main = "Histogram of absolute neighbor differences",
         xlab = "Abs difference", col = "gray", border = NA)
    abline(v = threshold, col = "red", lwd = 2, lty = 2)
  }
  
  threshold
}

# Build an edge mask for a grayscale matrix (TRUE where an edge is detected).
edge_mask_from_threshold <- function(gray_mat, threshold) {
  nr <- nrow(gray_mat)
  nc <- ncol(gray_mat)
  
  mask <- matrix(FALSE, nr, nc)
  
  dv <- abs(gray_mat[-1, , drop = FALSE] - gray_mat[-nr, , drop = FALSE]) > threshold
  dh <- abs(gray_mat[, -1, drop = FALSE] - gray_mat[, -nc, drop = FALSE]) > threshold
  
  mask[-nr, ] <- mask[-nr, ] | dv
  mask[, -nc] <- mask[, -nc] | dh
  
  mask
}

# Create a grob: original JPG + red dot edges based on a denoised grayscale matrix.
edge_overlay_grob <- function(original_jpg_path, denoised_gray_mat, threshold) {
  base_rgb <- readJPEG(original_jpg_path) # [row, col, 3]
  if (length(dim(base_rgb)) == 2L) {
    base_rgb <- array(base_rgb, dim = c(nrow(base_rgb), ncol(base_rgb), 3L))
  }
  
  mask <- edge_mask_from_threshold(denoised_gray_mat, threshold)
  coords <- which(mask, arr.ind = TRUE)
  
  xs <- (coords[,2] - 0.5) / ncol(mask)
  ys <- 1 - (coords[,1] - 0.5) / nrow(mask)
  
  raster <- rasterGrob(base_rgb, interpolate = FALSE)
  points <- pointsGrob(
    x = xs, y = ys,
    pch = 20, size = unit(0.5, "pt"),
    gp = gpar(col = "red")
  )
  
  grobTree(raster, points)
}

# Save side-by-side PDF: original | overlay
save_edge_overlay_pdf <- function(original_jpg_path, denoised_png_path, out_pdf_path, threshold = NULL) {
  original_rgb <- readJPEG(original_jpg_path)
  denoised_gray <- read_png_gray(denoised_png_path)
  
  if (is.null(threshold)) threshold <- auto_threshold_from_gap(denoised_gray)
  
  original_grob <- rasterGrob(original_rgb, interpolate = FALSE)
  overlay_grob  <- edge_overlay_grob(original_jpg_path, denoised_gray, threshold)
  
  height_in <- dim(original_rgb)[1] / 72
  width_in  <- dim(original_rgb)[2] / 72
  
  pdf(out_pdf_path, width = 2 * width_in, height = height_in)
  grid.arrange(original_grob, overlay_grob, nrow = 1)
  dev.off()
  
  invisible(threshold)
}

###############################################################################
# SECTION 6: RUN (EDIT ONLY THIS SECTION FOR YOUR PROJECT)
###############################################################################

# Inputs
y_path <- "yes/Y3.jpg"
n_path <- "no/4 no.jpg"

y_matrix <- read_jpg_gray(y_path)
n_matrix <- read_jpg_gray(n_path)

# Customizable lambda:
lambda_spec <- function(img) {
  seq(0.0035, 0.0025, length.out = 3)
}

# Denoise + save PNGs
irfl_denoise_and_save(y_matrix, label = "y", lambda_spec = lambda_spec, num_iters = 5)
irfl_denoise_and_save(n_matrix, label = "n", lambda_spec = lambda_spec, num_iters = 5)

# Edge overlays (example uses the last file produced for lambda=0.0035, iter=5)
y_png <- "y_lambda_0.0035_iter_5.png"
n_png <- "n_lambda_0.0035_iter_5.png"

# Automatically choose thresholds from the denoised results and export PDFs
thresh_y <- save_edge_overlay_pdf(y_path, y_png, "y_edges_overlay.pdf", threshold = NULL)
thresh_n <- save_edge_overlay_pdf(n_path, n_png, "n_edges_overlay.pdf", threshold = NULL)

###############################################################################
