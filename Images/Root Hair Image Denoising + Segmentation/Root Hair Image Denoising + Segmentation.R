###############################################################################
# IRFL image denoising (weighted TV) + edge overlays
###############################################################################

########## ============== LIBRARIES ============== ##########
library(Matrix)
library(png)
library(jpeg)
library(grid)
library(gridExtra)
library(magick)

###############################################################################
# SECTION 1: IMAGE I/O (CONSISTENT ORIENTATION)
###############################################################################

read_tif_gray <- function(path) {
  img_magick <- image_read(path)
  img_gray <- image_convert(img_magick, colorspace = "gray")
  
  img_array <- image_data(img_gray, channels = "gray")
  mat <- as.numeric(img_array[1, , ]) / 255
  matrix(mat, nrow = dim(img_array)[2], ncol = dim(img_array)[3])
}

read_jpg_gray <- function(path) {
  img <- readJPEG(path)
  if (length(dim(img)) == 2L) return(img)
  0.299 * img[,,1] + 0.587 * img[,,2] + 0.114 * img[,,3]
}

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

read_png_gray <- function(path) {
  img <- readPNG(path)
  if (length(dim(img)) == 2L) return(img)
  0.299 * img[,,1] + 0.587 * img[,,2] + 0.114 * img[,,3]
}

###############################################################################
# SECTION 2: SPARSE DIFFERENCE OPERATORS
###############################################################################

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
# SECTION 3: WEIGHTED TV VIA ADMM (FIXED RHO)
###############################################################################

tv_admm <- function(
    y, Dx, Dy, wx, wy, lambda,
    rho = 1,
    maxit = 200,
    abs_tol = 1e-7,
    rel_tol = 1e-4,
    alpha = 1.5,
    x_init = NULL
) {
  n <- length(y)
  
  Dx_w <- Diagonal(x = wx) %*% Dx
  Dy_w <- Diagonal(x = wy) %*% Dy
  D <- rbind(Dx_w, Dy_w)
  m <- nrow(D)
  Dt <- t(D)
  
  AtA <- Diagonal(n) + rho * crossprod(D)
  AtA_chol <- Cholesky(AtA)
  
  x <- if (!is.null(x_init)) as.numeric(x_init) else as.numeric(y)
  z <- as.numeric(D %*% x)
  u <- rep(0, m)
  
  for (iter in 1:maxit) {
    rhs <- y + rho * as.numeric(Dt %*% (z - u))
    x <- as.numeric(solve(AtA_chol, rhs))
    
    Dx_x <- as.numeric(D %*% x)
    z_old <- z
    v <- alpha * Dx_x + (1 - alpha) * z_old + u
    z <- pmax(0, v - lambda / rho) - pmax(0, -v - lambda / rho)
    
    r <- Dx_x - z
    u <- u + r
    
    r_norm <- sqrt(sum(r * r))
    dual_vec <- rho * as.numeric(Dt %*% (z - z_old))
    s_norm <- sqrt(sum(dual_vec * dual_vec))
    
    Dx_norm <- sqrt(sum(Dx_x * Dx_x))
    z_norm  <- sqrt(sum(z * z))
    rho_Dt_u_norm <- sqrt(sum((as.numeric(rho * (Dt %*% u)))^2))
    
    eps_pri  <- sqrt(m) * abs_tol + rel_tol * max(Dx_norm, z_norm)
    eps_dual <- sqrt(n) * abs_tol + rel_tol * rho_Dt_u_norm
    
    if (r_norm < eps_pri && s_norm < eps_dual) break
  }
  
  x
}

###############################################################################
# SECTION 4: IRFL LOOP (REWEIGHTED TV)
###############################################################################

make_lambda <- function(lambda_spec, img_matrix) {
  if (is.function(lambda_spec)) return(lambda_spec(img_matrix))
  as.numeric(lambda_spec)
}

irfl_denoise_and_save <- function(
    img_matrix,
    label,
    lambda_spec,
    num_iters = 5,
    eps_w = 1e-3,
    admm = list(rho = 1, maxit = 200, abs_tol = 1e-7, rel_tol = 1e-4, alpha = 1.5)
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

neighbor_abs_diffs <- function(mat) {
  dv <- abs(mat[-1, , drop = FALSE] - mat[-nrow(mat), , drop = FALSE])
  dh <- abs(mat[, -1, drop = FALSE] - mat[, -ncol(mat), drop = FALSE])
  c(dv, dh)
}

auto_threshold_from_gap <- function(gray_mat, bins = 1001, plot = FALSE) {
  diffs <- neighbor_abs_diffs(gray_mat)
  diffs <- diffs[diffs > 1e-8]
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
    legend("topright", legend = paste("Threshold =", signif(threshold, 4)),
           col = "red", lty = 2, bty = "n")
  }
  
  threshold
}

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

edge_overlay_grob <- function(original_jpg_path, denoised_gray_mat, threshold) {
  base_rgb <- readJPEG(original_jpg_path)
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

save_edge_overlay_pdf <- function(original_jpg_path, denoised_png_path, out_pdf_path, threshold = NULL) {
  original_img <- readJPEG(original_jpg_path)
  denoised_gray <- read_png_gray(denoised_png_path)
  
  if (is.null(threshold)) {
    threshold <- auto_threshold_from_gap(denoised_gray, bins = 1001, plot = FALSE)
  }
  
  # Ensure original is 3-channel for consistent rendering
  if (length(dim(original_img)) == 2L) {
    original_img <- array(original_img, dim = c(nrow(original_img), ncol(original_img), 3L))
  }
  
  # Reflect across main diagonal (top-left to bottom-right)
  original_reflect <- aperm(original_img, c(2, 1, 3))
  denoised_reflect <- t(denoised_gray)
  
  # Edge mask on reflected denoised image
  mask <- edge_mask_from_threshold(denoised_reflect, threshold)
  coords <- which(mask, arr.ind = TRUE)
  
  xs <- (coords[, 2] - 0.5) / ncol(mask)
  ys <- 1 - (coords[, 1] - 0.5) / nrow(mask)
  
  original_grob <- rasterGrob(original_reflect, interpolate = FALSE)
  overlay_grob <- grobTree(
    rasterGrob(original_reflect, interpolate = FALSE),
    pointsGrob(
      x = xs, y = ys,
      pch = 20, size = unit(0.5, "pt"),
      gp = gpar(col = "red")
    )
  )
  
  height_in <- dim(original_reflect)[1] / 72
  width_in  <- dim(original_reflect)[2] / 72
  
  pdf(out_pdf_path, width = 2 * width_in, height = height_in)
  grid.arrange(original_grob, overlay_grob, nrow = 1)
  dev.off()
  
  invisible(threshold)
}

###############################################################################
# SECTION 6: RUN
###############################################################################

tif_path <- "DR-BT642-6-11012017-1_3 - EDF.tif"
img_matrix <- read_tif_gray(tif_path)

writeJPEG(img_matrix, "picture.jpg")

lambda_spec <- function(img) 0.003

irfl_denoise_and_save(
  img_matrix,
  label = "denoised",
  lambda_spec = lambda_spec,
  num_iters = 5,
  eps_w = 1e-3,
  admm = list(rho = 1, maxit = 200, abs_tol = 1e-7, rel_tol = 1e-4, alpha = 1.5)
)

denoised_png <- "denoised_lambda_0.003_iter_5.png"
original_jpg <- "picture.jpg"
out_pdf <- "root_hair_overlay.pdf"

thresh <- save_edge_overlay_pdf(original_jpg, denoised_png, out_pdf, threshold = NULL)
###############################################################################
