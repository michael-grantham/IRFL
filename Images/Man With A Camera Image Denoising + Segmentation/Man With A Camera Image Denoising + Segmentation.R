###############################################################################
# BW IRFL workflow (denoise grid + edge overlays)
###############################################################################

########## ============== LIBRARIES ============== ##########
library(imager)
library(Matrix)
library(igraph)
library(png)
library(grid)
library(gridExtra)
library(ggplot2)

###############################################################################
# SECTION 1: HELPERS
###############################################################################

read_png_gray <- function(path) {
  img <- readPNG(path)
  if (length(dim(img)) == 2L) return(img)
  if (length(dim(img)) == 3L && dim(img)[3] >= 3L) {
    return(0.299 * img[, , 1] + 0.587 * img[, , 2] + 0.114 * img[, , 3])
  }
  stop("Unsupported PNG format: ", path)
}

load_image_grob <- function(path) {
  rasterGrob(readPNG(path), interpolate = FALSE)
}

###############################################################################
# SECTION 2: SPARSE DIFFERENCE OPERATORS
###############################################################################

get_diff_operators <- function(nrow_img, ncol_img) {
  pixel_index <- function(i, j) (j - 1L) * nrow_img + i
  
  rows_i <- list()
  cols_j <- list()
  vals_x <- list()
  horiz_flag <- logical(0)
  counter <- 1L
  
  for (j in 1L:ncol_img) {
    for (i in 1L:nrow_img) {
      curr <- pixel_index(i, j)
      
      if (i < nrow_img) {
        down <- pixel_index(i + 1L, j)
        rows_i[[counter]] <- c(counter, counter)
        cols_j[[counter]] <- c(curr, down)
        vals_x[[counter]] <- c(-1, 1)
        horiz_flag[counter] <- FALSE
        counter <- counter + 1L
      }
      
      if (j < ncol_img) {
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
    dims = c(counter - 1L, nrow_img * ncol_img)
  )
  
  list(
    Dx = D_all[which(horiz_flag), , drop = FALSE],
    Dy = D_all[which(!horiz_flag), , drop = FALSE]
  )
}

###############################################################################
# SECTION 3: WEIGHTED TV VIA ADMM (NO ADAPTIVE RHO; WARM START)
###############################################################################

tv_admm <- function(
    y, Dx, Dy, wx, wy, lambda,
    rho = 1.0, maxit = 200,
    abs_tol = 1e-6, rel_tol = 1e-4,
    alpha = 1.2,
    min_iter = 25,
    x_init = NULL
) {
  n <- length(y)
  
  Dx_w <- Matrix::Diagonal(x = wx) %*% Dx
  Dy_w <- Matrix::Diagonal(x = wy) %*% Dy
  D <- rbind(Dx_w, Dy_w)
  m <- nrow(D)
  Dt <- Matrix::t(D)
  
  AtA <- Matrix::Diagonal(n) + rho * Matrix::crossprod(D)
  AtA_chol <- Matrix::Cholesky(AtA, LDL = FALSE)
  
  x <- if (is.null(x_init)) as.numeric(y) else as.numeric(x_init)
  z <- as.numeric(D %*% x)
  u <- numeric(m)
  
  soft <- function(v, k) sign(v) * pmax(abs(v) - k, 0)
  
  for (k in 1:maxit) {
    rhs <- y + rho * as.numeric(Dt %*% (z - u))
    x <- as.numeric(Matrix::solve(AtA_chol, rhs))
    
    Dx_x <- as.numeric(D %*% x)
    z_old <- z
    
    Dx_x_hat <- alpha * Dx_x + (1 - alpha) * z_old
    z <- soft(Dx_x_hat + u, lambda / rho)
    
    r <- Dx_x - z
    u <- u + r
    
    s <- rho * as.numeric(Dt %*% (z - z_old))
    r_norm <- sqrt(sum(r * r))
    s_norm <- sqrt(sum(s * s))
    
    Dx_x_norm <- sqrt(sum(Dx_x * Dx_x))
    z_norm    <- sqrt(sum(z * z))
    Du_norm   <- sqrt(sum((rho * as.numeric(Dt %*% u))^2))
    
    eps_pri  <- sqrt(m) * abs_tol + rel_tol * max(Dx_x_norm, z_norm)
    eps_dual <- sqrt(n) * abs_tol + rel_tol * Du_norm
    
    if (k >= min_iter && r_norm <= eps_pri && s_norm <= eps_dual) break
  }
  
  x
}

###############################################################################
# SECTION 4: IRFL RUNNER (SAVES PNGS; VERBOSITY ONLY)
###############################################################################

write_png_gray_norm <- function(mat, path) {
  lo <- min(mat)
  hi <- max(mat)
  if (hi - lo < 1e-12) {
    norm <- matrix(0, nrow(mat), ncol(mat))
  } else {
    norm <- (mat - lo) / (hi - lo)
  }
  writePNG(norm, target = path)
}

run_irfl <- function(
    noised_matrix,
    num_iters = 6,
    lambda_grid = c(0.005, 0.01, 0.015, 0.02, 0.025),
    eps_w = 1e-5,
    rho = 1.0,
    admm = list(maxit = 200, abs_tol = 1e-6, rel_tol = 1e-4, alpha = 1.2, min_iter = 25)
) {
  nrow_img <- nrow(noised_matrix)
  ncol_img <- ncol(noised_matrix)
  
  diff_ops <- get_diff_operators(nrow_img, ncol_img)
  Dx <- diff_ops$Dx
  Dy <- diff_ops$Dy
  
  y <- as.vector(noised_matrix)
  
  for (lambda in lambda_grid) {
    cat(sprintf("Starting λ = %.5g\n", lambda))
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
        rho = rho,
        maxit = admm$maxit,
        abs_tol = admm$abs_tol,
        rel_tol = admm$rel_tol,
        alpha = admm$alpha,
        min_iter = admm$min_iter,
        x_init = beta
      )
      
      out_png <- sprintf("lambda_%g_iter_%d.png", lambda, iter)
      write_png_gray_norm(matrix(beta, nrow = nrow_img, ncol = ncol_img), out_png)
      cat(sprintf("  Completed iter %d/%d for λ = %.5g (wrote %s)\n",
                  iter, num_iters, lambda, out_png))
    }
    
    cat(sprintf("Finished λ = %.5g\n", lambda))
  }
  
  cat("IRFL run completed.\n")
  invisible(TRUE)
}

###############################################################################
# SECTION 5: GRID PDF (ONLY USES FILES THAT EXIST)
###############################################################################

make_irfl_grid_pdf <- function(lambda_vals, num_iters, pdf_path,
                               original_path = "camera_original_128x128.png",
                               noised_path   = "camera_heavily_noised_128x128.png") {
  if (!file.exists(original_path)) stop("Missing file: ", original_path)
  if (!file.exists(noised_path))   stop("Missing file: ", noised_path)
  
  image_paths <- matrix(NA_character_, nrow = length(lambda_vals), ncol = num_iters + 2)
  colnames(image_paths) <- c("Original", "Noised", paste0("Iter ", 1:num_iters))
  rownames(image_paths) <- paste0("lambda = ", lambda_vals)
  
  for (i in seq_along(lambda_vals)) {
    lam <- lambda_vals[i]
    image_paths[i, 1] <- original_path
    image_paths[i, 2] <- noised_path
    for (j in 1:num_iters) {
      f <- sprintf("lambda_%g_iter_%d.png", lam, j)
      image_paths[i, j + 2] <- if (file.exists(f)) f else NA_character_
    }
  }
  
  image_grobs <- matrix(vector("list", length(image_paths)),
                        nrow = nrow(image_paths), ncol = ncol(image_paths))
  for (i in 1:nrow(image_paths)) {
    for (j in 1:ncol(image_paths)) {
      p <- image_paths[i, j]
      image_grobs[[i, j]] <- if (!is.na(p)) load_image_grob(p) else nullGrob()
    }
  }
  
  col_labels <- c(
    list(textGrob("Original", rot = 90, just = "center")),
    list(textGrob("Noised",   rot = 90, just = "center")),
    lapply(1:num_iters, function(j) textGrob(sprintf("Iter %d", j), rot = 90, just = "center"))
  )
  row_labels <- lapply(lambda_vals, function(lam) textGrob(bquote(lambda == .(lam)), rot = 0, just = "center"))
  
  plot_grid <- list()
  plot_grid[[1]] <- nullGrob()
  plot_grid <- c(plot_grid, col_labels)
  
  for (i in 1:nrow(image_paths)) {
    plot_grid[[length(plot_grid) + 1]] <- row_labels[[i]]
    for (j in 1:ncol(image_paths)) {
      plot_grid[[length(plot_grid) + 1]] <- image_grobs[[i, j]]
    }
  }
  
  layout_matrix <- matrix(seq_along(plot_grid), nrow = length(lambda_vals) + 1, byrow = TRUE)
  
  sample_img <- readPNG(original_path)
  aspect_ratio <- dim(sample_img)[1] / dim(sample_img)[2]
  
  pdf(pdf_path,
      width  = (num_iters + 2 + 1),
      height = (length(lambda_vals) + 1) * aspect_ratio,
      useDingbats = FALSE)
  grid.arrange(grobs = plot_grid, layout_matrix = layout_matrix)
  dev.off()
  
  cat(sprintf("Wrote grid PDF: %s\n", pdf_path))
  invisible(TRUE)
}

###############################################################################
# SECTION 6: EDGE OVERLAY (NOISED + EDGES FROM DENOISED)
###############################################################################

extract_edges_mask <- function(mat, threshold = 0.05) {
  nr <- nrow(mat); nc <- ncol(mat)
  mask <- matrix(FALSE, nr, nc)
  
  dv <- abs(mat[-1, , drop = FALSE] - mat[-nr, , drop = FALSE]) > threshold
  dh <- abs(mat[, -1, drop = FALSE] - mat[, -nc, drop = FALSE]) > threshold
  
  mask[-nr, ] <- mask[-nr, ] | dv
  mask[, -nc] <- mask[, -nc] | dh
  
  mask
}

edge_overlay_grob <- function(base_png_path, edge_mask) {
  base_img <- readPNG(base_png_path)
  if (length(dim(base_img)) == 2L) base_img <- array(base_img, dim = c(nrow(base_img), ncol(base_img), 3L))
  
  coords <- which(edge_mask, arr.ind = TRUE)
  if (nrow(coords) == 0L) {
    return(rasterGrob(base_img, interpolate = FALSE))
  }
  
  xs <- (coords[,2] - 0.5) / ncol(edge_mask)
  ys <- 1 - (coords[,1] - 0.5) / nrow(edge_mask)
  
  grobTree(
    rasterGrob(base_img, interpolate = FALSE),
    pointsGrob(x = xs, y = ys, pch = 20, size = unit(0.5, "pt"), gp = gpar(col = "red"))
  )
}

save_1x3_pdf_original_noised_edges <- function(original_path, noised_path, denoised_path,
                                               out_pdf = "bw_irfl_1x3_edges_overlay.pdf",
                                               edge_threshold = 0.05) {
  if (!file.exists(original_path)) stop("Missing file: ", original_path)
  if (!file.exists(noised_path))   stop("Missing file: ", noised_path)
  if (!file.exists(denoised_path)) stop("Missing file: ", denoised_path)
  
  den_mat <- read_png_gray(denoised_path)
  edge_mask <- extract_edges_mask(den_mat, threshold = edge_threshold)
  
  g1 <- load_image_grob(original_path)
  g2 <- load_image_grob(noised_path)
  g3 <- edge_overlay_grob(noised_path, edge_mask)
  
  img_dim <- dim(readPNG(original_path))
  aspect_ratio <- img_dim[1] / img_dim[2]
  cell_height_in <- 1.2
  cell_width_in  <- cell_height_in / aspect_ratio
  
  pdf(out_pdf, width = 3 * cell_width_in, height = cell_height_in, useDingbats = FALSE)
  grid.arrange(g1, g2, g3, nrow = 1)
  dev.off()
  
  cat(sprintf("Wrote overlay PDF: %s\n", out_pdf))
  invisible(TRUE)
}

###############################################################################
# SECTION 7: RUN (EDIT ONLY THIS SECTION)
###############################################################################

# Required inputs:
#   original_matrix.png
#   noised_matrix.png
if (!file.exists("camera_original_128x128.png")) 
  stop("Missing file: camera_original_128x128.png")
if (!file.exists("camera_heavily_noised_128x128.png"))   
  stop("Missing file: camera_heavily_noised_128x128.png")

cat("Reading original_matrix.png ...\n")
original_matrix <- read_png_gray("camera_original_128x128.png")

cat("Reading noised_matrix.png ...\n")
noised_matrix <- read_png_gray("camera_heavily_noised_128x128.png")

lambda_vals <- c(0.005, 0.01, 0.015, 0.02, 0.025)
num_iters_run  <- 6
num_iters_grid <- 4

cat("Starting IRFL denoising ...\n")
run_irfl(
  noised_matrix = noised_matrix,
  num_iters = num_iters_run,
  lambda_grid = lambda_vals,
  eps_w = 1e-5,
  rho = 1.0,
  admm = list(maxit = 200, abs_tol = 1e-6, rel_tol = 1e-4, alpha = 1.2, min_iter = 25)
)

cat("Building grid PDF ...\n")
make_irfl_grid_pdf(lambda_vals = lambda_vals, num_iters = num_iters_grid, pdf_path = "bw_irfl_grid_labeled.pdf")

# 1x3 overlay: original | noised | noised+edges, where edges come from a denoised result
denoised_for_edges <- "lambda_0.025_iter_6.png"
if (!file.exists(denoised_for_edges)) {
  cat(sprintf("Skipping edge overlay: missing file %s\n", denoised_for_edges))
} else {
  cat("Building edge overlay PDF ...\n")
  save_1x3_pdf_original_noised_edges(
    original_path = "camera_original_128x128.png",
    noised_path   = "camera_heavily_noised_128x128.png",
    denoised_path = denoised_for_edges,
    out_pdf = "bw_irfl_1x3_edges_overlay.pdf",
    edge_threshold = 0.05
  )
}
