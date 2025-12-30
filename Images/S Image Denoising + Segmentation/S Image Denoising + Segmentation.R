###############################################################################
# S-image IRFL workflow (denoise grid + edge overlay)
###############################################################################

########## ============== LIBRARIES ============== ##########
library(png)
library(grid)
library(gridExtra)
library(ggplot2)
library(igraph)
library(Matrix)

###############################################################################
# SECTION 1: COLOR PATH (RGB ANCHORS -> CONTINUOUS SCALAR IN [1,4])
###############################################################################

anchors_rgb <- matrix(c(
  102, 205,   0,
  0,   0, 255,
  238, 130, 238,
  255,   0,   0
), ncol = 3, byrow = TRUE)

compute_path_coords <- function(rgb_anchors) {
  dists <- sqrt(rowSums((rgb_anchors[-1, ] - rgb_anchors[-nrow(rgb_anchors), ])^2))
  cum_dist <- c(0, cumsum(dists))
  path_coords <- cum_dist / max(cum_dist)
  list(coords = path_coords, rgb = rgb_anchors)
}

anchor_path <- compute_path_coords(anchors_rgb)

rgb_from_scalar <- colorRamp(
  rgb(anchors_rgb[, 1], anchors_rgb[, 2], anchors_rgb[, 3], maxColorValue = 255),
  space = "rgb"
)

rgb_to_scalar_continuous <- function(rgb, path) {
  min_dist <- Inf
  best_scalar <- NA_real_
  
  for (i in 1:(nrow(path$rgb) - 1)) {
    p1 <- path$rgb[i, ]
    p2 <- path$rgb[i + 1, ]
    
    denom <- sum((p2 - p1)^2)
    if (denom == 0) next
    
    t <- sum((rgb - p1) * (p2 - p1)) / denom
    t_clamped <- min(max(t, 0), 1)
    
    proj <- p1 + t_clamped * (p2 - p1)
    dist <- sum((rgb - proj)^2)
    
    if (dist < min_dist) {
      min_dist <- dist
      s1 <- path$coords[i]
      s2 <- path$coords[i + 1]
      best_scalar <- s1 + t_clamped * (s2 - s1)
    }
  }
  
  1 + 3 * best_scalar
}

read_png_as_continuous_vector <- function(path, path_data) {
  img <- readPNG(path)
  
  if (length(dim(img)) != 3L || dim(img)[3] < 3L) {
    stop("Image must have 3 or 4 channels.")
  }
  
  img <- img[, , 1:3, drop = FALSE]
  img_rgb <- round(img * 255)
  
  h <- dim(img_rgb)[1]
  w <- dim(img_rgb)[2]
  
  pixel_vals <- numeric(h * w)
  k <- 1L
  
  for (col in 1:w) {
    for (row in 1:h) {
      rgb <- img_rgb[row, col, ]
      pixel_vals[k] <- rgb_to_scalar_continuous(rgb, path_data)
      k <- k + 1L
    }
  }
  
  list(scalar = pixel_vals, height = h, width = w)
}

###############################################################################
# SECTION 2: IMAGE HELPERS (FAST, SAFE RASTERS)
###############################################################################

load_image_grob <- function(path) {
  img <- readPNG(path)
  rasterGrob(img, interpolate = FALSE)
}

plot_scalar_matrix_as_color <- function(scalar_matrix, title = NULL) {
  height <- nrow(scalar_matrix)
  width  <- ncol(scalar_matrix)
  
  scalar_vec <- pmin(pmax(as.vector(scalar_matrix), 1), 4)
  norm_vals  <- (scalar_vec - 1) / 3
  
  rgb_array <- array(0, dim = c(height, width, 3))
  k <- 1L
  for (col in 1L:width) {
    for (row in 1L:height) {
      rgb_array[row, col, ] <- as.numeric(rgb_from_scalar(norm_vals[k])) / 255
      k <- k + 1L
    }
  }
  
  img_raster <- as.raster(rgb_array)
  
  gg <- ggplot() +
    annotation_raster(img_raster, xmin = 0.5, xmax = width + 0.5,
                      ymin = 0.5, ymax = height + 0.5) +
    coord_fixed(expand = FALSE) +
    theme_void()
  
  if (!is.null(title)) gg <- gg + labs(title = title)
  gg
}

save_image_plot <- function(scalar_matrix, out_path) {
  height <- nrow(scalar_matrix)
  width  <- ncol(scalar_matrix)
  
  scalar_vec <- pmin(pmax(as.vector(scalar_matrix), 1), 4)
  norm_vals  <- (scalar_vec - 1) / 3
  
  rgb_array <- array(0, dim = c(height, width, 3))
  k <- 1L
  for (col in 1L:width) {
    for (row in 1L:height) {
      rgb_array[row, col, ] <- as.numeric(rgb_from_scalar(norm_vals[k])) / 255
      k <- k + 1L
    }
  }
  
  writePNG(rgb_array, target = out_path)
}

###############################################################################
# SECTION 3: DIFFERENCE OPERATORS
###############################################################################

get_diff_operators <- function(nrow_img, ncol_img) {
  pixel_index <- function(i, j) (j - 1) * nrow_img + i
  
  rows_i <- list()
  cols_j <- list()
  vals_x <- list()
  horiz_flag <- c()
  counter <- 1
  
  for (j in 1:ncol_img) {
    for (i in 1:nrow_img) {
      curr <- pixel_index(i, j)
      
      if (i < nrow_img) {
        down <- pixel_index(i + 1, j)
        rows_i[[counter]] <- c(counter, counter)
        cols_j[[counter]] <- c(curr, down)
        vals_x[[counter]] <- c(-1, 1)
        horiz_flag[counter] <- FALSE
        counter <- counter + 1
      }
      
      if (j < ncol_img) {
        right <- pixel_index(i, j + 1)
        rows_i[[counter]] <- c(counter, counter)
        cols_j[[counter]] <- c(curr, right)
        vals_x[[counter]] <- c(-1, 1)
        horiz_flag[counter] <- TRUE
        counter <- counter + 1
      }
    }
  }
  
  D_all <- sparseMatrix(
    i = unlist(rows_i),
    j = unlist(cols_j),
    x = unlist(vals_x),
    dims = c(counter - 1, nrow_img * ncol_img)
  )
  
  Dx <- D_all[which(horiz_flag), ]
  Dy <- D_all[which(!horiz_flag), ]
  list(Dx = Dx, Dy = Dy)
}

###############################################################################
# SECTION 4: WEIGHTED TV VIA ADMM
###############################################################################

tv_admm <- function(
    y, Dx, Dy, wx, wy, lambda,
    rho = 1.0, maxit = 5000,
    eps_abs = 1e-5, eps_rel = 5e-4,
    alpha = 1.2,
    adapt_rho = FALSE, mu = 10, tau_inc = 2, tau_dec = 2,
    min_iter = 100,
    x_init = NULL, z_init = NULL, u_init = NULL
) {
  n <- length(y)
  stopifnot(
    ncol(Dx) == n, ncol(Dy) == n,
    length(wx) == nrow(Dx), length(wy) == nrow(Dy)
  )
  
  Dx_w <- Matrix::Diagonal(x = wx) %*% Dx
  Dy_w <- Matrix::Diagonal(x = wy) %*% Dy
  D <- rbind(Dx_w, Dy_w)
  m <- nrow(D)
  Dt <- Matrix::t(D)
  
  soft <- function(v, k) sign(v) * pmax(abs(v) - k, 0)
  
  make_factor <- function(rho_val) {
    A <- Matrix::Diagonal(n) + rho_val * (Dt %*% D)
    Matrix::Cholesky(A, LDL = FALSE)
  }
  cholA <- make_factor(rho)
  
  x <- if (is.null(x_init)) as.numeric(y) else as.numeric(x_init)
  z <- if (is.null(z_init)) as.numeric(D %*% x) else as.numeric(z_init)
  u <- if (is.null(u_init)) numeric(m) else as.numeric(u_init)
  
  for (k in 1:maxit) {
    rhs <- y + rho * as.numeric(Dt %*% (z - u))
    x <- as.numeric(Matrix::solve(cholA, rhs))
    
    Dx_x <- as.numeric(D %*% x)
    z_old <- z
    
    Dx_x_hat <- alpha * Dx_x + (1 - alpha) * z
    z <- soft(Dx_x_hat + u, lambda / rho)
    
    r <- Dx_x - z
    u <- u + r
    
    s <- rho * as.numeric(Dt %*% (z - z_old))
    r_norm <- sqrt(sum(r * r))
    s_norm <- sqrt(sum(s * s))
    
    Dx_x_norm <- sqrt(sum(Dx_x * Dx_x))
    z_norm    <- sqrt(sum(z * z))
    Du_norm   <- sqrt(sum((rho * as.numeric(Dt %*% u))^2))
    
    eps_pri  <- sqrt(m) * eps_abs + eps_rel * max(Dx_x_norm, z_norm)
    eps_dual <- sqrt(n) * eps_abs + eps_rel * Du_norm
    
    if (k >= min_iter && r_norm <= eps_pri && s_norm <= eps_dual) break
    
    if (adapt_rho) {
      rho_new <- rho
      if (r_norm >  mu * s_norm) rho_new <- rho * tau_inc
      else if (s_norm > mu * r_norm) rho_new <- rho / tau_dec
      
      if (rho_new != rho) {
        u <- (rho / rho_new) * u
        rho <- rho_new
        cholA <- make_factor(rho)
      }
    }
  }
  
  list(x = x, z = z, u = u)
}

###############################################################################
# SECTION 5: IRFL RUNNER
###############################################################################

run_irfl <- function(
    original_matrix, noised_matrix,
    num_iters = 6,
    lambda_grid = c(0.1, 0.5, 1, 1.5, 2),
    eps = 5e-3,
    w_min = 0.1, w_max = 1e3,
    eta = 0.5,
    rho = 1.0
) {
  nrow_img <- nrow(noised_matrix)
  ncol_img <- ncol(noised_matrix)
  
  diff_ops <- get_diff_operators(nrow_img, ncol_img)
  Dx <- diff_ops$Dx
  Dy <- diff_ops$Dy
  
  y_vec <- as.vector(noised_matrix)
  
  out <- list()
  
  for (lambda in lambda_grid) {
    cat(sprintf("Starting λ = %.3g\n", lambda))
    
    beta <- y_vec
    
    dx <- as.numeric(Dx %*% beta)
    dy <- as.numeric(Dy %*% beta)
    wx <- 1 / sqrt(dx * dx + eps * eps)
    wy <- 1 / sqrt(dy * dy + eps * eps)
    wx <- pmin(pmax(wx, w_min), w_max)
    wy <- pmin(pmax(wy, w_min), w_max)
    
    x_ws <- NULL
    z_ws <- NULL
    u_ws <- NULL
    
    for (iter in 1:num_iters) {
      dx <- as.numeric(Dx %*% beta)
      dy <- as.numeric(Dy %*% beta)
      
      wx_new <- 1 / sqrt(dx * dx + eps * eps)
      wy_new <- 1 / sqrt(dy * dy + eps * eps)
      wx_new <- pmin(pmax(wx_new, w_min), w_max)
      wy_new <- pmin(pmax(wy_new, w_min), w_max)
      
      wx <- (1 - eta) * wx + eta * wx_new
      wy <- (1 - eta) * wy + eta * wy_new
      
      for (m in c(1.5, 1.2, 1.0)) {
        fit <- tv_admm(
          y = y_vec, Dx = Dx, Dy = Dy, wx = wx, wy = wy,
          lambda = lambda * m,
          rho = rho,
          alpha = 1.2,
          eps_abs = 1e-5, eps_rel = 5e-4,
          min_iter = 100,
          adapt_rho = FALSE,
          x_init = x_ws, z_init = z_ws, u_init = u_ws
        )
        x_ws <- fit$x
        z_ws <- fit$z
        u_ws <- fit$u
      }
      
      beta <- x_ws
      
      beta_mat <- matrix(beta, nrow = nrow_img, ncol = ncol_img)
      out_png <- sprintf("lambda_%g_iter_%d.png", lambda, iter)
      save_image_plot(beta_mat, out_png)
      
      cat(sprintf("  Completed iter %d/%d for λ = %.3g (wrote %s)\n",
                  iter, num_iters, lambda, out_png))
    }
    
    out[[as.character(lambda)]] <- matrix(beta, nrow = nrow_img, ncol = ncol_img)
    cat(sprintf("Finished λ = %.3g\n", lambda))
  }
  
  cat("IRFL run completed.\n")
  out
}

###############################################################################
# SECTION 6: GRID PDF
###############################################################################

make_irfl_grid_pdf <- function(lambda_vals, num_iters, pdf_path) {
  image_paths <- matrix(NA_character_, nrow = length(lambda_vals), ncol = num_iters + 2)
  colnames(image_paths) <- c("Original", "Noised", paste0("Iter ", 1:num_iters))
  rownames(image_paths) <- paste0("lambda = ", lambda_vals)
  
  for (i in seq_along(lambda_vals)) {
    lam <- lambda_vals[i]
    image_paths[i, 1] <- "original_matrix.png"
    image_paths[i, 2] <- "noised_matrix.png"
    for (j in 1:num_iters) {
      image_paths[i, j + 2] <- sprintf("lambda_%g_iter_%d.png", lam, j)
    }
  }
  
  image_grobs <- matrix(vector("list", length(image_paths)),
                        nrow = nrow(image_paths), ncol = ncol(image_paths))
  for (i in 1:nrow(image_paths)) {
    for (j in 1:ncol(image_paths)) {
      image_grobs[[i, j]] <- load_image_grob(image_paths[i, j])
    }
  }
  
  col_labels <- c(
    list(textGrob("Original", rot = 90, just = "center")),
    list(textGrob("Noised",   rot = 90, just = "center")),
    lapply(1:num_iters, function(j) textGrob(sprintf("Iter %d", j), rot = 90, just = "center"))
  )
  
  row_labels <- lapply(lambda_vals, function(lam) {
    textGrob(bquote(lambda == .(lam)), rot = 0, just = "center")
  })
  
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
  
  sample_img <- readPNG("original_matrix.png")
  aspect_ratio <- dim(sample_img)[1] / dim(sample_img)[2]
  
  pdf(pdf_path,
      width  = (num_iters + 2 + 1),
      height = (length(lambda_vals) + 1) * aspect_ratio,
      useDingbats = FALSE)
  
  grid.arrange(grobs = plot_grid, layout_matrix = layout_matrix,
               heights = unit(rep(1, length(lambda_vals) + 1), "null"),
               widths  = unit(rep(1, num_iters + 2 + 1), "null"))
  dev.off()
  
  cat(sprintf("Wrote grid PDF: %s\n", pdf_path))
  invisible(TRUE)
}

###############################################################################
# SECTION 7: EDGE OVERLAY (FOR A SINGLE RESULT)
###############################################################################

extract_edges <- function(beta_mat, nrow_img, ncol_img, threshold = 1e-4) {
  horiz_diff <- abs(beta_mat[, -1] - beta_mat[, -ncol_img]) > threshold
  vert_diff  <- abs(beta_mat[-1, ] - beta_mat[-nrow_img, ]) > threshold
  
  horiz_idx <- which(horiz_diff, arr.ind = TRUE)
  vert_idx  <- which(vert_diff, arr.ind = TRUE)
  
  horiz_segments <- data.frame(
    x = horiz_idx[, 2] - 0.5,
    y = nrow_img - horiz_idx[, 1] + 1,
    xend = horiz_idx[, 2] + 0.5,
    yend = nrow_img - horiz_idx[, 1] + 1
  )
  
  vert_segments <- data.frame(
    x = vert_idx[, 2],
    y = nrow_img - vert_idx[, 1] + 1 + 0.5,
    xend = vert_idx[, 2],
    yend = nrow_img - vert_idx[, 1] + 1 - 0.5
  )
  
  rbind(horiz_segments, vert_segments)
}

plot_with_edges <- function(scalar_matrix, edge_df) {
  height <- nrow(scalar_matrix)
  width  <- ncol(scalar_matrix)
  
  scalar_vec <- pmin(pmax(as.vector(scalar_matrix), 1), 4)
  norm_vals  <- (scalar_vec - 1) / 3
  
  rgb_array <- array(0, dim = c(height, width, 3))
  k <- 1L
  for (col in 1L:width) {
    for (row in 1L:height) {
      rgb_array[row, col, ] <- as.numeric(rgb_from_scalar(norm_vals[k])) / 255
      k <- k + 1L
    }
  }
  
  img_raster <- as.raster(rgb_array)
  
  ggplot() +
    annotation_raster(img_raster, xmin = 0.5, xmax = width + 0.5,
                      ymin = 0.5, ymax = height + 0.5) +
    geom_segment(data = edge_df,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 color = "black", linewidth = 0.15, lineend = "round") +
    coord_fixed(expand = FALSE) +
    theme_void()
}

###############################################################################
# SECTION 8: RUN
###############################################################################

if (!file.exists("original_matrix.png")) stop("Missing file: original_matrix.png")
if (!file.exists("noised_matrix.png"))   stop("Missing file: noised_matrix.png")

cat("Reading original_matrix.png ...\n")
orig_info <- read_png_as_continuous_vector("original_matrix.png", anchor_path)
original_matrix <- matrix(orig_info$scalar, nrow = orig_info$height, ncol = orig_info$width)

cat("Reading noised_matrix.png ...\n")
noised_info <- read_png_as_continuous_vector("noised_matrix.png", anchor_path)
noised_matrix <- matrix(noised_info$scalar, nrow = noised_info$height, ncol = noised_info$width)

lambda_vals <- c(0.1, 0.5, 1, 1.5, 2)
num_iters  <- 4

run_irfl(
  original_matrix = original_matrix,
  noised_matrix   = noised_matrix,
  num_iters       = 6,
  lambda_grid     = lambda_vals,
  eps             = 5e-3,
  w_min           = 0.1,
  w_max           = 1e3,
  eta             = 0.5,
  rho             = 1.0
)

make_irfl_grid_pdf(lambda_vals = lambda_vals, num_iters = num_iters, pdf_path = "irfl_grid_final.pdf")

# Edge overlay and 1x3 PDF (only if the expected denoised PNG exists)
denoised_path <- "lambda_1_iter_4.png"
if (!file.exists(denoised_path)) {
  cat(sprintf("Skipping edge overlay: missing file %s\n", denoised_path))
} else {
  cat(sprintf("Building edge overlay from %s ...\n", denoised_path))
  
  denoised_info <- read_png_as_continuous_vector(denoised_path, anchor_path)
  lambda1_matrix <- matrix(denoised_info$scalar, nrow = denoised_info$height, ncol = denoised_info$width)
  
  edges_df <- extract_edges(lambda1_matrix, nrow(lambda1_matrix), ncol(lambda1_matrix))
  
  p_edges <- plot_with_edges(noised_matrix, edges_df) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
  
  img_width  <- ncol(original_matrix)
  img_height <- nrow(original_matrix)
  
  ggsave(
    filename = "noised_edges.png",
    plot = p_edges,
    width = img_width, height = img_height, units = "px",
    dpi = 300, bg = "transparent", limitsize = FALSE
  )
  cat("Wrote noised_edges.png\n")
  
  g1 <- rasterGrob(readPNG("original_matrix.png"), interpolate = FALSE)
  g2 <- rasterGrob(readPNG("noised_matrix.png"),   interpolate = FALSE)
  g3 <- rasterGrob(readPNG("noised_edges.png"),    interpolate = FALSE)
  
  dpi <- 300
  pdf_width  <- 3 * img_width  / dpi
  pdf_height <-     img_height / dpi
  
  pdf("lambda1_grid_images_only.pdf", width = pdf_width, height = pdf_height, useDingbats = FALSE)
  grid.arrange(g1, g2, g3, nrow = 1, widths = unit(c(1, 1, 1), "null"))
  dev.off()
  
  cat("Wrote lambda1_grid_images_only.pdf\n")
}
