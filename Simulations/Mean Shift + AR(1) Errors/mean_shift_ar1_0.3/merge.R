#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lpSolve)
library(patchwork)

parent_path <- getwd()

## ---- Distance metric (with NA/empty-string robustness) ----
cpt.dist <- function(C1, C2, N) {
  # If C2 is NA or empty, return count of numbers in C1
  if (length(C2) == 0 || any(is.na(C2))) {
    return(length(C1))
  }
  
  m <- length(C1)
  k <- length(C2)
  
  # Guard: if C1 empty, distance is just |m-k| = k (since m=0)
  if (m == 0) return(k)
  
  pair <- expand.grid(C1, C2)
  if (m == k) {
    cost.mat <- matrix(abs(pair[, 1] - pair[, 2]), nrow = m, ncol = k, byrow = TRUE)
  } else if (m > k) {
    cost.mat <- cbind(
      matrix(abs(pair[, 1] - pair[, 2]), nrow = m, ncol = k, byrow = TRUE),
      matrix(0, nrow = m, ncol = (m - k), byrow = TRUE)
    )
  } else {
    cost.mat <- rbind(
      matrix(abs(pair[, 1] - pair[, 2]), nrow = m, ncol = k, byrow = FALSE),
      matrix(0, nrow = (k - m), ncol = k, byrow = TRUE)
    )
  }
  
  cpt.asgn <- lp.assign(cost.mat, direction = "min")
  cpt.asgn$objval / N + abs(m - k)
}

## Parse a CP string like "10, 50, 120" into numeric vector, treating "" as NA/empty
.parse_cp <- function(x) {
  if (is.na(x) || length(x) == 0) return(numeric(0))
  parts <- unlist(strsplit(as.character(x), ","))
  parts <- trimws(parts)
  # Treat empty tokens as NA; drop NAs
  parts[parts == ""] <- NA_character_
  as.numeric(parts[!is.na(parts)])
}

## ---- Apply distance (now accepts a normalization N; defaults to 1000) ----
apply_distance <- function(df, col1, col2, N_norm = 1000) {
  distances <- numeric(nrow(df))
  for (i in seq_len(nrow(df))) {
    v1 <- .parse_cp(df[[col1]][i])
    v2 <- .parse_cp(df[[col2]][i])
    
    if (length(v2) == 0) {
      distances[i] <- length(v1)
    } else {
      distances[i] <- cpt.dist(v1, v2, N_norm)
    }
  }
  distances
}

## ---- Collect all subfolders and read CSVs (one per folder is fine; robust anyway) ----
subfolders <- list.dirs(parent_path, recursive = TRUE)[-1]
folder_paths <- file.path(subfolders)

out.bind <- NULL
for (folder in folder_paths) {
  files <- list.files(folder, pattern = "[.]csv$", full.names = TRUE)
  if (length(files) >= 1) {
    # You said there's only one per folder; this still works if there are more.
    this <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
    out.bind <- rbind(out.bind, this)
  }
}

out <- out.bind
if (ncol(out) > 1) {
  first_col <- out[[1]]
  # Check if the first column is sequential integers 1:nrow
  if (is.numeric(first_col) && all(first_col == seq_len(nrow(out)))) {
    out <- out[, -1, drop = FALSE]
  }
}

# Save outbind
write.csv(out, "outbind.csv", row.names = FALSE)

## ---- Selections for min-BIC iteration among flasso variants ----
sel  <- c("bic_adaflasso", paste0("bic_iter", 2:10, "flasso"))
sel2 <- c("adaflasso",    paste0("iter", 2:10, "flasso"))
sel3 <- c("phi_yw_adaflasso", paste0("phi_yw_iter", 2:10, "flasso"))
sel4 <- c("phi_adaflasso",    paste0("phi_iter", 2:10, "flasso"))
sel5 <- c("nknots_adaflasso", paste0("nknots_iter", 2:10, "flasso"))

min_bic_iter <- out |>
  dplyr::select(dplyr::all_of(sel)) |>
  round(5) |>
  apply(1, which.min)

min_bic <- out |>
  dplyr::select(dplyr::all_of(sel)) |>
  round(5) |>
  apply(1, min)

flasso         <- out |> dplyr::select(dplyr::all_of(sel2))
phi_yw_flasso  <- out |> dplyr::select(dplyr::all_of(sel3))
phi_flasso     <- out |> dplyr::select(dplyr::all_of(sel4))
nknots_flasso  <- out |> dplyr::select(dplyr::all_of(sel5))

min_bic_iter_flasso <- numeric(nrow(flasso))
best_iter_phi       <- numeric(nrow(flasso))
best_iter_phi_yw    <- numeric(nrow(flasso))
best_iter_nknots    <- numeric(nrow(flasso))

for (k in seq_len(nrow(flasso))) {
  idx <- min_bic_iter[k]
  min_bic_iter_flasso[k] <- flasso[k, idx, drop = TRUE]
  best_iter_phi[k]       <- phi_flasso[k, idx, drop = TRUE]
  best_iter_phi_yw[k]    <- phi_yw_flasso[k, idx, drop = TRUE]
  best_iter_nknots[k]    <- nknots_flasso[k, idx, drop = TRUE]
}


N_norm <- 1000

## ---- Assemble semifinal with same columns (and order later) ----
semifinal <- out |>
  dplyr::select(c(
    "true_model", "true_phi",
    "true_nknots", "bic_true_model",
    "adaflasso", "phi_adaflasso",
    "phi_yw_adaflasso", "nknots_adaflasso",
    "bic_adaflasso",
    "knots_WCM",
    "nknots_WCM",
    "bic_WCM",
    "phi_WCM",
    "knots_AR1seg",
    "nknots_AR1seg",
    "bic_AR1seg",
    "phi_AR1seg"
  )) |>
  dplyr::mutate(
    best_iter          = min_bic_iter,
    changepoints       = min_bic_iter_flasso,
    nknots_best_iter   = best_iter_nknots,
    phi_best_iter      = best_iter_phi,
    phi_yw_best_iter   = best_iter_phi_yw,
    min_bic_best_iter  = min_bic,
    distance_adaflasso = apply_distance(out, "true_model", "adaflasso",  N_norm),
    # keep these two (no-op) to preserve your pipeline; harmless
    phi_adaflasso      = phi_adaflasso,
    phi_yw_adaflasso   = phi_yw_adaflasso,
    distance_WCM       = apply_distance(out, "true_model", "knots_WCM",  N_norm),
    distance_AR1seg    = apply_distance(out, "true_model", "knots_AR1seg", N_norm)
  )

final <- semifinal |>
  dplyr::mutate(
    distance_best_iter_flasso = apply_distance(semifinal, "true_model", "changepoints", N_norm)
  )

## ---- Write in the exact order you specified ----
final <- final |>
  dplyr::select(
    "true_model","true_phi","true_nknots","bic_true_model",
    "adaflasso","phi_adaflasso","phi_yw_adaflasso","nknots_adaflasso","bic_adaflasso",
    "knots_WCM","nknots_WCM","bic_WCM","phi_WCM",
    "knots_AR1seg","nknots_AR1seg","bic_AR1seg","phi_AR1seg",
    "best_iter","changepoints","nknots_best_iter","phi_best_iter","phi_yw_best_iter","min_bic_best_iter",
    "distance_adaflasso","distance_WCM","distance_AR1seg","distance_best_iter_flasso"
  )

write.csv(final, "final_out_0.3.csv", row.names = FALSE)


# Produce Graphs (optional)

data = final

phi = data$true_phi[1]


# Plot of slope
phi_plot = data |> 
  mutate(ada_flasso = phi_adaflasso,
         YW_ada_flasso = phi_yw_adaflasso,
         best_iter = phi_best_iter,
         YW_best_iter = phi_yw_best_iter,
         AR1seg = phi_AR1seg,
         WCM = phi_WCM) |> 
  pivot_longer(cols = c(18,28:32),names_to = "method",
               values_to = "phihat") |>
  mutate(method = factor(method, levels = c("ada_flasso","YW_ada_flasso",
                                            "best_iter","YW_best_iter","WCM",
                                            "AR1seg")))|>
  ggplot(aes(x = method, 
             y = phihat-true_phi, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) + 
  scale_fill_manual(
    values = c("ada_flasso" = "lightblue",
               "YW_ada_flasso" = "blue",
               "best_iter" = "lightcoral",
               "YW_best_iter" = "coral",
               "AR1seg" = "purple",
               "WCM" = "lightgreen"),
    labels = c("best_iter" = "IRFL",
               "YW_best_iter" = "IRFL (YW)",
               "ada_flasso" = "Ada FL",
               "YW_ada_flasso" = "Ada FL (YW)",
               "AR1seg" = "AR1Seg",
               "WCM" = "WCM")
  ) +
  geom_hline(yintercept = 0, color = "red", size = .75) +
  labs(x = NULL, y = expression(hat(phi) - phi), fill = "Method") +
  theme(axis.text.x = element_blank(),
        #axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = "none"
  ) 

# Plot of Number of Knots
nknots_plot = data |> 
  mutate(`Ada FL` = nknots_adaflasso,
         `IRFL` = nknots_best_iter,
         `AR1Seg` = nknots_AR1seg,
         `WCM` = nknots_WCM) |> 
  pivot_longer(cols = c(28:31),names_to = "method",
               values_to = "mhat") |>
  mutate(method = factor(method,
                         levels = c("Ada FL",
                                    "IRFL",
                                    "WCM",
                                    "AR1Seg")))|>
  ggplot(aes(x = method, 
             y = mhat-true_nknots, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) + 
  scale_fill_manual(
    values = c("Ada FL" = "lightblue",
               "IRFL" = "lightcoral",
               "AR1Seg" = "purple",
               "WCM" = "lightgreen"),
    labels = c("IRFL" = "IRFL",
               "Ada FL" = "Ada FL",
               "AR1Seg" = "AR1Seg",
               "WCM" = "WCM")
  ) +
  geom_hline(yintercept = 0, color = "red", size = .75) +
  labs(x=NULL, y = expression(hat(m)-m), fill = "Method") +
  theme(axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = "none") 

# Plot of Distance
distance_plot = data |> 
  pivot_longer(cols = 25:27,names_to = "method",
               values_to = "distance") |>
  mutate(method = factor(method,levels=c("distance_best_iter_flasso",
                                         "distance_WCM",
                                         "distance_AR1seg"))) |>
  ggplot(aes(x = method, 
             y = distance, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) + 
  scale_fill_manual(
    values = c("distance_best_iter_flasso" = "lightcoral",
               "distance_AR1seg" = "purple",
               "distance_WCM" = "lightgreen"),
    labels = c("distance_best_iter_flasso" = "IRFL",
               "distance_AR1seg" = "AR1Seg",
               "distance_WCM" = "WCM") 
  ) +
  geom_hline(yintercept = 0, color = "red", size = .75) +
  labs(x = NULL, y = "Distance", fill = "Method") +
  theme(axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = "none")+
  scale_x_discrete(
    labels=c("distance_best_iter_flasso" = "IRFL",
             "distance_AR1seg" = "AR1Seg",
             "distance_WCM" = "WCM")
  )


# Simulation graphs

# Set up phis
n = 1000
time = 1:n
phi = c(-.7,-.3,0,.3,.7)
true_model = rep(c(0,2,0,2),each=250)

y_sim = lapply(seq_along(phi),function(i){
  set.seed(134)
  y = numeric(n)
  y = rnorm(1)
  for(k in 2:n){
    y[k] = phi[i]*y[k-1]+rnorm(1)
  }
  y = y + true_model
  return(y)
})


sim1 = ggplot(NULL, aes(x = time)) +
  geom_line(aes(y = y_sim[[1]]), alpha = 0.25) + 
  geom_line(aes(y = true_model),linewidth=.75,color="red")+
  # geom_hline(yintercept = 0, color = "red", size = .75) + 
  labs(x = "t",
       y = expression(y[t])) + 
  ggtitle(expression(phi == -0.7)) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_y_continuous(limits = c(min(y_sim[[1]],y_sim[[2]],y_sim[[3]],y_sim[[4]]),
                                max(y_sim[[1]],y_sim[[2]],y_sim[[3]],y_sim[[4]])))

sim2 = ggplot(NULL, aes(x = time)) +
  geom_line(aes(y = y_sim[[2]]), alpha = 0.25) + 
  geom_line(aes(y = true_model),linewidth=.75,color="red")+
  # geom_hline(yintercept = 0, color = "red", size = .75) + 
  labs(x = "t",
       y = expression(y[t])) + 
  ggtitle(expression(phi == -0.3)) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_y_continuous(limits = c(min(y_sim[[1]],y_sim[[2]],y_sim[[3]],y_sim[[4]]),
                                max(y_sim[[1]],y_sim[[2]],y_sim[[3]],y_sim[[4]])))

sim3 = ggplot(NULL, aes(x = time)) +
  geom_line(aes(y = y_sim[[3]]), alpha = 0.25) + 
  geom_line(aes(y = true_model),linewidth=.75,color="red")+
  # geom_hline(yintercept = 0, color = "red", size = .75) + 
  labs(x = "t",
       y = expression(y[t])) + 
  ggtitle(expression(phi == 0)) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_y_continuous(limits = c(min(y_sim[[1]],y_sim[[2]],y_sim[[3]],y_sim[[4]]),
                                max(y_sim[[1]],y_sim[[2]],y_sim[[3]],y_sim[[4]])))

sim4 = ggplot(NULL, aes(x = time)) +
  geom_line(aes(y = y_sim[[4]]), alpha = 0.25) + 
  geom_line(aes(y = true_model),linewidth=.75,color="red")+
  # geom_hline(yintercept = 0, color = "red", size = .75) + 
  labs(x = "t",
       y = expression(y[t])) + 
  ggtitle(expression(phi == 0.3)) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_y_continuous(limits = c(min(y_sim[[1]],y_sim[[2]],y_sim[[3]],y_sim[[4]]),
                                max(y_sim[[1]],y_sim[[2]],y_sim[[3]],y_sim[[4]])))

sim5 = ggplot(NULL, aes(x = time)) +
  geom_line(aes(y = y_sim[[5]]), alpha = 0.25) + 
  geom_line(aes(y = true_model),linewidth=.75,color="red")+
  # geom_hline(yintercept = 0, color = "red", size = .75) + 
  labs(x = "t",
       y = expression(y[t])) + 
  ggtitle(expression(phi == 0.7)) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_y_continuous(limits = c(min(y_sim[[1]],y_sim[[2]],y_sim[[3]],y_sim[[4]]),
                                max(y_sim[[1]],y_sim[[2]],y_sim[[3]],y_sim[[4]])))


# Set up legend plot
custom_legend = ggplot() +
  # Adding colored points
  geom_point(aes(x = .75, y = 6, color = "lightblue"), size = 5, shape=15) +
  geom_point(aes(x = .75, y = 5, color = "blue"), size = 5, shape=15) +
  geom_point(aes(x = .75, y = 4, color = "lightcoral"), size = 5, shape=15) +
  geom_point(aes(x = .75, y = 3, color = "coral"), size = 5, shape=15) +
  geom_point(aes(x = .75, y = 2, color = "lightgreen"), size = 5, shape=15) +
  geom_point(aes(x = .75, y = 1, color = "purple"), size = 5, shape=15) +
  scale_color_manual(values = c("blue","coral","lightblue","lightcoral","lightgreen","purple")) +
  
  # Adding labels next to the points
  geom_text(aes(x = .8, y = 6, label = "Ada FL"), hjust = 0, size = 5) +
  geom_text(aes(x = .8, y = 5, label = "Ada FL (YW)"), hjust = 0, size = 5) +
  geom_text(aes(x = .8, y = 4, label = "IRFL"), hjust = 0, size = 5) +
  geom_text(aes(x = .8, y = 3, label = "IRFL (YW)"), hjust = 0, size = 5) +
  geom_text(aes(x = .8, y = 2, label = "WCM"), hjust = 0, size = 5) +
  geom_text(aes(x = .8, y = 1, label = "AR1Seg"), hjust = 0, size = 5) +
  
  scale_x_continuous(limits = c(0, 2)) + 
  scale_y_continuous(limits = c(0, 7)) +
  theme_void() +
  theme(legend.position = "none") +  # No default legend
  theme(legend.title = element_text(size = 14))

# Save plots
sig = sim2 + sim1 + sim4 + sim5 + plot_layout(ncol=2)

compare = phi_plot + nknots_plot + 
  distance_plot + custom_legend + plot_layout(ncol=2) +
  plot_annotation(title = expression(phi==0.3),
                  theme = theme(
                    plot.title = element_text(hjust = 0.5) # Center the title
                  ))


ggsave("Fixed Trend Simulations.pdf", 
       plot = sig, device = "pdf", width = 14, height = 6, dpi = 300)
ggsave("Fixed Trend Comparisons (phi=0.3).pdf", 
       plot = compare, device = "pdf", width = 14, height = 6, dpi = 300)





