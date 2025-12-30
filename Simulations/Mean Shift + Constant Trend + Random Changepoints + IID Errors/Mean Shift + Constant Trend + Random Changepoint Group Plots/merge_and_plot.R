
library(data.table)  # fread(), rbindlist(), fwrite()
library(dplyr)       # mutate(), %>%
library(tidyr)       # pivot_longer()
library(ggplot2)     # plotting
library(patchwork)   # plot composition
library(grid)        # unit() for legend sizing

cwd    <- normalizePath(getwd(), mustWork = TRUE)
parent <- dirname(cwd)

files <- unlist(lapply(list.dirs(parent, recursive = FALSE, full.names = TRUE), function(d)
  if (d != cwd) list.files(d, "^final_out(\\.csv)?$", recursive = TRUE, full.names = TRUE)
))

out <- if (length(files)) rbindlist(lapply(files, fread), fill = TRUE) else NULL

data = out

# Print out the merged result before performing calculations
write.csv(data,"combined_FTWI_alpha_results.csv",row.names=FALSE)

# Plot of slope
slope_plot = data |> 
  mutate(ada_flasso = slope_OLS_adaflasso,
         best_iter = slope_best_iter) |>
  pivot_longer(cols = c(6,12),names_to = "method",
               values_to = "alphahat") |>
  ggplot(aes(x = factor(true_slope, levels = c(-0.01, -0.005, 0.005, 0.01)), 
             y = alphahat-true_slope, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) + 
  scale_fill_manual(
    values = c("slope_OLS_adaflasso" = "lightblue", 
               "slope_best_iter" = "lightcoral"),
    labels = c("slope_best_iter" = "IRFL", 
               "slope_OLS_adaflasso" = "Ada FL") 
  ) +
  geom_hline(yintercept = 0, color = "red", size = .75) +
  labs(x = expression(alpha), y = expression(hat(alpha) - alpha), fill = "Method") +
  theme(axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = "none")

# Plot of Number of Knots
nknots_plot = data |> 
  mutate(ada_flasso = nknots_adaflasso,
         best_iter = nknots_best_iter) |>
  pivot_longer(cols = c(7,11),names_to = "method",
               values_to = "mhat") %>%
  ggplot(aes(x = factor(true_slope, levels = c(-0.01, -0.005, 0.005, 0.01)), 
             y = mhat-true_nknots, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) + 
  scale_fill_manual(
    values = c("nknots_adaflasso" = "lightblue", 
               "nknots_best_iter" = "lightcoral"),
    labels = c("nknots_best_iter" = "IRFL", 
               "nknots_adaflasso" = "Ada FL") 
  ) +
  geom_hline(yintercept = 0, color = "red", size = .75) +
  labs(x = expression(alpha), y = expression(hat(m)-m), fill = "Method") +
  theme(axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = "none")

# Plot of Distance
distance_plot = data |> 
  pivot_longer(cols = c(14,16),names_to = "method",
               values_to = "distance") %>% 
  ggplot(aes(x = factor(true_slope, levels = c(-0.01, -0.005, 0.005, 0.01)), 
             y = distance, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) + 
  scale_fill_manual(
    values = c("distance_adaflasso" = "lightblue", 
               "distance_best_iter_flasso" = "lightcoral"),
    labels = c("distance_best_iter_flasso" = "IRFL", 
               "distance_adaflasso" = "Ada FL") 
  ) +
  geom_hline(yintercept = 0, color = "red", size = .75) +
  labs(x = expression(alpha), y = "Distance", fill = "Method") +
  theme(axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = "none")


# Simulation graphs

# Function to generate random cps
generate_cp_vector = function(N, m, d) {
  
  # Check to see if possible
  if(N-(m+1)*d<m){
    cat("Not possible if N-(m+1)d<m. Try decreasing m or d.")
    break;
  }
  
  # Step 1: Set R
  R = N - (m + 1) * d
  
  # Step 2 & 3: Generate m+1 random parts that sum to R
  random_parts = function(total, n) {
    parts = runif(n)                # Generate random numbers
    parts = parts / sum(parts)       # Normalize to sum to 1
    parts = round(parts * total)     # Scale to sum to `total` and round
    diff = total - sum(parts)        # Adjust to ensure the sum is exactly `total`
    
    # Set up a variable to use in while loop to enforce minimum distance = d
    proceed = FALSE
    # Distribute any rounding differences
    if (diff != 0) { 
      counter = 1
      while(proceed==FALSE){
        test = parts
        for (i in 1:abs(diff)) {
          idx = sample(1:n, 1)
          test[idx] = test[idx] + sign(diff)
        }
        proceed = !any(parts<0)
        counter = counter + 1
        if(counter==1000){
          return(rep(NA,m+1))
          break
        }
      }
      parts[idx] = parts[idx] + sign(diff)
    }
    return(parts)
  }
  
  # Generate m+1 random parts that sum to R
  r_values = random_parts(R, m + 1)
  
  # Step 4: Calculate the `out` vector by summing adjusted values
  semi = r_values + d
  out = cumsum(semi)
  out = out[1:(length(out) - 1)]  # Remove the last element
  out = c(0,out,N)
  
  return(out)
}

# Function to generate mean vector
generate_means = function(cp.vec){
  # Set up difference and mean vec
  v = c(0)
  means = c()
  
  for (i in 2:(length(cp.vec)-1)) {
    step = runif(1, min = 1.5, max = 2.5)  # Generate a random step
    direction =sample(c(-1, 1), 1)  # Randomly choose -1 or 1
    v[i] = v[i-1] + direction * step  # Update the vector
  }
  
  for(k in 2:length(cp.vec)){
    means = c(means,rep(v[k-1],as.numeric(cp.vec[k]-cp.vec[k-1])))
  }
  return(means)
}


# Set up alphas
n = 1000
time = 1:n
alpha = c(-.01,-.005,0.005,.01)

set.seed(134)
# Simulate 
m = sample(1:7,1)
g = generate_cp_vector(n, m, d=100)
cpts = g[-c(1,length(g))]+1
true_model = generate_means(g)
y1 = rnorm(n,true_model) + alpha[1]*time
y2 = rnorm(n,true_model) + alpha[2]*time
y3 = rnorm(n,true_model) + alpha[3]*time
y4 = rnorm(n,true_model) + alpha[4]*time

sim1 = ggplot(NULL, aes(x = time)) +
  geom_line(aes(y = y1), alpha = 0.25) + 
  geom_line(aes(y = true_model+alpha[1]*time),linewidth=.75,color="red")+
  # geom_hline(yintercept = 0, color = "red", size = .75) + 
  labs(x = "t",
       y = expression(y[t])) + 
  ggtitle(expression(alpha == -0.01)) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_y_continuous(limits = c(min(y1,y2,y3,y4),max(y1,y2,y3,y4)))

sim2 = ggplot(NULL, aes(x = time)) +
  geom_line(aes(y = y2), alpha = 0.25) + 
  geom_line(aes(y = true_model+alpha[2]*time),linewidth=.75,color="red")+
  # geom_hline(yintercept = 0, color = "red", size = .75) + 
  labs(x = "t",
       y = expression(y[t])) + 
  ggtitle(expression(alpha == -0.005)) + 
  theme(plot.title = element_text(hjust = 0.5))+ 
  scale_y_continuous(limits = c(min(y1,y2,y3,y4),max(y1,y2,y3,y4)))

sim3 = ggplot(NULL, aes(x = time)) +
  geom_line(aes(y = y3), alpha = 0.25) + 
  geom_line(aes(y = true_model+alpha[3]*time),linewidth=.75,color="red")+
  # geom_hline(yintercept = 0, color = "red", size = .75) + 
  labs(x = "t",
       y = expression(y[t])) + 
  ggtitle(expression(alpha == 0.005)) + 
  theme(plot.title = element_text(hjust = 0.5))+ 
  scale_y_continuous(limits = c(min(y1,y2,y3,y4),max(y1,y2,y3,y4)))

sim4 = ggplot(NULL, aes(x = time)) +
  geom_line(aes(y = y4), alpha = 0.25) + 
  geom_line(aes(y = true_model+alpha[4]*time),linewidth=.75,color="red")+
  # geom_hline(yintercept = 0, color = "red", size = .75) + 
  labs(x = "t",
       y = expression(y[t])) + 
  ggtitle(expression(alpha == 0.01)) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_y_continuous(limits = c(min(y1,y2,y3,y4),max(y1,y2,y3,y4)))


# Save plots

sig = sim1 + sim2 + sim3 + sim4 + plot_layout(ncol=2)

compare = slope_plot + nknots_plot + 
  distance_plot + plot_layout(ncol=2) + 
  guides(color = guide_legend(keyheight = 4, keywidth = 4)) +
  theme(legend.position = c(1.6, 0.5),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(1, "cm"))

ggsave("Slope Simulations.pdf", 
       plot = sig, device = "pdf", width = 14*.5, height = 6*.5, dpi = 300)
ggsave("Slope Comparisons.pdf", 
       plot = compare, device = "pdf", width = 14*.5, height = 6*.5, dpi = 300)
