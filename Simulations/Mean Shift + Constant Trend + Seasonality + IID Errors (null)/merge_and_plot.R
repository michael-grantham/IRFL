
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lpSolve)
library(patchwork)
library(tidyverse)
library(cowplot)  # Used for stable legend extraction

parent_path = getwd() 

## Bring in distance metric
cpt.dist = function(C1, C2, N){
  #C1 = c(26,51)
  #C2 = c(25,78,99)
  #C2 = c(26,51)
  #C1 = c(25,78,99)
  m = length(C1)
  k = length(C2)
  
  ##Generate Cost Matrix via all paired distance
  pair = expand.grid(C1, C2)
  if(m==k){ 
    cost.mat = matrix(abs(pair[,1]-pair[,2]), 
                      nrow=m,ncol=k,byrow=T )
  }else if(m > k){ #C1 has more changepoints than C2
    cost.mat = cbind(matrix(abs(pair[,1]-pair[,2]), 
                            nrow=m,ncol=k,byrow=T ), 
                     matrix(0, nrow=m, ncol=(m-k), byrow=T))
  }else{ #C1 has less changepoints than C2, nrow < ncol
    cost.mat = rbind(matrix(abs(pair[,1]-pair[,2]), 
                            nrow=m,ncol=k,byrow=F ), 
                     matrix(0, nrow=(k-m), ncol=k, byrow= T ))
  }
  cpt.asgn = lp.assign(cost.mat, direction = "min")
  cpt.asgn$objval/N + abs(m-k)
}

## Bring in function to apply distance
apply_distance = function(df, col1, col2) {
  # Initialize an empty vector to store distances
  distances <- numeric(nrow(df))
  
  # Loop over each row in the dataframe
  for (i in seq_len(nrow(df))) {
    # Split the strings into character vectors for the given row
    vec1 <- unlist(strsplit(as.character(df[[col1]][i]), ","))
    vec2 <- unlist(strsplit(as.character(df[[col2]][i]), ","))
    
    # Convert the character vectors to numeric
    num_vec1 <- as.numeric(vec1)
    num_vec2 <- as.numeric(vec2)
    
    # Compute the distance using cpt.dist() and store in the distances vector
    distances[i] <- cpt.dist(num_vec1, num_vec2,1200)
  }
  
  return(distances)
}

## Get all subfolders where the results are stored to
subfolders = list.dirs(parent_path, recursive=T)[-1] 
folder.path = file.path(subfolders)

## Initialize empty data frame for storing merged results
#MCPT.bind = data.frame(matrix(0, nrow=5, ncol=22, byrow=T))

out.bind = NULL

for (folder in folder.path){
  # obtain file.path of csv's
  if(length(list.files(folder,pattern="*.csv"))>=1){
    file.path = paste(folder, list.files(folder, pattern="*.csv"), sep="/")
    out.ind = read.csv(file.path)
    out.bind = rbind(out.bind, out.ind)
  }
}

# Print out the merged result before performing calculations
out = out.bind[,-1]
write.csv(out,"outbind.csv",row.names=FALSE)

###### Produce nknots Plots
out = read.csv("outbind.csv")

sel = c("bic_adaflasso", paste("bic_iter", 2:6, "flasso", sep=""))
sel2 = c("adaflasso", paste("iter", 2:6, "flasso", sep=""))
sel3 = c("nknots_adaflasso", paste("nknots_iter", 2:6, "flasso", sep=""))

min_bic_iter = out |> select(all_of(sel)) |> round(5) |> apply(1, which.min)
min_bic = out |> select(all_of(sel)) |> round(5) |> apply(1, min)
flasso = out |> select(all_of(sel3))
min_bic_iter_flasso = numeric(nrow(flasso))

for (k in 1:nrow(flasso)) {
  min_bic_iter_flasso[k] = flasso[k, min_bic_iter[k]]
}

final = out |> 
  select(c("true_model", "bic_true_model", "nknots_flasso", "bic_flasso", "nknots_adaflasso", "bic_adaflasso")) |> 
  mutate(best_iter = min_bic_iter,
         nknots_best_iter = min_bic_iter_flasso,
         min_bic_flasso = min_bic)

df2 = final %>% 
  pivot_longer(cols = c(3, 5, 8), names_to = "Method", values_to = "Mean Shifts") %>%
  mutate(`Mean Shifts` = `Mean Shifts`,
         Label = factor(Method,
                        levels = c("nknots_flasso", "nknots_adaflasso", "nknots_best_iter"),
                        labels = c("FL", "Ada FL", "IRFL")))

knots = ggplot(df2, aes(x = Label, y = `Mean Shifts`, fill = Label)) +
  geom_boxplot() +
  scale_fill_manual(values = c("FL" = "yellow", "Ada FL" = "lightblue", "IRFL" = "lightcoral"),
                    breaks = c("FL", "Ada FL", "IRFL")) +
  scale_y_continuous(breaks = c(0, 1)) +  # <--- Only show ticks at 0 and 1
  labs(y = expression(hat(m) - m), x = NULL, fill = "Method") +
  geom_hline(yintercept = 0, color = "red", linewidth = .75)


###### Produce Simulation Plot
n = 1200
x = time = 1:n
s = 12
alpha = 0.005
mu = rep(0, n)

set.seed(123456)
Z = rnorm(s)
W = cumsum(Z)
t_grid = (1:s) / s
B = W - t_grid * W[s]  
B = B - mean(B)
seasonal_effect = rep(B, 100)
y = mu + alpha * 1:n + seasonal_effect + rnorm(n)

per = ggplot(NULL, aes(x = 1:36)) + 
  geom_line(aes(y = rep(B, 3))) + 
  labs(x = "t", y = expression(y[t])) + 
  scale_x_continuous(limits = c(0, 36), breaks = seq(0, 36, by = 12))

sim = ggplot(NULL, aes(x = time)) +
  geom_line(aes(y = y), alpha = 0.25) + 
  geom_line(aes(y = mu + alpha * 1:n), linewidth = .75, color = "red") +
  labs(x = "t", y = expression(y[t])) + 
  scale_x_continuous(limits = c(0, 1200), breaks = seq(0, 1200, by = 300))

#### Season plot
sel3 = c("season_adaflasso", paste("season_iter", 2:6, "flasso", sep=""))
sel4 = c("alpha_adaflasso", paste("alpha_iter", 2:6, "flasso", sep=""))
seasons = out |> select(all_of(sel3))
slopes = out |> select(all_of(sel4))

season_best_iter = numeric(nrow(out))
slope_best_iter = numeric(nrow(out))
for (k in 1:nrow(out)) {
  season_best_iter[k] = seasons[k, min_bic_iter[k]]
  slope_best_iter[k] = slopes[k, min_bic_iter[k]]
}

out2 = out |> select(c("true.season", "season_flasso", "season_adaflasso")) |> mutate(season_best_iter = season_best_iter)
true_season = strsplit(out2$true.season, ",") |> lapply(as.numeric)
flasso_season = strsplit(out2$season_flasso, ",") |> lapply(as.numeric)
adaflasso_season = strsplit(out2$season_adaflasso, ",") |> lapply(as.numeric)
best_iter_season = strsplit(out2$season_best_iter, ",") |> lapply(as.numeric)

flasso_df = matrix(NA, 11, 1201) |> as.data.frame()
ada_df = matrix(NA, 11, 1201) |> as.data.frame()
best_df = matrix(NA, 11, 1201) |> as.data.frame()

flasso_df[,1] = rep("FL", 11)
ada_df[,1] = rep("Ada FL", 11)
best_df[,1] = rep("IRFL", 11)

for (k in 1:1200) {
  flasso_df[, k+1] = true_season[[k]][1:11] - flasso_season[[k]]
  ada_df[, k+1] = true_season[[k]][1:11] - adaflasso_season[[k]]
  best_df[, k+1] = true_season[[k]][1:11] - best_iter_season[[k]]
}

cnames = c("Label", paste0("Rep", 1:1200))
df = rbind(flasso_df, ada_df, best_df)
colnames(df) = cnames
df = df |> mutate(s = factor(0:32 %% 11 + 1)) |>
  pivot_longer(cols = 2:1201, names_to = "Rep", values_to = "Value")

seasonality = ggplot(df, aes(x = s, y = Value, fill = Label)) +
  geom_boxplot(position = position_dodge(width = 0.8)) + 
  scale_fill_manual(values = c("FL" = "yellow", "Ada FL" = "lightblue", "IRFL" = "lightcoral"),
                    breaks = c("FL", "Ada FL", "IRFL")) +
  labs(y = expression(hat(s)[i] - s[i]), x = "i", fill = "Method") +
  geom_hline(yintercept = 0, color = "red", size = .75)

#### Trend Plot
df = out |> 
  select(c("true.alpha", "alpha_flasso", "alpha_adaflasso")) |> 
  mutate(alpha_IRFL = slope_best_iter)

dat = df %>%
  pivot_longer(cols = 2:4, names_to = "method", values_to = "slope") %>%
  mutate(Label = factor(method,
                        levels = c("alpha_flasso", "alpha_adaflasso", "alpha_IRFL"),
                        labels = c("FL", "Ada FL", "IRFL")))

slope = ggplot(dat, aes(x = Label, y = slope - 0.005, fill = Label)) +
  geom_boxplot() +
  scale_fill_manual(values = c("FL" = "yellow", "Ada FL" = "lightblue", "IRFL" = "lightcoral"),
                    breaks = c("FL", "Ada FL", "IRFL")) +
  labs(y = expression(hat(alpha) - alpha), x = NULL, fill = "Method") +
  geom_hline(yintercept = 0, color = "red", size = .75)

#### Final Plots with Correct Legend Placement
sim_clean = sim + theme(legend.position = "none")
per_clean = per + theme(legend.position = "none")
knots_clean = knots + theme(legend.position = "none")
seasonality_clean = seasonality + theme(legend.position = "none")
slope_clean = slope + theme(legend.position = "none")

legend_only = cowplot::get_legend(
  knots + theme(legend.position = "right",
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 12),
                legend.key.size = unit(.75, "cm"))
)
legend_plot = patchwork::wrap_elements(legend_only)

simplot = sim_clean + per_clean + plot_layout(ncol = 2)
grouppanel = (knots_clean + seasonality_clean) / (slope_clean + legend_plot) + plot_layout(guides = "collect")
combined = sim_clean + per_clean + knots_clean + seasonality_clean + slope_clean + legend_plot + plot_layout(ncol = 2)

ggsave("Seasonality Null Simulation Plot.pdf", plot = simplot, width = 14*.5, height = 3*.5, dpi = 300)
ggsave("Seasonality Null Group Plot.pdf", plot = grouppanel, width = 14*.5, height = 6*.5, dpi = 300)
ggsave("Seasonality Null Combined Plot.pdf", plot = combined, width = 14*.5, height = 9*.5, dpi = 300)


