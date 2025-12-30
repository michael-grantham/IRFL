library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)


###### Produce Distance Plots
final = read.csv("final_out.csv")

df = final %>% 
  pivot_longer(cols=16:21,names_to="method",values_to="distance") %>%
  mutate(method = factor(method,
                         levels = c("distance_flasso",
                                    "distance_adaflasso",
                                    "distance_best_iter_flasso",
                                    "distance_bs",
                                    "distance_wbs",
                                    "distance_pelt")))

dist = ggplot(df, aes(x = method, y = distance)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("distance_flasso"="FL",
                              "distance_adaflasso"="Ada FL",
                              "distance_best_iter_flasso"="IRFL",
                              "distance_bs"="BS",
                              "distance_wbs"="WBS",
                              "distance_pelt"="PELT")) +
  labs(y = "Distance", x = NULL)  +
  geom_hline(yintercept = 0, color = "red", size = .75) 
 # theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))

###### Produce nknots Plots
out = read.csv("outbind.csv")

sel = c("bic_adaflasso",paste("bic_iter",2:10,"flasso",sep=""))
sel2= c("adaflasso",paste("iter",2:10,"flasso",sep=""))
sel3= c("nknots_adaflasso",paste("nknots_iter",2:10,"flasso",sep=""))

min_bic_iter = out |> select(all_of(sel)) |> round(5) |> apply(1,which.min)
min_bic = out |> select(all_of(sel)) |> round(5) |> apply(1,min)
flasso = out |> select(all_of(sel3))
min_bic_iter_flasso = c()

for(k in 1:nrow(flasso)){
  min_bic_iter_flasso[k] = flasso[k,min_bic_iter[k]]
}

final = out |> 
  select(c("true_model","bic_true_model",
           "nknots_PELT","bic_PELT",
           "nknots_BS","bic_BS",
           "nknots_WBS","bic_WBS",
           "nknots_flasso","bic_flasso",
           "nknots_adaflasso","bic_adaflasso")) |>
  mutate(best_iter = min_bic_iter,
         nknots_best_iter = min_bic_iter_flasso,
         min_bic_flasso = min_bic) 

df2 = final %>% 
  pivot_longer(cols=c(3,5,7,9,11,14),names_to="Method",values_to="Mean Shifts") %>%
  mutate(Method = factor(Method,
                         levels = c("nknots_flasso",
                                    "nknots_adaflasso",
                                    "nknots_best_iter",
                                    "nknots_BS",
                                    "nknots_WBS",
                                    "nknots_PELT")))

knots = ggplot(df2, aes(x = Method, y = `Mean Shifts`)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("nknots_flasso"="FL",
                              "nknots_adaflasso"="Ada FL",
                              "nknots_best_iter"="IRFL",
                              "nknots_BS"="BS",
                              "nknots_WBS"="WBS",
                              "nknots_PELT"="PELT")) +
  labs(y = expression(hat(m)-m), x = NULL) +
  scale_y_continuous(breaks = c(3, seq(5, 20, by = 5))) +
  geom_hline(yintercept = 3, color = "red", linewidth = .75)

####### Produce Simulation Plot
n = 1000
x = time = 1:n
true_model = rep(c(0,2,0,2),each=250) + 0.005*x
y = rnorm(n,rep(c(0,2,0,2),each=250)) + 0.005*x

sim = ggplot(NULL, aes(x = time)) +
  geom_line(aes(y = y), alpha = 0.25) + 
  geom_line(aes(y = true_model),linewidth=.75,color="red")+
  # geom_hline(yintercept = 0, color = "red", size = .75) + 
  labs(x = "t",
       y = expression(y[t]))

### Produce Slope Plot
out = read.csv("outbind.csv")
sel = c("bic_adaflasso",paste("bic_iter",2:10,"flasso",sep=""))
sel3=c("slope_adaflasso",paste("slope_iter",2:10,"flasso",sep=""))

min_bic_iter = out |> select(all_of(sel)) |> round(5) |> apply(1,which.min)

slopes = out |> select(all_of(sel3))
min_bic_slope_iter_flasso = c()

for(k in 1:nrow(slopes)){
  min_bic_slope_iter_flasso[k] = slopes[k,min_bic_iter[k]]
}

slopedf = out |> select("slope_flasso","slope_adaflasso") %>%
  mutate(slope_flasso = slope_flasso-.005,
         slope_adaflasso = slope_adaflasso-.005) %>%
  mutate(slope_best_iter = min_bic_slope_iter_flasso-.005)

dat = slopedf %>% 
  pivot_longer(cols=1:3,names_to="method",values_to="slope") %>%
  mutate(method = factor(method,
                         levels = c("slope_flasso",
                                    "slope_adaflasso",
                                    "slope_best_iter")))

slope = ggplot(dat, aes(x = method, y = slope)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("slope_flasso"="FL",
                              "slope_adaflasso"="Ada FL",
                              "slope_best_iter"="IRFL")) +
  labs(y = expression(hat(alpha)-alpha), x = NULL) +
  geom_hline(yintercept = 0, color = "red", size = .75)

### Bring all plots together using patchwork
g = sim + slope + knots + dist + plot_layout(ncol=2)

ggsave("Mean Shift + Slope Group Plot.pdf", 
       plot = g, device = "pdf", width = 14, height = 6, dpi = 300)
