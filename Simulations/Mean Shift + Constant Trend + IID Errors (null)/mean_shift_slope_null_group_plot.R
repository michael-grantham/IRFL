library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(scales)


###### Produce Distance Plots
final = read.csv("final_out.csv")

df = final %>% 
  pivot_longer(cols=c(3,5,7,9,11,14),names_to="Method",values_to="Mean Shifts") %>%
  mutate(Method = factor(Method,
                         levels = c("nknots_flasso",
                                    "nknots_adaflasso",
                                    "nknots_best_iter",
                                    "nknots_BS",
                                    "nknots_WBS",
                                    "nknots_PELT")))

dist = ggplot(df, aes(x = Method, y = `Mean Shifts`)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("nknots_flasso"="FL",
                              "nknots_adaflasso"="Ada FL",
                              "nknots_best_iter"="IRFL",
                              "nknots_BS"="BS",
                              "nknots_WBS"="WBS",
                              "nknots_PELT"="PELT")) +
  labs(y = expression(hat(m)-m), x = NULL) +
  geom_hline(yintercept = 0, color = "red", size = .75)

####### Produce Simulation Plot
n = 1000
x = time = 1:n
true_model = 0.005*x
y = rnorm(n) + 0.005*x

sim = ggplot(NULL, aes(x = time)) +
  geom_line(aes(y = y), alpha = 0.25) + 
  geom_line(aes(y = true_model),linewidth=.75,color="red")+
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
  scale_y_continuous(labels = label_number(accuracy = 0.0001)) +
  geom_hline(yintercept = 0, color = "red", size = .75) 

### Bring all plots together using patchwork
g = sim + dist + slope + plot_layout(ncol=2)

ggsave("Mean Shift + Slope (Null) Group Plot.pdf", 
       plot = g, device = "pdf", width = 14*.5, height = 6*.5, dpi = 300)
