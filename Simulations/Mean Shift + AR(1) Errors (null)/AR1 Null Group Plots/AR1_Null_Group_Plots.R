
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

out.bind = NULL

file_names = c("final_out_null_neg0.7.csv",
               "final_out_null_neg0.3.csv",
               "final_out_null_0.0.csv",
               "final_out_null_0.3.csv",
               "final_out_null_0.7.csv")

for (file in file_names){
    out.ind = read.csv(file)
    out.bind = rbind(out.bind, out.ind)
}

out = out.bind

write.csv(out, "data_group_plot.csv", row.names=FALSE)

data = read.csv("data_group_plot.csv")

names(data)


# Plot of BIC
bic_plot = data |> 
  mutate(`Ada FL` = bic_adaflasso,
         `IRFL` = min_bic_best_iter,
         WCM = bic_WCM,
         AR1Seg = bic_AR1seg) |> 
  pivot_longer(cols = 28:31,names_to = "Method",
               values_to = "bic") |>
  mutate(Method = factor(Method, levels = c(
    "Ada FL","IRFL",
    "WCM","AR1Seg"
  ))) |>
  ggplot(aes(x = Method, 
             y = bic-bic_true_model, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) + 
  scale_fill_manual(
    values = c("Ada FL" = "lightblue", 
               "IRFL" = "lightcoral",
               WCM = "lightgreen",
               AR1Seg = "purple")) +
  geom_hline(yintercept = 0, color = "red", size = .25) +
  labs(x = NULL, 
       y = "BIC(Estimated) - BIC(True Model)", 
       fill = "Method")+
  theme(
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
    legend.position="none"
  ) +
  facet_wrap(~ true_phi, ncol = 3)


# Plot of slope
phi_plot = data |> 
  mutate(`Ada FL` = phi_adaflasso,
         #`Ada Flasso (YW)` = phi_yw_best_iter,
         `IRFL` = phi_best_iter,
         # `IRFL (YW)` = phi_yw_best_iter,
         WCM = phi_WCM,
         AR1Seg = phi_AR1seg) |> 
  pivot_longer(cols = 28:31,names_to = "Method",
               values_to = "phi_hat") |>
  mutate(Method = factor(Method, levels = c(
    "Ada FL","IRFL",
    "WCM","AR1Seg"
  ))) |>
  ggplot(aes(x = Method, 
             y = phi_hat-true_phi, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) + 
  scale_fill_manual(
    values = c("Ada FL" = "lightblue", 
               #"Ada Flasso (YW)" = "blue",
               "IRFL" = "lightcoral",
               #"IRFL (YW)" = "coral",
               WCM = "lightgreen",
               AR1Seg = "purple")#,
    # labels = c("slope_best_iter" = "IRFL", 
    #            "slope_OLS_adaflasso" = "Adaptive Fused Lasso") 
  ) +
  geom_hline(yintercept = 0, color = "red", size = .25) +
  labs(x = NULL, 
       y = expression(hat(phi) - phi), 
       fill = "Method")+
  theme(
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
    legend.position="none"
  ) +
  facet_wrap(~ true_phi, ncol = 3)


# Plot of nknots
nknots_plot = data |> 
  mutate(`Ada FL` = nknots_adaflasso,
         `IRFL` = nknots_best_iter,
         WCM = nknots_WCM,
         AR1Seg = nknots_AR1seg) |>
  pivot_longer(cols = 28:31,names_to = "Method",
               values_to = "nknots") |>
  mutate(Method = factor(Method, levels = c(
    "Ada FL","IRFL",
    "WCM","AR1Seg"
  ))) |>
  ggplot(aes(x = Method, 
             y = nknots, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) + 
  scale_fill_manual(
    values = c("Ada FL" = "lightblue", 
               "IRFL" = "lightcoral",
               WCM = "lightgreen",
               AR1Seg = "purple")) +
  geom_hline(yintercept = 0, color = "red", size = .25) +
  labs(x = NULL, 
       y = expression(hat(m) - m), 
       fill = "Method")+
  theme(
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
    legend.position="none"
  ) +
  facet_wrap(~ true_phi, ncol = 3)

ggsave("AR1 Group BIC Plot null.pdf", 
       plot = bic_plot, device = "pdf", width = 14*.5, height = 6*.5, dpi = 300)
ggsave("AR1 Group Phi Plot null.pdf", 
       plot = phi_plot, device = "pdf", width = 14*.5, height = 6*.5, dpi = 300)
ggsave("AR1 Group nknot Plot null.pdf", 
       plot = nknots_plot, device = "pdf", width = 14*.5, height = 6*.5, dpi = 300)
