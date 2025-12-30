

simulate = function(){
  alpha = 0.005
  n = 1000
  x = 1:n
  tau = c(0,251,501,751,1000)
y = alpha*x + rnorm(n,rep(c(0,2,0,2),each=250))

# Create design matrix X
dm = matrix(NA,nrow=1000,ncol=length(tau)-1)
for(j in 2:length(tau)){
  dm[,j-1] = c(rep(0,tau[j-1]),
               rep(1,tau[j]-tau[j-1]),
               rep(0,1000-tau[j]))
}
X_OLS = cbind(x, dm)

# Initial estimation of beta
beta_hat_OLS = solve(t(X_OLS)%*%X_OLS)%*%t(X_OLS)%*%y
slope_OLS = beta_hat_OLS[1]

# plot(y,type="l",
#      xlab = "t",ylab=expression(y[t]),
#      main = paste0("Simulation: ",i,"; Iteration: ",k,collapse=""))
# legend("topright",
#        legend = c("OLS", "FL"),
#        col=c("red","blue"),lty=c(1,1))

# Randomly increase tau by one changepoint and re-estimate
remaining_numbers = setdiff(0:n, tau)

random_sample = sample(remaining_numbers, size = 1)

# Re-assign tau and re-estimate slope
tau2 = sort(c(random_sample,tau))

# Create design matrix X
dm = matrix(NA,nrow=1000,ncol=length(tau2)-1)
for(j in 2:length(tau2)){
  dm[,j-1] = c(rep(0,tau2[j-1]),
               rep(1,tau2[j]-tau2[j-1]),
               rep(0,1000-tau2[j]))
}
X_OLS = cbind(x, dm)

# Initial estimation of beta
beta_hat_OLS_extra = solve(t(X_OLS)%*%X_OLS)%*%t(X_OLS)%*%y
slope_OLS_extra = beta_hat_OLS_extra[1]

# Remove a changepoint and do it again
ind = sample(1:length(tau), size = 1)

# Re-assign tau and re-estimate slope
tau3 = tau[-ind]

# Create design matrix X
dm = matrix(NA,nrow=1000,ncol=length(tau3)-1)
for(j in 2:length(tau3)){
  dm[,j-1] = c(rep(0,tau3[j-1]),
               rep(1,tau3[j]-tau3[j-1]),
               rep(0,1000-tau3[j]))
}
X_OLS = cbind(x, dm)

# Initial estimation of beta
beta_hat_OLS_less = solve(t(X_OLS)%*%X_OLS)%*%t(X_OLS)%*%y
slope_OLS_less = beta_hat_OLS_less[1]

return(c(slope_OLS,slope_OLS_extra,slope_OLS_less))
}

mat = replicate(1000,simulate())

df = as.data.frame(t(mat))
colnames(df) = c(
  "Truth",
  "Extra Changepoint",
  "Missing Changepoint"
)

library(tidyr)
library(ggplot2)
library(dplyr)

p = df |>
  pivot_longer(cols = 1:3,
               names_to = "Method",
               values_to = "Slope") |>
  mutate(Method = factor(Method,
                         levels = c(
                           "Truth",
                           "Extra Changepoint",
                           "Missing Changepoint"
                         ))) |>
  ggplot(aes(x=Method,y=Slope-.005)) + 
  geom_boxplot() + 
  labs(x = NULL,
       y = expression(hat(alpha)-alpha),
       title = expression("OLS Estimation of "*alpha*" with True and Faulty " * tau))+
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = 0, color = "red", size = .55) +
  scale_x_discrete(
    labels = c(
      "Truth" = expression("True " * tau),
      "Extra Changepoint" = expression(tau ~ "+ Random Changepoint"),
      "Missing Changepoint" = expression(tau ~ "- Random Changepoint")
    )
  ) 

ggsave("OLS Noised FTWI Estimates.pdf", 
       plot = p, device = "pdf", width = 7, height = 3, dpi = 300)
