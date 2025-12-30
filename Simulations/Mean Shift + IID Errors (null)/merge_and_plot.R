
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lpSolve)
library(patchwork)


parent_path = getwd() 

## Bring in distance metric
cpt.dist = function(C1, C2, N){
  m = length(C1); k = length(C2)
  
  # --- ADD: handle empty cases up front ---
  if (m == 0L && k == 0L) return(0)       # both empty → distance 0
  if (m == 0L || k == 0L) return(abs(m-k))# one empty  → penalty only
  
  ##Generate Cost Matrix via all paired distance
  pair = expand.grid(C1, C2)
  if(m==k){ 
    cost.mat = matrix(abs(pair[,1]-pair[,2]), nrow=m, ncol=k, byrow=TRUE)
  } else if(m > k){
    cost.mat = cbind(matrix(abs(pair[,1]-pair[,2]), nrow=m, ncol=k, byrow=TRUE),
                     matrix(0, nrow=m, ncol=(m-k)))
  } else {
    cost.mat = rbind(matrix(abs(pair[,1]-pair[,2]), nrow=m, ncol=k, byrow=FALSE),
                     matrix(0, nrow=(k-m), ncol=k))
  }
  cpt.asgn = lp.assign(cost.mat, direction = "min")
  cpt.asgn$objval / N + abs(m - k)
}

## Bring in function to apply distance
apply_distance = function(df, col1, col2) {
  distances <- numeric(nrow(df))
  for (i in seq_len(nrow(df))) {
    # split, then DROP empty strings
    v1 <- unlist(strsplit(as.character(df[[col1]][i]), ",", fixed = TRUE))
    v2 <- unlist(strsplit(as.character(df[[col2]][i]), ",", fixed = TRUE))
    v1 <- v1[nzchar(v1)]
    v2 <- v2[nzchar(v2)]
    
    # convert to numeric, drop any stray NAs
    num_vec1 <- suppressWarnings(as.numeric(v1)); num_vec1 <- num_vec1[!is.na(num_vec1)]
    num_vec2 <- suppressWarnings(as.numeric(v2)); num_vec2 <- num_vec2[!is.na(num_vec2)]
    
    distances[i] <- cpt.dist(num_vec1, num_vec2, 1000)
  }
  distances
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
out = out.bind[1:1000,]
write.csv(out,"outbind.csv",row.names=FALSE)

sel = c("bic_adaflasso",paste("bic_iter",2:10,"flasso",sep=""))
sel2= c("adaflasso",paste("iter",2:10,"flasso",sep=""))

min_bic_iter = out |> select(all_of(sel)) |> round(5) |> apply(1,which.min)
min_bic = out |> select(all_of(sel)) |> round(5) |> apply(1,min)
flasso = out |> select(all_of(sel2))
min_bic_iter_flasso = c()

for(k in 1:nrow(flasso)){
  min_bic_iter_flasso[k] = flasso[k,min_bic_iter[k]]
}

semifinal = out |> 
  select(c("true_model","bic_true_model",
           "PELT","bic_PELT",
           "knots_BS","bic_BS",
           "knots_WBS","bic_WBS",
           "flasso","bic_flasso",
           "adaflasso","bic_adaflasso")) |>
  mutate(best_iter = min_bic_iter,
         changepoints = min_bic_iter_flasso,
         best_iter = min_bic_iter,
         changepoints = min_bic_iter_flasso,
         min_bic_flasso = min_bic,
         distance_pelt = apply_distance(out,"true_model","PELT"),
         distance_wbs = apply_distance(out,"true_model","knots_WBS"),
         distance_bs = apply_distance(out,"true_model","knots_BS"),
         distance_flasso = apply_distance(out,"true_model","flasso"),
         distance_adaflasso = apply_distance(out,"true_model","adaflasso")
  ) 

final = semifinal |>
  mutate(distance_best_iter_flasso = apply_distance(semifinal,"true_model","changepoints"))

write.csv(final,"final_out.csv",row.names=FALSE)

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
  scale_y_continuous(breaks = seq(1, 2, by = 1)) +
  geom_hline(yintercept = 0, color = "red", linewidth = .75)

####### Produce Simulation Plot
n = 1000
x = time = 1:n
true_model = rep(0,n)
y = rnorm(n,0) 

sim = ggplot(NULL, aes(x = time)) +
  geom_line(aes(y = y), alpha = 0.25) + 
  geom_line(aes(y = true_model),linewidth=.75,color="red")+
  # geom_hline(yintercept = 0, color = "red", size = .75) + 
  labs(x = "t",
       y = expression(y[t]))

### Bring all plots together using patchwork
g = sim +  knots + plot_layout(ncol=2)

ggsave("Mean Shift (Null) Group Plot.pdf", 
       plot = g, device = "pdf", width = 14*.5, height = 6*.25, dpi = 300)
