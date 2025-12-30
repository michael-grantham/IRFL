library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lpSolve)
library(patchwork)

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
    distances[i] <- cpt.dist(num_vec1, num_vec2,1000)
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
out = out = out.bind[1:1000,]
write.csv(out,"outbind.csv",row.names=FALSE)

sel = c("bic_adaflasso",paste("bic_iter",2:10,"flasso",sep=""))
sel2= c("adaflasso",paste("iter",2:10,"flasso",sep=""))
sel3 = c("nknots_adaflasso",paste("nknots_iter",2:10,"flasso",sep=""))

min_bic_iter = out |> select(all_of(sel)) |> round(5) |> apply(1,which.min)
min_bic = out |> select(all_of(sel)) |> round(5) |> apply(1,min)
flasso = out |> select(all_of(sel2))
nknots_flasso = out |> select(all_of(sel3))
min_bic_iter_flasso = c()
nknots_iter_flasso = c()

for(k in 1:nrow(flasso)){
  min_bic_iter_flasso[k] = flasso[k,min_bic_iter[k]]
  nknots_iter_flasso[k] = nknots_flasso[k,min_bic_iter[k]]
}

semifinal = out |> 
  select(c("true_model","bic_true_model",
           # "PELT","bic_PELT",
           # "knots_BS","bic_BS",
           # "knots_WBS","bic_WBS",
           
           "cpop","bic_cpop",
           
           "flasso","bic_flasso",
           "adaflasso","bic_adaflasso")) |>
  mutate(best_iter = min_bic_iter,
         changepoints = min_bic_iter_flasso,
         best_iter = min_bic_iter,
         changepoints = min_bic_iter_flasso,
         min_bic_flasso = min_bic,
         distance_cpop = apply_distance(out,"true_model","cpop"),
         # distance_pelt = apply_distance(out,"true_model","PELT"),
         # distance_wbs = apply_distance(out,"true_model","knots_WBS"),
         # distance_bs = apply_distance(out,"true_model","knots_BS"),
         distance_flasso = apply_distance(out,"true_model","flasso"),
         distance_adaflasso = apply_distance(out,"true_model","adaflasso")
  ) 

final = semifinal |>
  mutate(distance_best_iter_flasso = apply_distance(semifinal,"true_model","changepoints"))

write.csv(final,"final_out.csv",row.names=FALSE)

## Produce Simulation plot
n = 1000
time = 1:n
true_model = c(1:250,251:2,1:250,251:2)/100
y = rnorm(n,true_model)

p1 = ggplot(NULL, aes(x = time)) +
  geom_line(aes(y = y), alpha = 0.25) + 
  geom_line(aes(y = true_model),linewidth=.75,color="red")+
  # geom_hline(yintercept = 0, color = "red", size = .75) + 
  labs(x = "t",
       y = expression(y[t]))

## Produce boxplot by group

df = final %>% 
  pivot_longer(cols=12:15,names_to="method",values_to="distance") %>%
  mutate(method = factor(method,
                         levels = c("distance_flasso",
                                    "distance_adaflasso",
                                    "distance_best_iter_flasso",
                                    # "distance_bs",
                                    # "distance_wbs",
                                    # "distance_pelt"
                                    "distance_cpop")))

p2 = ggplot(df, aes(x = method, y = distance)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("distance_flasso"="FL",
                              "distance_adaflasso"="Ada FL",
                              "distance_best_iter_flasso"="IRFL",
                              # "distance_bs" = "BS",
                              # "distance_wbs"="WBS",
                              # "distance_pelt"="PELT",
                              "distance_cpop"="CPOP")) +
  labs(x = NULL, y = "Distance") +
  geom_hline(yintercept = 0, color = "red", size = .75) 

## Produce nknots plot
knots = out |> select(nknots_flasso, nknots_adaflasso, nknots_cpop) |>
  mutate(nknots_best_iter = nknots_iter_flasso) |>
  pivot_longer(cols=1:4,names_to="method",values_to="mhat") %>%
  mutate(method = factor(method,
                         levels = c("nknots_flasso",
                                    "nknots_adaflasso",
                                    "nknots_best_iter",
                                    "nknots_cpop")))

p3 = ggplot(knots, aes(x = method, y = mhat-3)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("nknots_flasso"="FL",
                              "nknots_adaflasso"="Ada FL",
                              "nknots_best_iter"="IRFL",
                              # "distance_bs" = "BS",
                              # "distance_wbs"="WBS",
                              # "distance_pelt"="PELT",
                              "nknots_cpop"="CPOP")) +
  labs(x = NULL, y = expression(hat(m)-m)) +
  geom_hline(yintercept = 0, color = "red", size = .75) 

## Produce combined plot
combined = p1 + p2 + p3 + plot_layout(ncol=2)


ggsave("Trend Shift by Group.pdf", 
       plot = combined, device = "pdf", width = 7, height = 3, dpi = 300)



