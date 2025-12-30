## For General Purpose ##
## This file is designed to merge results from Cluster
## Author: Xueheng Shi ##
## Date: 11/01/2018 ##
## Version: V0
## Updated on 10/13/24 by Michael Grantham

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lpSolve)

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
out = out.bind[,-1]
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

## Produce boxplot by group

df = final %>% 
  pivot_longer(cols=16:21,names_to="method",values_to="distance") %>%
  mutate(method = factor(method,
                         levels = c("distance_flasso",
                                    "distance_adaflasso",
                                    "distance_best_iter_flasso",
                                    "distance_bs",
                                    "distance_wbs",
                                    "distance_pelt")))

p = ggplot(df, aes(x = method, y = distance)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("distance_flasso"="Flasso",
                              "distance_adaflasso"="Adaflasso",
                              "distance_best_iter_flasso"="Best Iteration",
                              "distance_bs"="BS",
                              "distance_wbs"="WBS",
                              "distance_pelt"="PELT")) +
  labs(y = "Distance", x = NULL)

ggsave("FTWI Distance by Group.pdf", 
       plot = p, device = "pdf", width = 7, height = 3, dpi = 300)




