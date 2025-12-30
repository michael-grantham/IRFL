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
out = out[1:1000,]
write.csv(out,"outbind.csv",row.names=FALSE)

sel = c("bic_adaflasso",paste("bic_iter",2:10,"flasso",sep=""))
sel2= c("adaflasso",paste("iter",2:10,"flasso",sep=""))
sel3= c("slope_OLS_adaflasso",paste("slope_OLS_iter",2:10,"flasso",sep=""))
sel4= c("nknots_adaflasso",paste("nknots_iter",2:10,"flasso",sep=""))

min_bic_iter = out |> select(all_of(sel)) |> round(5) |> apply(1,which.min)
min_bic = out |> select(all_of(sel)) |> round(5) |> apply(1,min)
flasso = out |> select(all_of(sel2))
slope_flasso = out |> select(all_of(sel3))
nknots_flasso = out |> select(all_of(sel4))
min_bic_iter_flasso = c()
best_iter_slope_OLS = c()
best_iter_nknots = c()

for(k in 1:nrow(flasso)){
  min_bic_iter_flasso[k] = flasso[k,min_bic_iter[k]]
  best_iter_slope_OLS[k] = slope_flasso[k,min_bic_iter[k]]
  best_iter_nknots[k] = nknots_flasso[k,min_bic_iter[k]]
}

semifinal = out |> 
  select(c("true_model","true_slope","true_nknots","bic_true_model",
           "adaflasso","slope_OLS_adaflasso","nknots_adaflasso","bic_adaflasso")) |>
  mutate(best_iter = min_bic_iter,
         changepoints = min_bic_iter_flasso,
         best_iter = min_bic_iter,
         nknots_best_iter = best_iter_nknots,
         slope_best_iter = best_iter_slope_OLS,
         min_bic_best_iter = min_bic,
         distance_adaflasso = apply_distance(out,"true_model","adaflasso"),
         slope_adaflasso = slope_OLS_adaflasso,
         slope_best_iter = best_iter_slope_OLS
  ) 

final = semifinal |>
  mutate(distance_best_iter_flasso = 
           apply_distance(semifinal,"true_model","changepoints"))

write.csv(final,"final_out.csv",row.names=FALSE)





