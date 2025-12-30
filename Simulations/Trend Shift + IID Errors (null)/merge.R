
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lpSolve)
library(patchwork)
library(stringr)

parent_path = getwd() 

## Bring in functions for random signal
# Define functions to generate signal
generate_signal = function(changepoints, slopes, n) {
  # Ensure inputs are valid
  if (length(slopes) != length(changepoints) + 1) {
    stop("The number of slopes must be one more than the number of changepoints.")
  }
  
  # Calculate intercepts to ensure continuity
  intercepts <- numeric(length(slopes))
  for (i in 2:length(slopes)) {
    intercepts[i] <- intercepts[i - 1] + slopes[i - 1] * (changepoints[i - 1] - ifelse(i == 2, 0, changepoints[i - 2]))
  }
  
  # Define the piecewise linear function
  piecewise_linear <- function(x) {
    # Determine which segment the input x belongs to
    segment <- findInterval(x, changepoints, left.open = TRUE)
    
    # Calculate the y-value for the given segment
    x_start <- ifelse(segment == 0, 0, changepoints[segment])
    y <- intercepts[segment + 1] + slopes[segment + 1] * (x - x_start)
    return(y)
  }
  
  # Generate y values for 1:n
  y_values <- sapply(1:n, piecewise_linear)
  return(y_values)
}
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
  # Initialize empty logical vectors to store TRUE/FALSE for empty strings
  distances = numeric(nrow(df))
  empty1 <- logical(nrow(df))
  empty2 <- logical(nrow(df))
  
  # Loop over each row in the dataframe
  for (i in seq_len(nrow(df))) {
    # Split the strings into character vectors for the given row
    vec1 <- unlist(strsplit(as.character(df[[col1]][i]), ",", fixed = TRUE))
    vec2 <- unlist(strsplit(as.character(df[[col2]][i]), ",", fixed = TRUE))
    
    # Convert the character vectors to numeric
    num_vec1 <- as.numeric(vec1)
    num_vec2 <- as.numeric(vec2)
    
    # Check if the values in the columns are empty strings
    empty1[i] <- df[[col1]][i] == ""
    empty2[i] <- df[[col2]][i] == ""
    
    # Check for empty vectors and handle accordingly
    if (empty1[i]+empty2[i]==2) {
      distances[i] <- 0
    } else if (empty1[i]+empty2[i]==1) {
      distances[i] <- max(length(num_vec1),length(num_vec2))
    } else {
      # Compute the distance using cpt.dist() and store in the distances vector
      distances[i] <- cpt.dist(num_vec1, num_vec2, 1000)
    }
  }
  return(distances)
}


## Get all subfolders where the results are stored to
subfolders = list.dirs(parent_path, recursive=T)[-1] 
folder.path = file.path(subfolders)

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

out[is.na(out)] = ""

sel  = c("bic_adaflasso",paste("bic_iter",2:10,"flasso",sep=""))
sel2 = c("adaflasso",paste("iter",2:10,"flasso",sep=""))
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


## Produce Simulation plot
n = N = 1000
time = 1:n
d = 50

# Plot signal
true_model = .005*time
y = rnorm(n,true_model)

p1 = ggplot(NULL, aes(x = time)) +
  geom_line(aes(y = y), alpha = 0.25) + 
  geom_line(aes(y = true_model),linewidth=.75,color="red")+
  labs(x = "t",
       y = expression(y[t]))

## Produce nknots plot
knots = out |> select(nknots_flasso, nknots_adaflasso, 
                      nknots_cpop, true_nknots) |>
  mutate(nknots_best_iter = nknots_iter_flasso) |>
  pivot_longer(cols=c(1:3,5),names_to="method",values_to="mhat") %>%
  mutate(method = factor(method,
                         levels = c("nknots_flasso",
                                    "nknots_adaflasso",
                                    "nknots_best_iter",
                                    "nknots_cpop")))

p2 = ggplot(knots, aes(x = method, y = mhat)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("nknots_flasso"="FL",
                              "nknots_adaflasso"="Ada FL",
                              "nknots_best_iter"="IRFL",
                              "nknots_cpop"="CPOP")) +
  labs(x = NULL, y = expression(hat(m)-m)) +
  geom_hline(yintercept = 0, color = "red", size = .75) +
  scale_y_continuous(breaks = seq(0, 2, by = 1)) 

## Produce combined plot
combined = p1 + p2 + plot_layout(ncol=2)


ggsave("Trend Shift (null) by Group.pdf", 
       plot = combined, device = "pdf", width = 7, height = 3*.5, dpi = 300)



