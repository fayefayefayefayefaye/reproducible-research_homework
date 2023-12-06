#install.packages("ggplot2")
#install.packages("gridExtra")

library(ggplot2)
library(gridExtra)

# This is code to make the random_walk function
random_walk  <- function (n_steps) {
  
  df <- data.frame(x = rep(NA, n_steps), y = rep(NA, n_steps), time = 1:n_steps)
  
  df[1,] <- c(0,0,1)
  
  for (i in 2:n_steps) {
    
    h <- 0.25
    
    angle <- runif(1, min = 0, max = 2*pi)
    
    df[i,1] <- df[i-1,1] + cos(angle)*h
    
    df[i,2] <- df[i-1,2] + sin(angle)*h
    
    df[i,3] <- i
    
  }
  
  return(df)
  
}

# Applying the random_walk function (for 500 steps)
set.seed(888) # use the set.seed() function to set the random seed
data1 <- random_walk(500)

# Plot data1 in 2D
plot1 <- ggplot(aes(x = x, y = y), data = data1) +
  
  geom_path(aes(colour = time)) +
  
  theme_bw() +
  
  xlab("x-coordinate") +
  
  ylab("y-coordinate")

# Applying the random_walk function again
set.seed(888) # this has to be done twice because the .random.seed changes after
  # random_walk function is run
data2 <- random_walk(500) 

# Plot data2 in 2D
plot2 <- ggplot(aes(x = x, y = y), data = data2) +
  
  geom_path(aes(colour = time)) +
  
  theme_bw() +
  
  xlab("x-coordinate") +
  
  ylab("y-coordinate")

# Combine the plots for data1 and data2 on the same figure
grid.arrange(plot1, plot2, ncol=2)
