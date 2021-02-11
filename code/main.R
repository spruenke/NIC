### Main Script for Bootstrap Analysis ###

pkgs = c("Rcpp", "RcppArmadillo", "microbenchmark", "ggplot2")
if(pkgs %in% installed.packages() == F) install.packages(pkgs = pkgs, dependencies = T)
lapply(as.list(pkgs), library, quietly = T, character.only = T)

############################################
############# Source Functions #############
############################################

# These load the respective scripts containing the functions. Note,
# that these scrips already load the cpp-functions into the environment.
# However, only use them carefully. There are R-wrappers making use of them
# to embed them properly.

source("summary.R") # Bootstrap Summary aka Parameter Estimation
source("t_test.R") # Bootstrap t-Test
source("regression.R") # Bootstrap Linear Regression

nboots = c(10, 100, 500, 1000, 10000) # define different bootstrap iterations
n      = c(5, 10, 15, 20, 50, 100, 200) # define different sample sizes
dat.R  = expand.grid(nboots, n)
colnames(dat.R) = c("nboots", "n")
dat.R  = data.frame(dat.R, MSE = NA, MAE = NA)
dat.Cpp  = dat.R

mae.fun = function(x, value){
    return(mean(abs(x - value)))
}
mse.fun = function(x, value){
    return(mean((x - value)^2))
}

##### Summary statistics

# True values (arbitrary)
q.025 = qnorm(0.25)
q.05  = 0
x.m   = 0
q.075 = qnorm(0.75)
x.sd  = 1

for(i in 1:nrow(dat.R)){
    x.sample = rnorm(n = dat.R$n, mean = x.m, sd = x.sd)
    
}

######## dependent on language

######## dependent on nboot and n
