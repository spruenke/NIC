### Main Script for Bootstrap Analysis ###

pkgs = c("Rcpp", "RcppArmadillo", "microbenchmark", "ggplot2")
if(pkgs %in% installed.packages() == F) install.packages(pkgs = pkgs, dependencies = T)
lapply(pkgs, library, quietly = T)

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
dat.R    = expand.grid(nboots, n)
colnames(dat) = c("nboots", "n")
dat.Cpp  = dat.R

##### Precision of NP and Wild dependent on n, nboot and language

# True values
q.025 = qnorm(0.25)
q.05  = 0
x.m   = 0
q.075 = qnorm(0.75)
x.sd  = 1

######## dependent on language

######## dependent on nboot and n
######## The following double loop is rather inefficient but convenient to read and write
for(i in 1:length(nboots)){
    for(j in 1:length(n)){
        
    }
}