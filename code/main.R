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

################# Summary Function
##### A Matter of Speed, Boxplots for Accuracy (e.g. via replicate)

################# t-Test Function
##### Speed, Type-1-Error, some Type-2-Errors

################# Linear Regression Function
##### Definitely a matter of Speed, Accuracies of Beta, bootstrap std vs sample std of betaHat

