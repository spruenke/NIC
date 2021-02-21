### Main Script for Bootstrap Analysis ###

pkgs = c("Rcpp", "RcppArmadillo", "microbenchmark", "ggplot2", "gridExtra", "xtable")
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
n      = c(5, 10, 50, 100, 200) # define different sample sizes
base_grid  = expand.grid(nboots, n)
colnames(base_grid) = c("nboots", "n")
dat.Cpp  = dat.R = base_grid

mae.fun = function(x, value){
    return(mean(abs(x - value)))
}
mse.fun = function(x, value){
    return(mean((x - value)^2))
}

mae.fun.comp = function(x, values){
    return(mean(rowMeans(abs(x -  values))))
}

mse.fun.comp = function(x, values){
    return(mean(rowMeans((x -  values)^2)))
}

##### Summary statistics

# True values Normal
q.025 = qnorm(0.25)
q.05  = 0
x.m   = 0
q.075 = qnorm(0.75)
x.sd  = 1
vals = c(q.025,q.05,x.m,q.075,x.sd)

prec.r.n = paste0(rep(c("MSE_", "MAE_"), each = 10), paste0(rep(c("np_", "w_"), each = 5), c("0.25 Q", "0.5 Q", "Mean", "0.75 Q", "Sd")))
prec.r = matrix(NA, ncol = length(prec.r.n), nrow = nrow(dat.R))
colnames(prec.r) = prec.r.n
prec.c = prec.r
time.r = time.c = matrix(NA, nrow = nrow(dat.R), ncol = 6)

for(i in 1:nrow(dat.R)){
    A.temp.r.np = replicate(1000, summary.boot(rnorm(n = dat.R$n[i], mean = x.m, sd = x.sd), nboot = dat.R$nboots[i]))[-c(1,6),]
    A.temp.r.w = replicate(1000, summary.boot(rnorm(n = dat.R$n[i], mean = x.m, sd = x.sd), nboot = dat.R$nboots[i], boot.type = "wild"))[-c(1,6),]
    A.temp.c.np = replicate(1000, summary.boot.c(rnorm(n = dat.R$n[i], mean = x.m, sd = x.sd), nboot = dat.R$nboots[i]))[-c(1,6),]
    A.temp.c.w = replicate(1000, summary.boot.c(rnorm(n = dat.R$n[i], mean = x.m, sd = x.sd), nboot = dat.R$nboots[i], boot.type = "wild"))[-c(1,6),]
    for(j in 1:length(vals)){
        prec.r[i,j] = mse.fun(A.temp.r.np[j,], vals[j])
        prec.r[i,j + 5] = mse.fun(A.temp.r.w[j,], vals[j])
        prec.r[i,j + 10] = mae.fun(A.temp.r.np[j,], vals[j])
        prec.r[i,j + 15] = mae.fun(A.temp.r.w[j,], vals[j])
        
        prec.c[i,j] = mse.fun(A.temp.c.np[j,], vals[j])
        prec.c[i,j + 5] = mse.fun(A.temp.c.w[j,], vals[j])
        prec.c[i,j + 10] = mae.fun(A.temp.c.np[j,], vals[j])
        prec.c[i,j + 15] = mae.fun(A.temp.c.w[j,], vals[j])
    }
    
    x.sample = rnorm(n = dat.R$n[i], mean = x.m, sd = x.sd)
    time.r[i,] = summary(microbenchmark(summary.boot(x.sample, nboot = dat.R$nboots[i]))$time / 1000) # microseconds
    time.c[i,] = summary(microbenchmark(summary.boot.c(x.sample, nboot = dat.R$nboots[i]))$time / 1000)
}

time.r = cbind(base_grid, time.r)
time.c = cbind(base_grid, time.c)

dat.R = cbind(dat.R, prec.r)
dat.Cpp = cbind(dat.Cpp, prec.c)

l.p = list()
for(i in seq(1, nrow(dat.R), 5)){
        
        l.p[[i]] = ggplot()+
            geom_line(mapping = aes(x = dat.R[i + c(0:4), "nboots"], y = dat.R[i + c(0:4), "MSE_np_Mean"]), color = "red") + 
            geom_line(mapping = aes(x = dat.R[i + c(0:4), "nboots"], y = dat.R[i + c(0:4), "MSE_w_Mean"]), color = "blue") +
            labs(x = "Bootstrap Iterations", y = "MSE")
        
        l.p[[i + 1]] = ggplot()+
            geom_line(mapping = aes(x = dat.R[i + c(0:4), "nboots"], y = dat.R[i + c(0:4), "MAE_np_Mean"]), color = "red") + 
            geom_line(mapping = aes(x = dat.R[i + c(0:4), "nboots"], y = dat.R[i + c(0:4), "MAE_w_Mean"]), color = "blue") +
            labs(x = "Bootstrap Iterations", y = "MAE")
        
        l.p[[i + 2]] = ggplot()+
            geom_line(mapping = aes(x = dat.Cpp[i + c(0:4), "nboots"], y = dat.Cpp[i + c(0:4), "MSE_np_Mean"]), color = "red") + 
            geom_line(mapping = aes(x = dat.Cpp[i + c(0:4), "nboots"], y = dat.Cpp[i + c(0:4), "MSE_w_Mean"]), color = "blue") +
            labs(x = "Bootstrap Iterations", y = "MSE")
        
        l.p[[i + 3]] = ggplot()+
            geom_line(mapping = aes(x = dat.Cpp[i + c(0:4), "nboots"], y = dat.Cpp[i + c(0:4), "MAE_np_Mean"]), color = "red") + 
            geom_line(mapping = aes(x = dat.Cpp[i + c(0:4), "nboots"], y = dat.Cpp[i + c(0:4), "MAE_w_Mean"]), color = "blue") +
            labs(x = "Bootstrap Iterations", y = "MAE")
        
        l.p[[i + 4]] = ggplot()+
            geom_line(mapping = aes(x = dat.R[i + c(0:4), "nboots"], y = time.r[i + c(0:4), 3]), color = "red")+
            geom_line(mapping = aes(x = dat.Cpp[i + c(0:4), "nboots"], y = time.c[i + c(0:4), 3]), color = "blue") +
            labs(x = "Bootstrap Iterations", y = "Microseconds")
        
}

summary_plot = do.call("grid.arrange", c(l.p, nrow = 5, ncol = 5))
ggsave("./resuls/plot_summary.pdf", plot = summary_plot, dpi = 300, width = 16, height = 9)

save(file = "./results/summary_data.RData", summary_plot, dat.R, dat.Cpp, time.r, time.c)
summary_xtable_R = xtable(t(dat.R))
summary_xtable_Cpp = xtable(t(dat.Cpp))

summary_time_xtab_R = xtable(time.r)
summary_time_xtab_C = xtable(time.c)
print(summary_xtable_R, file = "./results/tab_summary_r.tex")
print(summary_xtable_Cpp, file = "./results/tab_summary_cpp.tex")
print(summary_time_xtab_R, file = "./results/tab_summaryTime_r.tex")
print(summary_time_xtab_C, file = "./results/tab_summaryTime_Cpp.tex")

######### Linear Regression

    #### Specify a true data generating model (to assess accuracy of estimates)
    set.seed(314)
    beta = runif(4, -20, 20)
    prec = matrix(NA, ncol = 4, nrow = nrow(base_grid))
    colnames(prec) = c("MSE_np", "MSE_w", "MAE_np", "MAE_w")
    dat.Cpp  = dat.R = cbind(base_grid, prec)
    
    time.r = time.c = matrix(NA, nrow = nrow(dat.R), ncol = 6)
    
    for(i in 4:nrow(dat.R)){
        # A.temp.r.np = A.temp.r.w = A.temp.c.np = A.temp.c.w = matrix(NA, nrow = length(beta), ncol = 1000)
        # for(j in 1:1000){
        #         # structure of X is: one variable coming (discrete) uniformly from 30 to 100, one (continuously) from 50 to 230
        #         # one is log(normal) and one is normally distributed around 100 with sd 7.5
        #         x = cbind(sample(c(30:100), size = dat.R$n[i], replace = T), runif(dat.R$n[i], 50, 230), log(rnorm(dat.R$n[i], 100, 7.5)), rnorm(dat.R$n[i], 100, 7.5))
        #         y = x%*%beta + rnorm(dat.R$n[i])
        #         A.temp.r.np[,j] = tryCatch(reg.boot(y, x, nboot = dat.R$nboots[i], intercept = F)$estimate, error=function(cond){return(NA)}) # maybe trycatch
        #         A.temp.r.w[,j] = reg.boot(y, x, nboot = dat.R$nboots[i], boot.type = "wild", intercept = F)$estimate
        #         A.temp.c.np[,j] = tryCatch(reg.boot.c(y, x, nboot = dat.R$nboots[i], intercept = F)$estimate, error=function(cond){return(NA)})
        #         A.temp.c.w[,j] = reg.boot.c(y, x, nboot = dat.R$nboots[i], boot.type = "wild", intercept = F)$estimate
        # }
        # dat.R$MSE_np = mse.fun.comp(A.temp.r.np, beta)
        # dat.R$MSE_w  = mse.fun.comp(A.temp.r.w, beta)
        # dat.R$MAE_np = mae.fun.comp(A.temp.r.np, beta)
        # dat.R$MAE_w  = mae.fun.comp(A.temp.r.w, beta)
        # 
        # dat.Cpp$MSE_np = mse.fun.comp(A.temp.c.np, beta)
        # dat.Cpp$MSE_w  = mse.fun.comp(A.temp.c.w, beta)
        # dat.Cpp$MAE_np = mae.fun.comp(A.temp.c.np, beta)
        # dat.Cpp$MAE_w  = mae.fun.comp(A.temp.c.w, beta)
        
        x = cbind(sample(c(30:100), size = dat.R$n[i], replace = T), runif(dat.R$n[i], 50, 230), log(rnorm(dat.R$n[i], 100, 7.5)), rnorm(dat.R$n[i], 100, 7.5))
        y = x%*%beta + rnorm(dat.R$n[i])
        
        time.r[i,] = summary(microbenchmark(reg.boot(y, x, nboot = dat.R$nboots[i], boot.type = "wild"))$time / 1000) # microseconds
        time.c[i,] = summary(microbenchmark(reg.boot.c(y, x, nboot = dat.R$nboots[i], boot.type = "wild"))$time / 1000)
        print(i)
    }

    time.r = cbind(base_grid, time.r)
    time.c = cbind(base_grid, time.c)
    
    #save(file = "./results/reg_data.RData", dat.R, dat.Cpp, time.r, time.c)
    
    l.p = list()
    #for(i in seq(1, nrow(dat.R), 5)){
    for(i in 1:5){
        
        # l.p[[i]] = ggplot()+
        #     geom_line(mapping = aes(x = dat.R[i + c(0:4), "nboots"], y = dat.R[i + c(0:4), "MSE_np"]), color = "red") + 
        #     geom_line(mapping = aes(x = dat.R[i + c(0:4), "nboots"], y = dat.R[i + c(0:4), "MSE_w"]), color = "blue") +
        #     labs(x = "Bootstrap Iterations", y = "MSE")
        # 
        # l.p[[i + 1]] = ggplot()+
        #     geom_line(mapping = aes(x = dat.R[i + c(0:4), "nboots"], y = dat.R[i + c(0:4), "MAE_np"]), color = "red") + 
        #     geom_line(mapping = aes(x = dat.R[i + c(0:4), "nboots"], y = dat.R[i + c(0:4), "MAE_w"]), color = "blue") +
        #     labs(x = "Bootstrap Iterations", y = "MAE")
        # 
        # l.p[[i + 2]] = ggplot()+
        #     geom_line(mapping = aes(x = dat.Cpp[i + c(0:4), "nboots"], y = dat.Cpp[i + c(0:4), "MSE_np"]), color = "red") + 
        #     geom_line(mapping = aes(x = dat.Cpp[i + c(0:4), "nboots"], y = dat.Cpp[i + c(0:4), "MSE_w"]), color = "blue") +
        #     labs(x = "Bootstrap Iterations", y = "MSE")
        # 
        # l.p[[i + 3]] = ggplot()+
        #     geom_line(mapping = aes(x = dat.Cpp[i + c(0:4), "nboots"], y = dat.Cpp[i + c(0:4), "MAE_np"]), color = "red") + 
        #     geom_line(mapping = aes(x = dat.Cpp[i + c(0:4), "nboots"], y = dat.Cpp[i + c(0:4), "MAE_w"]), color = "blue") +
        #     labs(x = "Bootstrap Iterations", y = "MAE")
        
        l.p[[i ]] = ggplot()+
            geom_line(mapping = aes(x = dat.R[i + c(0:4), "nboots"], y = log(time.r[i + c(0:4), 3])), color = "red")+
            geom_line(mapping = aes(x = dat.Cpp[i + c(0:4), "nboots"], y = log(time.c[i + c(0:4), 3])), color = "blue") +
            labs(x = "Bootstrap Iterations", y = "Log-Microseconds")
        
    }
    
    reg_plot = do.call("grid.arrange", c(l.p, nrow = 5))
    ggsave("./results/plot_regression.pdf", plot = reg_plot, dpi = 300, width = 16, height = 9)
    