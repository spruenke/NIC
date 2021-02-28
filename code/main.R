### Main Script for Bootstrap Analysis ###
start.time = Sys.time()
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

t.power.fun = function(x, alpha = 0.05){
        return(x < 0.05)
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

for(i in 1:nrow(base_grid)){
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
colnames(time.r) = colnames(time.c) = c(colnames(base_grid), names(summary(rnorm(10))))

dat.R = cbind(dat.R, prec.r)
dat.Cpp = cbind(dat.Cpp, prec.c)

l.p = list()
for(i in seq(1, nrow(base_grid), 5)){
        
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
            labs(x = "Bootstrap Iterations", y = "Microseconds") # on Median
        
}

summary_plot = do.call("grid.arrange", c(l.p, nrow = 5, ncol = 5))
ggsave("./results/plot_summary.pdf", plot = summary_plot, dpi = 300, width = 16, height = 9)

save(file = "./results/summary_data.RData", summary_plot, dat.R, dat.Cpp, time.r, time.c)
summary_xtable_R = xtable(t(dat.R))
summary_xtable_Cpp = xtable(t(dat.Cpp))

summary_time_xtab_R = xtable(time.r)
summary_time_xtab_C = xtable(time.c)
print(summary_xtable_R, file = "./results/tab_summary_r.tex")
print(summary_xtable_Cpp, file = "./results/tab_summary_cpp.tex")
print(summary_time_xtab_R, file = "./results/tab_summaryTime_r.tex")
print(summary_time_xtab_C, file = "./results/tab_summaryTime_Cpp.tex")
    
######## t-Test study
    
    ### type I Error for one sample t-Test
    
    time.r = time.c = matrix(NA, nrow = nrow(dat.R), ncol = 6)
    dat.R.np = dat.R.w = dat.C.np = dat.C.w = cbind(base_grid, "normal" = NA, "poisson" = NA, "exponential" =  NA, "chisq" = NA)
    
    for(i in 1:nrow(base_grid)){
        A.temp.r.np = A.temp.r.w = A.temp.c.np = A.temp.c.w = matrix(NA, nrow = 4, ncol = 1000) # X from: Normal, Poisson, Exponential, ChiSquare
        
        
        A.temp.r.np[1,] = replicate(1000, t.testBoot(rnorm(n = dat.R$n[i]), nboot = dat.R$nboots[i])$reject)
        A.temp.r.w[1,] = replicate(1000, t.testBoot(rnorm(n = dat.R$n[i]), nboot = dat.R$nboots[i], boot.type = "wild")$reject)
        A.temp.c.np[1,] = replicate(1000, t.testBoot.cpp(rnorm(n = dat.R$n[i]), nboot = dat.R$nboots[i])$reject)
        A.temp.c.w[1,] = replicate(1000, t.testBoot.cpp(rnorm(n = dat.R$n[i]), nboot = dat.R$nboots[i], boot.type = "wild")$reject)
        
        A.temp.r.np[2,] = replicate(1000, t.testBoot(rpois(n = dat.R$n[i], lambda = 5), mu.0 = 5, nboot = dat.R$nboots[i])$reject)
        A.temp.r.w[2,] = replicate(1000, t.testBoot(rpois(n = dat.R$n[i], lambda = 5), mu.0 = 5, nboot = dat.R$nboots[i], boot.type = "wild")$reject)
        A.temp.c.np[2,] = replicate(1000, t.testBoot.cpp(rpois(n = dat.R$n[i], lambda = 5), mu.0 = 5, nboot = dat.R$nboots[i])$reject)
        A.temp.c.w[2,] = replicate(1000, t.testBoot.cpp(rpois(n = dat.R$n[i], lambda = 5), mu.0 = 5, nboot = dat.R$nboots[i], boot.type = "wild")$reject)
        
        A.temp.r.np[3,] = replicate(1000, t.testBoot(rexp(n = dat.R$n[i], rate = 3), mu.0 = 1/3, nboot = dat.R$nboots[i])$reject)
        A.temp.r.w[3,] = replicate(1000, t.testBoot(rexp(n = dat.R$n[i], rate = 3), mu.0 = 1/3, nboot = dat.R$nboots[i], boot.type = "wild")$reject)
        A.temp.c.np[3,] = replicate(1000, t.testBoot.cpp(rexp(n = dat.R$n[i], rate = 3), mu.0 = 1/3, nboot = dat.R$nboots[i])$reject)
        A.temp.c.w[3,] = replicate(1000, t.testBoot.cpp(rexp(n = dat.R$n[i], rate = 3), mu.0 = 1/3, nboot = dat.R$nboots[i], boot.type = "wild")$reject)
        
        A.temp.r.np[4,] = replicate(1000, t.testBoot(rchisq(n = dat.R$n[i], df = 2), mu.0 = 2, nboot = dat.R$nboots[i])$reject)
        A.temp.r.w[4,] = replicate(1000, t.testBoot(rchisq(n = dat.R$n[i], df = 2), mu.0 = 2, nboot = dat.R$nboots[i], boot.type = "wild")$reject)
        A.temp.c.np[4,] = replicate(1000, t.testBoot.cpp(rchisq(n = dat.R$n[i], df = 2), mu.0 = 2, nboot = dat.R$nboots[i])$reject)
        A.temp.c.w[4,] = replicate(1000, t.testBoot.cpp(rchisq(n = dat.R$n[i], df = 2), mu.0 = 2, nboot = dat.R$nboots[i], boot.type = "wild")$reject)
        
        dat.R.np[i, c(3:6)] = rowMeans(A.temp.r.np)
        dat.R.w[i, c(3:6)] = rowMeans(A.temp.r.w)
        dat.C.np[i, c(3:6)] = rowMeans(A.temp.c.np)
        dat.C.w[i, c(3:6)] = rowMeans(A.temp.c.w)
        
        x.sample = rchisq(n = dat.R$n[i], df = 2)
        
        time.r[i,] = summary(microbenchmark(t.testBoot(x.sample, mu.0 = 2, nboot = dat.R$nboots[i]))$time / 1000) # microseconds
        time.c[i,] = summary(microbenchmark(t.testBoot.cpp(x.sample, mu.0 = 2, nboot = dat.R$nboots[i]))$time / 1000)
    }
    
    
    time.r = cbind(base_grid, time.r)
    time.c = cbind(base_grid, time.c)
    
    colnames(time.r) = colnames(time.c) = c(colnames(base_grid), names(summary(rnorm(10))))
    
    l.p = list()
    ind = seq(1, nrow(base_grid), 5)
    for(i in seq(1, nrow(base_grid), 5)){
        
        l.p[[i]] = ggplot(data = dat.R.np[i + c(0:4),])+
            geom_line(mapping = aes(x = nboots, y = normal), color = "black") + 
            geom_line(mapping = aes(x = nboots, y = poisson), color = "blue") +
            geom_line(mapping = aes(x = nboots, y = exponential), color = "red") + 
            geom_line(mapping = aes(x = nboots, y = chisq), color = "green") +
            labs(x = "Bootstrap Iterations", y = expression(paste("Type-I-Error (", alpha, ")")))
        
        l.p[[i + 1]] = ggplot(data = dat.C.np[i + c(0:4),])+
            geom_line(mapping = aes(x = nboots, y = normal), color = "black") + 
            geom_line(mapping = aes(x = nboots, y = poisson), color = "blue") +
            geom_line(mapping = aes(x = nboots, y = exponential), color = "red") + 
            geom_line(mapping = aes(x = nboots, y = chisq), color = "green") +
            labs(x = "Bootstrap Iterations", y = expression(paste("Type-I-Error (", alpha, ")")))
        
        l.p[[i + 2]] = ggplot(data = dat.R.w[i + c(0:4),])+
            geom_line(mapping = aes(x = nboots, y = normal), color = "black") + 
            geom_line(mapping = aes(x = nboots, y = poisson), color = "blue") +
            geom_line(mapping = aes(x = nboots, y = exponential), color = "red") + 
            geom_line(mapping = aes(x = nboots, y = chisq), color = "green") +
            labs(x = "Bootstrap Iterations", y = expression(paste("Type-I-Error (", alpha, ")")))
        
        l.p[[i + 3]] = ggplot(data = dat.C.w[i + c(0:4),])+
            geom_line(mapping = aes(x = nboots, y = normal), color = "black") + 
            geom_line(mapping = aes(x = nboots, y = poisson), color = "blue") +
            geom_line(mapping = aes(x = nboots, y = exponential), color = "red") + 
            geom_line(mapping = aes(x = nboots, y = chisq), color = "green") +
            labs(x = "Bootstrap Iterations", y = expression(paste("Type-I-Error (", alpha, ")")))
        
        
        l.p[[i + 4]] = ggplot()+
            geom_line(mapping = aes(x = time.r[i + c(0:4), "nboots"], y = time.r[i + c(0:4), 3]), color = "red") +
            geom_line(mapping = aes(x = time.c[i + c(0:4), "nboots"], y = time.c[i + c(0:4), 3]), color = "blue") +
            labs(x = "Bootstrap Iterations", y = "Microseconds") # on Median
        
    }
    
    type1_plot = do.call("grid.arrange", c(l.p, nrow = 5, ncol = 5))
    ggsave("./results/plot_t1s1.pdf", plot = type1_plot, dpi = 300, width = 16, height = 9)
    
    
    ts1_time_xtab_R = xtable(time.r)
    ts1_time_xtab_C = xtable(time.c)
    
    ts1_xtab_R_np = xtable(dat.R.np)
    ts1_xtab_R_w = xtable(dat.R.w)
    ts1_xtab_C_np = xtable(dat.C.np)
    ts1_xtab_C_w = xtable(dat.C.w)
    
    print(ts1_time_xtab_R, file = "./results/tab_ts1Time_r.tex")
    print(ts1_time_xtab_C, file = "./results/tab_ts1Time_Cpp.tex")
    
    ### Power Study
    pow.true.1 = pow.true.2 = matrix(NA, nrow = 25, ncol = 5)
    
    dat.1 = dat.2 = dat.3 = dat.4 = dat.5 = cbind(base_grid, "R_np" = NA, "R_w" = NA, "C_np" =  NA, "C_w" = NA)
    
    
    for(i in 1:nrow(base_grid)){
        A.temp.1 = A.temp.2 = A.temp.3 = A.temp.4 = A.temp.5 = matrix(NA, nrow = 4, ncol = 1000)
        
        A.temp.1[1,] = replicate(1000, t.testBoot(rnorm(base_grid$n[i], 0.5, 1), nboot = base_grid$nboots[i])$reject)
        A.temp.1[2,] = replicate(1000, t.testBoot(rnorm(base_grid$n[i], 0.5, 1), nboot = base_grid$nboots[i], boot.type = "wild")$reject)
        A.temp.1[3,] = replicate(1000, t.testBoot.cpp(rnorm(base_grid$n[i], 0.5, 1), nboot = base_grid$nboots[i])$reject)
        A.temp.1[4,] = replicate(1000, t.testBoot.cpp(rnorm(base_grid$n[i], 0.5, 1), nboot = base_grid$nboots[i], boot.type = "wild")$reject)
        
        A.temp.2[1,] = replicate(1000, t.testBoot(rnorm(base_grid$n[i], 0.75, 1), nboot = base_grid$nboots[i])$reject)
        A.temp.2[2,] = replicate(1000, t.testBoot(rnorm(base_grid$n[i], 0.75, 1), nboot = base_grid$nboots[i], boot.type = "wild")$reject)
        A.temp.2[3,] = replicate(1000, t.testBoot.cpp(rnorm(base_grid$n[i], 0.75, 1), nboot = base_grid$nboots[i])$reject)
        A.temp.2[4,] = replicate(1000, t.testBoot.cpp(rnorm(base_grid$n[i], 0.75, 1), nboot = base_grid$nboots[i], boot.type = "wild")$reject)
        
        A.temp.3[1,] = replicate(1000, t.testBoot(rnorm(base_grid$n[i], 1, 1), nboot = base_grid$nboots[i])$reject)
        A.temp.3[2,] = replicate(1000, t.testBoot(rnorm(base_grid$n[i], 1, 1), nboot = base_grid$nboots[i], boot.type = "wild")$reject)
        A.temp.3[3,] = replicate(1000, t.testBoot.cpp(rnorm(base_grid$n[i], 1, 1), nboot = base_grid$nboots[i])$reject)
        A.temp.3[4,] = replicate(1000, t.testBoot.cpp(rnorm(base_grid$n[i], 1, 1), nboot = base_grid$nboots[i], boot.type = "wild")$reject)
        
        A.temp.4[1,] = replicate(1000, t.testBoot(rnorm(base_grid$n[i], 1.25, 1), nboot = base_grid$nboots[i])$reject)
        A.temp.4[2,] = replicate(1000, t.testBoot(rnorm(base_grid$n[i], 1.25, 1), nboot = base_grid$nboots[i], boot.type = "wild")$reject)
        A.temp.4[3,] = replicate(1000, t.testBoot.cpp(rnorm(base_grid$n[i], 1.25, 1), nboot = base_grid$nboots[i])$reject)
        A.temp.4[4,] = replicate(1000, t.testBoot.cpp(rnorm(base_grid$n[i], 1.25, 1), nboot = base_grid$nboots[i], boot.type = "wild")$reject)
        
        A.temp.5[1,] = replicate(1000, t.testBoot(rnorm(base_grid$n[i], 1.5, 1), nboot = base_grid$nboots[i])$reject)
        A.temp.5[2,] = replicate(1000, t.testBoot(rnorm(base_grid$n[i], 1.5, 1), nboot = base_grid$nboots[i], boot.type = "wild")$reject)
        A.temp.5[3,] = replicate(1000, t.testBoot.cpp(rnorm(base_grid$n[i], 1.5, 1), nboot = base_grid$nboots[i])$reject)
        A.temp.5[4,] = replicate(1000, t.testBoot.cpp(rnorm(base_grid$n[i], 1.5, 1), nboot = base_grid$nboots[i], boot.type = "wild")$reject)
        
        dat.1[i, c(3:6)] = rowMeans(A.temp.1)
        dat.2[i, c(3:6)] = rowMeans(A.temp.2)
        dat.3[i, c(3:6)] = rowMeans(A.temp.3)
        dat.4[i, c(3:6)] = rowMeans(A.temp.4)
        dat.5[i, c(3:6)] = rowMeans(A.temp.5)
    }
    
    for(i in seq(1, nrow(base_grid), 5)){#1:nrow(base_grid)){
        pow.true.1[i + c(0:4),1] = mean(t.power.fun(replicate(1000, t.test(rnorm(base_grid$n[i], 0.5, 1))$p.value)))
        pow.true.1[i + c(0:4),2] = mean(t.power.fun(replicate(1000, t.test(rnorm(base_grid$n[i], 0.75, 1))$p.value)))
        pow.true.1[i + c(0:4),3] = mean(t.power.fun(replicate(1000, t.test(rnorm(base_grid$n[i], 1, 1))$p.value)))
        pow.true.1[i + c(0:4),4] = mean(t.power.fun(replicate(1000, t.test(rnorm(base_grid$n[i], 1.25, 1))$p.value)))
        pow.true.1[i + c(0:4),5] = mean(t.power.fun(replicate(1000, t.test(rnorm(base_grid$n[i], 1.5, 1))$p.value)))
        
        pow.true.2[i + c(0:4),1] = mean(t.power.fun(replicate(1000, t.test(rnorm(base_grid$n[i], 0.5, 1), rnorm(base_grid$n[i]))$p.value)))
        pow.true.2[i + c(0:4),2] = mean(t.power.fun(replicate(1000, t.test(rnorm(base_grid$n[i], 0.75, 1), rnorm(base_grid$n[i]))$p.value)))
        pow.true.2[i + c(0:4),3] = mean(t.power.fun(replicate(1000, t.test(rnorm(base_grid$n[i], 1, 1), rnorm(base_grid$n[i]))$p.value)))
        pow.true.2[i + c(0:4),4] = mean(t.power.fun(replicate(1000, t.test(rnorm(base_grid$n[i], 1.25, 1), rnorm(base_grid$n[i]))$p.value)))
        pow.true.2[i + c(0:4),5] = mean(t.power.fun(replicate(1000, t.test(rnorm(base_grid$n[i], 1.5, 1), rnorm(base_grid$n[i]))$p.value)))
        
    }
    
    dat.1$true = pow.true.1[,1]
    dat.2$true = pow.true.1[,2]
    dat.3$true = pow.true.1[,3]
    dat.4$true = pow.true.1[,4]
    dat.5$true = pow.true.1[,5]
    
    l.p = list()
    for(i in seq(1, nrow(base_grid), 5)){
        
        l.p[[i]] = ggplot(data = dat.1[i + c(0:4),])+
            geom_line(mapping = aes(x = nboots, y = R_np), color = "black") + 
            geom_line(mapping = aes(x = nboots, y = R_w), color = "blue") +
            geom_line(mapping = aes(x = nboots, y = C_np), color = "red") + 
            geom_line(mapping = aes(x = nboots, y = C_w), color = "green") +
            geom_line(mapping = aes(x = nboots, y = true), linetype = "dashed", color = "chocolate3") + 
            labs(x = "Bootstrap Iterations", y = expression(paste("Power (1 - ", beta, ")")))
        
        l.p[[i + 1]] = ggplot(data = dat.2[i + c(0:4),])+
            geom_line(mapping = aes(x = nboots, y = R_np), color = "black") + 
            geom_line(mapping = aes(x = nboots, y = R_w), color = "blue") +
            geom_line(mapping = aes(x = nboots, y = C_np), color = "red") + 
            geom_line(mapping = aes(x = nboots, y = C_w), color = "green") +
            geom_line(mapping = aes(x = nboots, y = true), linetype = "dashed", color = "chocolate3") + 
            labs(x = "Bootstrap Iterations", y = expression(paste("Power (1 - ", beta, ")")))
        
        l.p[[i + 2]] = ggplot(data = dat.3[i + c(0:4),])+
            geom_line(mapping = aes(x = nboots, y = R_np), color = "black") + 
            geom_line(mapping = aes(x = nboots, y = R_w), color = "blue") +
            geom_line(mapping = aes(x = nboots, y = C_np), color = "red") + 
            geom_line(mapping = aes(x = nboots, y = C_w), color = "green") +
            geom_line(mapping = aes(x = nboots, y = true), linetype = "dashed", color = "chocolate3") + 
            labs(x = "Bootstrap Iterations", y = expression(paste("Power (1 - ", beta, ")")))
        
        l.p[[i + 3]] = ggplot(data = dat.4[i + c(0:4),])+
            geom_line(mapping = aes(x = nboots, y = R_np), color = "black") + 
            geom_line(mapping = aes(x = nboots, y = R_w), color = "blue") +
            geom_line(mapping = aes(x = nboots, y = C_np), color = "red") + 
            geom_line(mapping = aes(x = nboots, y = C_w), color = "green") +
            geom_line(mapping = aes(x = nboots, y = true), linetype = "dashed", color = "chocolate3") + 
            labs(x = "Bootstrap Iterations", y = expression(paste("Power (1 - ", beta, ")")))
        
        l.p[[i + 4]] = ggplot(data = dat.5[i + c(0:4),])+
            geom_line(mapping = aes(x = nboots, y = R_np), color = "black") + 
            geom_line(mapping = aes(x = nboots, y = R_w), color = "blue") +
            geom_line(mapping = aes(x = nboots, y = C_np), color = "red") + 
            geom_line(mapping = aes(x = nboots, y = C_w), color = "green") +
            geom_line(mapping = aes(x = nboots, y = true), linetype = "dashed", color = "chocolate3") + 
            labs(x = "Bootstrap Iterations", y = expression(paste("Power (1 - ", beta, ")")))
        
    }
    
    
    
    type2_plot = do.call("grid.arrange", c(l.p, nrow = 5, ncol = 5))
    ggsave("./results/plot_t2s1.pdf", plot = type2_plot, dpi = 300, width = 16, height = 9)
    
    save(file = "./results/ts1_data.RData", type1_plot, dat.R.np, dat.C.np, dat.R.w, dat.C.w, time.r, time.c, dat.1, dat.2, dat.3, dat.4, dat.5, type2_plot)
    
    
    ### type I Error for two sample t-Test
    time.r = time.c = matrix(NA, nrow = nrow(dat.R), ncol = 6)
    dat.R.np = dat.R.w = dat.C.np = dat.C.w = dat.R.npg = dat.C.npg = cbind(base_grid, "normal" = NA, "poisson" = NA, "exponential" =  NA, "chisq" = NA)
    
    for(i in 1:nrow(base_grid)){
        A.temp.r.np = A.temp.r.w = A.temp.c.np = A.temp.c.w =  matrix(NA, nrow = 4, ncol = 1000) # X from: Normal, Poisson, Exponential, ChiSquare
        
        
        A.temp.r.np[1,] = replicate(1000, t.testBoot(rnorm(n = dat.R$n[i]), rnorm(n = dat.R$n[i]), nboot = dat.R$nboots[i])$reject)
        A.temp.r.w[1,] = replicate(1000, t.testBoot(rnorm(n = dat.R$n[i]), rnorm(n = dat.R$n[i]), nboot = dat.R$nboots[i], boot.type = "wild")$reject)
        A.temp.c.np[1,] = replicate(1000, t.testBoot.cpp(rnorm(n = dat.R$n[i]), rnorm(n = dat.R$n[i]), nboot = dat.R$nboots[i])$reject)
        A.temp.c.w[1,] = replicate(1000, t.testBoot.cpp(rnorm(n = dat.R$n[i]), rnorm(n = dat.R$n[i]), nboot = dat.R$nboots[i], boot.type = "wild")$reject)
        # A.temp.r.npg[1,] = replicate(1000, t.testBoot(rnorm(n = dat.R$n[i]), rnorm(n = dat.R$n[i]), nboot = dat.R$nboots[i], boot.type = "npg")$reject)
        # A.temp.c.npg[1,] = replicate(1000, t.testBoot.cpp(rnorm(n = dat.R$n[i]), rnorm(n = dat.R$n[i]), nboot = dat.R$nboots[i], boot.type = "npg")$reject)
        
        A.temp.r.np[2,] = replicate(1000, t.testBoot(rpois(n = dat.R$n[i], lambda = 5), rpois(n = dat.R$n[i], lambda = 5), mu.0 = 5, nboot = dat.R$nboots[i])$reject)
        A.temp.r.w[2,] = replicate(1000, t.testBoot(rpois(n = dat.R$n[i], lambda = 5), rpois(n = dat.R$n[i], lambda = 5), mu.0 = 5, nboot = dat.R$nboots[i], boot.type = "wild")$reject)
        A.temp.c.np[2,] = replicate(1000, t.testBoot.cpp(rpois(n = dat.R$n[i], lambda = 5), rpois(n = dat.R$n[i], lambda = 5), mu.0 = 5, nboot = dat.R$nboots[i])$reject)
        A.temp.c.w[2,] = replicate(1000, t.testBoot.cpp(rpois(n = dat.R$n[i], lambda = 5), rpois(n = dat.R$n[i], lambda = 5), mu.0 = 5, nboot = dat.R$nboots[i], boot.type = "wild")$reject)
        # A.temp.r.npg[2,] = replicate(1000, t.testBoot(rpois(n = dat.R$n[i], lambda = 5), rpois(n = dat.R$n[i], lambda = 5), mu.0 = 5, nboot = dat.R$nboots[i], boot.type = "npg")$reject)
        # A.temp.c.npg[2,] = replicate(1000, t.testBoot.cpp(rpois(n = dat.R$n[i], lambda = 5), rpois(n = dat.R$n[i], lambda = 5), mu.0 = 5, nboot = dat.R$nboots[i], boot.type = "npg")$reject)
        
        A.temp.r.np[3,] = replicate(1000, t.testBoot(rexp(n = dat.R$n[i], rate = 3), rexp(n = dat.R$n[i], rate = 3), mu.0 = 1/3, nboot = dat.R$nboots[i])$reject)
        A.temp.r.w[3,] = replicate(1000, t.testBoot(rexp(n = dat.R$n[i], rate = 3), rexp(n = dat.R$n[i], rate = 3), mu.0 = 1/3, nboot = dat.R$nboots[i], boot.type = "wild")$reject)
        A.temp.c.np[3,] = replicate(1000, t.testBoot.cpp(rexp(n = dat.R$n[i], rate = 3), rexp(n = dat.R$n[i], rate = 3), mu.0 = 1/3, nboot = dat.R$nboots[i])$reject)
        A.temp.c.w[3,] = replicate(1000, t.testBoot.cpp(rexp(n = dat.R$n[i], rate = 3), rexp(n = dat.R$n[i], rate = 3), mu.0 = 1/3, nboot = dat.R$nboots[i], boot.type = "wild")$reject)
        # A.temp.r.npg[3,] = replicate(1000, t.testBoot(rexp(n = dat.R$n[i], rate = 3), rexp(n = dat.R$n[i], rate = 3), mu.0 = 1/3, nboot = dat.R$nboots[i], boot.type = "npg")$reject)
        # A.temp.c.npg[3,] = replicate(1000, t.testBoot.cpp(rexp(n = dat.R$n[i], rate = 3), rexp(n = dat.R$n[i], rate = 3), mu.0 = 1/3, nboot = dat.R$nboots[i], boot.type = "npg")$reject)
        
        A.temp.r.np[4,] = replicate(1000, t.testBoot(rchisq(n = dat.R$n[i], df = 2), rchisq(n = dat.R$n[i], df = 2), mu.0 = 2, nboot = dat.R$nboots[i])$reject)
        A.temp.r.w[4,] = replicate(1000, t.testBoot(rchisq(n = dat.R$n[i], df = 2), rchisq(n = dat.R$n[i], df = 2), mu.0 = 2, nboot = dat.R$nboots[i], boot.type = "wild")$reject)
        A.temp.c.np[4,] = replicate(1000, t.testBoot.cpp(rchisq(n = dat.R$n[i], df = 2), rchisq(n = dat.R$n[i], df = 2), mu.0 = 2, nboot = dat.R$nboots[i])$reject)
        A.temp.c.w[4,] = replicate(1000, t.testBoot.cpp(rchisq(n = dat.R$n[i], df = 2), rchisq(n = dat.R$n[i], df = 2), mu.0 = 2, nboot = dat.R$nboots[i], boot.type = "wild")$reject)
        # A.temp.r.npg[4,] = replicate(1000, t.testBoot(rchisq(n = dat.R$n[i], df = 2), rchisq(n = dat.R$n[i], df = 2), mu.0 = 2, nboot = dat.R$nboots[i], boot.type = "npg")$reject)
        # A.temp.c.npg[4,] = replicate(1000, t.testBoot.cpp(rchisq(n = dat.R$n[i], df = 2), rchisq(n = dat.R$n[i], df = 2), mu.0 = 2, nboot = dat.R$nboots[i], boot.type = "npg")$reject)
        
        dat.R.np[i, c(3:6)] = rowMeans(A.temp.r.np)
        dat.R.w[i, c(3:6)] = rowMeans(A.temp.r.w)
        dat.C.np[i, c(3:6)] = rowMeans(A.temp.c.np)
        dat.C.w[i, c(3:6)] = rowMeans(A.temp.c.w)
        # dat.R.npg[i, c(3:6)] = rowMeans(A.temp.r.npg)
        # dat.C.npg[i, c(3:6)] = rowMeans(A.temp.c.npg)
        
        x.sample.1 = rchisq(n = base_grid$n[i], df = 2)
        x.sample.2 = rchisq(n = base_grid$n[i], df = 2)
        
        time.r[i,] = summary(microbenchmark(t.testBoot(x.sample.1, x.sample.2, mu.0 = 2, nboot = dat.R$nboots[i]))$time / 1000) # microseconds
        time.c[i,] = summary(microbenchmark(t.testBoot.cpp(x.sample.1, x.sample.2, mu.0 = 2, nboot = dat.R$nboots[i]))$time / 1000)
    }
    
    
    time.r = cbind(base_grid, time.r)
    time.c = cbind(base_grid, time.c)
    
    colnames(time.r) = colnames(time.c) = c(colnames(base_grid), names(summary(rnorm(10))))
    
    l.p = list()
    # ind = seq(1, 35, 7)
    ind = ind.2 = seq(1, nrow(base_grid), 5)
    for(i in c(1:5)){
        
        l.p[[ind[i]]] = ggplot(data = dat.R.np[ind.2[i] + c(0:4),])+
            geom_line(mapping = aes(x = nboots, y = normal), color = "black") + 
            geom_line(mapping = aes(x = nboots, y = poisson), color = "blue") +
            geom_line(mapping = aes(x = nboots, y = exponential), color = "red") + 
            geom_line(mapping = aes(x = nboots, y = chisq), color = "green") +
            labs(x = "Bootstrap Iterations", y = expression(paste("Type-I-Error (", alpha, ")")))
        
        l.p[[ind[i] + 1]] = ggplot(data = dat.C.np[ind.2[i] + c(0:4),])+
            geom_line(mapping = aes(x = nboots, y = normal), color = "black") + 
            geom_line(mapping = aes(x = nboots, y = poisson), color = "blue") +
            geom_line(mapping = aes(x = nboots, y = exponential), color = "red") + 
            geom_line(mapping = aes(x = nboots, y = chisq), color = "green") +
            labs(x = "Bootstrap Iterations", y = expression(paste("Type-I-Error (", alpha, ")")))
        
        l.p[[ind[i] + 2]] = ggplot(data = dat.R.w[ind.2[i] + c(0:4),])+
            geom_line(mapping = aes(x = nboots, y = normal), color = "black") + 
            geom_line(mapping = aes(x = nboots, y = poisson), color = "blue") +
            geom_line(mapping = aes(x = nboots, y = exponential), color = "red") + 
            geom_line(mapping = aes(x = nboots, y = chisq), color = "green") +
            labs(x = "Bootstrap Iterations", y = expression(paste("Type-I-Error (", alpha, ")")))
        
        l.p[[ind[i] + 3]] = ggplot(data = dat.C.w[ind.2[i] + c(0:4),])+
            geom_line(mapping = aes(x = nboots, y = normal), color = "black") + 
            geom_line(mapping = aes(x = nboots, y = poisson), color = "blue") +
            geom_line(mapping = aes(x = nboots, y = exponential), color = "red") + 
            geom_line(mapping = aes(x = nboots, y = chisq), color = "green") +
            labs(x = "Bootstrap Iterations", y = expression(paste("Type-I-Error (", alpha, ")")))
        
        # l.p[[ind[i] + 4]] = ggplot(data = dat.R.npg[ind.2[i] + c(0:4),])+
        #     geom_line(mapping = aes(x = nboots, y = normal), color = "black") + 
        #     geom_line(mapping = aes(x = nboots, y = poisson), color = "blue") +
        #     geom_line(mapping = aes(x = nboots, y = exponential), color = "red") + 
        #     geom_line(mapping = aes(x = nboots, y = chisq), color = "green") +
        #     labs(x = "Bootstrap Iterations", y = paste0("Type-I-Error (", expression(alpha), ")"))
        # 
        # l.p[[ind[i] + 5]] = ggplot(data = dat.C.npg[ind.2[i] + c(0:4),])+
        #     geom_line(mapping = aes(x = nboots, y = normal), color = "black") + 
        #     geom_line(mapping = aes(x = nboots, y = poisson), color = "blue") +
        #     geom_line(mapping = aes(x = nboots, y = exponential), color = "red") + 
        #     geom_line(mapping = aes(x = nboots, y = chisq), color = "green") +
        #     labs(x = "Bootstrap Iterations", y = paste0("Type-I-Error (", expression(alpha), ")"))
        
        l.p[[ind[i] + 4]] = ggplot()+
            geom_line(mapping = aes(x = time.r[ind.2[i] + c(0:4), "nboots"], y = time.r[ind.2[i] + c(0:4), 3]), color = "red") +
            geom_line(mapping = aes(x = time.c[ind.2[i] + c(0:4), "nboots"], y = time.c[ind.2[i] + c(0:4), 3]), color = "blue") +
            labs(x = "Bootstrap Iterations", y = "Microseconds") # on Median
        
    }
    
    type1_plot = do.call("grid.arrange", c(l.p, nrow = 5, ncol = 5))
    ggsave("./results/plot_t1s2.pdf", plot = type1_plot, dpi = 300, width = 16, height = 9)
    
    ts2_time_xtab_R = xtable(time.r)
    ts2_time_xtab_C = xtable(time.c)
    
    print(ts2_time_xtab_R, file = "./results/tab_ts2Time_r.tex")
    print(ts2_time_xtab_C, file = "./results/tab_ts2Time_Cpp.tex")
    
    ### Power Study
    
    dat.1 = dat.2 = dat.3 = dat.4 = dat.5 = cbind(base_grid, "R_np" = NA, "R_w" = NA, "C_np" =  NA, "C_w" = NA)
    
    for(i in 1:nrow(base_grid)){
        A.temp.1 = A.temp.2 = A.temp.3 = A.temp.4 = A.temp.5 = matrix(NA, nrow = 4, ncol = 1000)
        
        A.temp.1[1,] = replicate(1000, t.testBoot(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 0.5, 1), nboot = base_grid$nboots[i])$reject)
        A.temp.1[2,] = replicate(1000, t.testBoot(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 0.5, 1), nboot = base_grid$nboots[i], boot.type = "wild")$reject)
        A.temp.1[3,] = replicate(1000, t.testBoot.cpp(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 0.5, 1), nboot = base_grid$nboots[i])$reject)
        A.temp.1[4,] = replicate(1000, t.testBoot.cpp(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 0.5, 1), nboot = base_grid$nboots[i], boot.type = "wild")$reject)
        # A.temp.1[5,] = replicate(1000, t.testBoot(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 0.5, 1), nboot = base_grid$nboots[i], boot.type = "npg")$reject)
        # A.temp.1[6,] = replicate(1000, t.testBoot.cpp(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 0.5, 1), nboot = base_grid$nboots[i], boot.type = "npg")$reject)
        # 
        A.temp.2[1,] = replicate(1000, t.testBoot(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 0.75, 1), nboot = base_grid$nboots[i])$reject)
        A.temp.2[2,] = replicate(1000, t.testBoot(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 0.75, 1), nboot = base_grid$nboots[i], boot.type = "wild")$reject)
        A.temp.2[3,] = replicate(1000, t.testBoot.cpp(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 0.75, 1), nboot = base_grid$nboots[i])$reject)
        A.temp.2[4,] = replicate(1000, t.testBoot.cpp(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 0.75, 1), nboot = base_grid$nboots[i], boot.type = "wild")$reject)
        # A.temp.2[5,] = replicate(1000, t.testBoot(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 0.75, 1), nboot = base_grid$nboots[i], boot.type = "npg")$reject)
        # A.temp.2[6,] = replicate(1000, t.testBoot.cpp(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 0.75, 1), nboot = base_grid$nboots[i], boot.type = "npg")$reject)
        # 
        A.temp.3[1,] = replicate(1000, t.testBoot(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 1, 1), nboot = base_grid$nboots[i])$reject)
        A.temp.3[2,] = replicate(1000, t.testBoot(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 1, 1), nboot = base_grid$nboots[i], boot.type = "wild")$reject)
        A.temp.3[3,] = replicate(1000, t.testBoot.cpp(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 1, 1), nboot = base_grid$nboots[i])$reject)
        A.temp.3[4,] = replicate(1000, t.testBoot.cpp(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 1, 1), nboot = base_grid$nboots[i], boot.type = "wild")$reject)
        # A.temp.3[5,] = replicate(1000, t.testBoot(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 1, 1), nboot = base_grid$nboots[i], boot.type = "npg")$reject)
        # A.temp.3[6,] = replicate(1000, t.testBoot.cpp(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 1, 1), nboot = base_grid$nboots[i], boot.type = "npg")$reject)
        # 
        A.temp.4[1,] = replicate(1000, t.testBoot(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 1.25, 1), nboot = base_grid$nboots[i])$reject)
        A.temp.4[2,] = replicate(1000, t.testBoot(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 1.25, 1), nboot = base_grid$nboots[i], boot.type = "wild")$reject)
        A.temp.4[3,] = replicate(1000, t.testBoot.cpp(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 1.25, 1), nboot = base_grid$nboots[i])$reject)
        A.temp.4[4,] = replicate(1000, t.testBoot.cpp(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 1.25, 1), nboot = base_grid$nboots[i], boot.type = "wild")$reject)
        # A.temp.4[5,] = replicate(1000, t.testBoot(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 1.25, 1), nboot = base_grid$nboots[i], boot.type = "npg")$reject)
        # A.temp.4[6,] = replicate(1000, t.testBoot.cpp(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 1.25, 1), nboot = base_grid$nboots[i], boot.type = "npg")$reject)
        # 
        A.temp.5[1,] = replicate(1000, t.testBoot(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 1.5, 1), nboot = base_grid$nboots[i])$reject)
        A.temp.5[2,] = replicate(1000, t.testBoot(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 1.5, 1), nboot = base_grid$nboots[i], boot.type = "wild")$reject)
        A.temp.5[3,] = replicate(1000, t.testBoot.cpp(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 1.5, 1), nboot = base_grid$nboots[i])$reject)
        A.temp.5[4,] = replicate(1000, t.testBoot.cpp(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 1.5, 1), nboot = base_grid$nboots[i], boot.type = "wild")$reject)
        # A.temp.5[5,] = replicate(1000, t.testBoot(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 1.5, 1), nboot = base_grid$nboots[i], boot.type = "npg")$reject)
        # A.temp.5[6,] = replicate(1000, t.testBoot.cpp(rnorm(base_grid$n[i], 0, 1), rnorm(base_grid$n[i], 1.5, 1), nboot = base_grid$nboots[i], boot.type = "npg")$reject)
        # 
        dat.1[i, c(3:6)] = rowMeans(A.temp.1)
        dat.2[i, c(3:6)] = rowMeans(A.temp.2)
        dat.3[i, c(3:6)] = rowMeans(A.temp.3)
        dat.4[i, c(3:6)] = rowMeans(A.temp.4)
        dat.5[i, c(3:6)] = rowMeans(A.temp.5)
    }
    
    dat.1$true = pow.true.2[,1]
    dat.2$true = pow.true.2[,2]
    dat.3$true = pow.true.2[,3]
    dat.4$true = pow.true.2[,4]
    dat.5$true = pow.true.2[,5]
    
    l.p = list()
    for(i in seq(1, nrow(base_grid), 5)){
        
        l.p[[i]] = ggplot(data = dat.1[i + c(0:4),])+
            geom_line(mapping = aes(x = nboots, y = R_np), color = "black") + 
            geom_line(mapping = aes(x = nboots, y = R_w), color = "blue") +
            geom_line(mapping = aes(x = nboots, y = C_np), color = "red") + 
            geom_line(mapping = aes(x = nboots, y = C_w), color = "green") +
            # geom_line(mapping = aes(x = nboots, y = C_npg), color = "orange") + 
            # geom_line(mapping = aes(x = nboots, y = R_npg), color = "purple") +
            geom_line(mapping = aes(x = nboots, y = true), linetype = "dashed", color = "chocolate3") + 
            labs(x = "Bootstrap Iterations", y = expression(paste("Power (1 - ", beta, ")")))
        
        l.p[[i + 1]] = ggplot(data = dat.2[i + c(0:4),])+
            geom_line(mapping = aes(x = nboots, y = R_np), color = "black") + 
            geom_line(mapping = aes(x = nboots, y = R_w), color = "blue") +
            geom_line(mapping = aes(x = nboots, y = C_np), color = "red") + 
            geom_line(mapping = aes(x = nboots, y = C_w), color = "green") +
            # geom_line(mapping = aes(x = nboots, y = C_npg), color = "orange") + 
            # geom_line(mapping = aes(x = nboots, y = R_npg), color = "purple") +
            geom_line(mapping = aes(x = nboots, y = true), linetype = "dashed", color = "chocolate3") + 
            labs(x = "Bootstrap Iterations", y = expression(paste("Power (1 - ", beta, ")")))
        
        l.p[[i + 2]] = ggplot(data = dat.3[i + c(0:4),])+
            geom_line(mapping = aes(x = nboots, y = R_np), color = "black") + 
            geom_line(mapping = aes(x = nboots, y = R_w), color = "blue") +
            geom_line(mapping = aes(x = nboots, y = C_np), color = "red") + 
            geom_line(mapping = aes(x = nboots, y = C_w), color = "green") +
            # geom_line(mapping = aes(x = nboots, y = C_npg), color = "orange") + 
            # geom_line(mapping = aes(x = nboots, y = R_npg), color = "purple") +
            geom_line(mapping = aes(x = nboots, y = true), linetype = "dashed", color = "chocolate3") + 
            labs(x = "Bootstrap Iterations", y = expression(paste("Power (1 - ", beta, ")")))
        
        l.p[[i + 3]] = ggplot(data = dat.4[i + c(0:4),])+
            geom_line(mapping = aes(x = nboots, y = R_np), color = "black") + 
            geom_line(mapping = aes(x = nboots, y = R_w), color = "blue") +
            geom_line(mapping = aes(x = nboots, y = C_np), color = "red") + 
            geom_line(mapping = aes(x = nboots, y = C_w), color = "green") +
            # geom_line(mapping = aes(x = nboots, y = C_npg), color = "orange") + 
            # geom_line(mapping = aes(x = nboots, y = R_npg), color = "purple") +
            geom_line(mapping = aes(x = nboots, y = true), linetype = "dashed", color = "chocolate3") + 
            labs(x = "Bootstrap Iterations", y = expression(paste("Power (1 - ", beta, ")")))
        
        l.p[[i + 4]] = ggplot(data = dat.5[i + c(0:4),])+
            geom_line(mapping = aes(x = nboots, y = R_np), color = "black") + 
            geom_line(mapping = aes(x = nboots, y = R_w), color = "blue") +
            geom_line(mapping = aes(x = nboots, y = C_np), color = "red") + 
            geom_line(mapping = aes(x = nboots, y = C_w), color = "green") +
            # geom_line(mapping = aes(x = nboots, y = C_npg), color = "orange") + 
            # geom_line(mapping = aes(x = nboots, y = R_npg), color = "purple") +
            geom_line(mapping = aes(x = nboots, y = true), linetype = "dashed", color = "chocolate3") + 
            labs(x = "Bootstrap Iterations", y = expression(paste("Power (1 - ", beta, ")")))
        
    }
    
    type2_plot = do.call("grid.arrange", c(l.p, nrow = 5, ncol = 5))
    ggsave("./results/plot_t2s2.pdf", plot = type2_plot, dpi = 300, width = 16, height = 9)
    
    save(file = "./results/ts2_data.RData", type1_plot, dat.R.np, dat.C.np, dat.R.npg, dat.C.npg, dat.R.w, dat.C.w, time.r, time.c, dat.1, dat.2, dat.3, dat.4, dat.5, type2_plot)
    
######### Linear Regression
    
    #### Specify a true data generating model (to assess accuracy of estimates)
    set.seed(314)
    beta = runif(4, -20, 20)
    prec = matrix(NA, ncol = 4, nrow = nrow(base_grid))
    colnames(prec) = c("MSE_np", "MSE_w", "MAE_np", "MAE_w")
    dat.Cpp  = dat.R = cbind(base_grid, prec)
    
    time.r = time.c = matrix(NA, nrow = nrow(dat.R), ncol = 6)
    
    for(i in 1:nrow(base_grid)){
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
    colnames(time.r) = colnames(time.c) = c(colnames(base_grid), names(summary(rnorm(10))))
    
    ## Exemplary accuracy
    i = 14
    x = cbind(sample(c(30:100), size = dat.R$n[i], replace = T), runif(dat.R$n[i], 50, 230), log(rnorm(dat.R$n[i], 100, 7.5)), rnorm(dat.R$n[i], 100, 7.5))
    y = x%*%beta + rnorm(dat.R$n[i])
    temp.r = replicate(500, reg.boot(y, x, nboot = dat.R$nboots[i], intercept = F)$estimate)
    
    temp.c = replicate(500, reg.boot(y, x, nboot = dat.R$nboots[i], intercept = F)$estimate)
    
    temp.r.w = replicate(500, reg.boot(y, x, nboot = dat.R$nboots[i], intercept = F, boot.type = "wild")$estimate)
    
    temp.c.w = replicate(500, reg.boot(y, x, nboot = dat.R$nboots[i], intercept = F, boot.type = "wild")$estimate)
    
    err.r = c(mse.fun.comp(temp.r, beta), mae.fun.comp(temp.r, beta))
    err.c = c(mse.fun.comp(temp.c, beta), mae.fun.comp(temp.c, beta))
    err.r.w = c(mse.fun.comp(temp.r.w, beta), mae.fun.comp(temp.r.w, beta))
    err.c.w = c(mse.fun.comp(temp.c.w, beta), mae.fun.comp(temp.c.w, beta))
    
    
    errors = rbind(err.r, err.r.w, err.c, err.c.w)
    colnames(errors) = c("Compound MSE", "Compound MAE")
    rownames(errors) = c("R - NP", "R - W", "C++ - NP", "C++ - W")
    
    beta.names = paste0("b", c(1:length(beta)))
    rownames(temp.r) = rownames(temp.c) = rownames(temp.r.w) = rownames(temp.c.w) = beta.names
    
    temp.final = cbind(temp.r, temp.c, temp.r.w, temp.c.w)
    temp.final = as.data.frame(t(temp.final))
    temp.final$id = as.factor(rep(c("R - NP", "C++ - NP", "R - W", "C++ - W"), each = 500))
    #colnames(temp.final) = c(beta.names, "id")
    
    save(file = "./results/reg_data.RData", time.r, time.c, temp.r, temp.c, temp.r.w, temp.c.w, errors)
    
    l.p = list()
    #for(i in seq(1, nrow(dat.R), 5)){
    for(i in 1:5){
        
        l.p[[i ]] = ggplot()+
            geom_line(mapping = aes(x = dat.R[i + c(0:4), "nboots"], y = log(time.r[i + c(0:4), 3])), color = "red")+
            geom_line(mapping = aes(x = dat.Cpp[i + c(0:4), "nboots"], y = log(time.c[i + c(0:4), 3])), color = "blue") +
            labs(x = "Bootstrap Iterations", y = "Log-Microseconds") # on Median
        
    }
    
    p1 = ggplot(temp.final, aes(id, b1)) +
        geom_boxplot() + 
        xlab("Type") +
        ylab(expression(hat(beta[1])))
    
    p2 = ggplot(temp.final, aes(id, b2)) +
        geom_boxplot() + 
        xlab("Type") +
        ylab(expression(hat(beta[2])))
    
    p3 = ggplot(temp.final, aes(id, b3)) +
        geom_boxplot() + 
        xlab("Type") +
        ylab(expression(hat(beta[3])))
    
    p4 = ggplot(temp.final, aes(id, b4)) +
        geom_boxplot() + 
        xlab("Type") +
        ylab(expression(hat(beta[4])))
    
    
    reg_plot = do.call("grid.arrange", c(l.p, nrow = 5))
    ggsave("./results/plot_regression.pdf", plot = reg_plot, dpi = 300, width = 16, height = 9)
    
    reg_err_plot = grid.arrange(p1, p2, p3, p4, nrow = 4)
    ggsave("./results/plot_regression_error.pdf", plot = reg_err_plot, dpi = 300, width = 16, height = 9)
    
    
    reg_xtable = xtable(errors)
    print(reg_xtable, file = "./results/tab_reg.tex")
    
    reg_time_xtab_R = xtable(time.r)
    reg_time_xtab_C = xtable(time.c)
    
    print(reg_time_xtab_R, file = "./results/tab_regTime_r.tex")
    print(reg_time_xtab_C, file = "./results/tab_regTime_Cpp.tex")
end.time = Sys.time()
time.diff.total = end.time - start.time
    