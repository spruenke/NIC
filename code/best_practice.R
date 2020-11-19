## Note: The following two functions are only for demonstration purposes; Of course, the actual functions would contain more code
##        and would not return mean(T.x);
##        However, the body of these two is the part of the t-Tests which requires most computation!

best.practice = function(x, mu.0, nboot){
  T.x = numeric(nboot)
  n = length(x)
  t.true = (mean(x) - mu.0) / sqrt(var(x)) * sqrt(n)
  
  x.boot = matrix(sample(x, size = n*nboot, replace = T), ncol = nboot)
  x.mean = colMeans(x.boot)
  x.var  = (colSums(x.boot^2) - n*x.mean^2)/(n-1)
  T.x    = sqrt(n) * (x.mean - mean(x)) / sqrt(x.var)
  
  return(mean(T.x))
}

not.best.practice = function(x, mu.0, nboot){
  T.x    = c()
  t.true = (mean(x) - mu.0) / sqrt(var(x)) * sqrt(length(x))
  
  for(i in 1:nboot){
    x.boot = sample(x, size = length(x), replace = T)
    t.boot = (mean(x.boot) - mean(x))/sqrt(var(x.boot))*sqrt(length(x))
    T.x = c(T.x, t.boot)
  }
  return(mean(T.x))
}

n.boot = c(20, 50, 100, 200, 500, 1000, 2000, 5000, 10000)
comp.time = list()
z      = c(10, 20, 50, 100)
mu.0   = 0

for(k in 1:length(z)){
    t.not  = numeric(length(n.boot))
    t.bes  = numeric(length(n.boot))
    x      = rnorm(z[k])
    for(j in 1:length(n.boot)){
        t.st.1 = Sys.time()
        not.best.practice(x, mu.0, n.boot[j])
        t.en.1 = Sys.time()
        
        t.st.2 = Sys.time()
        best.practice(x, mu.0, n.boot[j])
        t.en.2 = Sys.time()
        
        t.not[j] = t.en.1 - t.st.1
        t.bes[j] = t.en.2 - t.st.2
        
    }
    comp.time[[k]] = cbind(t.not, t.bes)
}

pdf(file = "time1.pdf", width = 16, height = 9)
plot(x = n.boot, y = comp.time[[1]][,1], type = "b", lwd = 2, col = "red", ylab = "Time (in sec)", xlab = "Bootstrap Iterations", main = "t-Test (Nonparametric Bootstrap) for Sample Size n=10")
lines(x = n.boot, y = comp.time[[1]][,2], type = "b", lwd = 2, col = "darkgreen")
dev.off()