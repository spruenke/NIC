# # # # # # REGRESSION # # # # # #
sourceCpp("rcpp_beta.cpp")
# see Fox and Weisberg (2019)
reg.boot = function(y, x, nboot = 100, boot.type = "np", intercept = T){
  if(length(names(Filter(is.factor, x))) != 0){warning("There are factors in your data! Please convert them to appropriate dummies before running this function!")}

  var.names = colnames(x)
    if(intercept == T){
        x = cbind(rep(1, nrow(x)), x)
        var.names = c("Intercept", var.names)
    }
    beta.boot = matrix(NA, ncol = nboot, nrow = ncol(x))
    data = cbind(y,x)
    data = as.matrix(data)

    beta.true           = as.vector(solve(t(data[,-1])%*%data[,-1])%*%t(data[,-1])%*%data[,1])
    names(beta.true)    = var.names
    y.hat               = data[,-1]%*%beta.true

    if(boot.type == "np"){
        for(i in 1:nboot){ # nice, if X is random
            index.boot    = sample(c(1:nrow(data)), size = nrow(data), replace = T)
            data.boot     = data[index.boot, ]
            y.boot        = data.boot[, 1]
            x.boot        = data.boot[, -1]
            beta.boot[,i] = solve(t(x.boot)%*%x.boot)%*%t(x.boot)%*%y.boot
        }
    }
    if(boot.type == "wild"){ # if X is deterministic/fixed

        resid.hat   = y - y.hat
        for(i in 1:nboot){
            w          = sample(c(-1,1), replace = T, size = length(y))
            resid.boot = w * resid.hat#sample(resid.fit, size = length(y), replace = T)
            y.boot     = y.hat + resid.boot
            beta.boot[,i] = solve(t(data[,-1])%*%data[,-1])%*%t(data[,-1])%*%y.boot
        }
    }

    beta                = rowMeans(beta.boot)
    names(beta)         = var.names
    beta.var            = (rowSums(beta.boot^2) - nboot * rowMeans(beta.boot)^2) / (nboot - 1)
    beta.sd             = sqrt(beta.var)
    names(beta.sd)      = var.names

    conf.int            = apply(beta.boot, 1, quantile, probs = c(0.025, 0.975))
    colnames(conf.int)  = var.names
    return(list("true"  = beta.true, "estimate" = beta, "conf.int" = conf.int, "Std.D." = beta.sd))

}

reg.boot.c = function(y, x, nboot = 100, boot.type = "np", intercept = T){
  if(length(names(Filter(is.factor, x))) != 0){warning("There are factors in your data! Please convert them to appropriate dummies before running this function!")}
  
  var.names = colnames(x)
  if(intercept == T){
    x = cbind(rep(1, nrow(x)), x)
    var.names = c("Intercept", var.names)
  }
  
  data = cbind(y,x)
  data = as.matrix(data)
  
switch(boot.type,
       np = {
          beta.hat = betaBoot(y, data[,-1], nboot)
       },
       wild = {
          beta.hat = BetaBootWild(y, data[,-1], nboot)
       }
       )
  
  beta.true           = as.vector(solve(t(data[,-1])%*%data[,-1])%*%t(data[,-1])%*%data[,1])
  names(beta.true)    = var.names
  y.hat               = data[,-1]%*%beta.true
  
  beta                = rowMeans(beta.boot)
  names(beta)         = var.names
  beta.var            = (rowSums(beta.boot^2) - nboot * rowMeans(beta.boot)^2) / (nboot - 1)
  beta.sd             = sqrt(beta.var)
  names(beta.sd)      = var.names
  
  conf.int            = apply(beta.boot, 1, quantile, probs = c(0.025, 0.975))
  colnames(conf.int)  = var.names
  return(list("true"  = beta.true, "estimate" = beta, "conf.int" = conf.int, "Std.D." = beta.sd))
}

microbenchmark(reg.boot(y.test, x.test, nboot = 1000, intercept = F, boot.type = "wild"),
               BetaBootWild(y.test, x.test, 1000))

microbenchmark(reg.boot(y.test, x.test, nboot = 1000, intercept = F),
               betaBoot(y.test, x.test, 1000))

A = matrix(rnorm(10000), nrow = 1000)