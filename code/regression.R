# # # # # # REGRESSION # # # # # #

reg.boot = function(y, x, nboot = 100, boot.type = "np", intercept = T){
    var.names = colnames(x)
    if(intercept == T){
        x = cbind(rep(1, nrow(x)), x)
        var.names = c("Intercept", var.names)
    }
    beta.boot = matrix(NA, ncol = nboot, nrow = ncol(x))
    data = cbind(y,x)
    data = as.matrix(data)
    

    for(i in 1:nboot){
        index.boot    = sample(c(1:nrow(data)), size = nrow(data), replace = T)
        data.boot     = data[index.boot, ]
        y.boot        = data.boot[, 1]
        x.boot        = data.boot[, -1]
        beta.boot[,i] = solve(t(x.boot)%*%x.boot)%*%t(x.boot)%*%y.boot
    }
    
    
    beta                = rowMeans(beta.boot)
    names(beta)         = var.names
    beta.var            = (rowSums(beta.boot^2) - nboot * rowMeans(beta.boot)^2) / (nboot - 1)
    beta.sd             = sqrt(beta.var)
    names(beta.sd)      = var.names
    beta.true           = as.vector(solve(t(data[,-1])%*%data[,-1])%*%t(data[,-1])%*%data[,1])
    names(beta.true)    = var.names
    conf.int            = apply(beta.boot, 1, quantile, probs = c(0.025, 0.975))
    colnames(conf.int)  = var.names
    return(list("true"  = beta.true, "estimate" = beta, "conf.int" = conf.int, "Std.D." = beta.sd))
    
}
