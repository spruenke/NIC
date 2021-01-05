# t tests w.r.t. Toutenbourg

# # # # # # ONE SAMPLE t-TESTS # # # # # #

  # Not best practice
    t.one.1 = function(x, mu.0 = 0, alpha = 0.05, alternative = "two.sided", nboot = 100, wild = T, ){
        T.x    = c()
        t.true = (mean(x) - mu.0) / sqrt(var(x)) * sqrt(length(x))
        
        for(i in 1:nboot){
            x.boot = sample(x, size = length(x), replace = T)
            if(wild == T){
                x.boot = (x - mean(x)) * sample(c(-1,1), size = length(x), replace = T)
                t.boot = (mean(x.boot))/sqrt(var(x.boot))*sqrt(length(x))
            } else {
                t.boot = (mean(x.boot) - mean(x))/sqrt(var(x.boot))*sqrt(length(x))
            }
            T.x = c(T.x, t.boot)
        }
        
        crit = quantile(T.x, probs = c(alpha/2, alpha, 1-alpha, 1-alpha/2))
        
        if(alternative == "two.sided"){
            ret  = (t.true < crit[1] || t.true > crit[4])
            
        } else if(alternative == "greater") {
            ret  = t.true > crit[3]
            
        } else if(alternative == "less") {
            ret  = t.true < crit[2]
            
        }
        
        p.v  = mean(T.x > t.true)
        p.v = min(2 * p.v, 2 - 2 * p.v)
        return(c(ret, p.v))
    }
    
  # R best practice - Why: Vectorization
    t.one.2 = function(x, mu.0 = 0, alpha = 0.05, alternative = "two.sided", nboot = 100, boot.type = "np"){
        T.x = numeric(nboot)
        n = length(x)
        t.true = (mean(x) - mu.0) / sqrt(var(x)) * sqrt(n)
        
        switch(boot.type, # switch between nonparametric and wild bootstrap
               np = {
                 x.boot = matrix(sample(x, size = n*nboot, replace = T), ncol = nboot)
                 x.mean = colMeans(x.boot)
                 x.var  = (colSums(x.boot^2) - n*x.mean^2)/(n-1)
                 T.x    = sqrt(n) * (x.mean - mean(x)) / sqrt(x.var)
                 est    = "Nonparametric Bootstrap"
               },
               wild = {
                 z = x - mean(x)
                 w.boot = matrix(sample(c(-1,1), size = n*nboot, replace = T), ncol = nboot) # Rademacher Weights
                 x.boot = w.boot * z
                 x.mean = colMeans(x.boot)
                 x.var  = (colSums(x.boot^2) - n*x.mean^2)/(n-1)
                 T.x    = sqrt(n) * x.mean / sqrt(x.var)
                 est    = "Wild Bootstrap"
               }
          )
        
        crit   = quantile(T.x, probs = c(alpha/2, alpha, 1-alpha, 1-alpha/2)) # compute critical values for one- and twosided
        switch(alternative, # switch between alternative hypotheses and choose whether to reject or not reject
               two.sided = {
                  ret    = (t.true < crit[1] || t.true > crit[4])
                  
               },
               greater   = {
                  ret  = t.true > crit[3]
                  
               },
               less      = {
                  ret  = t.true < crit[2]
                  
               })
        
        # compute p-value
        p.v  = mean(T.x > t.true)
        p.v = min(2 * p.v, 2 - 2 * p.v)
        
        # Prepare Results #######################################
          names(mu.0) = "mean"
          pop.var = mean(x.mean)
          names(pop.var) = "Mean"
          pop.stat = t.true
          names(pop.stat) = "t"
          pop.par = length(x) - 1
          names(pop.par) = "df"
          conf.int = quantile(x.mean, probs = c(alpha/2, 1-alpha/2))
          attr(conf.int, "conf.level") = 1 - alpha
        ##########################################################
          
        t.list = list("reject" = ret, "null.value" = mu.0, "alternative" = alternative, "method" = "t.test", "estimate" = pop.var, "data.name" = deparse(substitute(x)), "statistic" = pop.stat, "parameters" = pop.par, "p.value" = p.v, "estimation.method" = est, "sample.size" = length(x), "conf.int" = conf.int)
        class(t.list) = "htest"
        return(t.list)
    }
    
    t.one.c = function(x, mu.0 = 0, alpha = 0.05, alternative = "two.sided", nboot = 100, boot.type = "np"){
          y   = matrix(sample(x, size = nboot * length(x), replace = T), ncol = nboot, nrow = length(x))
          t.c = .Call("t_stat", t(y), as.integer(length(x)), 
                    as.integer(nboot), mean(x))
          return(t.c)
    }
    AB = replicate(1000, t.one.2(rnorm(20))$estimate)
        
    
# # # # # # TWO SAMPLE t-TESTS # # # # # #
    
    
    t.two  = function(x.1, x.2, alpha = 0.05, alternative = "two.sided", nboot = 100, boot.type = "np"){
      T.x = numeric(nboot)
      n.1 = length(x.1)
      n.2 = length(x.2)
      n   = n.1 + n.2
      
      t.true = (mean(x.1) - mean(x.2)) / (sqrt((var(x.1)/n.1) + (var(x.2)/n.2))) # see Toutenburg p. 145
      
      switch(boot.type, # switch between (groupwise) nonparametric and wild bootstrap
             npg = { # groupwise nonparametric bootstrap
                 x.1.boot = matrix(sample(x.1, size = n.1*nboot, replace = T), ncol = nboot)
                 x.2.boot = matrix(sample(x.2, size = n.2*nboot, replace = T), ncol = nboot)
                 x.1.mean = colMeans(x.1.boot)
                 x.2.mean = colMeans(x.2.boot)
                 x.1.var  = (colSums(x.1.boot^2) - n.1*x.1.mean^2)/(n.1-1)
                 x.2.var  = (colSums(x.2.boot^2) - n.2*x.2.mean^2)/(n.2-1)
                 T.x      = (x.1.mean - x.2.mean) / (sqrt((x.1.var / n.1) + (x.2.var / n.2)))
                 est      = "Groupwise Nonparametric Bootstrap"
             },
             np  = { # Nonparametric Bootstrap
                 x.boot   = matrix(sample(c(x.1,x.2), size = n*nboot, replace = T), ncol = nboot)
                 x.1.boot = x.boot[1:n.1, ]
                 x.2.boot = x.boot[(n.1+1):n, ]
                 x.1.mean = colMeans(x.1.boot)
                 x.2.mean = colMeans(x.2.boot)
                 x.1.var  = (colSums(x.1.boot^2) - n.1*x.1.mean^2)/(n.1-1)
                 x.2.var  = (colSums(x.2.boot^2) - n.2*x.2.mean^2)/(n.2-1)
                 T.x      = (x.1.mean - x.2.mean) / (sqrt((x.1.var / n.1) + (x.2.var / n.2)))
                 est      = "Nonparametric Bootstrap"
                
             },
             wild = { # Wild (Rademacher) Bootstrap
                 z.1      = x.1 - mean(x.1)
                 z.2      = x.2 - mean(x.2)
                 w.boot   = matrix(sample(c(-1,1), size = n*nboot, replace = T), ncol = nboot)
                 x.1.boot = z.1 * w.boot[1:n.1, ]
                 x.2.boot = z.2 * w.boot[(n.1+1):n, ]
                 x.1.mean = colMeans(x.1.boot)
                 x.2.mean = colMeans(x.2.boot)
                 x.1.var  = (colSums(x.1.boot^2) - n.1*x.1.mean^2)/(n.1-1)
                 x.2.var  = (colSums(x.2.boot^2) - n.2*x.2.mean^2)/(n.2-1)
                 T.x      = (x.1.mean - x.2.mean) / (sqrt((x.1.var / n.1) + (x.2.var / n.2)))
                 est      = "Wild Bootstrap"
             }
      )
      
      
      crit   = quantile(T.x, probs = c(alpha/2, alpha, 1-alpha, 1-alpha/2)) # compute critical values for one- and twosided
      switch(alternative, # switch between alternative hypotheses and choose whether to reject or not reject
             two.sided = {
               ret    = (t.true < crit[1] || t.true > crit[4])
               
             },
             greater   = {
               ret  = t.true > crit[3]
               
             },
             less      = {
               ret  = t.true < crit[2]
               
             })
      
      # compute p-value
      p.v  = mean(T.x > t.true)
      p.v  = min(p.v, 1 - p.v)
      
      # Prepare Results #######################################
      dat.nam         = c(deparse(substitute(x.1)), deparse(substitute(x.2)))
      dat.nam.2       = paste(dat.nam[1], "and", dat.nam[2])
      h.0             = 0
      names(h.0)      = "difference in means"
      pop.var         = c(mean(x.1.mean), mean(x.2.mean))
      names(pop.var)  = paste("Mean of", dat.nam)
      pop.stat        = t.true
      names(pop.stat) = "t"
      pop.par         = n - 2
      names(pop.par)  = "df"
      conf.int        = quantile((x.1.mean - x.2.mean), probs = c(alpha/2, 1-alpha/2))
      
      attr(conf.int, "conf.level") = 1 - alpha
      ##########################################################
      
      t.list = list("reject" = ret ,"null.value" = h.0, "alternative" = alternative, "method" = "Welch Two Sample t-Test", "estimate" = pop.var, "data.name" = dat.nam.2, "statistic" = pop.stat, "parameters" = pop.par, "p.value" = p.v, "estimation.method" = est, "sample.size" = n, "conf.int" = conf.int)
      class(t.list) = "htest"
      return(t.list)
    }
    
    
    AB.2 = replicate(1000, t.two(rnorm(20,0,1), rnorm(20,0,3), boot.type = "np")$reject)
    mean(AB.2)
    