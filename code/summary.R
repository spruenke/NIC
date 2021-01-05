# # # # # # SUMMARY STATISTICS # # # # # #

  summary.boot = function(x, nboot = 100, boot.type = "np"){
    n = length(x)
    results.summary = numeric(7) # min, l-quartile, median, mean, u-quartile, max, sd
    names(results.summary) = c("Min", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max", "Std.Dev.")
    

    switch(boot.type,
             np = {
                x.boot  = matrix(sample(x, size = n*nboot, replace = T), ncol = nboot)
                
                x.quant = apply(x.boot, 2, quantile, probs = c(0.25, 0.5, 0.75))
                x.mean  = colMeans(x.boot)
                x.var   = (colSums(x.boot^2) - n*x.mean^2)/(n-1)
                
                results.summary[c(2,3,5)] = rowMeans(x.quant)
                results.summary[4]        = mean(x.mean)
                results.summary[7]        = sqrt(mean(x.var))
                results.summary[c(1,6)]   = c(min(x), max(x))
                
             },
             wild = {
                z       = x - mean(x)
                w.boot  = matrix(sample(c(-1,1), size = n*nboot, replace = T), ncol = nboot)
                x.boot  = z * w
                
                x.quant = apply(x.boot, 2, quantile, probs = c(0.25, 0.5, 0.75))
                x.mean  = colMeans(x.boot)
                x.var   = (colSums(x.boot^2) - n*x.mean^2)/(n-1)
                
                results.summary[c(2,3,5)] = rowMeans(x.quant)
                results.summary[4]        = mean(x.mean)
                results.summary[7]        = sqrt(mean(x.var))
                results.summary[c(1,6)]   = c(min(x), max(x))
               
             }
             )
        results.summary = round(results.summary, 4)
        return(results.summary)
  }
  
  summary.boot.c = function(x, nboot = 100, boot.type = "np"){
    results.summary = numeric(7) # min, l-quartile, median, mean, u-quartile, max, sd
    names(results.summary) = c("Min", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max", "Std.Dev.")
    results.summary[c(1,6)] = c(min(x), max(x))
    
    switch(boot.type,
          np = {
                results.summary[c(2,3,4,5,7)] = summaryBoot(x, nboot)
          },
          wild = {
            print("hi")
          }
    )
    results.summary = round(results.summary, 4)
    return(results.summary)
  }
  