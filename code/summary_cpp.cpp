
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector summaryBoot(NumericVector x, const int nboot) {
    int n = x.length();
    NumericVector summaryBoot(5);
    IntegerVector v = seq(0, (n-1));
    NumericVector lQ(nboot);
    NumericVector med(nboot);
    NumericVector men(nboot);
    NumericVector uQ(nboot);
    NumericVector var(nboot);
    for(int i = 0; i < nboot; i++){
        NumericVector xBoot = Rcpp::sample(x, n, true);
        xBoot.sort();
        
        double ind_1 = 1 + (n - 1) * 0.25;
        int lo_1 = floor(ind_1);
        int hi_1 = ceil(ind_1);
        double qs_1 = xBoot[lo_1 - 1];
        double h_1 = ind_1 - lo_1;
        lQ[i] = ((1 - h_1) * qs_1) + (h_1 * xBoot[hi_1 - 1]);
        
        double ind_2 = 1 + (n - 1) * 0.5;
        int lo_2 = floor(ind_2);
        int hi_2 = ceil(ind_2);
        double qs_2 = xBoot[lo_2 - 1];
        double h_2 = ind_2 - lo_2;
        med[i] = ((1 - h_2) * qs_2) + (h_2 * xBoot[hi_2-1]);
        
        double ind = 1 + (n - 1) * 0.75;
        int lo = floor(ind);
        int hi = ceil(ind);
        double qs = xBoot[lo-1];
        double h = ind - lo;
        uQ[i] = ((1 - h) * qs) + (h * xBoot[hi-1]);
        
        men[i] = mean(xBoot);
        
        double xSq = 0;
        for(int j = 0; j < n; j++){
            xSq += pow(xBoot[j] - mean(xBoot), 2.0);
        }
        var[i] = std::sqrt(xSq / (n-1));
    }
    
    summaryBoot[0] = mean(lQ);
    summaryBoot[1] = mean(med);
    summaryBoot[2] = mean(men);
    summaryBoot[3] = mean(uQ);
    summaryBoot[4] = mean(var);
    
    return summaryBoot;
}


// [[Rcpp::export]]
NumericVector summaryBootWild(NumericVector x, const int nboot) {
    int n = x.length();
    NumericVector summaryBoot(5);
    NumericVector lQ(nboot);
    NumericVector med(nboot);
    NumericVector men(nboot);
    NumericVector uQ(nboot);
    NumericVector var(nboot);
    double xMean = mean(x);
    NumericVector xSample = x - xMean;
    NumericVector w(2);
    w[0] = -1;
    w[1] = 1;
    for(int i = 0; i < nboot; i++){
        NumericVector z = Rcpp::sample(w, n, true);
        NumericVector xBoot = z * xSample;
        xBoot.sort();
        
        double ind_1 = 1 + (n - 1) * 0.25;
        int lo_1 = floor(ind_1);
        int hi_1 = ceil(ind_1);
        double qs_1 = xBoot[lo_1 - 1];
        double h_1 = ind_1 - lo_1;
        lQ[i] = ((1 - h_1) * qs_1) + (h_1 * xBoot[hi_1 - 1]);
        
        double ind_2 = 1 + (n - 1) * 0.5;
        int lo_2 = floor(ind_2);
        int hi_2 = ceil(ind_2);
        double qs_2 = xBoot[lo_2 - 1];
        double h_2 = ind_2 - lo_2;
        med[i] = ((1 - h_2) * qs_2) + (h_2 * xBoot[hi_2-1]);
        
        double ind = 1 + (n - 1) * 0.75;
        int lo = floor(ind);
        int hi = ceil(ind);
        double qs = xBoot[lo-1];
        double h = ind - lo;
        uQ[i] = ((1 - h) * qs) + (h * xBoot[hi-1]);
        
        men[i] = mean(xBoot);
        
        double xSq = 0;
        for(int j = 0; j < n; j++){
            xSq += pow(xBoot[j] - mean(xBoot), 2.0);
        }
        var[i] = std::sqrt(xSq / (n-1));
    }
    
    summaryBoot[0] = mean(lQ) + xMean;
    summaryBoot[1] = mean(med) + xMean;
    summaryBoot[2] = mean(men) + xMean;
    summaryBoot[3] = mean(uQ) + xMean;
    summaryBoot[4] = mean(var);
    
    return summaryBoot;
}