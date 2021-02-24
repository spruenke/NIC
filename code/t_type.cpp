#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector tBootNp(NumericVector x, const int nboot) {
  int n = x.length();
  double xM = mean(x);
  NumericVector t_stat(nboot);
  for(int i = 0; i < nboot; i++){
    NumericVector x_boot = sample(x, n, true);
    double xMean = mean(x_boot);
    double xSq = 0;
    for(int j = 0; j < n; j++){
        xSq += std::pow(x_boot[j] - xMean, 2.0);
    }
    double xVar = xSq / (n-1);
    t_stat[i] = std::sqrt(n) * (xMean - xM) / std::sqrt(xVar);
  }
  
  return(t_stat);
}

// [[Rcpp::export]]
NumericVector tBootW(NumericVector x, const int nboot){
    int n = x.length();
    NumericVector ind(2);
    ind[0] = -1;
    ind[1] = 1;
    double xM = mean(x);
    NumericVector t_stat(nboot);
    for(int i = 0; i < nboot; i++){
        NumericVector w = sample(ind, n, true);
        NumericVector z = x - xM;
        NumericVector xBoot = z * w;
        double xMean = mean(xBoot);
        double xSq = 0;
        for(int j = 0; j < n; j++){
            xSq += std::pow(xBoot[j] - xMean, 2.0);
        }
        double xVar = xSq / (n-1);
        t_stat[i] = std::sqrt(n) * xMean / std::sqrt(xVar);
    }
      return(t_stat);
}

// [[Rcpp::export]]
NumericVector tBoot2Npg(const NumericVector x, const NumericVector y, const int nboot){
  const int n_1 = x.length();
  const int n_2 = y.length();
  NumericVector t_stat(nboot);
  for(int i = 0; i < nboot; i++){
      NumericVector x_boot = sample(x, n_1, true);
      NumericVector y_boot = sample(y, n_2, true);
      const double x_mean = mean(x_boot);
      const double y_mean = mean(y_boot);
      double xSq = 0;
      double ySq = 0;
      for(int j = 0; j < n_1; j++){
          xSq += std::pow(x_boot[j] - x_mean, 2.0);
      }
      for(int a = 0; a < n_2; a++){
          ySq += std::pow(y_boot[a] - y_mean, 2.0);
      }
      const double xVar = xSq / (n_1 - 1);
      const double yVar = ySq / (n_2 - 1);
      t_stat[i] = (x_mean - y_mean) / (std::sqrt((xVar / n_1) + (yVar / n_2)));
  }
    return(t_stat);
}

// [[Rcpp::export]]
NumericVector tBoot2Np(const NumericVector x, const NumericVector y, const int nboot){
  const int n_1 = x.length();
  const int n_2 = y.length();
  NumericVector t_stat(nboot);
  for(int i = 0; i < nboot; i++){
    NumericVector xy(n_1 + n_2);
    
    std::copy(x.begin(), x.end(), xy.begin());
    std::copy(y.begin(), y.end(), xy.begin() + x.size());
    IntegerVector seq_1 = seq(0, (n_1 - 1));
    IntegerVector seq_2 = seq(n_1, (n_1 + n_2 - 1));
    NumericVector xy_boot = sample(xy, (n_1 + n_2), true);
    NumericVector x_boot = xy_boot[seq_1];
    NumericVector y_boot = xy_boot[seq_2];
    
    const double x_mean = mean(x_boot);
    const double y_mean = mean(y_boot);
    double xSq = 0;
    double ySq = 0;
    for(int j = 0; j < n_1; j++){
      xSq += std::pow(x_boot[j] - x_mean, 2.0);
    }
    for(int a = 0; a < n_2; a++){
      ySq += std::pow(y_boot[a] - y_mean, 2.0);
    }
    const double xVar = xSq / (n_1 - 1);
    const double yVar = ySq / (n_2 - 1);
    t_stat[i] = (x_mean - y_mean) / (std::sqrt((xVar / n_1) + (yVar / n_2)));
  }
  return(t_stat);
}

// [[Rcpp::export]]
NumericVector tBoot2W(NumericVector x, NumericVector y, const int nboot){
  const int n_1 = x.length();
  const int n_2 = y.length();
  const double xM = mean(x);
  const double yM = mean(y);
   NumericVector z_x = x - xM;
   NumericVector z_y = y - yM;
  NumericVector ind(2);
  ind[0] = -1;
  ind[1] = 1;
  NumericVector t_stat(nboot);
  const IntegerVector seq_1 = seq(0, (n_1 - 1));
  const IntegerVector seq_2 = seq(n_1, (n_1 + n_2 - 1));
  for(int i = 0; i < nboot; i++){
    NumericVector w = sample(ind, (n_1 + n_2), true);
    NumericVector w_1 = w[seq_1];
    NumericVector w_2 = w[seq_2];
    NumericVector xBoot = z_x * w_1;
    NumericVector yBoot = z_y * w_2;
    double xMean = mean(xBoot);
    double yMean = mean(yBoot);
    double xSq = 0;
    double ySq = 0;
    for(int j = 0; j < n_1; j++){
      xSq += std::pow(xBoot[j] - xMean, 2.0);
    }
    for(int a = 0; a < n_2; a++){
        ySq += std::pow(yBoot[a] - yMean, 2.0);
    }
    double xVar = xSq / (n_1-1);
    double yVar = ySq / (n_2 -1);
    t_stat[i] = (xMean - yMean) / (std::sqrt((xVar / n_1) + (yVar / n_2)));
  }
  return(t_stat);
}