#include "FPOP.h"

#include "math.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

//' @title approx_fpop
//'
//' @description Detection changepoints using the Functional Pruning Optimal Partitioning method (FPOP) in p-variate time series in a p-variable time series of length n.
//'
//' @param data is a matrix of data(p-rows x n-columns).
//' @param penalty is a value of penalty (a non-negative real number).
//' @param type_approx is a value defining the  type of geometry for FPOP-pruning: type=1: ("intersection" of sets), approximation - rectangle; type=2:("intersection" of sets)"minus"("union" of sets), approximation - rectangle; type=3: (last disk)"minus"("union" of sets), approximation - disk.
//'
//' @return a list of  elements  = (changepoints, means, globalCost).
//'
//' \describe{
//' \item{\code{chpts}}{is the changepoint vector that gives the last index of each segment for the p-variate time series.}
//' \item{\code{means}}{is the list of successive means for the p-variate time series.}
//' \item{\code{globalCost}}{is a number equal to the global cost.}
//' }
//'
//' @examples approx_fpop(data = chpt_rnorm(p = 3, n = 100, chpts = 50, means = matrix(c (1,2,3,4, 5, 7), nrow = 3), noise = 1), penalty = 2*log(100), type_approx = 2)
//' data2 =  chpt_rnorm(p = 2, n = 20, chpts = 10, means = matrix(c(0,1,1,10), nrow = 2), noise = 1)
//'approx_fpop(data = data2, penalty = 2*log(20), type_approx = 2)
//'approx_fpop(data = chpt_rnorm(p = 2, n = 20, chpts = 50, means = matrix(c (1,2,7,9), nrow = 2), noise = 1), penalty = 2*log(100), type_approx = 2)
// [[Rcpp::export]]
List approx_fpop(Rcpp::NumericMatrix data, double penalty, int type_approx) {
  //----------stop--------------------------------------------------------------
  if(penalty < 0) {throw std::range_error("penalty should be a non-negative number");}
  if(type_approx < 1 || type_approx > 8)
  {throw std::range_error("type_approx must be one of: 1, 2, 3, 4, 5, 6, 7 or 8.");}
  //----------------------------------------------------------------------------
  List res;
  bool test;
  test = false;
  if (type_approx == 1){
    //  test = true;
    FPOP<Candidate_sphere_sphere_1> X = FPOP<Candidate_sphere_sphere_1>(data, penalty);
    X.algoFPOP(data, type_approx, test);
    res["chpts"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["globalCost"] = X.GetGlobalCost();
  }
  if (type_approx == 2){
    // test = true;
    FPOP<Candidate_Iall_Eall_2> X = FPOP<Candidate_Iall_Eall_2>(data, penalty);
    X.algoFPOP(data, type_approx, test);
    res["chpts"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["globalCost"] = X.GetGlobalCost();
  }
  if (type_approx == 3){
    // test = true;
    FPOP<Candidate_Iall_Eempty_3> X = FPOP<Candidate_Iall_Eempty_3>(data, penalty);
    X.algoFPOP(data, type_approx, test);
    res["chpts"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["globalCost"] = X.GetGlobalCost();
    }
  if (type_approx == 4){
    // test = true;
    FPOP<Candidate_Iempty_Eall_4> X = FPOP<Candidate_Iempty_Eall_4>(data, penalty);
    X.algoFPOP(data, type_approx, test);
    res["chpts"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["globalCost"] = X.GetGlobalCost();
  }
  if (type_approx == 5){
    // test = true;
    FPOP<Candidate_Ilast_Eall_5> X = FPOP<Candidate_Ilast_Eall_5>(data, penalty);
    X.algoFPOP(data, type_approx, test);
    res["chpts"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["globalCost"] = X.GetGlobalCost();
  }
  if (type_approx == 7){
    // test = true;
    FPOP<Candidate_Iall_Erandom_7> X = FPOP<Candidate_Iall_Erandom_7>(data, penalty);
    X.algoFPOP(data, type_approx, test);
    res["chpts"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["globalCost"] = X.GetGlobalCost();
  }
  return res;
}
