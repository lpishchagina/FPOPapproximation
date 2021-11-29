#include "FPOP.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

//' @title approx_fpop
//'
//' @description FPOP method (using the rectangle approximation of the sets) for  the multiple changepoint detection.
//' @param data is a matrix of data (p-rows x n-columns).
//' @param penalty is a value of penalty (a non-negative real number).
//' @param intersection is the type of intersection : 'empty', 'all', 'last', 'random' or 'sphere'.
//' @param exclusion is the type of intersection : 'empty', 'all', 'random' or 'sphere'.
//' The following parameter combinations are implemented:
//' (intersection ='sphere', exclusion ='sphere')
//' (intersection ='all', exclusion ='all')
//' (intersection ='all', exclusion ='empty')
//' (intersection ='empty', exclusion ='all')
//' (intersection ='last', exclusion ='all')
//' (intersection ='last', exclusion ='random')
//' (intersection ='all', exclusion ='random')
//' (intersection ='random', exclusion ='random')
//' (intersection ='empty', exclusion ='empty')
//' @param test_nb_cands is the logical parameter (if test_nb_cands = TRUE, than the file "Nb_cands.txt" contains the number of change candidates for each operation.
//' @param test_nb_exclus is the logical parameter (if test_nb_exclus = TRUE, than the file "Nb_cands.txt" contains the number of exclusion for change candidates for each operation. (!!Change!!)
//'
//' @return a list of  elements  = (changes, means, globalCost).
//'
//' \describe{
//' \item{\code{changes}}{is the changepoint vector that gives the last index of each segment for the p-variate time series.}
//' \item{\code{means}}{is the list of successive means for the p-variate time series.}
//' \item{\code{globalCost}}{is a number equal to the global cost.}
//' }
//'
//' @examples approx_fpop(data = chpt_rnorm(p = 3, n = 100, changes = 50, means = matrix(c (1,2,3,4, 5, 7), nrow = 3), noise = 1), penalty = 2*log(100), type_approx = 2)
//' N <- 10
//' Chpt <-5
//' Means <-  matrix(c(0,1,1,10), nrow = 2)
//' Noise <- 1
//' Dim <- 2
//' Penality <- 2*Dim*log(N)
//'time_series <- changes_rnorm(p = Dim, n = N, changes = Chpt, means = Means, noise = Noise)
//'Approx <- list()
//'Approx[[1]] <- approx_fpop(data = time_series, penalty = Penality, intersection = 'sphere', exclusion = 'sphere', test_nb_cands = FALSE, test_nb_exclus = FALSE)
//'Approx[[2]] <- approx_fpop(data = time_series, penalty = Penality, intersection = 'all', exclusion = 'all', test_nb_cands = FALSE, test_nb_exclus = FALSE)
//'Approx[[3]] <-approx_fpop(data = time_series, penalty = Penality, intersection = 'all', exclusion = 'empty', test_nb_cands = FALSE, test_nb_exclus = FALSE)
//'Approx[[4]] <-approx_fpop(data = time_series, penalty = Penality, intersection = 'empty', exclusion = 'all', test_nb_cands = FALSE, test_nb_exclus = FALSE)
//'Approx[[5]] <-approx_fpop(data = time_series, penalty = Penality, intersection = 'last', exclusion = 'all', test_nb_cands = FALSE, test_nb_exclus = FALSE)
//'Approx[[6]] <-approx_fpop(data = time_series, penalty = Penality, intersection = 'last', exclusion = 'random', test_nb_cands = FALSE, test_nb_exclus = FALSE)
//'Approx[[7]] <-approx_fpop(data = time_series, penalty = Penality, intersection = 'all', exclusion = 'random', test_nb_cands = FALSE, test_nb_exclus = FALSE)
//'Approx[[8]] <-approx_fpop(data = time_series, penalty = Penality, intersection = 'random', exclusion = 'random', test_nb_cands = FALSE, test_nb_exclus = FALSE)
//'Approx[[9]] <-approx_fpop(data = time_series, penalty = Penality, intersection = 'empty', exclusion = 'empty', test_nb_cands = FALSE, test_nb_exclus = FALSE)
//'Approx
// [[Rcpp::export]]
List approx_fpop(Rcpp::NumericMatrix data, double penalty, std::string intersection = "all",  std::string exclusion = "all", bool test_nb_cands = false, bool test_nb_exclus = false) {
  int type_approx = 0;
  // 2 = (Iall, Eall), 3 = (Iall,Eempty), 7 = (Iall,Erandom)
  if (intersection == "all" ) {
    if (exclusion == "all") {
      type_approx = 2;
    } else if (exclusion == "empty") {
      type_approx = 3;
      if (test_nb_exclus) {
        test_nb_exclus = false;
        Rcpp::Rcout << "The test of the exclusion number  is not available for the combination ('all', 'empty')."<< endl;
      }
    } else if (exclusion == "random") {
      type_approx = 7;
      if (test_nb_exclus) {
        test_nb_exclus = false;
        Rcpp::Rcout << "The exclusion number for each change candidate is 1 for the combination ('all', 'random')."<< endl;
      }
    }
  }
  //4 = (Iempty, Eall), 9 = (Iempty,Eempty) //<=PELT//
  if (intersection == "empty" ) {
    if (exclusion == "all") {
      type_approx = 4;
    } else if (exclusion == "empty") {
      type_approx = 9;
      if (test_nb_exclus) {
        test_nb_exclus = false;
        Rcpp::Rcout << "The test of the exclusion number  is not available for the combination ('empty', 'empty')."<< endl;
      }
    }
  }
  //5 = (Ilast,Eall), 6 = (Ilast, Erandom)
  if (intersection == "last" ) {
    if (exclusion == "all") {
      type_approx = 5;
    } else if (exclusion == "random") {
      type_approx = 6;
      if (test_nb_exclus) {
        test_nb_exclus = false;
        Rcpp::Rcout << "The exclusion number for each change candidate  is 1 for the combination ('last', 'random')."<< endl;
      }
    }
  }
  //8 = (Irandom, Erandom)
  if ((intersection == "random") && (exclusion == "random"))  {
    type_approx = 8;
    if (test_nb_exclus) {
      test_nb_exclus = false;
      Rcpp::Rcout << "The exclusion number for each change candidate  is 1 for the combination ('random', 'random')."<< endl;
    }
  }
  if ((intersection == "sphere") && (exclusion == "sphere"))  {
    type_approx = 1;
  }
  //----------stop--------------------------------------------------------------
  if (penalty < 0) {throw std::range_error("Penalty should be a non-negative number!");}
  if(type_approx < 1 || type_approx > 9)
  {throw std::range_error("This combination of parameters 'intersection' and 'exclusion' is not available. ");}
  //----------------------------------------------------------------------------
  List res;
  bool test;
  test = false;
  if (type_approx == 1) {
    //  test = true;
    FPOP<Candidate_sphere_sphere_1> X = FPOP<Candidate_sphere_sphere_1>(data, penalty);
    X.algoFPOP(data, type_approx, test_nb_cands, test_nb_exclus);
    res["changes"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["globalCost"] = X.GetGlobalCost();
  }
  if (type_approx == 2) {
    // test = true;
    FPOP<Candidate_Iall_Eall_2> X = FPOP<Candidate_Iall_Eall_2>(data, penalty);
    X.algoFPOP(data, type_approx, test_nb_cands, test_nb_exclus);
    res["changes"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["globalCost"] = X.GetGlobalCost();
  }
  if (type_approx == 3) {
    // test = true;
    FPOP<Candidate_Iall_Eempty_3> X = FPOP<Candidate_Iall_Eempty_3>(data, penalty);
    X.algoFPOP(data, type_approx, test_nb_cands, test_nb_exclus);
    res["changes"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["globalCost"] = X.GetGlobalCost();
  }
  if (type_approx == 4) {
    // test = true;
    FPOP<Candidate_Iempty_Eall_4> X = FPOP<Candidate_Iempty_Eall_4>(data, penalty);
    X.algoFPOP(data, type_approx, test_nb_cands, test_nb_exclus);
    res["changes"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["globalCost"] = X.GetGlobalCost();
  }
  if (type_approx == 5) {
    // test = true;
    FPOP<Candidate_Ilast_Eall_5> X = FPOP<Candidate_Ilast_Eall_5>(data, penalty);
    X.algoFPOP(data, type_approx, test_nb_cands, test_nb_exclus);
    res["changes"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["globalCost"] = X.GetGlobalCost();
  }
  if (type_approx == 6) {
    // test = true;
    FPOP<Candidate_Ilast_Erandom_6> X = FPOP<Candidate_Ilast_Erandom_6>(data, penalty);
    X.algoFPOP(data, type_approx, test_nb_cands, test_nb_exclus);
    res["changes"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["globalCost"] = X.GetGlobalCost();
  }
  if (type_approx == 7) {
    // test = true;
    FPOP<Candidate_Iall_Erandom_7> X = FPOP<Candidate_Iall_Erandom_7>(data, penalty);
    X.algoFPOP(data, type_approx, test_nb_cands, test_nb_exclus);
    res["changes"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["globalCost"] = X.GetGlobalCost();
  }
  if (type_approx == 8) {
    // test = true;
    FPOP<Candidate_Irandom_Erandom_8> X = FPOP<Candidate_Irandom_Erandom_8>(data, penalty);
    X.algoFPOP(data, type_approx,test_nb_cands, test_nb_exclus);
    res["changes"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["globalCost"] = X.GetGlobalCost();
  }
  if (type_approx == 9) {
    // test = true;
    FPOP<Candidate_Iempty_Eempty_9> X = FPOP<Candidate_Iempty_Eempty_9>(data, penalty);
    X.algoFPOP(data, type_approx, test_nb_cands, test_nb_exclus);
    res["changes"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["globalCost"] = X.GetGlobalCost();
  }
  return res;
}
