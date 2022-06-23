#include "FPOP.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

//converting parameters("approximation", "intersection", "exclusion") to a numeric value.
int NmbOfapproxFPOP(std::string approximation, std::string intersection,  std::string exclusion) {
  int type_approx = 0;
  if (approximation == "rectangle") {
    if ((intersection == "empty") && (exclusion == "empty")) { type_approx = 1; }
    if (intersection == "all" ) {
      if (exclusion == "empty") { type_approx = 10; }
      if (exclusion == "all") { type_approx = 11; }
      if (exclusion == "random") { type_approx = 21; }
    }
    if (intersection == "last" ) {
      if (exclusion == "all")  { type_approx = 13; }
      if (exclusion == "random")  { type_approx = 17; }

    }
    if (intersection == "random")  {
      if (exclusion == "random") { type_approx = 19; }
      if (exclusion == "all")  { type_approx = 15; }
    }
  }
  if (approximation == "sphere") {
    if  (intersection == "last")  {
      if (exclusion == "all")  { type_approx = 2; }
    }
  }
  return type_approx;
}

//Comparison of two FPOP-methods
bool TestOfComparisonTwoFPOP(Rcpp::NumericMatrix data, double penalty, unsigned int type_approx2, double UnpenalizedCost1, unsigned int* LastChpt1) {
  bool res = false;
  if (type_approx2 == 1) {
    FPOP<Rec_Empty_Empty> X = FPOP<Rec_Empty_Empty>(data, penalty);
    X.algoFPOP(data, type_approx2, false );
    res = X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else  if (type_approx2 == 2) {
    FPOP<Sph_Last_lAll> X = FPOP<Sph_Last_lAll>(data, penalty);
    X.algoFPOP(data, type_approx2, false );
    res = X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else if (type_approx2 == 10) {
    FPOP<Rec_All_Empty> X = FPOP<Rec_All_Empty>(data, penalty);
    X.algoFPOP(data, type_approx2, false );
    res = X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else if (type_approx2 == 11) {
    FPOP<Rec_All_lAll> X = FPOP<Rec_All_lAll>(data, penalty);
    X.algoFPOP(data, type_approx2, false );
    res =  X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else if (type_approx2 == 13) {
    FPOP<Rec_Last_lAll> X = FPOP<Rec_Last_lAll>(data, penalty);
    X.algoFPOP(data, type_approx2, false );
    res =  X.TestFPOP(UnpenalizedCost1, LastChpt1);
  }  else if (type_approx2 == 15) {
    FPOP<Rec_Rand_lAll> X = FPOP<Rec_Rand_lAll>(data, penalty);
    X.algoFPOP(data, type_approx2, false);
    res =  X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else if (type_approx2 == 17) {
    FPOP<Rec_Last_lRand> X = FPOP<Rec_Last_lRand>(data, penalty);
    X.algoFPOP(data, type_approx2, false);
    res =  X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else if (type_approx2 == 19) {
    FPOP<Rec_Rand_lRand> X = FPOP<Rec_Rand_lRand>(data, penalty);
    X.algoFPOP(data, type_approx2, false);
    res =  X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else if (type_approx2 == 21) {
    FPOP<Rec_All_lRand> X = FPOP<Rec_All_lRand>(data, penalty);
    X.algoFPOP(data, type_approx2, false);
    res =  X.TestFPOP(UnpenalizedCost1, LastChpt1);
  }

  return res;
}


//' @title approxFpop
//'
//' @description FPOP method (using the rectangle approximation of the sets) for  the multiple changepoint detection.
//' @param data is a matrix of data (p-rows x n-columns).
//' @param penalty is a value of penalty (a non-negative real number).
//' @param approximation is the type of approximation : 'rectangle' or 'sphere' (by default, 'rectangle').
//' @param intersection is the type of intersection : 'empty', 'all', 'last' or 'random' (by default, 'last').
//' @param exclusion is the type of intersection : 'empty', 'all'or 'random'(by default, 'all').
//' @param NbOfCands is the logical parameter (if NbOfCands = TRUE, than the file "NbOfCands.txt" contains the number of change candidates for each iteration.
//'
//' @return a list of  elements  = (changes, means, UnpenalizedCost, NumberOfCandidates).
//'
//' \describe{
//' \item{\code{changes}}{is the changepoint vector that gives the last index of each segment for the p-variate time series.}
//' \item{\code{means}}{is the list of successive means for the p-variate time series.}
//' \item{\code{UnpenalizedCost}}{is a number equal to the global cost.}
//' \item{\code{NumberOfCandidates}}{is a number of candidates at each iteration (vector).}
//' }
//'
//' @examples
//' N <- 100
//' Chpt <-50
//' Means <-  matrix(c(0,1,1,10), nrow = 2)
//' Noise <- 1
//' Dim <- 2
//' Penalty <- 2*Dim*log(N)
//' time_series <- rnormChanges(p = Dim, n = N, changes = Chpt, means = Means, noise = Noise)
//' time_series <- rnormChanges(p = 2, n = N, changes = NULL, means = matrix(0, ncol = 1, nrow = 2), noise = 1)
//' approxFpop(data = time_series, penalty = Penalty, approximation = 'rectangle', intersection = 'last', exclusion = 'all', NbOfCands = TRUE)
//' approxFpop(data = time_series, penalty = Penalty, approximation = 'rectangle', intersection = 'last', exclusion = 'random',  NbOfCands = FALSE)
//' approxFpop(data = time_series, penalty = Penalty, approximation = 'sphere', intersection = 'last', exclusion = 'all', NbOfCands = FALSE)

// [[Rcpp::export]]
List approxFpop(Rcpp::NumericMatrix data, double penalty, std::string approximation = "rectangle", std::string intersection = "all",  std::string exclusion = "all", bool NbOfCands = false) {
  List res;
  int type_approx = NmbOfapproxFPOP(approximation, intersection, exclusion);
  //----------stop--------------------------------------------------------------
  if (penalty < 0) {throw std::range_error("Penalty should be a non-negative number!");}
  if(type_approx == 0){throw std::range_error("This combination of parameters 'intersection' and 'exclusion' is not available. ");}
  //----------------------------------------------------------------------------
  if (type_approx == 1) {
   FPOP<Rec_Empty_Empty> X = FPOP<Rec_Empty_Empty>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands);
    res = X.ResAlgoFPOP();
  } else  if (type_approx == 2) {
    FPOP<Sph_Last_lAll> X = FPOP<Sph_Last_lAll>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands);
    res = X.ResAlgoFPOP();
  } else if (type_approx == 10) {
    FPOP<Rec_All_Empty> X = FPOP<Rec_All_Empty>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands);
    res = X.ResAlgoFPOP();
  } else if (type_approx == 11) {
    FPOP<Rec_All_lAll> X = FPOP<Rec_All_lAll>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands);
    res = X.ResAlgoFPOP();
  } else if (type_approx == 13) {
    FPOP<Rec_Last_lAll> X = FPOP<Rec_Last_lAll>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands);
    res = X.ResAlgoFPOP();
  } else if (type_approx == 15) {
    FPOP<Rec_Rand_lAll> X = FPOP<Rec_Rand_lAll>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands);
    res = X.ResAlgoFPOP();
  } else if (type_approx == 17) {
    FPOP<Rec_Last_lRand> X = FPOP<Rec_Last_lRand>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands);
    res = X.ResAlgoFPOP();
  } else if (type_approx == 19) {
    FPOP<Rec_Rand_lRand> X = FPOP<Rec_Rand_lRand>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands);
    res = X.ResAlgoFPOP();
  } else if (type_approx == 21) {
    FPOP<Rec_All_lRand> X = FPOP<Rec_All_lRand>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands);
    res = X.ResAlgoFPOP();
  }
  return res;
}

//'@title TestTwoApproxFpop
//'
//' @description Сomparing the parameters ("UnpenalizedCost", "LastChpt") for two different FPOP methods (using the rectangle approximation of the sets) .
//' @param data is a matrix of data (p-rows x n-columns).
//' @param penalty is a value of penalty (a non-negative real number).
//' @param approximation1 is the type of approximation : 'rectangle' or 'sphere' (by default, 'rectangle').
//' @param intersection1 is the type of intersection : 'empty', 'all', 'last', 'random' or 'sphere'(by default, 'last').
//' @param exclusion1 is the type of intersection : 'empty', 'all', 'random' or 'sphere'(by default, 'all').
//' @param approximation2 is the type of approximation : 'rectangle' or 'sphere' (by default, 'rectangle').
//' @param intersection2 is the type of intersection : 'empty', 'all', 'last', 'random' or 'sphere'(by default, 'empty').
//' @param exclusion2 is the type of intersection : 'empty', 'all', 'random' or 'sphere'(by default, 'empty').
//'
//' @return TRUE or FALSE
//'
//' \describe{
//' \item{\code{TRUE}}{ 'TRUE' if parameters are the same.}
//' \item{\code{FALSE}}{'TRUE' if parameters are different.}
//' }
//'
//' @examples
//' N <- 10
//' Chpt <-5
//' Means <-  matrix(c(0,1,1,10), nrow = 2)
//' Noise <- 1
//' Dim <- 2
//' Penality <- 2*Dim*log(N)
//' time_series <- rnormChanges(p = Dim, n = N, changes = Chpt, means = Means, noise = Noise)
//' TestTwoApproxFpop(data = time_series, penalty = Penality, approximation1 = 'rectangle', intersection1 = 'all', exclusion1 = 'all')
//' TestTwoApproxFpop(data = time_series, penalty = Penality, approximation1 = 'rectangle', intersection1 = 'last', exclusion1 = 'all')
//' TestTwoApproxFpop(data = time_series, penalty = Penality, approximation1 = 'rectangle', intersection1 = 'random', exclusion1 = 'random', approximation2 = 'rectangle', intersection2 = 'last', exclusion2 = 'all')
// [[Rcpp::export]]
bool TestTwoApproxFpop(Rcpp::NumericMatrix data, double penalty, std::string approximation1 = "rectangle", std::string intersection1 = "last",  std::string exclusion1 = "all", std::string str1 = "v", std::string approximation2 = "rectangle", std::string intersection2 = "empty",  std::string exclusion2 = "empty",  std::string str2 = "v") {
  int type_approx1 = NmbOfapproxFPOP(approximation1, intersection1, exclusion1);
  int type_approx2 = NmbOfapproxFPOP(approximation2, intersection2, exclusion2);

  //----------stop--------------------------------------------------------------
  if (penalty < 0) {
    throw std::range_error("Penalty should be a non-negative number!");
    return false;
  }
  if((type_approx1 == 0)||(type_approx2 == 0)) {
    throw std::range_error("These combinations of parameters 'intersection','exclusion' and 'str' are  not available. ");
    return false;
  }
  if(type_approx1 == type_approx2) {
    throw std::range_error("These combinations have same parameters 'intersection', 'exclusion' and 'str'. ");
    return true;
  }
  //----------------------------------------------------------------------------
  bool res = false;
  if (type_approx1 == 1) {
    FPOP<Rec_Empty_Empty> X1 = FPOP<Rec_Empty_Empty>(data, penalty);
     X1.algoFPOP(data, type_approx1, false );
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else  if (type_approx1 == 2) {
    FPOP<Sph_Last_lAll> X1 = FPOP<Sph_Last_lAll>(data, penalty);
    X1.algoFPOP(data, type_approx1, false );
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  }  else if (type_approx1 == 10) {
    FPOP<Rec_All_Empty> X1 = FPOP<Rec_All_Empty>(data, penalty);
    X1.algoFPOP(data, type_approx1, false );
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else if (type_approx1 == 11) {
    FPOP<Rec_All_lAll> X1 = FPOP<Rec_All_lAll>(data, penalty);
    X1.algoFPOP(data, type_approx1, false );
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else if (type_approx1 == 13) {
    FPOP<Rec_Last_lAll> X1 = FPOP<Rec_Last_lAll>(data, penalty);
    X1.algoFPOP(data, type_approx1, false );
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else if (type_approx1 == 15) {
    FPOP<Rec_Rand_lAll> X1 = FPOP<Rec_Rand_lAll>(data, penalty);
    X1.algoFPOP(data, type_approx1, false );
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else if (type_approx1 == 17) {
    FPOP<Rec_Last_lRand> X1 = FPOP<Rec_Last_lRand>(data, penalty);
    X1.algoFPOP(data, type_approx1, false );
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else if (type_approx1 == 19) {
    FPOP<Rec_Rand_lRand> X1 = FPOP<Rec_Rand_lRand>(data, penalty);
    X1.algoFPOP(data, type_approx1, false );
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else if (type_approx1 == 21) {
    FPOP<Rec_All_lRand> X1 = FPOP<Rec_All_lRand>(data, penalty);
    X1.algoFPOP(data, type_approx1, false );
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  }
  return res;
}

