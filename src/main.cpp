#include "FPOP.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

//converting parameters("approximation", "intersection", "exclusion") to a numeric value.
int NmbOfapproxFPOP(std::string approximation, std::string intersection,  std::string exclusion, bool NbOfExclus = false) {
  int type_approx = 0;
  if (approximation == "rectangle") {
    // 2 = (Iall, Eall), 3 = (Iall,Eempty), 7 = (Iall,Erandom)
    if (intersection == "all" ) {
      if (exclusion == "all") {
        type_approx = 2;
      } else if (exclusion == "empty") {
        type_approx = 3;
        if (NbOfExclus) {
          NbOfExclus = false;
          Rcpp::Rcout << "The exclusion number  is 0 for the combination ('all', 'empty')."<< endl;
        }
      } else if (exclusion == "random") {
        type_approx = 7;
        if (NbOfExclus) {
          Rcpp::Rcout << "The maximum number of exclusion for each change candidate is 1 for the combination ('all', 'random')."<< endl;
        }
      }
    }
    //4 = (Iempty, Eall), 9 = (Iempty,Eempty) //<=PELT//
    if (intersection == "empty" ) {
      if (exclusion == "all") {
        type_approx = 4;
      } else if (exclusion == "empty") {
        type_approx = 9;
        if (NbOfExclus) {
          NbOfExclus = false;
          Rcpp::Rcout << "The exclusion number is 0 for the combination ('empty', 'empty')."<< endl;
        }
      }
    }
    //5 = (Ilast,Eall), 6 = (Ilast, Erandom)
    if (intersection == "last" ) {
      if (exclusion == "all") {
        type_approx = 5;
      } else if (exclusion == "random") {
        type_approx = 6;
        if (NbOfExclus) {
          NbOfExclus = false;
          Rcpp::Rcout << "The maximum number of exclusions for each change candidate  is 1 for the combination ('last', 'random')."<< endl;
        }
      }
    }
    //8 = (Irandom, Erandom)
    if ((intersection == "random") && (exclusion == "random"))  {
      type_approx = 8;
      if (NbOfExclus) {
        Rcpp::Rcout << "The maximum number of exclusions for each change candidate  is 1 for the combination ('random', 'random')."<< endl;
      }
    }
    //10 = (Irandom,Eall)
    if (intersection == "random" ) {
      if (exclusion == "all") {
        type_approx = 10;
      }
    }
  }
  if (approximation == "sphere") {
    if (exclusion == "all")  {
      if (intersection == "last") {
        type_approx = 1;
      } else if (intersection == "random") {
        type_approx = 11;
      }
    }
    if (exclusion == "random")  {
      if (intersection == "last") {
        type_approx = 13;
      } else if (intersection == "random") {
        type_approx = 12;
      }
    }
  }
  return type_approx;
}

//Comparison of two FPOP-methods
bool TestOfComparisonTwoFPOP(Rcpp::NumericMatrix data, double penalty, unsigned int type_approx2, double UnpenalizedCost1, unsigned int* LastChpt1) {
  bool res = false;
  if (type_approx2 == 1) {
    FPOP<Candidate_Ilast_Eall_sphere_1> X = FPOP<Candidate_Ilast_Eall_sphere_1>(data, penalty);
    X.algoFPOP(data, type_approx2, false, false);
    res = X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else  if (type_approx2 == 2) {
    FPOP<Candidate_Iall_Eall_2> X = FPOP<Candidate_Iall_Eall_2>(data, penalty);
    X.algoFPOP(data, type_approx2, false, false);
    res = X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else if (type_approx2 == 3) {
    FPOP<Candidate_Iall_Eempty_3> X = FPOP<Candidate_Iall_Eempty_3>(data, penalty);
    X.algoFPOP(data, type_approx2, false, false);
    res = X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else if (type_approx2 == 4) {
    FPOP<Candidate_Iempty_Eall_4> X = FPOP<Candidate_Iempty_Eall_4>(data, penalty);
    X.algoFPOP(data, type_approx2, false, false);
    res =  X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else if (type_approx2 == 5) {
    FPOP<Candidate_Ilast_Eall_5> X = FPOP<Candidate_Ilast_Eall_5>(data, penalty);
    X.algoFPOP(data, type_approx2, false, false);
    res =  X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else if (type_approx2 == 6) {
    FPOP<Candidate_Ilast_Erandom_6> X = FPOP<Candidate_Ilast_Erandom_6>(data, penalty);
    X.algoFPOP(data, type_approx2, false, false);
    res =  X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else if (type_approx2 == 7) {
    FPOP<Candidate_Iall_Erandom_7> X = FPOP<Candidate_Iall_Erandom_7>(data, penalty);
    X.algoFPOP(data, type_approx2, false, false);
    res = X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else if (type_approx2 == 8) {
    FPOP<Candidate_Irandom_Erandom_8> X = FPOP<Candidate_Irandom_Erandom_8>(data, penalty);
    X.algoFPOP(data, type_approx2,false, false);
    res =  X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else if (type_approx2 == 9) {
    FPOP<Candidate_Iempty_Eempty_9> X = FPOP<Candidate_Iempty_Eempty_9>(data, penalty);
    X.algoFPOP(data, type_approx2,false, false);
    res =  X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else if (type_approx2 == 10) {
    FPOP<Candidate_Irandom_Eall_10> X = FPOP<Candidate_Irandom_Eall_10>(data, penalty);
    X.algoFPOP(data, type_approx2, false, false);
    res = X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else if (type_approx2 == 11) {
    FPOP<Candidate_Irandom_Eall_sphere_11> X = FPOP<Candidate_Irandom_Eall_sphere_11>(data, penalty);
    X.algoFPOP(data, type_approx2, false, false);
    res =  X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else if (type_approx2 == 12) {
    FPOP<Candidate_Irandom_Erandom_sphere_12> X = FPOP<Candidate_Irandom_Erandom_sphere_12>(data, penalty);
    X.algoFPOP(data, type_approx2, false, false);
    res =  X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else if (type_approx2 == 13) {
    FPOP<Candidate_Ilast_Erandom_sphere_13> X = FPOP<Candidate_Ilast_Erandom_sphere_13>(data, penalty);
    X.algoFPOP(data, type_approx2, false, false);
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
//' The following parameter combinations are implemented:
//' (approximation = 'rectangle', intersection ='empty', exclusion ='empty') (PELT-method)
//' (approximation = 'rectangle', intersection ='all', exclusion ='all')
//' (approximation = 'rectangle', intersection ='all', exclusion ='empty')
//' (approximation = 'rectangle', intersection ='empty', exclusion ='all')
//' (approximation = 'rectangle', intersection ='last', exclusion ='all')
//' (approximation = 'rectangle', intersection ='last', exclusion ='random')
//' (approximation = 'rectangle', intersection ='all', exclusion ='random')
//' (approximation = 'rectangle', intersection ='random', exclusion ='random')
//' (approximation = 'sphere', intersection ='last', exclusion ='all')
//' (approximation = 'sphere', intersection ='last', exclusion ='random')
//' (approximation = 'sphere', intersection ='random', exclusion ='random')
//' @param NbOfCands is the logical parameter (if NbOfCands = TRUE, than the file "NbOfCands.txt" contains the number of change candidates for each iteration.
//' @param NbOfExclus is the logical parameter (if NbOfExclus = TRUE, than the file "NbOfExclus.txt" contains the label of candidate and the number of exclusion for change candidates for each iteration.
//'
//' @return a list of  elements  = (changes, means, UnpenalizedCost).
//'
//' \describe{
//' \item{\code{changes}}{is the changepoint vector that gives the last index of each segment for the p-variate time series.}
//' \item{\code{means}}{is the list of successive means for the p-variate time series.}
//' \item{\code{UnpenalizedCost}}{is a number equal to the global cost.}
//' }
//'
//' @examples
//' N <- 100
//' Chpt <-50
//' Means <-  matrix(c(0,1,1,10), nrow = 2)
//' Noise <- 1
//' Dim <- 2
//' Penality <- 2*Dim*log(N)
//' time_series <- rnormChanges(p = Dim, n = N, changes = Chpt, means = Means, noise = Noise)
//' time_series <- rnormChanges(p = 2, n = N, changes = NULL, means = matrix(0, ncol = 1, nrow = 2), noise = 1)
//' Approx <- list()
//' Approx[[1]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'all', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE)
//' Approx[[2]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'all', exclusion = 'empty', NbOfCands = FALSE, NbOfExclus = FALSE)
//' Approx[[3]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'all', exclusion = 'random', NbOfCands = FALSE, NbOfExclus = FALSE)
//' Approx[[4]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'empty', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE)
//' Approx[[5]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'last', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE)
//' Approx[[6]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'last', exclusion = 'random', NbOfCands = FALSE, NbOfExclus = FALSE)
//' Approx[[7]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle',intersection = 'random', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE)
//' Approx[[8]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'random', exclusion = 'random', NbOfCands = FALSE, NbOfExclus = FALSE)
//' Approx[[9]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'sphere', intersection = 'last', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE)
//' Approx[[10]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'sphere', intersection = 'random', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE)
//' Approx[[11]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'sphere', intersection = 'random', exclusion = 'random', NbOfCands = FALSE, NbOfExclus = FALSE)
//' Approx[[12]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'sphere', intersection = 'last', exclusion = 'random', NbOfCands = FALSE, NbOfExclus = FALSE)
//'
//' Approx[[13]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'empty', exclusion = 'empty', NbOfCands = FALSE, NbOfExclus = FALSE)
//' Approx
// [[Rcpp::export]]
List approxFpop(Rcpp::NumericMatrix data, double penalty, std::string approximation = "rectangle", std::string intersection = "all",  std::string exclusion = "all", bool NbOfCands = false, bool NbOfExclus = false) {
  List res;
  int type_approx = NmbOfapproxFPOP(approximation, intersection, exclusion, NbOfExclus);
  //----------stop--------------------------------------------------------------
  if (penalty < 0) {throw std::range_error("Penalty should be a non-negative number!");}
  if(type_approx == 0){throw std::range_error("This combination of parameters 'intersection' and 'exclusion' is not available. ");}
  //----------------------------------------------------------------------------
  if (type_approx == 1) {
   FPOP<Candidate_Ilast_Eall_sphere_1> X = FPOP<Candidate_Ilast_Eall_sphere_1>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfExclus);
    res = X.ResAlgoFPOP();
  } else  if (type_approx == 2) {
    FPOP<Candidate_Iall_Eall_2> X = FPOP<Candidate_Iall_Eall_2>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfExclus);
    res = X.ResAlgoFPOP();
  } else if (type_approx == 3) {
    FPOP<Candidate_Iall_Eempty_3> X = FPOP<Candidate_Iall_Eempty_3>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfExclus);
    res = X.ResAlgoFPOP();
  } else if (type_approx == 4) {
    FPOP<Candidate_Iempty_Eall_4> X = FPOP<Candidate_Iempty_Eall_4>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfExclus);
    res = X.ResAlgoFPOP();
  } else if (type_approx == 5) {
    FPOP<Candidate_Ilast_Eall_5> X = FPOP<Candidate_Ilast_Eall_5>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfExclus);
    res = X.ResAlgoFPOP();
  } else if (type_approx == 6) {
    FPOP<Candidate_Ilast_Erandom_6> X = FPOP<Candidate_Ilast_Erandom_6>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfExclus);
    res = X.ResAlgoFPOP();
  } else if (type_approx == 7) {
    FPOP<Candidate_Iall_Erandom_7> X = FPOP<Candidate_Iall_Erandom_7>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfExclus);
    res = X.ResAlgoFPOP();
  } else if (type_approx == 8) {
    FPOP<Candidate_Irandom_Erandom_8> X = FPOP<Candidate_Irandom_Erandom_8>(data, penalty);
    X.algoFPOP(data, type_approx,NbOfCands, NbOfExclus);
    res = X.ResAlgoFPOP();
  } else if (type_approx == 9) {
    FPOP<Candidate_Iempty_Eempty_9> X = FPOP<Candidate_Iempty_Eempty_9>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfExclus);
    res = X.ResAlgoFPOP();
  } else if (type_approx == 10) {
    FPOP<Candidate_Irandom_Eall_10> X = FPOP<Candidate_Irandom_Eall_10>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfExclus);
    res = X.ResAlgoFPOP();
  } else if (type_approx == 11) {
    FPOP<Candidate_Irandom_Eall_sphere_11> X = FPOP<Candidate_Irandom_Eall_sphere_11>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfExclus);
    res = X.ResAlgoFPOP();
  } else if (type_approx == 12) {
    FPOP<Candidate_Irandom_Erandom_sphere_12> X = FPOP<Candidate_Irandom_Erandom_sphere_12>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfExclus);
    res = X.ResAlgoFPOP();
  } else if (type_approx == 13) {
    FPOP<Candidate_Ilast_Erandom_sphere_13> X = FPOP<Candidate_Ilast_Erandom_sphere_13>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfExclus);
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
//' N <- 100
//' Chpt <-50
//' Means <-  matrix(c(0,1,1,10), nrow = 2)
//' Noise <- 1
//' Dim <- 2
//' Penality <- 2*Dim*log(N)
//' time_series <- rnormChanges(p = Dim, n = N, changes = Chpt, means = Means, noise = Noise)
//' TestTwoApproxFpop(data = time_series, penalty = Penality, approximation1 = 'rectangle', intersection1 = 'all', exclusion1 = 'all')
//' TestTwoApproxFpop(data = time_series, penalty = Penality, approximation1 = 'rectangle', intersection1 = 'last', exclusion1 = 'all')
//' TestTwoApproxFpop(data = time_series, penalty = Penality, approximation1 = 'rectangle', intersection1 = 'random', exclusion1 = 'random', approximation2 = 'rectangle', intersection2 = 'last', exclusion2 = 'all')
// [[Rcpp::export]]
bool TestTwoApproxFpop(Rcpp::NumericMatrix data, double penalty, std::string approximation1 = "rectangle", std::string intersection1 = "last",  std::string exclusion1 = "all",std::string approximation2 = "rectangle", std::string intersection2 = "empty",  std::string exclusion2 = "empty") {
  int type_approx1 = NmbOfapproxFPOP(approximation1, intersection1, exclusion1, false);
  int type_approx2 = NmbOfapproxFPOP(approximation2, intersection2, exclusion2, false);

  if (type_approx1 == type_approx2){
    throw std::range_error("Combinations have the same parameters 'intersection' and 'exclusion'! ");
    return true;
  }
  //----------stop--------------------------------------------------------------
  if (penalty < 0) {
    throw std::range_error("Penalty should be a non-negative number!");
    return false;
  }
  if((type_approx1 == 0)||(type_approx2 == 0)) {
    throw std::range_error("These combinations of parameters 'intersection' and 'exclusion' are  not available. ");
    return false;
  }
  if(type_approx1 == type_approx2) {
    throw std::range_error("These combinations have same parameters 'intersection' and 'exclusion'. ");
    return true;
  }
  //----------------------------------------------------------------------------
  bool res = false;
  if (type_approx1 == 1) {
    FPOP<Candidate_Ilast_Eall_sphere_1> X1 = FPOP<Candidate_Ilast_Eall_sphere_1>(data, penalty);
     X1.algoFPOP(data, type_approx1, false, false);
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else  if (type_approx1 == 2) {
    FPOP<Candidate_Iall_Eall_2> X1 = FPOP<Candidate_Iall_Eall_2>(data, penalty);
    X1.algoFPOP(data, type_approx1, false, false);
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else if (type_approx1 == 3) {
    FPOP<Candidate_Iall_Eempty_3> X1 = FPOP<Candidate_Iall_Eempty_3>(data, penalty);
    X1.algoFPOP(data, type_approx1, false, false);
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else if (type_approx1 == 4) {
    FPOP<Candidate_Iempty_Eall_4> X1 = FPOP<Candidate_Iempty_Eall_4>(data, penalty);
    X1.algoFPOP(data, type_approx1, false, false);
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else if (type_approx1 == 5) {
    FPOP<Candidate_Ilast_Eall_5> X1 = FPOP<Candidate_Ilast_Eall_5>(data, penalty);
    X1.algoFPOP(data, type_approx1, false, false);
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else if (type_approx1 == 6) {
    FPOP<Candidate_Ilast_Erandom_6> X1 = FPOP<Candidate_Ilast_Erandom_6>(data, penalty);
    X1.algoFPOP(data, type_approx1, false, false);
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else if (type_approx1 == 7) {
    FPOP<Candidate_Iall_Erandom_7> X1 = FPOP<Candidate_Iall_Erandom_7>(data, penalty);
    X1.algoFPOP(data, type_approx1, false, false);
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else if (type_approx1 == 8) {
    FPOP<Candidate_Irandom_Erandom_8> X1 = FPOP<Candidate_Irandom_Erandom_8>(data, penalty);
    X1.algoFPOP(data, type_approx1, false, false);
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else if (type_approx1 == 9) {
    FPOP<Candidate_Iempty_Eempty_9> X1 = FPOP<Candidate_Iempty_Eempty_9>(data, penalty);
    X1.algoFPOP(data, type_approx1, false, false);
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else if (type_approx1 == 10) {
    FPOP<Candidate_Irandom_Eall_10> X1 = FPOP<Candidate_Irandom_Eall_10>(data, penalty);
    X1.algoFPOP(data, type_approx1, false, false);
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else if (type_approx1 == 11) {
    FPOP<Candidate_Irandom_Eall_sphere_11> X1 = FPOP<Candidate_Irandom_Eall_sphere_11>(data, penalty);
    X1.algoFPOP(data, type_approx1, false, false);
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else if (type_approx1 == 12) {
    FPOP<Candidate_Irandom_Erandom_sphere_12> X1 = FPOP<Candidate_Irandom_Erandom_sphere_12>(data, penalty);
    X1.algoFPOP(data, type_approx1, false, false);
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else if (type_approx1 == 13) {
    FPOP<Candidate_Ilast_Erandom_sphere_13> X1 = FPOP<Candidate_Ilast_Erandom_sphere_13>(data, penalty);
    X1.algoFPOP(data, type_approx1, false, false);
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  }
  return res;
}

