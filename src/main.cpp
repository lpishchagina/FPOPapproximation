#include "FPOP.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

//' @title approxFpop
//'
//' @description FPOP method (using the rectangle approximation of the sets) for  the multiple changepoint detection.
//' @param data is a matrix of data (p-rows x n-columns).
//' @param penalty is a value of penalty (a non-negative real number).
//' @param approximation is the type of approximation : 'rectangle' or 'sphere' (by default, 'rectangle').
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
//' @examples approxFpop(data = chpt_rnorm(p = 3, n = 100, changes = 50, means = matrix(c (1,2,3,4, 5, 7), nrow = 3), noise = 1), penalty = 2*log(100), type_approx = 2)
//' N <- 1000
//' Chpt <-NULL
//' Means <-  matrix(c(0,1,1,10), nrow = 2)
//' Noise <- 1
//' Dim <- 2
//' Penality <- 2*Dim*log(N)
//' time_series <- rnormChanges(p = Dim, n = N, changes = Chpt, means = Means, noise = Noise)
//' time_series <- rnormChanges(p = 2, n = N, changes = NULL, means = matrix(0, ncol = 1, nrow = 2), noise = 1)

//'Approx <- list()
//'Approx[[1]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'all', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE)
//'Approx[[2]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'all', exclusion = 'empty', NbOfCands = FALSE, NbOfExclus = FALSE)
//'Approx[[3]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'all', exclusion = 'random', NbOfCands = FALSE, NbOfExclus = FALSE)
//'
//'Approx[[4]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'empty', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE)
//'
//'Approx[[5]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'last', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE)
//'Approx[[6]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'last', exclusion = 'random', NbOfCands = FALSE, NbOfExclus = FALSE)
//'
//'Approx[[7]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle',intersection = 'random', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE)
//'Approx[[8]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'random', exclusion = 'random', NbOfCands = FALSE, NbOfExclus = FALSE)
//'
//'Approx[[9]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'sphere', intersection = 'last', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE)
//'Approx[[10]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'sphere', intersection = 'random', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE)
//'Approx[[11]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'sphere', intersection = 'random', exclusion = 'random', NbOfCands = FALSE, NbOfExclus = FALSE)
//'Approx[[12]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'sphere', intersection = 'last', exclusion = 'random', NbOfCands = FALSE, NbOfExclus = FALSE)
//'
//'Approx[[13]] <- approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'empty', exclusion = 'empty', NbOfCands = FALSE, NbOfExclus = FALSE)
//'Approx
//'
//'TimeApprox <- list()
//'TimeApprox[[1]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'all', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE))
//'TimeApprox[[2]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'all', exclusion = 'empty', NbOfCands = FALSE, NbOfExclus = FALSE))
//'TimeApprox[[3]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'all', exclusion = 'random', NbOfCands = FALSE, NbOfExclus = FALSE))
//'
//'TimeApprox[[4]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'empty', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE))
//'
//'TimeApprox[[5]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'last', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE))
//'TimeApprox[[6]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'last', exclusion = 'random', NbOfCands = FALSE, NbOfExclus = FALSE))
//'
//'TimeApprox[[7]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle',intersection = 'random', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE))
//'TimeApprox[[8]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'random', exclusion = 'random', NbOfCands = FALSE, NbOfExclus = FALSE))
//'
//'TimeApprox[[9]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'sphere', intersection = 'last', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE))
//'TimeApprox[[10]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'sphere', intersection = 'random', exclusion = 'all', NbOfCands = FALSE, NbOfExclus = FALSE))
//'TimeApprox[[11]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'sphere', intersection = 'random', exclusion = 'random', NbOfCands = FALSE, NbOfExclus = FALSE))
//'TimeApprox[[12]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'sphere', intersection = 'last', exclusion = 'random', NbOfCands = FALSE, NbOfExclus = FALSE))
//'
//'TimeApprox[[13]] <- system.time(approxFpop(data = time_series, penalty = Penality, approximation = 'rectangle', intersection = 'empty', exclusion = 'empty', NbOfCands = FALSE, NbOfExclus = FALSE))
//'TimeApprox
// [[Rcpp::export]]
List approxFpop(Rcpp::NumericMatrix data, double penalty, std::string approximation = "rectangle", std::string intersection = "all",  std::string exclusion = "all", bool NbOfCands = false, bool NbOfExclus = false) {
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

  //----------stop--------------------------------------------------------------
  if (penalty < 0) {throw std::range_error("Penalty should be a non-negative number!");}
  if(type_approx < 1 || type_approx > 13)
  {throw std::range_error("This combination of parameters 'intersection' and 'exclusion' is not available. ");}
  //----------------------------------------------------------------------------
  List res;
  if (type_approx == 1) {
    FPOP<Candidate_Ilast_Eall_sphere_1> X = FPOP<Candidate_Ilast_Eall_sphere_1>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfExclus);
    res["changes"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["UnpenalizedCost"] = X.GetUnpenalizedCost();
  }
  if (type_approx == 2) {
    // test = true;
    FPOP<Candidate_Iall_Eall_2> X = FPOP<Candidate_Iall_Eall_2>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfExclus);
    res["changes"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["UnpenalizedCost"] = X.GetUnpenalizedCost();
  }
  if (type_approx == 3) {
    // test = true;
    FPOP<Candidate_Iall_Eempty_3> X = FPOP<Candidate_Iall_Eempty_3>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfExclus);
    res["changes"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["UnpenalizedCost"] = X.GetUnpenalizedCost();
  }
  if (type_approx == 4) {
    // test = true;
    FPOP<Candidate_Iempty_Eall_4> X = FPOP<Candidate_Iempty_Eall_4>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfExclus);
    res["changes"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["UnpenalizedCost"] = X.GetUnpenalizedCost();
  }
  if (type_approx == 5) {
    // test = true;
    FPOP<Candidate_Ilast_Eall_5> X = FPOP<Candidate_Ilast_Eall_5>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfExclus);
    res["changes"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["UnpenalizedCost"] = X.GetUnpenalizedCost();
  }
  if (type_approx == 6) {
    // test = true;
    FPOP<Candidate_Ilast_Erandom_6> X = FPOP<Candidate_Ilast_Erandom_6>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfExclus);
    res["changes"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["UnpenalizedCost"] = X.GetUnpenalizedCost();
  }
  if (type_approx == 7) {
    // test = true;
    FPOP<Candidate_Iall_Erandom_7> X = FPOP<Candidate_Iall_Erandom_7>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfExclus);
    res["changes"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["UnpenalizedCost"] = X.GetUnpenalizedCost();
  }
  if (type_approx == 8) {
    // test = true;
    FPOP<Candidate_Irandom_Erandom_8> X = FPOP<Candidate_Irandom_Erandom_8>(data, penalty);
    X.algoFPOP(data, type_approx,NbOfCands, NbOfExclus);
    res["changes"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["UnpenalizedCost"] = X.GetUnpenalizedCost();
  }
  if (type_approx == 9) {
    // test = true;
    FPOP<Candidate_Iempty_Eempty_9> X = FPOP<Candidate_Iempty_Eempty_9>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfExclus);
    res["changes"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["UnpenalizedCost"] = X.GetUnpenalizedCost();
  }
  if (type_approx == 10) {
    // test = true;
    FPOP<Candidate_Irandom_Eall_10> X = FPOP<Candidate_Irandom_Eall_10>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfExclus);
    res["changes"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["UnpenalizedCost"] = X.GetUnpenalizedCost();
  }
  if (type_approx == 11) {
    // test = true;
    FPOP<Candidate_Irandom_Eall_sphere_11> X = FPOP<Candidate_Irandom_Eall_sphere_11>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfExclus);
    res["changes"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["UnpenalizedCost"] = X.GetUnpenalizedCost();
  }
  if (type_approx == 12) {
    // test = true;
    FPOP<Candidate_Irandom_Erandom_sphere_12> X = FPOP<Candidate_Irandom_Erandom_sphere_12>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfExclus);
    res["changes"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["UnpenalizedCost"] = X.GetUnpenalizedCost();
  }
  if (type_approx == 13) {
    // test = true;
    FPOP<Candidate_Ilast_Erandom_sphere_13> X = FPOP<Candidate_Ilast_Erandom_sphere_13>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfExclus);
    res["changes"] = X.GetChanges();
    res["means"] = X.GetSegmentMeans();
    res["UnpenalizedCost"] = X.GetUnpenalizedCost();
  }
  return res;
}
