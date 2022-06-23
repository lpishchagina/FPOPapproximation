// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// approxFpop
List approxFpop(Rcpp::NumericMatrix data, double penalty, std::string approximation, std::string intersection, std::string exclusion, bool NbOfCands);
RcppExport SEXP _FPOPapproximation_approxFpop(SEXP dataSEXP, SEXP penaltySEXP, SEXP approximationSEXP, SEXP intersectionSEXP, SEXP exclusionSEXP, SEXP NbOfCandsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type penalty(penaltySEXP);
    Rcpp::traits::input_parameter< std::string >::type approximation(approximationSEXP);
    Rcpp::traits::input_parameter< std::string >::type intersection(intersectionSEXP);
    Rcpp::traits::input_parameter< std::string >::type exclusion(exclusionSEXP);
    Rcpp::traits::input_parameter< bool >::type NbOfCands(NbOfCandsSEXP);
    rcpp_result_gen = Rcpp::wrap(approxFpop(data, penalty, approximation, intersection, exclusion, NbOfCands));
    return rcpp_result_gen;
END_RCPP
}
// TestTwoApproxFpop
bool TestTwoApproxFpop(Rcpp::NumericMatrix data, double penalty, std::string approximation1, std::string intersection1, std::string exclusion1, std::string str1, std::string approximation2, std::string intersection2, std::string exclusion2, std::string str2);
RcppExport SEXP _FPOPapproximation_TestTwoApproxFpop(SEXP dataSEXP, SEXP penaltySEXP, SEXP approximation1SEXP, SEXP intersection1SEXP, SEXP exclusion1SEXP, SEXP str1SEXP, SEXP approximation2SEXP, SEXP intersection2SEXP, SEXP exclusion2SEXP, SEXP str2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type penalty(penaltySEXP);
    Rcpp::traits::input_parameter< std::string >::type approximation1(approximation1SEXP);
    Rcpp::traits::input_parameter< std::string >::type intersection1(intersection1SEXP);
    Rcpp::traits::input_parameter< std::string >::type exclusion1(exclusion1SEXP);
    Rcpp::traits::input_parameter< std::string >::type str1(str1SEXP);
    Rcpp::traits::input_parameter< std::string >::type approximation2(approximation2SEXP);
    Rcpp::traits::input_parameter< std::string >::type intersection2(intersection2SEXP);
    Rcpp::traits::input_parameter< std::string >::type exclusion2(exclusion2SEXP);
    Rcpp::traits::input_parameter< std::string >::type str2(str2SEXP);
    rcpp_result_gen = Rcpp::wrap(TestTwoApproxFpop(data, penalty, approximation1, intersection1, exclusion1, str1, approximation2, intersection2, exclusion2, str2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_FPOPapproximation_approxFpop", (DL_FUNC) &_FPOPapproximation_approxFpop, 6},
    {"_FPOPapproximation_TestTwoApproxFpop", (DL_FUNC) &_FPOPapproximation_TestTwoApproxFpop, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_FPOPapproximation(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
