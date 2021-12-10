#ifndef CANDIDATE_IEMPTY_EEMPTY_9_H
#define CANDIDATE_IEMPTY_EEMPTY_9_H

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <iterator>
#include <stdio.h>
#include <math.h>

#include <Rcpp.h>

#include "Cost.h"

class Candidate_Iempty_Eempty_9{
private:
  unsigned int Dim;
  unsigned int Tau;
  double** CumSumData;
  double** CumSumData2;
  double* VectOfCosts;
  bool fl_empty;

public:
  Candidate_Iempty_Eempty_9(): Dim(0), Tau(0), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL),    fl_empty(false) { }
  Candidate_Iempty_Eempty_9(unsigned  int dim): Dim(dim), Tau(0),  CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL),  fl_empty(false) { }
  Candidate_Iempty_Eempty_9(unsigned int dim, unsigned int t): Dim(dim), Tau(t), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL),  fl_empty(false) { }
  Candidate_Iempty_Eempty_9(const Candidate_Iempty_Eempty_9 & candidate);
  ~Candidate_Iempty_Eempty_9();

  unsigned int GetTau()const;

  void CleanOfCandidate();
  bool EmptyOfCandidate();
  void InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts);
  void UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Candidate_Iempty_Eempty_9>::iterator> &vectlinktocands, unsigned int& RealNbExclus);
};
#endif //CANDIDATE_IEMPTY_EEMPTY_9_H
