#ifndef CANDIDATE_IRANDOM_EALL_10_H
#define CANDIDATE_IRANDOM_EALL_10_H

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <iterator>
#include <stdio.h>
#include <math.h>
#include <Rcpp.h>

#include "pRectangle.h"
#include "Cost.h"

class Candidate_Irandom_Eall_10{
private:
  unsigned int Dim;
  unsigned int Tau;
  pRectangle* Rect;
  double** CumSumData;
  double** CumSumData2;
  double* VectOfCosts;

public:
  Candidate_Irandom_Eall_10(): Dim(0), Tau(0), Rect(0), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL) { }
  Candidate_Irandom_Eall_10(unsigned  int dim): Dim(dim), Tau(0), Rect(new pRectangle(dim)), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL) { }
  Candidate_Irandom_Eall_10(unsigned int dim, unsigned int t): Dim(dim), Tau(t), Rect(new pRectangle(dim)), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL) { }
  Candidate_Irandom_Eall_10(const Candidate_Irandom_Eall_10& candidate);
  ~Candidate_Irandom_Eall_10();

  unsigned int GetTau()const;
  int get_Number(int N);

  void CleanOfCandidate();
  bool EmptyOfCandidate();
  void InitialOfCandidate(unsigned int tau, double** &cumsumdata,  double** &cumsumdata2, double* &vectofcosts);
  void UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Candidate_Irandom_Eall_10>::iterator> &vectlinktocands, unsigned int& RealNbExclus);
};
#endif // CANDIDATE_IRANDOM_EALL_10_H
