#ifndef CANDIDATE_IRANDOM_ERANDOM_8_H
#define CANDIDATE_IRANDOM_ERANDOM_8_H

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

class Candidate_Irandom_Erandom_8 {
private:
  unsigned int Dim;
  unsigned int Tau;
  pRectangle* Rect;
  double** CumSumData;
  double** CumSumData2;
  double* VectOfCosts;

public:
  Candidate_Irandom_Erandom_8(): Dim(0), Tau(0), Rect(0), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL) { }
  Candidate_Irandom_Erandom_8(unsigned  int dim): Dim(dim), Tau(0), Rect(new pRectangle(dim)), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL) { }
  Candidate_Irandom_Erandom_8(unsigned int dim, unsigned int t): Dim(dim), Tau(t), Rect(new pRectangle(dim)), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL) { }
  Candidate_Irandom_Erandom_8(const Candidate_Irandom_Erandom_8 & candidate);
  ~Candidate_Irandom_Erandom_8();

  unsigned int GetTau()const;
  int get_Number(int N);

  void CleanOfCandidate();
  bool EmptyOfCandidate();
  void InitialOfCandidate(unsigned int t, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts);
  void UpdateOfCandidate(unsigned int i, std::vector<std::list<Candidate_Irandom_Erandom_8>::iterator> &vectlinktocands, unsigned int& RealNbExclus);
};
#endif // CANDIDATE_IRANDOM_ERANDOM_8_H
