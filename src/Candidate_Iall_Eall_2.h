#ifndef CANDIDATE_IALL_EALL_2_H
#define CANDIDATE_IALL_EALL_2_H

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

class Candidate_Iall_Eall_2 {
private:
  unsigned int Dim;
  unsigned int Tau;
  pRectangle* Rect;
  double** CumSumData;
  double* VectOfCosts;

public:
  Candidate_Iall_Eall_2(): Dim(0), Tau(0), Rect(0), CumSumData(NULL), VectOfCosts(NULL){}
  Candidate_Iall_Eall_2(unsigned  int dim): Dim(dim), Tau(0), Rect(new pRectangle(dim)), CumSumData(NULL), VectOfCosts(NULL) {}
  Candidate_Iall_Eall_2(unsigned int dim, unsigned int t): Dim(dim), Tau(t), Rect(new pRectangle(dim)), CumSumData(NULL), VectOfCosts(NULL){}
  Candidate_Iall_Eall_2(const Candidate_Iall_Eall_2 & candidate);
  ~Candidate_Iall_Eall_2();

  unsigned int GetTau()const;
  void CleanOfCandidate();
  bool EmptyOfCandidate();
  void InitialOfCandidate(unsigned int t, double** &cumsumdata, double* &vectofcosts);
  void UpdateOfCandidate(unsigned int i, std::vector<std::list<Candidate_Iall_Eall_2>::iterator> &vectlinktocands, unsigned int& RealNbExclus);
};
#endif //CANDIDATE_IALL_EALL_2_H
