#ifndef CANDIDATE_ILAST_ERANDOM_6_H
#define CANDIDATE_ILAST_ERANDOM_6_H

#include <iostream>
#include <vector>
#include <list>

#include "pRectangle.h"
#include "Cost.h"


class Candidate_Ilast_Erandom_6{
private:
  unsigned int Dim;
  unsigned int Tau;
  pRectangle* Rect;
  double** CumSumData;
  double** CumSumData2;
  double* VectOfCosts;

public:
  Candidate_Ilast_Erandom_6(): Dim(0), Tau(0), Rect(0), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL) { }
  Candidate_Ilast_Erandom_6(unsigned  int dim): Dim(dim), Tau(0), Rect(new pRectangle(dim)), CumSumData(NULL), CumSumData2(NULL),  VectOfCosts(NULL) { }
  Candidate_Ilast_Erandom_6(unsigned int dim, unsigned int t): Dim(dim), Tau(t), Rect(new pRectangle(dim)), CumSumData(NULL), CumSumData2(NULL),  VectOfCosts(NULL) { }
  Candidate_Ilast_Erandom_6(const Candidate_Ilast_Erandom_6 & candidate);
  ~Candidate_Ilast_Erandom_6();

  unsigned int GetTau()const;
  int get_Number(int N);

  void CleanOfCandidate();
  bool EmptyOfCandidate();
  void InitialOfCandidate(unsigned int t, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts);
  void UpdateOfCandidate(unsigned int i, std::vector<std::list<Candidate_Ilast_Erandom_6>::iterator> &vectlinktocands, unsigned int& RealNbExclus);
};
#endif // CANDIDATE_ILAST_ERANDOM_6_H
