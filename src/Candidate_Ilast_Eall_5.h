#ifndef CANDIDATE_ILAST_EALL_5_H
#define CANDIDATE_ILAST_EALL_5_H

#include <iostream>
#include <vector>
#include <list>

#include "pRectangle.h"
#include "Cost.h"

class Candidate_Ilast_Eall_5 {
private:
  unsigned int Dim;
  unsigned int Tau;
  pRectangle* Rect;
  double** CumSumData;
  double* VectOfCosts;

public:
  Candidate_Ilast_Eall_5(): Dim(0), Tau(0), Rect(0), CumSumData(NULL), VectOfCosts(NULL) { }
  Candidate_Ilast_Eall_5(unsigned  int dim): Dim(dim), Tau(0), Rect(new pRectangle(dim)), CumSumData(NULL), VectOfCosts(NULL) { }
  Candidate_Ilast_Eall_5(unsigned int dim, unsigned int t): Dim(dim), Tau(t), Rect(new pRectangle(dim)), CumSumData(NULL), VectOfCosts(NULL) { }
  Candidate_Ilast_Eall_5(const Candidate_Ilast_Eall_5 & candidate);
  ~Candidate_Ilast_Eall_5();

  unsigned int GetTau()const;
  void CleanOfCandidate();
  bool EmptyOfCandidate();
  void InitialOfCandidate(unsigned int i, double** &cumsumdata, double* &vectofcosts);
  void UpdateOfCandidate(unsigned int t, std::vector<std::list<Candidate_Ilast_Eall_5>::iterator> &vectlinktocands);
};
#endif //CANDIDATE_ILAST_EALL_5_H
