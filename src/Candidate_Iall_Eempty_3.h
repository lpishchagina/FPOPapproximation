#ifndef CANDIDATE_Iall_Eempty_3_H
#define CANDIDATE_Iall_Eempty_3_H

#include <iostream>
#include <vector>
#include <list>

#include "pRectangle.h"
#include "Cost.h"


class Candidate_Iall_Eempty_3{
private:
  unsigned int Dim;
  unsigned int Tau;
  pRectangle* Rect;
  double** CumSumData;
  double* VectOfCosts;

public:
  Candidate_Iall_Eempty_3(): Dim(0), Tau(0), Rect(0), CumSumData(NULL), VectOfCosts(NULL){}
  Candidate_Iall_Eempty_3(unsigned  int dim): Dim(dim), Tau(0), Rect(new pRectangle(dim)), CumSumData(NULL), VectOfCosts(NULL) {}
  Candidate_Iall_Eempty_3(unsigned int dim, unsigned int t): Dim(dim), Tau(t), Rect(new pRectangle(dim)), CumSumData(NULL), VectOfCosts(NULL){}
  Candidate_Iall_Eempty_3(const Candidate_Iall_Eempty_3 & candidate);
  ~Candidate_Iall_Eempty_3();

  unsigned int GetTau()const;

  void CleanOfCandidate();
  bool EmptyOfCandidate();
  void InitialOfCandidate(unsigned int t, double** &cumsumdata, double* &vectofcosts);
  void UpdateOfCandidate(unsigned int t, std::vector<std::list<Candidate_Iall_Eempty_3>::reverse_iterator> &vectlinktocands);
};
#endif //CANDIDATE_Iall_Eempty_3_H