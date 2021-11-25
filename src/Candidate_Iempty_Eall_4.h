#ifndef CANDIDATE_IEMPTY_EALL_4_H
#define CANDIDATE_IEMPTY_EALL_4_H

#include <iostream>
#include <vector>
#include <list>

#include "pRectangle.h"
#include "Cost.h"

class Candidate_Iempty_Eall_4{
private:
  unsigned int Dim;
  unsigned int Tau;
  pRectangle* Rect;
  double** CumSumData;
  double* VectOfCosts;

public:
  Candidate_Iempty_Eall_4(): Dim(0), Tau(0), Rect(0), CumSumData(NULL), VectOfCosts(NULL){}
  Candidate_Iempty_Eall_4(unsigned  int dim): Dim(dim), Tau(0), Rect(new pRectangle(dim)), CumSumData(NULL), VectOfCosts(NULL) {}
  Candidate_Iempty_Eall_4(unsigned int dim, unsigned int t): Dim(dim), Tau(t), Rect(new pRectangle(dim)), CumSumData(NULL), VectOfCosts(NULL){}
  Candidate_Iempty_Eall_4(const Candidate_Iempty_Eall_4 & candidate);
  ~Candidate_Iempty_Eall_4();

  unsigned int GetTau()const;

  void CleanOfCandidate();
  bool EmptyOfCandidate();
  void InitialOfCandidate(unsigned int t, double** &cumsumdata, double* &vectofcosts);
  void UpdateOfCandidate(unsigned int t, std::vector<std::list<Candidate_Iempty_Eall_4>::reverse_iterator> &vectlinktocands);
};
#endif //CANDIDATE_IEMPTY_EALL_4_H
