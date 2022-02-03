#ifndef Rec_Empty_Empty_H
#define Rec_Empty_Empty_H

#include "pRectangle.h"
#include "Cost.h"
#include <vector>
#include <list>
#include <iterator>
#include <stdio.h>

class Rec_Empty_Empty{
private:
  unsigned int Dim;
  unsigned int Tau;
  double** CumSumData;
  double** CumSumData2;
  double* VectOfCosts;
  bool fl_empty;

public:
  Rec_Empty_Empty(): Dim(0), Tau(0), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL),    fl_empty(false) { }
  Rec_Empty_Empty(unsigned  int dim): Dim(dim), Tau(0),  CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL),  fl_empty(false) { }
  Rec_Empty_Empty(unsigned int dim, unsigned int t): Dim(dim), Tau(t), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL),  fl_empty(false) { }
  Rec_Empty_Empty(const Rec_Empty_Empty & candidate);
  ~Rec_Empty_Empty();

  unsigned int GetTau()const;

  void CleanOfCandidate();
  bool EmptyOfCandidate();
  void InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts);
  void UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Rec_Empty_Empty>::iterator> &vectlinktocands, unsigned int& RealNbExclus);
};
#endif //Rec_Empty_Empty_H
