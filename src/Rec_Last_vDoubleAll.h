#ifndef REC_LAST_VDOUBLEALL_H
#define REC_LAST_VDOUBLEALL_H

#include "pRectangle.h"
#include "Cost.h"
#include <vector>
#include <list>
#include <iterator>
#include <stdio.h>

class Rec_Last_vDoubleAll {
private:
  unsigned int Dim;
  unsigned int Tau;
  pRectangle* Rect;
  double** CumSumData;
  double** CumSumData2;
  double* VectOfCosts;

public:
  Rec_Last_vDoubleAll(): Dim(0), Tau(0), Rect(0), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL) { }
  Rec_Last_vDoubleAll(unsigned  int dim): Dim(dim), Tau(0), Rect(new pRectangle(dim)), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL) { }
  Rec_Last_vDoubleAll(unsigned int dim, unsigned int t): Dim(dim), Tau(t), Rect(new pRectangle(dim)), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL) { }
  Rec_Last_vDoubleAll(const Rec_Last_vDoubleAll& candidate);
  ~Rec_Last_vDoubleAll();

  unsigned int GetTau()const;

  void CleanOfCandidate();
  bool EmptyOfCandidate();
  void InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts);
  void UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Rec_Last_vDoubleAll>::iterator> &vectlinktocands, unsigned int& RealNbExclus);
};
#endif //REC_LAST_VDOUBLEALL_H
