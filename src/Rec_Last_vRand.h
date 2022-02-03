#ifndef RECT_LAST_VRAND_H
#define RECT_LAST_VRAND_H

#include "pRectangle.h"
#include "Cost.h"
#include <vector>
#include <list>
#include <iterator>
#include <stdio.h>

class Rec_Last_vRand{
private:
  unsigned int Dim;
  unsigned int Tau;
  pRectangle* Rect;
  double** CumSumData;
  double** CumSumData2;
  double* VectOfCosts;

public:
  Rec_Last_vRand(): Dim(0), Tau(0), Rect(0), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL) { }
  Rec_Last_vRand(unsigned  int dim): Dim(dim), Tau(0), Rect(new pRectangle(dim)), CumSumData(NULL), CumSumData2(NULL),  VectOfCosts(NULL) { }
  Rec_Last_vRand(unsigned int dim, unsigned int t): Dim(dim), Tau(t), Rect(new pRectangle(dim)), CumSumData(NULL), CumSumData2(NULL),  VectOfCosts(NULL) { }
  Rec_Last_vRand(const Rec_Last_vRand & candidate);
  ~Rec_Last_vRand();

  unsigned int GetTau()const;
  int get_Number(int N);

  void CleanOfCandidate();
  bool EmptyOfCandidate();
  void InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts);
  void UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Rec_Last_vRand>::iterator> &vectlinktocands, unsigned int& RealNbExclus);
};
#endif // RECT_LAST_VRAND_H
