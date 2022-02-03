#ifndef REC_RAND_LRAND_H
#define REC_RAND_LRAND_H

#include "pRectangle.h"
#include "Cost.h"
#include <vector>
#include <list>
#include <iterator>
#include <stdio.h>

class Rec_Rand_lRand{
private:
  unsigned int Dim;
  unsigned int Tau;
  pRectangle* Rect;
  double** CumSumData;
  double** CumSumData2;
  double* VectOfCosts;
  std::vector<unsigned int> IndexVectBefore;
  bool CreationFl;

public:
  Rec_Rand_lRand(): Dim(0), Tau(0), Rect(0), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL) , CreationFl (true) { }
  Rec_Rand_lRand(unsigned  int dim): Dim(dim), Tau(0), Rect(new pRectangle(dim)), CumSumData(NULL), CumSumData2(NULL),  VectOfCosts(NULL) , CreationFl (true) { }
  Rec_Rand_lRand(unsigned int dim, unsigned int t): Dim(dim), Tau(t), Rect(new pRectangle(dim)), CumSumData(NULL), CumSumData2(NULL),  VectOfCosts(NULL) , CreationFl (true) { }
  Rec_Rand_lRand(const Rec_Rand_lRand & candidate);
  ~Rec_Rand_lRand();

  unsigned int GetTau()const;
  int get_Number(int N);

  void CleanOfCandidate();
  bool EmptyOfCandidate();
  void InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts);
  void UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Rec_Rand_lRand>::iterator> &vectlinktocands, unsigned int& RealNbExclus);
};
#endif //REC_RAND_LRAND_H
