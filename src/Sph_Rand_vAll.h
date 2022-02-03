#ifndef SPH_RAND_VALL_H
#define SPH_RAND_VALL_H

#include "pSphere.h"
#include "Cost.h"
#include <vector>
#include <list>
#include <iterator>
#include <stdio.h>

class Sph_Rand_vAll {
private:
  unsigned int Dim;
  unsigned int Tau;
  double** CumSumData;
  double** CumSumData2;
  double* VectOfCosts;
  bool fl_empty;

public:
  Sph_Rand_vAll(): Dim(0), Tau(0), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL),  fl_empty(false) { }
  Sph_Rand_vAll(unsigned  int dim): Dim(dim), Tau(0), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL), fl_empty(false) { }
  Sph_Rand_vAll(unsigned int dim, unsigned int t): Dim(dim), Tau(t), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL), fl_empty(false) { }
  Sph_Rand_vAll(const Sph_Rand_vAll & candidate);
  ~Sph_Rand_vAll();

  unsigned int GetTau()const;
  int get_Number(int N);
  double Dist(double* a, double*b);
  void CleanOfCandidate();
  bool EmptyOfCandidate();
  void InitialOfCandidate(unsigned int tau, double** &cumsumdata,  double** &cumsumdata2, double* &vectofcosts);
  void UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Sph_Rand_vAll>::iterator> &vectlinktocands, unsigned int& RealNbExclus);
};
#endif //SPH_RAND_VALL_H
//------------------------------------------------------------------------------
