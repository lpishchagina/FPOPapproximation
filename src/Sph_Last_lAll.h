#ifndef SPH_LAST_LALL_H
#define SPH_LAST_LALL_H
#include "pSphere.h"
#include "Cost.h"
#include <vector>
#include <list>
#include <iterator>
#include <stdio.h>

class Sph_Last_lAll {
private:
  unsigned int Dim;
  unsigned int Tau;
  double** CumSumData;
  double** CumSumData2;
  double* VectOfCosts;
  bool fl_empty;
  std::list<pSphere> DiskListBefore;
  bool CreationFl;

public:
  Sph_Last_lAll(): Dim(0), Tau(0), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL),  fl_empty(false) { }
  Sph_Last_lAll(unsigned  int dim): Dim(dim), Tau(0), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL), fl_empty(false) { }
  Sph_Last_lAll(unsigned int dim, unsigned int t): Dim(dim), Tau(t), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL), fl_empty(false) { }
  Sph_Last_lAll(const Sph_Last_lAll & candidate);
  ~Sph_Last_lAll();

  unsigned int GetTau()const;
  double Dist(double* a, double*b);
  void CleanOfCandidate();
  bool EmptyOfCandidate();
  void InitialOfCandidate(unsigned int tau, double** &cumsumdata,  double** &cumsumdata2, double* &vectofcosts);
  void UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Sph_Last_lAll>::iterator> &vectlinktocands, unsigned int& RealNbExclus);
};
#endif //SPH_LAST_LALL_H
//------------------------------------------------------------------------------
