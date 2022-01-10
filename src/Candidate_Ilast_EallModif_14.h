#ifndef CANDIDATE_ILAST_EALLMODIF_14_H
#define CANDIDATE_ILAST_EALLMODIF_14_H

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <iterator>
#include <stdio.h>
#include <math.h>

#include <Rcpp.h>

#include "pRectangle.h"
#include "Cost.h"

class Candidate_Ilast_EallModif_14 {
private:
  unsigned int Dim;
  unsigned int Tau;
  pRectangle* Rect;
  double** CumSumData;
  double** CumSumData2;
  double* VectOfCosts;
  std::list<pSphere> disks_t_1;                            //list of disks(t-1)


public:
  Candidate_Ilast_EallModif_14(): Dim(0), Tau(0), Rect(0), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL) { }
  Candidate_Ilast_EallModif_14(unsigned  int dim): Dim(dim), Tau(0), Rect(new pRectangle(dim)), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL) { }
  Candidate_Ilast_EallModif_14(unsigned int dim, unsigned int t): Dim(dim), Tau(t), Rect(new pRectangle(dim)), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL) { }
  Candidate_Ilast_EallModif_14(const Candidate_Ilast_EallModif_14 & candidate);
  ~Candidate_Ilast_EallModif_14();

  unsigned int GetTau()const;
  std::list<pSphere> get_disks_t_1() const;
  void CleanOfCandidate();
  bool EmptyOfCandidate();
  void InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts);
  void UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Candidate_Ilast_EallModif_14>::iterator> &vectlinktocands, unsigned int& RealNbExclus);
};
#endif //CANDIDATE_ILAST_EALLMODIF_14_H
