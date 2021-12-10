#ifndef CANDIDATE_ILAST_EALL_SPHERE_1_H
#define CANDIDATE_ILAST_EALL_SPHERE_1_H
#include <vector>
#include <iostream>
#include <fstream>
#include <list>
#include <iterator>
#include <stdio.h>
#include <math.h>

#include <Rcpp.h>

#include "pSphere.h"
#include "Cost.h"

class Candidate_Ilast_Eall_sphere_1 {
private:
  unsigned int Dim;
  unsigned int Tau;
  double** CumSumData;
  double** CumSumData2;
  double* VectOfCosts;
  bool fl_empty;

public:
  Candidate_Ilast_Eall_sphere_1(): Dim(0), Tau(0), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL),  fl_empty(false) { }
  Candidate_Ilast_Eall_sphere_1(unsigned  int dim): Dim(dim), Tau(0), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL), fl_empty(false) { }
  Candidate_Ilast_Eall_sphere_1(unsigned int dim, unsigned int t): Dim(dim), Tau(t), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL), fl_empty(false) { }
  Candidate_Ilast_Eall_sphere_1(const Candidate_Ilast_Eall_sphere_1 & candidate);
  ~Candidate_Ilast_Eall_sphere_1();

  unsigned int GetTau()const;
  double Dist(double* a, double*b);
  void CleanOfCandidate();
  bool EmptyOfCandidate();
  void InitialOfCandidate(unsigned int tau, double** &cumsumdata,  double** &cumsumdata2, double* &vectofcosts);
  void UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Candidate_Ilast_Eall_sphere_1>::iterator> &vectlinktocands, unsigned int& RealNbExclus);
};
#endif //CANDIDATE_ILAST_EALL_SPHERE_1_H
//------------------------------------------------------------------------------
