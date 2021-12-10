#ifndef CANDIDATE_IRANDOM_ERANDOM_SPHERE_12_H
#define CANDIDATE_IRANDOM_ERANDOM_SPHERE_12_H
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

class Candidate_Irandom_Erandom_sphere_12 {
private:
  unsigned int Dim;
  unsigned int Tau;
  double** CumSumData;
  double** CumSumData2;
  double* VectOfCosts;
  bool fl_empty;

public:
  Candidate_Irandom_Erandom_sphere_12(): Dim(0), Tau(0), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL),  fl_empty(false) { }
  Candidate_Irandom_Erandom_sphere_12(unsigned  int dim): Dim(dim), Tau(0), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL), fl_empty(false) { }
  Candidate_Irandom_Erandom_sphere_12(unsigned int dim, unsigned int t): Dim(dim), Tau(t), CumSumData(NULL), CumSumData2(NULL), VectOfCosts(NULL), fl_empty(false) { }
  Candidate_Irandom_Erandom_sphere_12(const Candidate_Irandom_Erandom_sphere_12 & candidate);
  ~Candidate_Irandom_Erandom_sphere_12();

  unsigned int GetTau()const;
  int get_Number(int N);
  double Dist(double* a, double*b);
  void CleanOfCandidate();
  bool EmptyOfCandidate();
  void InitialOfCandidate(unsigned int tau, double** &cumsumdata,  double** &cumsumdata2, double* &vectofcosts);
  void UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Candidate_Irandom_Erandom_sphere_12>::iterator> &vectlinktocands, unsigned int& RealNbExclus);
};
#endif //CANDIDATE_IRANDOM_ERANDOM_SPHERE_12_H
//------------------------------------------------------------------------------
