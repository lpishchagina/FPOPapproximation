#ifndef CANDIDATE_SPHERE_SPHERE_1_H
#define CANDIDATE_SPHERE_SPHERE_1_H
#include <vector>
#include <iostream>
#include <list>
#include <iterator>

#include "pSphere.h"
#include "Cost.h"
class Candidate_sphere_sphere_1{
private:
  unsigned int Dim;
  unsigned int Tau;
  double** CumSumData;
  double* VectOfCosts;
  bool fl_empty;

public:
  Candidate_sphere_sphere_1(): Dim(0), Tau(0), CumSumData(NULL), VectOfCosts(NULL),  fl_empty(false){}
  Candidate_sphere_sphere_1(unsigned  int dim): Dim(dim), Tau(0), CumSumData(NULL), VectOfCosts(NULL), fl_empty(false){}
  Candidate_sphere_sphere_1(unsigned int dim, unsigned int t): Dim(dim), Tau(t), CumSumData(NULL), VectOfCosts(NULL), fl_empty(false){}
  Candidate_sphere_sphere_1(const Candidate_sphere_sphere_1 & candidate);
  ~Candidate_sphere_sphere_1();

  unsigned int GetTau()const;
  double Dist(double* a, double*b);

  void CleanOfCandidate();
  bool EmptyOfCandidate();
  void InitialOfCandidate(unsigned int t, double** &cumsumdata, double* &vectofcosts);
  void UpdateOfCandidate(unsigned int t, std::vector<std::list<Candidate_sphere_sphere_1>::iterator> &vectlinktocands);
};
#endif //CANDIDATE_SPHERE_SPHERE_1_H
//------------------------------------------------------------------------------
