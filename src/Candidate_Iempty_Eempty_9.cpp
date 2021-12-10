#include "pRectangle.h"
#include "Candidate_Iempty_Eempty_9.h"
#include "Cost.h"
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

Candidate_Iempty_Eempty_9::Candidate_Iempty_Eempty_9(const Candidate_Iempty_Eempty_9 & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData;
  VectOfCosts = candidate.VectOfCosts;
}

Candidate_Iempty_Eempty_9::~Candidate_Iempty_Eempty_9() { CumSumData = NULL; CumSumData2 = NULL;  VectOfCosts = NULL; }

unsigned int Candidate_Iempty_Eempty_9::GetTau()const { return Tau; }

void Candidate_Iempty_Eempty_9::CleanOfCandidate() { CumSumData = NULL;  CumSumData2 = NULL; VectOfCosts = NULL; fl_empty = false; }

bool Candidate_Iempty_Eempty_9::EmptyOfCandidate() { return fl_empty; }

void Candidate_Iempty_Eempty_9::InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts) {
  Tau = tau;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
  fl_empty = false;
}

void Candidate_Iempty_Eempty_9::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Candidate_Iempty_Eempty_9>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
  fl_empty = false;
  //pelt
  Cost cost = Cost(Dim);
  unsigned int LastT = vectlinktocands[vectlinktocands.size()-1] -> GetTau();
  cost.InitialCost(Dim, Tau, LastT, CumSumData, CumSumData2, VectOfCosts);
  double Radius2 = (VectOfCosts[LastT + 1] - VectOfCosts[Tau] - cost.get_coef_Var())/cost.get_coef();
  if (Radius2 < 0) {
    fl_empty = true;
  }
}
