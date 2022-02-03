#include "Rec_Empty_Empty.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

Rec_Empty_Empty::Rec_Empty_Empty(const Rec_Empty_Empty & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData;
  VectOfCosts = candidate.VectOfCosts;
}

Rec_Empty_Empty::~Rec_Empty_Empty() { CumSumData = NULL; CumSumData2 = NULL;  VectOfCosts = NULL; }

unsigned int Rec_Empty_Empty::GetTau()const { return Tau; }

void Rec_Empty_Empty::CleanOfCandidate() { CumSumData = NULL;  CumSumData2 = NULL; VectOfCosts = NULL; fl_empty = false; }

bool Rec_Empty_Empty::EmptyOfCandidate() { return fl_empty; }

void Rec_Empty_Empty::InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts) {
  Tau = tau;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
  fl_empty = false;
}

void Rec_Empty_Empty::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Rec_Empty_Empty>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
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
