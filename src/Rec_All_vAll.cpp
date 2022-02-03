#include "Rec_All_vAll.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

Rec_All_vAll::Rec_All_vAll(const Rec_All_vAll & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  Rect= new pRectangle(Dim);
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
}

Rec_All_vAll::~Rec_All_vAll() { delete Rect;  CumSumData = NULL; CumSumData2 = NULL; VectOfCosts = NULL; }

unsigned int Rec_All_vAll::GetTau()const { return Tau; }
void Rec_All_vAll::CleanOfCandidate() { CumSumData = NULL; CumSumData2 = NULL; VectOfCosts = NULL; }
bool Rec_All_vAll::EmptyOfCandidate() { return Rect -> IsEmpty_rect(); }

void Rec_All_vAll::InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts) {
  Tau = tau;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
}

void Rec_All_vAll::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Rec_All_vAll>::iterator> &vectlinktocands, unsigned int &RealNbExclus) {
  RealNbExclus = 0;
  double Radius2;
  pSphere Disk = pSphere(Dim);
  Cost cost = Cost(Dim);
  unsigned int j;
  //intersection :
  for (unsigned int i = IndexToLinkOfUpdCand; i < vectlinktocands.size(); i++) {
    j = vectlinktocands[i] -> GetTau();
    cost.InitialCost(Dim, Tau, j, CumSumData, CumSumData2, VectOfCosts);
    Radius2 = (VectOfCosts[j + 1] - VectOfCosts[Tau] - cost.get_coef_Var())/cost.get_coef();

    if (Radius2 < 0) { Rect -> DoEmpty_rect(); return; }//pelt

    Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));

    Rect -> Intersection_disk(Disk);
    if (Rect -> IsEmpty_rect()) {  return;}
  }
  //exclusion :
  if ((IndexToLinkOfUpdCand > 0) && (!Rect -> IsEmpty_rect())) {
    unsigned int NbEmptyInter = 0;
    for (unsigned int i = 0; i < IndexToLinkOfUpdCand; i++) {
      j = vectlinktocands[i] -> GetTau();
      cost.InitialCost(Dim, j, Tau-1, CumSumData, CumSumData2, VectOfCosts);
      Radius2 = (VectOfCosts[Tau] - VectOfCosts[j] - cost.get_coef_Var()) / cost.get_coef();
      Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
      if (Rect -> EmptyIntersection((Disk))) {
        NbEmptyInter++;
      } else {
        Rect -> Exclusion_disk(Disk);
        if (Rect -> IsEmpty_rect()) {
          RealNbExclus = i + 1 - NbEmptyInter;
          return;
        }
      }
    }
    RealNbExclus = IndexToLinkOfUpdCand - NbEmptyInter;
  }
}
