#include "Rec_Last_vDoubleAll.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

Rec_Last_vDoubleAll::Rec_Last_vDoubleAll(const Rec_Last_vDoubleAll & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  Rect= new pRectangle(Dim);
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
}

Rec_Last_vDoubleAll::~Rec_Last_vDoubleAll() { delete Rect;  CumSumData = NULL;  CumSumData2 = NULL;  VectOfCosts = NULL; }

unsigned int Rec_Last_vDoubleAll::GetTau()const { return Tau; }
void Rec_Last_vDoubleAll::CleanOfCandidate() { CumSumData = NULL;  CumSumData2 = NULL;  VectOfCosts = NULL; }
bool Rec_Last_vDoubleAll::EmptyOfCandidate() { return Rect -> IsEmpty_rect(); }

void Rec_Last_vDoubleAll::InitialOfCandidate(unsigned int tau, double** &cumsumdata,  double** &cumsumdata2, double* &vectofcosts) {
  Tau = tau;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
}

void Rec_Last_vDoubleAll::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Rec_Last_vDoubleAll>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
  RealNbExclus = 0;
  Cost cost = Cost(Dim);
  unsigned int j = vectlinktocands[vectlinktocands.size() - 1] -> GetTau();
  cost.InitialCost(Dim, Tau, j, CumSumData, CumSumData2, VectOfCosts);
  double Radius2 = (VectOfCosts[j + 1] - VectOfCosts[Tau] - cost.get_coef_Var())/cost.get_coef();

  if (Radius2 < 0) { Rect -> DoEmpty_rect(); return; }   //pelt

  //intersection :
  pSphere Disk = pSphere(Dim);
  Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
  Rect -> Intersection_disk(Disk);
  //exclusion :
  if ((IndexToLinkOfUpdCand > 0) && (!Rect -> IsEmpty_rect())) {
    unsigned int NbEmptyInter = 0;
    unsigned int NbLoops = 2;
    for (unsigned int k = 0; k < NbLoops; k++){
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
}
