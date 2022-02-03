#include "Rec_Rand_vAll.h"
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

Rec_Rand_vAll::Rec_Rand_vAll(const  Rec_Rand_vAll & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  Rect = new pRectangle(Dim);
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
}

Rec_Rand_vAll::~Rec_Rand_vAll() { delete Rect;  CumSumData = NULL;  CumSumData2 = NULL; VectOfCosts = NULL; }

unsigned int Rec_Rand_vAll::GetTau()const { return Tau; }
void  Rec_Rand_vAll::CleanOfCandidate() { CumSumData = NULL; CumSumData2 = NULL; VectOfCosts = NULL; }
bool  Rec_Rand_vAll::EmptyOfCandidate() { return Rect -> IsEmpty_rect(); }

int Rec_Rand_vAll::get_Number(int N) {
  srand(time(NULL));
  int res = rand() % N + 1;
  return res;
}

void Rec_Rand_vAll::InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts) {
  Tau = tau;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
}

void  Rec_Rand_vAll::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Rec_Rand_vAll>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
  RealNbExclus = 0;
  //random after Tau:
  unsigned int IndexRandAfterTau = get_Number (vectlinktocands.size() - IndexToLinkOfUpdCand) + IndexToLinkOfUpdCand - 1;
  unsigned int j = vectlinktocands[IndexRandAfterTau] -> GetTau();
  Cost cost = Cost(Dim);
  cost.InitialCost(Dim, Tau, j, CumSumData, CumSumData2, VectOfCosts);
  double Radius2 = (VectOfCosts[j+ 1] - VectOfCosts[Tau] - cost.get_coef_Var())/cost.get_coef();

  if (Radius2 < 0) { Rect -> DoEmpty_rect(); return; }   //pelt

  pSphere Disk = pSphere(Dim);
  Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
  //intersection :
  Rect -> Intersection_disk(Disk);
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
