// Rect^tau_t = Cube(S^Tau_LastT)\ (union_{j=1^Tau-1}S^j_(Tau-1))
#include "Candidate_Ilast_Eall_5.h"

using namespace Rcpp;
using namespace std;

Candidate_Ilast_Eall_5::Candidate_Ilast_Eall_5(const Candidate_Ilast_Eall_5 & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  Rect= new pRectangle(Dim);
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
}

Candidate_Ilast_Eall_5::~Candidate_Ilast_Eall_5() { delete Rect;  CumSumData = NULL;  CumSumData2 = NULL;  VectOfCosts = NULL; }

unsigned int Candidate_Ilast_Eall_5::GetTau()const { return Tau; }

void Candidate_Ilast_Eall_5::CleanOfCandidate() { CumSumData = NULL;  CumSumData2 = NULL;  VectOfCosts = NULL; }

bool Candidate_Ilast_Eall_5::EmptyOfCandidate() { return Rect -> IsEmpty_rect(); }

void Candidate_Ilast_Eall_5::InitialOfCandidate(unsigned int tau, double** &cumsumdata,  double** &cumsumdata2, double* &vectofcosts) {
  Tau = tau;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
}

void Candidate_Ilast_Eall_5::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Candidate_Ilast_Eall_5>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
  RealNbExclus = 0;
  //pelt
  Cost cost = Cost(Dim);
  unsigned int tLast = vectlinktocands[vectlinktocands.size() -1] -> GetTau();
  cost.InitialCost(Dim, Tau, tLast, CumSumData, CumSumData2, VectOfCosts);
  double Radius2 = (VectOfCosts[tLast + 1] - VectOfCosts[Tau] - cost.get_coef_Var())/cost.get_coef();
  if (Radius2 < 0) {
    Rect -> DoEmpty_rect();
    return;
  }
  //last : Rect^Tau_t = Cube(S^Tau_LastT)
  pSphere Disk = pSphere(Dim);
  Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
  Rect -> CubeApproximation(Disk);
  //exclusion : Rect^Tau_t= Rect^Tau_t \ (union_{j=1^Tau-1}S^j_(Tau-1))
  if (IndexToLinkOfUpdCand > 0) {
    unsigned int j;
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
