// Rect^Tau_t = RectApprox(S^Tau_Tau) \ (union_{j=1^Tau-1}S^j_(Tau-1))?????
#include "Candidate_Iempty_Eall_4.h"

using namespace Rcpp;
using namespace std;

Candidate_Iempty_Eall_4::Candidate_Iempty_Eall_4(const Candidate_Iempty_Eall_4 & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  Rect= new pRectangle(Dim);
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
}

Candidate_Iempty_Eall_4::~Candidate_Iempty_Eall_4() { delete Rect;  CumSumData = NULL; CumSumData2 = NULL; VectOfCosts = NULL; }

unsigned int Candidate_Iempty_Eall_4::GetTau()const { return Tau; }

void Candidate_Iempty_Eall_4::CleanOfCandidate() { CumSumData = NULL; CumSumData2 = NULL; VectOfCosts = NULL; }

bool Candidate_Iempty_Eall_4::EmptyOfCandidate() { return Rect -> IsEmpty_rect(); }

void Candidate_Iempty_Eall_4::InitialOfCandidate(unsigned int t, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts) {
  Tau = t;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
}

void Candidate_Iempty_Eall_4::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Candidate_Iempty_Eall_4>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
  RealNbExclus = 0;
  if (IndexToLinkOfUpdCand == vectlinktocands.size() - 1 ) {
    //pelt
    Cost cost = Cost(Dim);
    cost.InitialCost(Dim, Tau, Tau, CumSumData, CumSumData2, VectOfCosts);
    double Radius2 = (VectOfCosts[Tau + 1] - VectOfCosts[Tau] - cost.get_coef_Var())/cost.get_coef();
    if (Radius2 < 0) {
      Rect -> DoEmpty_rect();
      return;
    }
    //Rect^Tau_t = RectApprox(S^Tau_Tau)
    pSphere Disk = pSphere(Dim);
    Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
    Rect -> CubeApproximation(Disk);
    //exclusion : Rect^Tau_t = Rect^Tau_t \ (union_{j=1^Tau-1}S^j_(Tau-1))
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
}
