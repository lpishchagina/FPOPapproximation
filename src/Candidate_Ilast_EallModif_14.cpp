// Rect^tau_t = (Rect^tau-(LastT-1) interS^Tau_LastT)\ (union_{j=1^Tau-1}S^j_(Tau-1))
#include "Candidate_Ilast_EallModif_14.h"

using namespace Rcpp;
using namespace std;

Candidate_Ilast_EallModif_14::Candidate_Ilast_EallModif_14(const Candidate_Ilast_EallModif_14 & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  Rect= new pRectangle(Dim);
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
  disks_t_1.clear();
  disks_t_1 = candidate.disks_t_1;
  CreationFl = candidate.CreationFl;
}

Candidate_Ilast_EallModif_14::~Candidate_Ilast_EallModif_14() { delete Rect;  CumSumData = NULL;  CumSumData2 = NULL;  VectOfCosts = NULL; }

unsigned int Candidate_Ilast_EallModif_14::GetTau()const { return Tau; }

std::list<pSphere> Candidate_Ilast_EallModif_14::get_disks_t_1()const { return disks_t_1; }

void Candidate_Ilast_EallModif_14::CleanOfCandidate() { CumSumData = NULL;  CumSumData2 = NULL;  VectOfCosts = NULL; disks_t_1.clear(); }

bool Candidate_Ilast_EallModif_14::EmptyOfCandidate() { return Rect -> IsEmpty_rect(); }

void Candidate_Ilast_EallModif_14::InitialOfCandidate(unsigned int tau, double** &cumsumdata,  double** &cumsumdata2, double* &vectofcosts) {
  Tau = tau;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
  disks_t_1.clear();
  CreationFl = false;
}

void Candidate_Ilast_EallModif_14::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Candidate_Ilast_EallModif_14>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
  RealNbExclus = 0;
  //pelt
  Cost cost = Cost(Dim);
  unsigned int LastT = vectlinktocands[vectlinktocands.size() - 1] -> GetTau();
  cost.InitialCost(Dim, Tau, LastT, CumSumData, CumSumData2, VectOfCosts);
  double Radius2 = (VectOfCosts[LastT + 1] - VectOfCosts[Tau] - cost.get_coef_Var())/cost.get_coef();
  if (Radius2 < 0) {
    Rect -> DoEmpty_rect();
    return;
  }
  //
  pSphere Disk = pSphere(Dim);
  if (CreationFl) {
    CreationFl = false;
    if (IndexToLinkOfUpdCand > 0) {
      unsigned int j;
      for (unsigned int i = 0; i < IndexToLinkOfUpdCand; i++) {
        j = vectlinktocands[i] -> GetTau();
        cost.InitialCost(Dim, j, Tau-1, CumSumData, CumSumData2, VectOfCosts);
        Radius2 = (VectOfCosts[Tau] - VectOfCosts[j] - cost.get_coef_Var()) / cost.get_coef();
        Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
        disks_t_1.push_back(Disk);
      }
    }
  }
  //intersection : Rect^Tau_t = Rect^Tau_(t-1) inter S^Tau_t
  cost.InitialCost(Dim, Tau, LastT, CumSumData, CumSumData2, VectOfCosts);
  Radius2 = (VectOfCosts[LastT + 1] - VectOfCosts[Tau] - cost.get_coef_Var())/cost.get_coef();
  Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
  Rect -> Intersection_disk(Disk);
  //exclusion : Rect^Tau_t= Rect^Tau_t \ (union_{j=1^Tau-1}S^j_(Tau-1))
  if ((disks_t_1.size() > 0) && (!Rect -> IsEmpty_rect())) {
    std::list<pSphere>::iterator iter = disks_t_1.begin();
    while(iter != disks_t_1.end() && (!Rect -> IsEmpty_rect())){
      Rect -> Exclusion_disk(*iter);
      RealNbExclus++;
      ++iter;
    }
  }
}
