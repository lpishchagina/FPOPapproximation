#include "Rec_All_lRand.h"
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

Rec_All_lRand::Rec_All_lRand(const Rec_All_lRand & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  Rect = new pRectangle(Dim);
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
  IndexVectBefore.clear();
  IndexVectBefore = candidate.IndexVectBefore;
  CreationFl = candidate.CreationFl;
}

Rec_All_lRand::~Rec_All_lRand() { delete Rect;  CumSumData = NULL;  CumSumData2 = NULL;  VectOfCosts = NULL; }

int Rec_All_lRand::get_Number(int N) {
  srand(time(NULL));
  int res = rand() % N + 1;
  return res;
}

unsigned int Rec_All_lRand::GetTau()const { return Tau; }
void Rec_All_lRand::CleanOfCandidate() { CumSumData = NULL; CumSumData2 = NULL;  VectOfCosts = NULL; }
bool Rec_All_lRand::EmptyOfCandidate() { return Rect -> IsEmpty_rect(); }

void Rec_All_lRand::InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts) {
  Tau = tau;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
  CreationFl = true;
  IndexVectBefore.clear();
}

void Rec_All_lRand::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Rec_All_lRand>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
  RealNbExclus = 0;
  double Radius2;
  pSphere Disk = pSphere(Dim);
  Cost cost = Cost(Dim);
  unsigned int j;
  //intersection :
  for (unsigned int i = IndexToLinkOfUpdCand; i < vectlinktocands.size(); i++) {
    j = vectlinktocands[i] -> GetTau();
    cost.InitialCost(Dim, Tau, j, CumSumData, CumSumData, VectOfCosts);
    Radius2 = (VectOfCosts[j + 1] - VectOfCosts[Tau] - cost.get_coef_Var())/cost.get_coef();

    if (Radius2 < 0) { Rect -> DoEmpty_rect(); return; } //pelt

    Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
    Rect -> Intersection_disk(Disk);
    if (Rect -> IsEmpty_rect()) { return; }
  }

  //CreationFl = true =>1 iteration : Creation of IndexVectBefore
  if (CreationFl) {
    CreationFl = false;
    if (IndexToLinkOfUpdCand > 0) {
      for (unsigned int i = 0; i < IndexToLinkOfUpdCand; i++) {
        IndexVectBefore.push_back(vectlinktocands[i] -> GetTau());
      }
    }
  }
  //exclusion :
  if ((!Rect->IsEmpty_rect()) && (IndexVectBefore.size() > 0)) {
    unsigned int IndexRandBeforeTau = get_Number(IndexVectBefore.size()) - 1;
    j = IndexVectBefore[IndexRandBeforeTau];
    cost.InitialCost(Dim, j, Tau-1, CumSumData, CumSumData2, VectOfCosts);
    Radius2 = (VectOfCosts[Tau] - VectOfCosts[j] - cost.get_coef_Var()) / cost.get_coef();
    Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));

    if (!Rect -> EmptyIntersection(Disk)) {
      Rect -> Exclusion_disk(Disk);
    } else {
      if (IndexRandBeforeTau < (IndexVectBefore.size() -1 )) {
        IndexVectBefore[IndexRandBeforeTau] = IndexVectBefore.back();
      }
      IndexVectBefore.pop_back();
    }
    RealNbExclus = 1;
  }
}

