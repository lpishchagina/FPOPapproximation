#include "Rec_Last_lRand.h"
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

Rec_Last_lRand::Rec_Last_lRand(const Rec_Last_lRand & candidate) {
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

Rec_Last_lRand::~Rec_Last_lRand() { delete Rect;  CumSumData = NULL;  CumSumData2 = NULL;  VectOfCosts = NULL; }

int Rec_Last_lRand::get_Number(int N) {
  srand(time(NULL));
  int res = rand() % N + 1;
  return res;
}

unsigned int Rec_Last_lRand::GetTau()const { return Tau; }
void Rec_Last_lRand::CleanOfCandidate() { CumSumData = NULL; CumSumData2 = NULL;  VectOfCosts = NULL; }
bool Rec_Last_lRand::EmptyOfCandidate() { return Rect -> IsEmpty_rect(); }

void Rec_Last_lRand::InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts) {
  Tau = tau;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
  CreationFl = true;
  IndexVectBefore.clear();
}

void Rec_Last_lRand::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Rec_Last_lRand>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
  RealNbExclus = 0;
  Cost cost = Cost(Dim);
  unsigned int j = vectlinktocands[vectlinktocands.size() - 1] -> GetTau();
  cost.InitialCost(Dim, Tau, j, CumSumData, CumSumData2, VectOfCosts);
  double Radius2 = (VectOfCosts[j + 1] - VectOfCosts[Tau] - cost.get_coef_Var())/cost.get_coef();

  if (Radius2 < 0) { Rect -> DoEmpty_rect(); return; }//pelt

  //intersection :
  pSphere Disk = pSphere(Dim);
  Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
  Rect -> Intersection_disk(Disk);

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

