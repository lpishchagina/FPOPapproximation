#include "Rec_All_lAll.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

Rec_All_lAll::Rec_All_lAll(const Rec_All_lAll & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  Rect= new pRectangle(Dim);
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
  DiskListBefore.clear();
  DiskListBefore = candidate.DiskListBefore;
  CreationFl = candidate.CreationFl;
}

Rec_All_lAll::~Rec_All_lAll() { delete Rect;  CumSumData = NULL; CumSumData2 = NULL; VectOfCosts = NULL; }

std::list<pSphere> Rec_All_lAll::GetDiskListBefore()const { return DiskListBefore; }
unsigned int Rec_All_lAll::GetTau()const { return Tau; }

void Rec_All_lAll::CleanOfCandidate() { CumSumData = NULL; CumSumData2 = NULL; VectOfCosts = NULL;  DiskListBefore.clear(); }
bool Rec_All_lAll::EmptyOfCandidate() { return Rect -> IsEmpty_rect(); }

void Rec_All_lAll::InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts) {
  Tau = tau;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
  DiskListBefore.clear();
  CreationFl = true;
}

double Rec_All_lAll::Dist(double* a, double*b) {
  double dist = 0;
  for (unsigned int k = 0; k < Dim; k++) {
    dist = dist + (a[k] - b[k])*(a[k] - b[k]);
  }
  return sqrt(dist);
}

void Rec_All_lAll::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Rec_All_lAll>::iterator> &vectlinktocands, unsigned int &RealNbExclus) {
  RealNbExclus = 0;
  Cost cost = Cost(Dim);
  unsigned int j = vectlinktocands[vectlinktocands.size() - 1] -> GetTau();
  cost.InitialCost(Dim, Tau, j, CumSumData, CumSumData2, VectOfCosts);
  double Radius2 = (VectOfCosts[j + 1] - VectOfCosts[Tau] - cost.get_coef_Var())/cost.get_coef();

  if (Radius2 < 0) { Rect -> DoEmpty_rect(); return;}   //pelt

  pSphere Disk = pSphere(Dim);
  Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
  //CreationFl = true =>1 iteration : Creation of DiskListBefore
  if (CreationFl) {
    CreationFl = false;
    if (IndexToLinkOfUpdCand > 0) {
      double dist;
      pSphere DiskTau_1 = pSphere(Dim);
      for (unsigned int i = 0; i < IndexToLinkOfUpdCand; i++) {
        j = vectlinktocands[i] -> GetTau();
        cost.InitialCost(Dim, j, Tau-1, CumSumData, CumSumData2, VectOfCosts);
        Radius2 = (VectOfCosts[Tau] - VectOfCosts[j] - cost.get_coef_Var()) / cost.get_coef();
        DiskTau_1.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
        //check for 1 iteration: inclusion
        dist = Dist(Disk.get_center(), DiskTau_1.get_center());
        if (dist < (DiskTau_1.get_radius() + Disk.get_radius())) {
          if (dist <= (DiskTau_1.get_radius() - Disk.get_radius())) { Rect -> DoEmpty_rect(); return; }
        }
        DiskListBefore.push_back(DiskTau_1);
      }
    }
  }
  //intersection :
  for (unsigned int i = IndexToLinkOfUpdCand; i < vectlinktocands.size(); i++) {
    j = vectlinktocands[i] -> GetTau();
    cost.InitialCost(Dim, Tau, j, CumSumData, CumSumData2, VectOfCosts);
    Radius2 = (VectOfCosts[j + 1] - VectOfCosts[Tau] - cost.get_coef_Var())/cost.get_coef();

    if (Radius2 < 0) { Rect -> DoEmpty_rect(); return; }     //pelt

    Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));

    Rect -> Intersection_disk(Disk);
    if (Rect -> IsEmpty_rect()) { return; }
  }
  //exclusion :
  if ((DiskListBefore.size() > 0) && (!Rect -> IsEmpty_rect())) {
    std::list<pSphere>::iterator iter = DiskListBefore.begin();
    while ( (iter != DiskListBefore.end()) && (!Rect -> IsEmpty_rect())) {
      if (Rect -> EmptyIntersection(*iter)) {
        iter = DiskListBefore.erase(iter);
      } else {
        Rect -> Exclusion_disk(*iter);
        ++iter;
        ++RealNbExclus;
      }
    }
  }
}

