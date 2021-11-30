#include "Candidate_sphere_sphere_1.h"

using namespace Rcpp;
using namespace std;

Candidate_sphere_sphere_1::Candidate_sphere_sphere_1(const Candidate_sphere_sphere_1 & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  CumSumData = candidate.CumSumData;
  VectOfCosts = candidate.VectOfCosts;
}

Candidate_sphere_sphere_1::~Candidate_sphere_sphere_1(){ CumSumData = NULL;  VectOfCosts = NULL; }

unsigned int Candidate_sphere_sphere_1::GetTau()const { return Tau; }

double Candidate_sphere_sphere_1::Dist(double* a, double*b) {
  double dist = 0;
  for (unsigned int k = 0; k < Dim; k++) {
    dist = dist + (a[k] - b[k])*(a[k] - b[k]);
  }
  return sqrt(dist);
}

void Candidate_sphere_sphere_1::Clean_fl_empty(){
  fl_empty = false;
}

void Candidate_sphere_sphere_1::CleanOfCandidate() { CumSumData = NULL;  VectOfCosts = NULL; }

bool Candidate_sphere_sphere_1::EmptyOfCandidate() { return fl_empty; }

void Candidate_sphere_sphere_1::InitialOfCandidate(unsigned int t, double** &cumsumdata, double* &vectofcosts) {
  Tau = t;
  CumSumData = cumsumdata;
  VectOfCosts = vectofcosts;
  fl_empty = false;
}

void Candidate_sphere_sphere_1::UpdateOfCandidate(unsigned int i, std::vector<std::list<Candidate_sphere_sphere_1>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
  Clean_fl_empty();
  unsigned int N = vectlinktocands.size();
  unsigned int t = vectlinktocands[N-1] -> GetTau();
  unsigned int u;
  double r2;
  Cost cost = Cost(Dim);
  cost.InitialCost(Dim, Tau, t, CumSumData[Tau], CumSumData[t + 1], VectOfCosts[Tau]);
  r2 = (VectOfCosts[t + 1] - VectOfCosts[Tau] - cost.get_coef_Var()) / cost.get_coef();
  if (r2 < 0){
    fl_empty = true;
    return;
  }
  pSphere Disk_Taut = pSphere(Dim);
  Disk_Taut.InitialpSphere(Dim, cost.get_mu(), sqrt(r2));
  if (i > 0)  {
    pSphere Disk_ut1 = pSphere(Dim);
    double dist;
    for (unsigned int j = 0; j < i; j++) {
      u = vectlinktocands[j] -> GetTau();
      cost.InitialCost(Dim, u, t-1, CumSumData[u], CumSumData[t], VectOfCosts[u]);
      r2 = (VectOfCosts[t] - VectOfCosts[u] - cost.get_coef_Var()) / cost.get_coef();
      Disk_ut1.InitialpSphere(Dim, cost.get_mu(), sqrt(r2));
      dist = Dist(Disk_Taut.get_center(), Disk_ut1.get_center());
      if (dist < (Disk_ut1.get_radius() + Disk_Taut.get_radius())) {
        if (dist <= (Disk_ut1.get_radius() - Disk_Taut.get_radius())) {
          fl_empty = true;
          return;
        }
      }
    }
  }
}
