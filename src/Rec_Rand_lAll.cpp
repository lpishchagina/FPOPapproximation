#include "Rec_Rand_lAll.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

//constructor copy, destructor--------------------------------------------------
Rec_Rand_lAll::Rec_Rand_lAll(const Rec_Rand_lAll & candidate) {
  p = candidate.p;
  tau = candidate.tau;
  rectangle = new pRectangle(p);
  csY = candidate.csY;
  csY2 = candidate.csY2;
  locCosts = candidate.locCosts;
  spheresBefore.clear();
  spheresBefore = candidate.spheresBefore;
  flCreate = candidate.flCreate;
}

Rec_Rand_lAll::~Rec_Rand_lAll() { delete rectangle;  csY = NULL; csY2 = NULL; locCosts = NULL; }

//accessory---------------------------------------------------------------------
std::list<pSphere> Rec_Rand_lAll::get_spheresBefore() const { return spheresBefore; }
unsigned int Rec_Rand_lAll::get_tau() const { return tau; }

//tools-------------------------------------------------------------------------
double Rec_Rand_lAll::get_dist(double* pnt1, double* pnt2) {
  double res = 0;
  for (unsigned int k = 0; k < p; k++) {
    res = res + (pnt2[k] - pnt1[k]) * (pnt2[k] - pnt1[k]);
  }
  return sqrt(res);
}

int Rec_Rand_lAll::get_Number(int N) {
  srand(time(NULL));
  int res = rand() % N + 1;
  return res;
}

bool Rec_Rand_lAll::EmptyOfCandidate() { return rectangle -> IsEmptyRect(); }

void Rec_Rand_lAll::idCandidate(unsigned int dim, unsigned int t, double** &csy, double** &csy2, double* &loccosts) {
  p = dim;
  tau = t;
  csY = csy;
  csY2 = csy2;
  locCosts = loccosts;
  spheresBefore.clear();
  flCreate = true;
}


void Rec_Rand_lAll::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Rec_Rand_lAll>::iterator> &vectlinktocands, unsigned int &RealNbExclus) {
  std::list<pSphere>::iterator iter;
  pSphere sphere = pSphere(p);
  //exclusion set
  if (flCreate) {//flCreate = true =>1 iteration : Creation of spheresBefore
    flCreate = false;
    if (IndexToLinkOfUpdCand > 0) {
      for (unsigned int i = 0; i < IndexToLinkOfUpdCand; i++) {
        sphere.createSphere(p, vectlinktocands[i] -> get_tau(), tau-1, csY, csY2, locCosts);
        spheresBefore.push_back(sphere);
      }
    }
  }
  //random sphere from intersection set:
  unsigned int IndexRandAfterTau = get_Number (vectlinktocands.size() - IndexToLinkOfUpdCand) + IndexToLinkOfUpdCand - 1;
  sphere.createSphere(p, tau, vectlinktocands[IndexRandAfterTau] -> get_tau(), csY, csY2, locCosts);
  if (sphere.get_r() == 0) { rectangle -> DoEmptyRect(); return;}   //pelt
  //intersection approximation with random sphere from intersection set
  rectangle -> IntersectionSphere(sphere);
  if (rectangle -> IsEmptyRect()) { return; }
  //exclusion approximation:
  if ((spheresBefore.size() > 0) && (!rectangle -> IsEmptyRect())) {
    iter = spheresBefore.begin();
    while ( (iter != spheresBefore.end()) && (!rectangle -> IsEmptyRect())) {
      if (rectangle -> EmptyIntersection(*iter)) {
        iter = spheresBefore.erase(iter);
      } else {
        rectangle -> ExclusionSphere(*iter);
        ++iter;
      }
    }
  }
}
//----------------------------------------------------------------------------//
