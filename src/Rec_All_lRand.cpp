#include "Rec_All_lRand.h"
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

//constructor copy, destructor--------------------------------------------------
Rec_All_lRand::Rec_All_lRand(const Rec_All_lRand & candidate) {
  p = candidate.p;
  tau = candidate.tau;
  rectangle = new pRectangle(p);
  csY = candidate.csY;
  csY2 = candidate.csY2;
  locCosts = candidate.locCosts;
  indexSpheresBefore = candidate.indexSpheresBefore;
  flCreate = candidate.flCreate;
}

Rec_All_lRand::~Rec_All_lRand() { delete rectangle;  csY = NULL; csY2 = NULL; locCosts = NULL; }

//accessory---------------------------------------------------------------------
unsigned int Rec_All_lRand::get_tau() const { return tau; }

//tools-------------------------------------------------------------------------
double Rec_All_lRand::get_dist(double* pnt1, double* pnt2) {
  double res = 0;
  for (unsigned int k = 0; k < p; k++) {
    res = res + (pnt2[k] - pnt1[k]) * (pnt2[k] - pnt1[k]);
  }
  return sqrt(res);
}

int Rec_All_lRand::get_Number(int N) {
  srand(time(NULL));
  int res = rand() % N + 1;
  return res;
}

bool Rec_All_lRand::EmptyOfCandidate() { return rectangle -> IsEmptyRect(); }

void Rec_All_lRand::idCandidate(unsigned int dim, unsigned int t, double** &csy, double** &csy2, double* &loccosts) {
  p = dim;
  tau = t;
  csY = csy;
  csY2 = csy2;
  locCosts = loccosts;
  indexSpheresBefore.clear();
  flCreate = true;
}

void Rec_All_lRand::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Rec_All_lRand>::iterator> &vectlinktocands, unsigned int &RealNbExclus) {
  pSphere sphere = pSphere(p);
  //labels of the elements from exclusion set
  if (flCreate) {//flCreate = true =>1 iteration : create labels of the elements from exclusion set
    flCreate = false;
    if (IndexToLinkOfUpdCand > 0) {
      for (unsigned int i = 0; i < IndexToLinkOfUpdCand; i++) {
        indexSpheresBefore.push_back(vectlinktocands[i] -> get_tau());
      }
    }
  }
  //intersection approximation    //intersection set:
  for (int i = IndexToLinkOfUpdCand; i < vectlinktocands.size(); i++) {
    sphere.createSphere(p, tau, vectlinktocands[i] -> get_tau(), csY, csY2, locCosts);
    if (sphere.get_r()  == 0) {
      rectangle -> DoEmptyRect();
      return;
    }     //pelt
    rectangle -> IntersectionSphere(sphere);
    if (rectangle -> IsEmptyRect()) {
      return;
    }
  }
  //exclusion approximation with one random sphere from exclusion set :
  if ((!rectangle->IsEmptyRect()) && (indexSpheresBefore.size() > 0)) {
    unsigned int IndexRandBeforeTau = get_Number(indexSpheresBefore.size()) - 1;
    sphere.createSphere(p, indexSpheresBefore[IndexRandBeforeTau],  tau-1, csY, csY2, locCosts);
    if (!rectangle -> EmptyIntersection(sphere)) {
      rectangle -> ExclusionSphere(sphere);
    } else {
      if (IndexRandBeforeTau < (indexSpheresBefore.size() - 1 )) {
        indexSpheresBefore[IndexRandBeforeTau] = indexSpheresBefore.back();
      }
      indexSpheresBefore.pop_back();
    }
  }
}
//----------------------------------------------------------------------------//
