#include "Rec_Last_lRand.h"
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

//constructor copy, destructor--------------------------------------------------
Rec_Last_lRand::Rec_Last_lRand(const Rec_Last_lRand & candidate) {
  p = candidate.p;
  tau = candidate.tau;
  rectangle = new pRectangle(p);
  csY = candidate.csY;
  csY2 = candidate.csY2;
  locCosts = candidate.locCosts;

  indexSpheresBefore = candidate.indexSpheresBefore;
  flCreate = candidate.flCreate;
}

Rec_Last_lRand::~Rec_Last_lRand() { delete rectangle;  csY = NULL; csY2 = NULL; locCosts = NULL; }

//accessory---------------------------------------------------------------------
unsigned int Rec_Last_lRand::get_tau() const { return tau; }

//tools-------------------------------------------------------------------------
double Rec_Last_lRand::get_dist(double* pnt1, double* pnt2) {
  double res = 0;
  for (unsigned int k = 0; k < p; k++) {
    res = res + (pnt2[k] - pnt1[k]) * (pnt2[k] - pnt1[k]);
  }
  return sqrt(res);
}

int Rec_Last_lRand::get_Number(int N) {
  srand(time(NULL));
  int res = rand() % N + 1;
  return res;
}

bool Rec_Last_lRand::EmptyOfCandidate() { return rectangle -> IsEmptyRect(); }

void Rec_Last_lRand::idCandidate(unsigned int dim, unsigned int t, double** &csy, double** &csy2, double* &loccosts) {
  p = dim;
  tau = t;
  csY = csy;
  csY2 = csy2;
  locCosts = loccosts;
  indexSpheresBefore.clear();
  flCreate = true;
}

void Rec_Last_lRand::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Rec_Last_lRand>::iterator> &vectlinktocands, unsigned int &RealNbExclus) {
  std::list<pSphere> spheresAfter;
  typename std::list<pSphere>::reverse_iterator riter;
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
  //last sphere from intersection set:
  sphere.createSphere(p, tau, vectlinktocands[vectlinktocands.size() - 1] -> get_tau(), csY, csY2, locCosts);//last
  if (sphere.get_r() == 0) {
    rectangle -> DoEmptyRect();
    return;
  }   //pelt
  //intersection approximation with last sphere from intersection set
  rectangle -> IntersectionSphere(sphere);
  if (rectangle -> IsEmptyRect()) { return; }
  //exclusion approximation with one random sphere from exclusion set :
  if ((!rectangle->IsEmptyRect()) && (indexSpheresBefore.size() > 0)) {
    unsigned int IndexRandBeforeTau = get_Number(indexSpheresBefore.size()) - 1;
    sphere.createSphere(p, indexSpheresBefore[IndexRandBeforeTau],  tau-1, csY, csY2, locCosts);
    if (!rectangle -> EmptyIntersection(sphere)) {
      rectangle -> ExclusionSphere(sphere);
    } else {
      if (IndexRandBeforeTau < (indexSpheresBefore.size() -1 )) {
        indexSpheresBefore[IndexRandBeforeTau] = indexSpheresBefore.back();
      }
      indexSpheresBefore.pop_back();
    }
  }
}
//----------------------------------------------------------------------------//
