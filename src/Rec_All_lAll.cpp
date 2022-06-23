#include "Rec_All_lAll.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

//constructor copy, destructor--------------------------------------------------
Rec_All_lAll::Rec_All_lAll(const Rec_All_lAll & candidate) {
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

Rec_All_lAll::~Rec_All_lAll() { delete rectangle;  csY = NULL; csY2 = NULL; locCosts = NULL; }

//accessory---------------------------------------------------------------------
std::list<pSphere> Rec_All_lAll::get_spheresBefore() const { return spheresBefore; }
unsigned int Rec_All_lAll::get_tau() const { return tau; }

//tools-------------------------------------------------------------------------
double Rec_All_lAll::get_dist(double* pnt1, double* pnt2) {
  double res = 0;
  for (unsigned int k = 0; k < p; k++) {
    res = res + (pnt2[k] - pnt1[k]) * (pnt2[k] - pnt1[k]);
  }
  return sqrt(res);
}

bool Rec_All_lAll::EmptyOfCandidate() { return rectangle -> IsEmptyRect(); }

void Rec_All_lAll::idCandidate(unsigned int dim, unsigned int t, double** &csy, double** &csy2, double* &loccosts) {
  p = dim;
  tau = t;
  csY = csy;
  csY2 = csy2;
  locCosts = loccosts;
  spheresBefore.clear();
  flCreate = true;
}


void Rec_All_lAll::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Rec_All_lAll>::iterator> &vectlinktocands, unsigned int &RealNbExclus) {
  std::list<pSphere> spheresAfter;
  std::list<pSphere>::iterator iter;
  typename std::list<pSphere>::reverse_iterator riter;
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
  //intersection set:
  for (unsigned int i = IndexToLinkOfUpdCand; i < vectlinktocands.size(); i++) {
    sphere.createSphere(p, tau, vectlinktocands[i] -> get_tau(), csY, csY2, locCosts);
    if (sphere.get_r()  == 0) {
      rectangle -> DoEmptyRect();
      return;
    }     //pelt
    spheresAfter.push_back(sphere);
  }
  //intersection approximation
  riter = spheresAfter.rbegin();
  while (riter != (spheresAfter.rend())) {
    rectangle -> IntersectionSphere(*riter);
    if (rectangle -> IsEmptyRect()) {
      return;
    }
    ++riter;
  }
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
