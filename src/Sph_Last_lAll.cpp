#include "Sph_Last_lAll.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

//constructor copy, destructor--------------------------------------------------
Sph_Last_lAll::Sph_Last_lAll(const Sph_Last_lAll & candidate) {
  p = candidate.p;
  tau = candidate.tau;
  csY = candidate.csY;
  csY2 = candidate.csY2;
  locCosts = candidate.locCosts;
  fl_empty = candidate.fl_empty;
  spheresBefore.clear();
  spheresBefore = candidate.spheresBefore;
  flCreate = candidate.flCreate;
}

Sph_Last_lAll::~Sph_Last_lAll() { csY = NULL; csY2 = NULL;  locCosts = NULL; }

//accessory---------------------------------------------------------------------
unsigned int Sph_Last_lAll::get_tau() const { return tau; }

//tools-------------------------------------------------------------------------
double Sph_Last_lAll::get_dist(double* pnt1, double* pnt2) {
  double res = 0;
  for (unsigned int k = 0; k < p; k++) {
    res = res + (pnt2[k] - pnt1[k]) * (pnt2[k] - pnt1[k]);
  }
  return sqrt(res);
}

bool Sph_Last_lAll::EmptyOfCandidate() { return fl_empty; }

void Sph_Last_lAll::idCandidate(unsigned  int dim, unsigned int t, double** &csy, double** &csy2, double* &loccosts) {
  p = dim;
  tau = t;
  csY = csy;
  csY2 = csy2;
  locCosts = loccosts;
  fl_empty = false;
  flCreate = true;
  spheresBefore.clear();
}

void Sph_Last_lAll::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Sph_Last_lAll>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
  std::list<pSphere>::iterator iter;
  pSphere sphere = pSphere(p);
  pSphere testSphere = pSphere(p);
  //exclusion set
  if (flCreate) { //flCreate = true =>1 iteration : Creation of spheresBefore
    flCreate = false;
    if (IndexToLinkOfUpdCand > 0) {
      for (unsigned int i = 0; i < IndexToLinkOfUpdCand; i++) {
        testSphere.createSphere(p, vectlinktocands[i] -> get_tau(), tau-1, csY, csY2, locCosts);
        spheresBefore.push_back(testSphere);
      }
    }
  }
  sphere.createSphere(p, tau, vectlinktocands[vectlinktocands.size() - 1] -> get_tau(), csY, csY2, locCosts);
  if (sphere.get_r() == 0) {
    fl_empty = true;
    return;
  }//pelt
  //spheres (check exclusion)
  if (spheresBefore.size() > 0)   {
    iter = spheresBefore.begin();
    while (iter != spheresBefore.end()) {
      if (sphere.isInclusion(*iter)) {
        fl_empty = true;
        return;
      } else if (sphere.isnotIntersection(*iter)) {
        iter = spheresBefore.erase(iter);
      } else {
        ++iter;
      }
    }
  }
  //spheres (check intersection)
  unsigned int end = vectlinktocands.size() - 1;
  for (int i = IndexToLinkOfUpdCand; i < end; i++) {
      testSphere.createSphere(p, tau, vectlinktocands[i] -> get_tau(), csY, csY2, locCosts);
      if ((testSphere.get_r()  == 0) || (sphere.isnotIntersection(testSphere))) {
        rectangle -> DoEmptyRect();
        return;
      }
  }
}
//----------------------------------------------------------------------------//
