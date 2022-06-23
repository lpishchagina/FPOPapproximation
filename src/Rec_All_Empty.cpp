#include "Rec_All_Empty.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

//constructor copy, destructor--------------------------------------------------
Rec_All_Empty::Rec_All_Empty(const Rec_All_Empty & candidate) {
  p = candidate.p;
  tau = candidate.tau;
  rectangle = new pRectangle(p);
  csY = candidate.csY;
  csY2 = candidate.csY2;
  locCosts = candidate.locCosts;
}

Rec_All_Empty::~Rec_All_Empty() { delete rectangle;  csY = NULL; csY2 = NULL; locCosts = NULL; }

//accessory---------------------------------------------------------------------
unsigned int Rec_All_Empty::get_tau() const { return tau; }

//tools-------------------------------------------------------------------------
bool Rec_All_Empty::EmptyOfCandidate() { return rectangle -> IsEmptyRect(); }

void Rec_All_Empty::idCandidate(unsigned int dim, unsigned int t, double** &csy,  double** &csy2, double* &loccosts) {
  p = dim;
  tau = t;
  csY = csy;
  csY2 = csy2;
  locCosts = loccosts;
}

void Rec_All_Empty::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Rec_All_Empty>::iterator> &vectlinktocands, unsigned int &RealNbExclus) {
  pSphere sphere = pSphere(p);
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
}
//----------------------------------------------------------------------------//
