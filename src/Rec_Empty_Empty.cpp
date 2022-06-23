#include "Rec_Empty_Empty.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

//constructor copy, destructor--------------------------------------------------
Rec_Empty_Empty::Rec_Empty_Empty(const Rec_Empty_Empty & candidate) {
  p = candidate.p;
  tau = candidate.tau;
  csY = candidate.csY;
  csY2 = candidate.csY2;
  locCosts = candidate.locCosts;
}

Rec_Empty_Empty::~Rec_Empty_Empty() { csY = NULL; csY2 = NULL; locCosts = NULL; }

//accessory---------------------------------------------------------------------
unsigned int Rec_Empty_Empty::get_tau() const { return tau; }
pRectangle* Rec_Empty_Empty::getRectangle() const {return rectangle;}

//tools-------------------------------------------------------------------------
bool Rec_Empty_Empty::EmptyOfCandidate() { return fl_empty; }

void Rec_Empty_Empty::idCandidate(unsigned int dim, unsigned int t, double** &csy, double** &csy2, double* &loccosts) {
  p = dim;
  tau = t;
  csY = csy;
  csY2 = csy2;
  locCosts = loccosts;
  fl_empty = false;
}

void Rec_Empty_Empty::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Rec_Empty_Empty>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
  fl_empty = false;
  pCost cost = pCost(p);
  //last t : vectlinktocands[vectlinktocands.size()-1] -> get_tau();
  cost.idpCost(p, tau, vectlinktocands[vectlinktocands.size() - 1] -> get_tau(), csY, csY2, locCosts);
  double r2 = (locCosts[(vectlinktocands[vectlinktocands.size() - 1] -> get_tau()) + 1] - locCosts[tau] - cost.get_kVYit())/cost.get_k();
  if (r2 < 0) {
    fl_empty = true;
  }
}
//----------------------------------------------------------------------------//
