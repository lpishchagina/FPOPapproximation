#ifndef SPH_LAST_LALL_H
#define SPH_LAST_LALL_H

#include "pRectangle.h"
#include <vector>
#include <list>
#include <iterator>
#include <stdio.h>

/*+++
 Class Sph_Last_lAll
 -------------------------------------------------------------------------------
 Description:
 Sphere (S-type) approximation

 Parameters:
 "p" - dimension;
 "tau" - change-point candidate;
 "csY" - cumsum of data;
 "csY2" - cumsum data^2;
 "locCosts" - min value of costs
 "spheresBefore" - exclusion set;
 "flCreate" - indicator, if "true" => new change-point candidate => create exclusion set;
 "fl_empty" - indicator, if "true" => approximation zone is empty;
 "rectangle" - empty (not used, added for template)
 -------------------------------------------------------------------------------
 */

class Sph_Last_lAll {
private:
  unsigned int p;
  unsigned int tau;
  double** csY;
  double** csY2;
  double* locCosts;
  std::list<pSphere> spheresBefore;
  bool flCreate;
  bool fl_empty;
  pRectangle* rectangle;//empty

public:
  Sph_Last_lAll(): p(0), tau(0), csY(NULL), csY2(NULL), locCosts(NULL), flCreate(true), fl_empty(false), rectangle(new pRectangle(p)) { }
  Sph_Last_lAll(unsigned  int dim): p(dim), tau(0), csY(NULL), csY2(NULL), locCosts(NULL), fl_empty(false) { }
  Sph_Last_lAll(unsigned int dim, unsigned int t): p(dim), tau(t), csY(NULL), csY2(NULL), locCosts(NULL), fl_empty(false) { }
  Sph_Last_lAll(const Sph_Last_lAll & candidate);
  ~Sph_Last_lAll();

  unsigned int get_tau()const;

  double get_dist(double* pnt1, double* pnt2);
  bool EmptyOfCandidate();
  void idCandidate(unsigned  int dim, unsigned int t, double** &csy,  double** &csy2, double* &loccosts);
  void UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Sph_Last_lAll>::iterator> &vectlinktocands, unsigned int& RealNbExclus);
};
#endif //SPH_LAST_LALL_H

