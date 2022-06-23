#ifndef REC_ALL_LALL_H
#define REC_ALL_LALL_H

#include "pRectangle.h"
#include <vector>
#include <list>
#include <iterator>
#include <stdio.h>

/*+++
 Class Rec_All_lAll
 -------------------------------------------------------------------------------
 Description:
 Rectangular approximation (Zone = AI(all) + AE(all))

 Parameters:
 "p" - dimension;
 "tau" - change-point candidate;
 "csY" - cumsum of data;
 "csY2" - cumsum data^2;
 "locCosts" - min value of costs
 "rectangle" - approximation zone
 "spheresBefore" - exclusion set;
 "flCreate" - indicator, if "true" => new change-point candidate => create exclusion set;
 -------------------------------------------------------------------------------
 */

class Rec_All_lAll {
private:
  unsigned int p;
  unsigned int tau;
  pRectangle* rectangle;
  double** csY;
  double** csY2;
  double* locCosts;
  std::list<pSphere> spheresBefore;
  bool flCreate;

public:
  Rec_All_lAll(): p(0), tau(0), rectangle(0), csY(NULL), csY2(NULL), locCosts(NULL), flCreate(true) {}
  Rec_All_lAll(unsigned  int dim): p(dim), tau(0), rectangle(new pRectangle(dim)),  csY(NULL),csY2(NULL), locCosts(NULL), flCreate(true) {}
  Rec_All_lAll(unsigned int dim, unsigned int t): p(dim), tau(t), rectangle(new pRectangle(dim)), csY(NULL), csY2(NULL), locCosts(NULL), flCreate(true) {}
  Rec_All_lAll(const Rec_All_lAll & candidate);
  ~Rec_All_lAll();

  unsigned int get_tau()const;
  std::list<pSphere> get_spheresBefore() const;

  double get_dist(double* pnt1, double* pnt2);
  bool EmptyOfCandidate();
  void idCandidate(unsigned int dim, unsigned int t, double** &csy, double** &csy2, double* &loccosts);
  void UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Rec_All_lAll>::iterator> &vectlinktocands, unsigned int& RealNbExclus);
};
#endif //REC_ALL_LALL_H
