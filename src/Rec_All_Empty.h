#ifndef REC_ALL_EMPTY_H
#define REC_ALL_EMPTY_H

#include "pRectangle.h"
#include <vector>
#include <list>
#include <iterator>
#include <stdio.h>

/*+++
 Class Rec_All_Empty
 -------------------------------------------------------------------------------
 Description:
 Rectangular approximation without operator "exclusion approximation" (Zone = AI(all))

 Parameters:
 "p" - dimension;
 "tau" - change-point candidate;
 "csY" - cumsum of data;
 "csY2" - cumsum data^2;
 "locCosts" - min value of costs
 "rectangle" - approximation zone
 -------------------------------------------------------------------------------
 */

class Rec_All_Empty {
private:
  unsigned int p;
  unsigned int tau;
  pRectangle* rectangle;
  double** csY;
  double** csY2;
  double* locCosts;

public:
  Rec_All_Empty(): p(0), tau(0), rectangle(0), csY(NULL), csY2(NULL), locCosts(NULL) { }
  Rec_All_Empty(unsigned  int dim): p(dim), tau(0), rectangle(new pRectangle(dim)), csY(NULL), csY2(NULL), locCosts(NULL) { }
  Rec_All_Empty(unsigned int dim, unsigned int t): p(dim), tau(t), rectangle(new pRectangle(dim)), csY(NULL), csY2(NULL), locCosts(NULL){ }
  Rec_All_Empty(const Rec_All_Empty & candidate);
  ~Rec_All_Empty();

  unsigned int get_tau()const;

  bool EmptyOfCandidate();
  void idCandidate(unsigned int dim, unsigned int t, double** &csy, double** &csy2, double* &loccosts);
  void UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Rec_All_Empty>::iterator> &vectlinktocands, unsigned int& RealNbExclus);
};
#endif //REC_ALL_EMPTY_H
