#ifndef REC_EMPTY_EMPTY_H
#define REC_EMPTY_EMPTY_H

#include "pRectangle.h"
#include "pCost.h"
#include <vector>
#include <list>
#include <iterator>
#include <stdio.h>

/*+++
 Class Rec_Empty_Empty
 -------------------------------------------------------------------------------
 Description:
 //PELT//

 Parameters:
 "p" - dimension;
 "tau" - change-point candidate;
 "csY" - cumsum of data;
 "csY2" - cumsum data^2;
 "locCosts" - min value of costs
 "fl_empty" - indicator, if "true" => approximation zone is empty;
 "rectangle" - empty (not used, added for template)
 -------------------------------------------------------------------------------
 */

class Rec_Empty_Empty {
private:
  unsigned int p;
  unsigned int tau;
  double** csY;
  double** csY2;
  double* locCosts;
  bool fl_empty;
  pRectangle* rectangle;// empty

public:
  Rec_Empty_Empty(): p(0), tau(0), csY(NULL), csY2(NULL), locCosts(NULL),    fl_empty(false), rectangle(NULL) { }
  Rec_Empty_Empty(unsigned  int dim): p(dim), tau(0),  csY(NULL), csY2(NULL), locCosts(NULL),  fl_empty(false), rectangle(NULL) { }
  Rec_Empty_Empty(unsigned int dim, unsigned int t): p(dim), tau(t), csY(NULL), csY2(NULL), locCosts(NULL),  fl_empty(false), rectangle(NULL) { }
  Rec_Empty_Empty(const Rec_Empty_Empty & candidate);
  ~Rec_Empty_Empty();

  unsigned int get_tau() const;
  pRectangle* getRectangle() const;

  bool EmptyOfCandidate();
  void idCandidate(unsigned int dim, unsigned int t, double** &csy, double** &csy2, double* &loccosts);
  void UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Rec_Empty_Empty>::iterator> &vectlinktocands, unsigned int& RealNbExclus);
};
#endif //REC_EMPTY_EMPTY_H
