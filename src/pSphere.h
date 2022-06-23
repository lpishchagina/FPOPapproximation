#ifndef PSPHERE_H
#define PSPHERE_H

#include "pCost.h"

/*+++
 Class pSphere
 -------------------------------------------------------------------------------
 Description:
 pSphere in p-dimension.

 Parameters:
 "c" - vector  of  the pSphere  center coordinates;
 "r" - value of the pSphere radius;
 "p" - dimension.
 -------------------------------------------------------------------------------
 */
class pSphere{
private:
  unsigned int p;
  double r;
  double* c;

public:
  pSphere() {};
  pSphere(unsigned int dim): p(dim),  r(0), c(new double[dim]) {}
  pSphere(unsigned int dim, double* center, double radius);
  pSphere(const pSphere &sphere);
  ~pSphere();

  unsigned int  get_p()const;
  double get_r() const;
  double* get_c()const;

  void createSphere(unsigned int dim, unsigned int i, unsigned int t, double** &csdY, double** &csdY2,  double* &lCosts);
  void idpSphere(unsigned int dim, double* center, double radius);
  double get_dist(double* pnt1, double* pnt2);
  bool isIntersection(const pSphere & sphere);
  bool isnotIntersection(const pSphere & sphere);
  bool isInclusion(const pSphere & sphere);
};
#endif //PSPHERE_H
//------------------------------------------------------------------------------
