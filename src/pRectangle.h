#ifndef PRECTANGLE_H
#define PRECTANGLE_H

#include "math.h"
#include "pSphere.h"

/*+++
 Class pRectangle
 -------------------------------------------------------------------------------
 Description:
 Rectangle in p-dimension.

 Parameters:
 "borders" - the values of two constraints for each axis;
 "p" - dimension.
 -------------------------------------------------------------------------------
 */

class pRectangle {
private:
  unsigned int p;
  double** borders;//matrix(px2) of constraints for x ,each xi =(xi0,xi1)  i = 0, p-1

public:
  pRectangle(): p(0),   borders(NULL) {}
  pRectangle(unsigned int dim);
  pRectangle(unsigned int dim, double** coords);
  pRectangle(const pRectangle &rect);
  ~pRectangle();

  double** get_borders()const;
  unsigned int get_p()const;

  double min_ab(double a, double b);
  double max_ab(double a, double b);
  double get_dist(double* pnt1, double* pnt2);

  void DoEmptyRect();
  bool IsEmptyRect();
  bool EmptyIntersection(const pSphere &disk);
  void SphereApproximation(const pSphere &sphere);
  void IntersectionSphere(const pSphere &sphere);
  void ExclusionSphere(const pSphere &sphere);


};

#endif //PRECTANGLE_H




