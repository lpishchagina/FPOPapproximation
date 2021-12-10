#ifndef PRECTANGLE_H
#define PRECTANGLE_H

#include "math.h"
#include "pSphere.h"

/*
 Class pRectangle
 -------------------------------------------------------------------------------
 Description:
 Rectangle in p-dimension.

 Parameters:
 "coordinates" - the values of two constraints for each axis;
 "p" - dimension.
 -------------------------------------------------------------------------------
 */

class pRectangle{
private:
  unsigned int p;
  double** coordinates;//matrix(px2) of constraints for x ,each xi =(xi0,xi1)  i = 0, p-1

public:
  pRectangle(): p(0), coordinates(NULL){}
  pRectangle(unsigned int dim);
  pRectangle(unsigned int dim, double** coords);
  pRectangle(const pRectangle &rect);
  ~pRectangle();

  double** get_coordinates()const;
  unsigned int get_p()const;

  double min_ab(double a, double b);
  double max_ab(double a, double b);
  void Clean_rect();

  bool EmptyIntersection(const pSphere &disk);

  bool IsEmpty_rect();
  void DoEmpty_rect();
  void Exclusion_disk(const pSphere &disk);
  void Intersection_disk(const pSphere &disk);
  void CubeApproximation(const pSphere &disk);
};

#endif //PRECTANGLE_H
