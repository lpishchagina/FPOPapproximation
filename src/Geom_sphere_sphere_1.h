#ifndef GEOM_SPHERE_SPHERE_1_H
#define GEOM_SPHERE_SPHERE_1_H

#include <iostream>
#include <list>
#include <iterator>

#include "pSphere.h"
#include "Cost.h"
/*
 Class Geom_sphere_sphere_1
 --------------------------------------------------------------------------------
 Description of geometry "Geom_sphere_sphere_1":
 Geometry for FPOP-Algorithm in p-dimension.

 Geometry parameters:
 "p" - dimension;
 "label_t" - moment in time;
 "disks_t_1" - list of active disks at moment (t-1);
 "fl_empty"  - "false" if geometry exists, otherwise "true".

 The updated geometry is the result of exclusion from the disk at time t of disks at time t-1.

 Check for emptiness - the distance between the centers of the disks.
 --------------------------------------------------------------------------------
 */
class Geom_sphere_sphere_1{
private:
  unsigned int p;
  unsigned int label_t;
  std::list<pSphere> disks_t_1;
  bool fl_empty;

public:
  Geom_sphere_sphere_1(){};
  Geom_sphere_sphere_1(unsigned  int dim): p(dim), label_t(0), fl_empty(false){}
  Geom_sphere_sphere_1(unsigned int dim, unsigned int t): p(dim), label_t(t), fl_empty(false){}
  Geom_sphere_sphere_1(const Geom_sphere_sphere_1 & geom3);

  unsigned int get_p() const;
  unsigned int get_label_t()const;
  bool get_fl_empty()const;
  std::list<pSphere> get_disks_t_1()const;

  double Dist(double* a, double* b);

  void CleanGeometry();
  bool EmptyGeometry();
  void InitialGeometry(unsigned int dim, unsigned int t, const std::list<pSphere> &disks);
  void UpdateGeometry(const pSphere &disk_t);
};
#endif //GEOM_SPHERE_SPHERE_1_H
//------------------------------------------------------------------------------
