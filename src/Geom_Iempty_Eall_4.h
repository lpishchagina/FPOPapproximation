#ifndef GEOM_IEMPTY_EALL_4_H
#define GEOM_IEMPTY_EALL_4_H

#include <iostream>
#include <vector>
#include <list>

#include "pRectangle.h"
#include "Cost.h"

/*
 Class Geom_Iempty_Eall_4
 --------------------------------------------------------------------------------
 Description of geometry "Geom_Iempty_Eall_4":
 Geometry for FPOP-Algorithm in p-dimension.

 Geometry parameters:
 "p" - dimension;
 "label_t" - moment in time;
 "rect_t" -  pointer to rectangle (approximated set);
 "disks_t_1" - list of active disks at moment (t-1);

 The updated geometry is a rectangle that approximates (the intersection of the rectangle and disk at the moment t) minus (list of disks(t-1)) .

 Check for emptiness - the correctness of the  rectangle coordinates.
 --------------------------------------------------------------------------------
 */

class Geom_Iempty_Eall_4{
private:
  unsigned int p;
  unsigned int label_t;
  pRectangle* rect_t;
  std::list<pSphere> disks_t_1;

public:
  Geom_Iempty_Eall_4(): p(0), label_t(0), rect_t(0){}
  Geom_Iempty_Eall_4(unsigned  int dim): p(dim), label_t(0), rect_t(new pRectangle(dim)){}
  Geom_Iempty_Eall_4(unsigned int dim, unsigned int t): p(dim), label_t(t), rect_t(new pRectangle(dim)){}
  Geom_Iempty_Eall_4(const Geom_Iempty_Eall_4 & geom2);
  ~Geom_Iempty_Eall_4();

  unsigned int get_p()const;
  unsigned int get_label_t()const;
  std::list<pSphere> get_disks_t_1()const;

  void CleanGeometry();
  bool EmptyGeometry();
  void InitialGeometry(unsigned int dim, unsigned int t, const std::list<pSphere> &disks);
  void UpdateGeometry(const pSphere &disk_t);
};
#endif //GEOM_IEMPTY_EALL_4_H
//------------------------------------------------------------------------------
