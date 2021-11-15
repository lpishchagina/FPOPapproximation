#ifndef GEOM_IALL_EEMPTY_3_H
#define GEOM_IALL_EEMPTY_3_H

#include <iostream>
#include <vector>
#include <list>

#include "pRectangle.h"
#include "Cost.h"

/*
 Class Geom_Iall_Eempty_3
 --------------------------------------------------------------------------------
 Description of geometry "Geom_Iall_Eempty_3":
 Geometry for FPOP-Algorithm in p-dimension.

 Geometry parameters:
 "p" - dimension;
 "label_t" - moment in time;
 "rect_t" -  pointer to rectangle (approximated set);

 The updated geometry  is a rectangle that approximates the intersection of the rectangle and disk.

 Check for emptiness - the correctness of the  rectangle coordinates.
 --------------------------------------------------------------------------------
 */
class Geom_Iall_Eempty_3{
private:
  unsigned int p;
  unsigned int label_t;
  pRectangle* rect_t;


public:
  Geom_Iall_Eempty_3(): p(0), label_t(0), rect_t(0){}
  Geom_Iall_Eempty_3(unsigned  int dim): p(dim), label_t(0), rect_t(new pRectangle(dim)){}
  Geom_Iall_Eempty_3(unsigned int dim, unsigned int t): p(dim), label_t(t), rect_t(new pRectangle(dim)){}
  Geom_Iall_Eempty_3(const Geom_Iall_Eempty_3 & geom1);
  ~Geom_Iall_Eempty_3();

  unsigned int get_p() const;
  unsigned int get_label_t() const;
  std::list<pSphere> get_disks_t_1() const;

  void CleanGeometry(){};
  bool EmptyGeometry();
  void InitialGeometry(unsigned int dim, unsigned int t, const std::list<pSphere> &disks);
  void UpdateGeometry(const pSphere &disk_t);
};
#endif //GEOM_IALL_EEMPTY_3_H
//------------------------------------------------------------------------------
