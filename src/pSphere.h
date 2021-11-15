#ifndef PSPHERE_H
#define PSPHERE_H

/*
 Class pSphere
 -------------------------------------------------------------------------------
 Description:
 pSphere in p-dimension.

 Parameters:
 "center" - vector  of  the pSphere  center coordinates;
 "radius" - value of the pSphere radius;
 "p" - dimension.
 -------------------------------------------------------------------------------
 */
class pSphere{
private:
  unsigned int p;
  double radius;
  double* center;

public:
  pSphere(){};
  pSphere(unsigned int dim);
  pSphere(unsigned int dim, double* c, double r);
  pSphere(const pSphere &sphere);
  ~pSphere();

  void InitialpSphere(unsigned int dim, double* c, double r);

  unsigned int  get_p()const;
  double get_radius() const;
  double* get_center()const;
};
#endif //PSPHERE_H
//------------------------------------------------------------------------------
