#ifndef PCOST_H
#define PCOST_H
/* +++
 Class pCost
 -------------------------------------------------------------------------------
 Description:
 The Gaussian cost for the interval (i,t) in p-dimension.

 Parameters:
 "p" - dimension;
 "k" = (t - i + 1);
 "EYit" - vector of means for the interval (i,t);
 "mi1beta" - sum of the value of the optimal cost at moment (i-1) and penalty;
 "kVYit" = (t - i + 1) * Var(Yit).

 Remark :  min{The Gaussian cost} = kVYit+ mi1beta
 -------------------------------------------------------------------------------
 */
class pCost{
private:
  unsigned int p;
  unsigned int k;
  double kVYit;
  double mi1beta;
  double* EYit;

public:
  pCost(){};
  pCost(unsigned int dim);
  pCost(unsigned int dim, unsigned int i, unsigned int t, double** &csdY, double** &csdY2, double* &lCosts);
  pCost(const pCost &cost);
  ~pCost();

  unsigned int get_p() const;
  unsigned int get_k() const;
  double get_kVYit() const;
  double get_mi1beta() const;
  double* get_EYit();

  void idpCost(unsigned int dim, unsigned int i, unsigned int t, double** &csdY, double** &csdY2, double* &lCosts);
  double get_min();
};

#endif // PCOST
