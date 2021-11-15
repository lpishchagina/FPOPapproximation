#ifndef OPDP_H
#define OPDP_H
#include "Geom_sphere_sphere_1.h"
#include "Geom_Iall_Eall_2.h"
#include "Geom_Iall_Eempty_3.h"
#include "Geom_Iall_Erandom_4.h"
#include "Geom_last1_Eall_5.h"



#include <vector>
#include <list>
#include <iterator>

#include <string>
#include <fstream>
#include <iostream>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

/*
 Class OPDp
 -------------------------------------------------------------------------------
 Description:
 Template for the realization of FPOP-Algorithm in p-dimension.

 Parameters:
 "p" - dimension;
 "n" - data length;
 "penalty" - value of penalty;

 "sx12" - matrix(n+1x2*p) of sum:x1:xp, x1^2:xp^2;

 "chpts" - vector of changepoints;
 "means" - means matrix for  vector of changepoints;
 "globalCost" - value of global cost.
 --------------------------------------------------------------------------------
 */
template <class GeomX>
class OPDp{
private:
  unsigned int n;
  unsigned int p;
  double penalty;
  double** sx12;

  std::vector <unsigned int> chpts;
  std::vector <std::vector <double>> means;
  double globalCost;

public:
  //constructor*****************************************************************
  OPDp<GeomX>(){}

  OPDp<GeomX> (Rcpp::NumericMatrix x, double beta){
    p  = (unsigned int)x.nrow();
    n = (unsigned int)x.ncol();
    penalty = beta;
    sx12 = new double*[n + 1];
    for(unsigned int i = 0; i < n + 1; i++) {sx12[i] = new double[(2*p)];}
  }
  //constructor copy************************************************************
  OPDp<GeomX> (const OPDp<GeomX> &geomX){
    p  = geomX.p;
    n = geomX.n;
    penalty = geomX.penalty;
    chpts = geomX.chpts;
    means = geomX.means;
    globalCost = geomX.globalCost;

    sx12 = new double*[n + 1];
    for(unsigned int i = 0; i < n + 1; i++) {
      sx12[i] = new double[(2*p)];
      for (unsigned int k = 0; k < p; k++){
        sx12[i][k] = geomX.sx12[i][k];
        sx12[i][p + k] = geomX.sx12[i][p + k];
      }
    }
  }
  //destructor******************************************************************
  ~OPDp<GeomX>(){
    for(unsigned int i = 0; i < n + 1; i++) {delete(sx12[i]);}
    delete [] sx12;
    sx12 = NULL;
  }
  //accessory*******************************************************************
  std::vector <unsigned int> get_chpts() const {return chpts;}

  std::vector <std::vector <double>> get_means() const{return means;}

  double get_globalCost() const {return globalCost;}

  unsigned int get_n() const {return n;}

  unsigned int get_p() const {return p;}

  double get_penalty() const {return penalty;}

  //preprocessing***************************************************************
  double** vect_sx12(Rcpp::NumericMatrix x) {
    for (unsigned int k = 0; k < p; k++){ sx12[0][k] = 0; sx12[0][p + k] = 0;}
    for (unsigned int j = 1; j < (n + 1); j++){
      for (unsigned int k = 0; k < p; k++){
        sx12[j][k] = sx12[j - 1][k] + x(k, j-1);
        sx12[j][p + k] = sx12[j - 1][p + k] + x(k, j-1)*x(k, j-1);
      }
    }
    return(sx12);
  }

  //algorithm FPOP**************************************************************
  void algoFPOP(Rcpp::NumericMatrix x, int type_approx, bool test_mode){
    //preprocessing-------------------------------------------------------------
    double* m = new double[n + 1];                       //globalCost = m[n] - chpts.size()*penalty
    double* mus = new double[p];                         //values of temporary means
    unsigned int* last_chpt = new unsigned int[n];       //vector of the best last changepoints
    double** last_mean = new double*[n];                 //matrix (nxp) of means for the best last changepoints
    for(unsigned int i = 0; i < n; i++) {last_mean[i] = new double[p];}

    m[0] = 0;
    sx12 = vect_sx12(x);

    std::ofstream test_file;                                  // candidates test
    if (test_mode == true){test_file.open("test.txt");}
    GeomX geom = GeomX(p);
    pSphere disk = pSphere(p);
    Cost cost = Cost(p);

    std::list<GeomX> list_geom;                     // list of active geometries
    std::list<pSphere> list_disk;                     //list of active disks(t-1)
    double min_val;
    double r2;
    unsigned int lbl;
    unsigned int u;
    //Algorithm-----------------------------------------------------------------
    for (unsigned int t = 0; t < n; t++){
      cost.InitialCost(p, t, t, sx12[t], sx12[t+1], m[t]);
      min_val = cost.get_min();                       //min value of cost
      lbl = t;                                 //best last position
      for (unsigned j = 0; j < p; j++){mus[j] = cost.get_mu()[j];}

      //First run: searching min------------------------------------------------
      typename std::list<GeomX>::reverse_iterator rit_geom = list_geom.rbegin();
      while(rit_geom != list_geom.rend()){
        u = rit_geom -> get_label_t();
        // Searching: min
        cost.InitialCost(p, u, t, sx12[u], sx12[t + 1], m[u]);
        if( min_val >= cost.get_min()){
          for (unsigned j = 0; j < p; j++){mus[j] = cost.get_mu()[j];}
          min_val = cost.get_min();
          lbl = u;
        }
        //list of active disks(t-1)
        cost.InitialCost(p, u, t-1, sx12[u], sx12[t], m[u]);
        r2 = (m[t] - m[u] - cost.get_coef_Var())/cost.get_coef();
        disk.InitialpSphere(p, cost.get_mu(), sqrt(r2));
        list_disk.push_back(disk);

        ++rit_geom;
      }//First run: end
      //new min, best last changepoint and means--------------------------------
      for (unsigned int j = 0; j < p; j++){last_mean[t][j] = mus[j];}
      m[t + 1] = min_val + penalty;
      last_chpt[t] = lbl;

      //Initialisation of geometry----------------------------------------------
      geom.CleanGeometry();                  //if necessary, we clear the memory
      geom.InitialGeometry(p,t,list_disk);
      list_disk.clear();    //we clear the disk list(t-1) for the next iteration
      list_geom.push_back(geom);

      //Second run: Update list of geometry-------------------------------------
      typename std::list<GeomX>::iterator it_geom = list_geom.begin();
      while (it_geom != list_geom.end()){
        lbl = it_geom->get_label_t();
        cost.InitialCost(p, lbl, t, sx12[lbl], sx12[t + 1], m[lbl]);
        r2 = (m[t + 1] - m[lbl] - cost.get_coef_Var())/cost.get_coef();
        //Pruning "PELT"
        if (r2 <= 0){it_geom = list_geom.erase(it_geom); --it_geom;}
        //Pruning "FPOP"
        if (r2 > 0){
          disk.InitialpSphere(p, cost.get_mu(), sqrt(r2));
          it_geom -> UpdateGeometry(disk);
          if (it_geom -> EmptyGeometry()){it_geom = list_geom.erase(it_geom);--it_geom;}
          else {
            if (test_mode == true){//test: the number of candidates and exclusions
              test_file << it_geom ->get_label_t() << " ";
              if (type_approx == 2 || type_approx == 3){test_file << it_geom ->get_disks_t_1().size() << " ";}
            }
          }
        }
        ++it_geom;
      }//Second run: end
      if (test_mode == true){test_file << "\n";}
    }
    if (test_mode == true){test_file.close();}

    //Result vectors------------------------------------------------------------
    std::vector<double> means_chp;

    unsigned int chp = n;
    while (chp > 0){
      chpts.push_back(chp);
      means_chp.clear();
      for (unsigned int i = 0; i < p; i++){means_chp.push_back(last_mean[chp-1][i]);}
      means.push_back(means_chp);
      chp = last_chpt[chp-1];
    }
    reverse(chpts.begin(), chpts.end());
    globalCost = m[n] - penalty * chpts.size();
    chpts.pop_back();
    reverse(means.begin(), means.end());


    //memory--------------------------------------------------------------------
    for(unsigned int i = 0; i < n; i++) {delete(last_mean[i]);}
    delete [] last_mean;
    delete [] last_chpt;
    delete [] mus;
    delete [] m;
    m = NULL;
    last_mean = NULL;
    last_chpt = NULL;
    mus = NULL;
  }
};

#endif //OPDP_H
