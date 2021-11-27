#ifndef FPOP_H
#define FPOP_H

#include <vector>
#include <list>
#include <iterator>

#include <string>
#include <fstream>
#include <iostream>
#include <Rcpp.h>
#include "Candidate_sphere_sphere_1.h"
#include "Candidate_Iall_Eall_2.h"
#include "Candidate_Iall_Eempty_3.h"
#include "Candidate_Iempty_Eall_4.h"
#include "Candidate_Ilast_Eall_5.h"
#include "Candidate_Ilast_Erandom_6.h"
#include "Candidate_Iall_Erandom_7.h"
#include "Candidate_Irandom_Erandom_8.h"


using namespace Rcpp;
using namespace std;

template <class CandidateOfChange>
class FPOP{
private:
  unsigned int N;
  unsigned int Dim;
  double Penality;
  double** CumSumData;

  std::vector <unsigned int> Changes;
  std::vector <std::vector <double>> SegmentMeans;
  double GlobalCost;

public:
  //constructor*****************************************************************
  FPOP<CandidateOfChange>(){}

  FPOP<CandidateOfChange> (Rcpp::NumericMatrix data, double penality){
    Dim  = (unsigned int)data.nrow();
    N = (unsigned int)data.ncol();
    Penality = penality;
    CumSumData = new double*[N + 1];
    for(unsigned int i = 0; i < (N + 1); i++) {CumSumData[i] = new double[(2*Dim)];}
  }
  //constructor copy************************************************************
  FPOP<CandidateOfChange> (const FPOP<CandidateOfChange> &candidate){
    Dim  = candidate.Dim;
    N = candidate.N;
    Penality = candidate.Penality;
    Changes = candidate.Changes;
    SegmentMeans = candidate.SegmentMeans;
    GlobalCost = candidate.GlobalCost;

    CumSumData = new double*[N + 1];
    for(unsigned int i = 0; i < N + 1; i++) {
      CumSumData[i] = new double[(2*Dim)];
      for (unsigned int k = 0; k < Dim; k++){
        CumSumData[i][k] = candidate.CumSumData[i][k];
        CumSumData[i][Dim + k] = candidate.CumSumData[i][Dim + k];
      }
    }
  }
  //destructor******************************************************************
  ~FPOP<CandidateOfChange>(){
    for(unsigned int i = 0; i < N + 1; i++) {delete(CumSumData[i]);}
    delete [] CumSumData;
    CumSumData = NULL;
  }
  //accessory*******************************************************************
  std::vector <unsigned int> GetChanges() const {return Changes;}
  std::vector <std::vector <double>> GetSegmentMeans() const{return SegmentMeans;}
  double GetGlobalCost() const {return GlobalCost;}
  unsigned int GetN() const {return N;}
  unsigned int GetDim() const {return Dim;}
  double GetPenality() const {return Penality;}

  //preprocessing***************************************************************
  double** CalcCumSumData(Rcpp::NumericMatrix data) {
    for (unsigned int k = 0; k < Dim; k++){ CumSumData[0][k] = 0; CumSumData[0][Dim + k] = 0;}
    for (unsigned int j = 1; j < (N + 1); j++){
      for (unsigned int k = 0; k < Dim; k++){
        CumSumData[j][k] = CumSumData[j - 1][k] + data(k, j-1);
        CumSumData[j][Dim + k] = CumSumData[j - 1][Dim + k] + data(k, j-1)*data(k, j-1);
      }
    }
    return(CumSumData);
  }

  //algorithm FPOP**************************************************************
  void algoFPOP(Rcpp::NumericMatrix data, int type_approx, bool test_mode){
    //preprocessing-------------------------------------------------------------
    double* VectOfCosts = new double[N + 1];                       //GlobalCost = VectOfCosts[n] - Changes.size()*Penality
    double* TempMean = new double[Dim];                         //values of temporary Means
    unsigned int* LastChpt = new unsigned int[N];       //vector of the best last changepoints
    double** LastMean = new double*[N];                 //matrix (nxp) of SegmentMeans for the best last changepoints
    for(unsigned int i = 0; i < N; i++) {LastMean[i] = new double[Dim];}

    VectOfCosts[0] = 0;
    CumSumData = CalcCumSumData(data);

    std::ofstream test_file;                                  // candidates test
    if (test_mode == true){test_file.open("test.txt");}
    CandidateOfChange candidate = CandidateOfChange(Dim);
    pSphere disk = pSphere(Dim);
    Cost cost = Cost(Dim);

    std::list<CandidateOfChange> ListOfCandidates;                     // list of active geometries
    std::vector< typename std::list<CandidateOfChange>::iterator> VectLinkToCandidates ;

  //  std::vector<unsigned int> VectTauOfCandidates;                     //list of active disks(t-1)
    double min_val;
    unsigned int tau;
    unsigned int u;
    //Algorithm-----------------------------------------------------------------
    for (unsigned int t = 0; t < N; t++){
      Rcpp::Rcout << "t ="<< t <<endl;
      cost.InitialCost(Dim, t, t, CumSumData[t], CumSumData[t+1], VectOfCosts[t]);
      min_val = cost.get_min();                       //min value of cost
      tau = t;                                 //best last position
      for (unsigned j = 0; j < Dim; j++){TempMean[j] = cost.get_mu()[j];}
      VectLinkToCandidates.clear();
      //First run: searching min------------------------------------------------
      typename std::list<CandidateOfChange>::reverse_iterator rit_candidate = ListOfCandidates.rbegin();

      while(rit_candidate != ListOfCandidates.rend()){
        u = rit_candidate -> GetTau();
        // Searching: min
        cost.InitialCost(Dim, u, t, CumSumData[u], CumSumData[t + 1], VectOfCosts[u]);
        if( min_val >= cost.get_min()){
          for (unsigned j = 0; j < Dim; j++){TempMean[j] = cost.get_mu()[j];}
          min_val = cost.get_min();
          tau = u;
        }
        ++rit_candidate;
      }//First run: end

      //new min, best last changepoint and SegmentMeans--------------------------------
      VectOfCosts[t + 1] = min_val + Penality;
      LastChpt[t] = tau;
      Rcpp::Rcout<<"last chpt = "<<tau<<endl;
      for (unsigned int j = 0; j < Dim; j++){LastMean[t][j] = TempMean[j];}

      //Candidate of Change.Initialisation.----------------------------------------------
      candidate.CleanOfCandidate();                  //if necessary, we clear the memory
      candidate.InitialOfCandidate(t, CumSumData, VectOfCosts);

      ListOfCandidates.push_back(candidate);

      //generate vector
      typename std::list<CandidateOfChange>::iterator VecIt = ListOfCandidates.begin();
      while(VecIt != ListOfCandidates.end()){
        VectLinkToCandidates.push_back(VecIt);
        VecIt++;
      }

      //Second run: Update list of geometry-------------------------------------
      unsigned int SizeVectLink = VectLinkToCandidates.size();
      for (unsigned int i = 0; i< SizeVectLink; i++){VectLinkToCandidates[i]->UpdateOfCandidate(i,VectLinkToCandidates);}

      typename std::list<CandidateOfChange>::iterator it_candidate = ListOfCandidates.begin();

      while (it_candidate != ListOfCandidates.end()){
        if (it_candidate -> EmptyOfCandidate()){it_candidate = ListOfCandidates.erase(it_candidate);--it_candidate;}
        ++it_candidate;
      }
    }
    //Result vectors------------------------------------------------------------
    std::vector<double> SegmentMeans_chp;
    unsigned int chp = N;
    while (chp > 0){
      Changes.push_back(chp);
      SegmentMeans_chp.clear();
      for (unsigned int i = 0; i < Dim; i++){SegmentMeans_chp.push_back(LastMean[chp-1][i]);}
      SegmentMeans.push_back(SegmentMeans_chp);
      chp = LastChpt[chp-1];
    }
    reverse(Changes.begin(), Changes.end());
    GlobalCost = VectOfCosts[N] - Penality * (Changes.size()-1);//verifier!!!
    Changes.pop_back();
    reverse(SegmentMeans.begin(), SegmentMeans.end());
    //memory--------------------------------------------------------------------
    for(unsigned int i = 0; i < N; i++) {delete(LastMean[i]);}
    delete [] LastMean;
    delete [] LastChpt;
    delete [] TempMean;
    delete [] VectOfCosts;
    VectOfCosts = NULL;
    LastMean = NULL;
    LastChpt = NULL;
    TempMean = NULL;
  }
};

#endif //FPOP_H
