#ifndef FPOP_H
#define FPOP_H

#include "Candidate_Ilast_Eall_sphere_1.h"
#include "Candidate_Iall_Eall_2.h"
#include "Candidate_Iall_Eempty_3.h"
#include "Candidate_Iempty_Eall_4.h"
#include "Candidate_Ilast_Eall_5.h"
#include "Candidate_Ilast_Erandom_6.h"
#include "Candidate_Iall_Erandom_7.h"
#include "Candidate_Irandom_Erandom_8.h"
#include "Candidate_Iempty_Eempty_9.h"
#include "Candidate_Irandom_Eall_10.h"
#include "Candidate_Irandom_Eall_sphere_11.h"
#include "Candidate_Irandom_Erandom_sphere_12.h"
#include "Candidate_Ilast_Erandom_sphere_13.h"
#include "Candidate_Ilast_EallModif_14.h"
#include "Candidate_Ilast_EallModif2_15.h"

#include <Rcpp.h>
#include "math.h"

using namespace Rcpp;
using namespace std;

template <class CandidateOfChange>
class FPOP {
private:
  unsigned int N;
  unsigned int Dim;
  double Penality;
  double** CumSumData;
  double** CumSumData2;
  std::vector <unsigned int> Changes;
  std::vector <std::vector <double>> SegmentMeans;
  double UnpenalizedCost;

  double* VectOfCosts;                    //UnpenalizedCost = VectOfCosts[n] - Changes.size()*Penality
  unsigned int* LastChpt;       //vector of the best last changepoints
  //
  std::vector <unsigned int> NbOfCandidats;  //vector of the candidats number

public:
  FPOP<CandidateOfChange>() { }

  FPOP<CandidateOfChange> (Rcpp::NumericMatrix data, double penality) {
    Dim  = (unsigned int)data.nrow();
    N = (unsigned int)data.ncol();
    Penality = penality;

    VectOfCosts = new double[N + 1];
    LastChpt = new unsigned int[N];

    CumSumData = new double*[N + 1];
    CumSumData2 = new double*[N + 1];
    for (unsigned int i = 0; i < (N + 1); i++) {
      CumSumData[i] = new double[Dim];
      CumSumData2[i] = new double[Dim];
    }
  }

  FPOP<CandidateOfChange> (const FPOP<CandidateOfChange> &candidate) {
    Dim  = candidate.Dim;
    N = candidate.N;
    Penality = candidate.Penality;
    Changes = candidate.Changes;
    SegmentMeans = candidate.SegmentMeans;
    UnpenalizedCost = candidate.UnpenalizedCost;
    //
    NbOfCandidats = candidate.NbOfCandidats;

    VectOfCosts = new double[N + 1];
    LastChpt = new unsigned int[N];

    CumSumData = new double*[N + 1];
    CumSumData2 = new double*[N + 1];

    for (unsigned int i = 0; i < N + 1; i++) {

      VectOfCosts[i] = candidate.VectOfCosts[i];

      CumSumData[i] = new double[Dim];
      CumSumData2[i] = new double[Dim];
      for (unsigned int k = 0; k < Dim; k++) {
        CumSumData[i][k] = candidate.CumSumData[i][k];
        CumSumData2[i][k] = candidate.CumSumData2[i][k];
      }
    }
    for (unsigned int i = 0; i < N; i++) {
      LastChpt[i] = candidate.LastChpt[i];
    }
  }

  ~FPOP<CandidateOfChange>() {
    for (unsigned int i = 0; i < N + 1; i++) {
      delete(CumSumData[i]);
      delete(CumSumData2[i]);
    }
    delete [] CumSumData;
    delete [] CumSumData2;
    delete [] VectOfCosts;
    delete [] LastChpt;

    CumSumData = NULL;
    CumSumData2 = NULL;
    VectOfCosts = NULL;
    LastChpt = NULL;
  }

  std::vector <unsigned int> GetChanges() const { return Changes; }
  std::vector <std::vector <double>> GetSegmentMeans() const { return SegmentMeans; }
  double GetUnpenalizedCost() const { return UnpenalizedCost; }
  unsigned int GetN() const { return N; }
  unsigned int GetDim() const { return Dim; }
  double GetPenality() const { return Penality; }
  double* GetVectOfCosts() const { return VectOfCosts; }
  unsigned int* GetLastChpt()  { return LastChpt; }
  //
  std::vector <unsigned int> GetNbOfCandidats()  { return NbOfCandidats; }

  double** CalcCumSumData(Rcpp::NumericMatrix data) {
    for (unsigned int k = 0; k < Dim; k++) {
      CumSumData[0][k] = 0;
    }
    for (unsigned int j = 1; j < (N + 1); j++) {
      for (unsigned int k = 0; k < Dim; k++) {
        CumSumData[j][k] = CumSumData[j - 1][k] + data(k, j-1);
      }
    }
    return(CumSumData);
  }

  double** CalcCumSumData2(Rcpp::NumericMatrix data) {
    for (unsigned int k = 0; k < Dim; k++) {
      CumSumData2[0][k] = 0;
    }
    for (unsigned int j = 1; j < (N + 1); j++) {
      for (unsigned int k = 0; k < Dim; k++) {
        CumSumData2[j][k] = CumSumData2[j - 1][k] + data(k, j-1) * data(k, j - 1);
      }
    }
    return(CumSumData2);
  }

  void algoFPOP(Rcpp::NumericMatrix data, int type_approx, bool NbOfCands, bool NbOfExclus){
    /*
    //file initialisation for tests
    std::ofstream FileNbCands;
    std::ofstream FileNbExclus;
     */
    unsigned int RealNbExclus = 0;
    //
    if (!NbOfCands) {
      NbOfCandidats.push_back(NULL);
    }
    /*
    if (NbOfCands == true) {
      FileNbCands.open("NbOfCands.txt", ios_base::trunc);
    }
    if (NbOfExclus == true) {
      FileNbExclus.open("NbOfExclus.txt", ios_base::trunc);
    }
    */
    VectOfCosts[0] = 0;
    CumSumData = CalcCumSumData(data);
    CumSumData2 = CalcCumSumData2(data);
    CandidateOfChange candidate = CandidateOfChange(Dim);
    pSphere disk = pSphere(Dim);
    Cost cost = Cost(Dim);
    std::list<CandidateOfChange> ListOfCandidates;                     // list of active geometries
    std::vector< typename std::list<CandidateOfChange>::iterator> VectLinkToCandidates ;
    double min_val;
    unsigned int tau;
    unsigned int u;
    //Algorithm-----------------------------------------------------------------
    for (unsigned int t = 0; t < N; t++) {
      cost.InitialCost(Dim, t, t, CumSumData, CumSumData2, VectOfCosts); // Guillem : a link to CumSumData ? check in cost object (cost.h)
      min_val = cost.get_min();                       //min value of cost
      tau = t;                                 //best last position
      //First run: searching min
      typename std::list<CandidateOfChange>::reverse_iterator rit_candidate = ListOfCandidates.rbegin();
      while (rit_candidate != ListOfCandidates.rend()) {
        u = rit_candidate -> GetTau();
        cost.InitialCost(Dim, u, t, CumSumData, CumSumData2, VectOfCosts);
        if (min_val >= cost.get_min()) {
          min_val = cost.get_min();
          tau = u;
        }
        ++rit_candidate;
      }
      //new min, best last changepoint and SegmentMeans--------------------------------
      VectOfCosts[t + 1] = min_val + Penality;
      LastChpt[t] = tau;
      //Candidate of Change.Initialisation.
      candidate.InitialOfCandidate(t, CumSumData, CumSumData2, VectOfCosts);
      ListOfCandidates.push_back(candidate);

      //Generate vector of link
      VectLinkToCandidates.clear();
      typename std::list<CandidateOfChange>::iterator VecIt = ListOfCandidates.begin();
      while (VecIt != ListOfCandidates.end()) {
        VectLinkToCandidates.push_back(VecIt);
        VecIt++;
      }
      //Second run:
      //Update ListOfCandidates
      unsigned int SizeVectLink = VectLinkToCandidates.size();
      for (unsigned int IndexOfCandVectLink = 0; IndexOfCandVectLink < SizeVectLink; IndexOfCandVectLink++) {
        VectLinkToCandidates[IndexOfCandVectLink] -> UpdateOfCandidate(IndexOfCandVectLink,VectLinkToCandidates, RealNbExclus);
/*        if (NbOfExclus) {
          FileNbExclus << VectLinkToCandidates[IndexOfCandVectLink] -> GetTau() << " " << RealNbExclus << " ";
        }
 */
      }

/*
      if (NbOfExclus) {
        FileNbExclus << "\n";
      }
 */
      //Remove empty candidates
      typename std::list<CandidateOfChange>::iterator it_candidate = ListOfCandidates.begin();
      while (it_candidate != ListOfCandidates.end()) {
        if (it_candidate -> EmptyOfCandidate()) {
          it_candidate = ListOfCandidates.erase(it_candidate);
          --it_candidate;
        }
        ++it_candidate;
      }
      //fixe nb_cands
      /*
      if (NbOfCands) {
        FileNbCands << VectLinkToCandidates.size() << " ";
      }
       */
      if (NbOfCands) {
        NbOfCandidats.push_back(VectLinkToCandidates.size());
      }
    }
    /*
    //close file
    if (NbOfCands) {
      FileNbCands << "\n";
      FileNbCands.close();
    }
    if (NbOfExclus) {
      FileNbExclus.close();
    }
     */
    //Result
    //vector of Changes
    unsigned int chp = N;
    while (chp > 0) {
      Changes.push_back(chp);
      chp = LastChpt[chp-1];
    }
    Changes.push_back(0);

    ///
    unsigned int j = 1;
    std::vector<double> MeanOneSegment;
    chp = N - 1;
    while (chp > 0) {
      MeanOneSegment.clear();
      for (unsigned int k = 0; k < Dim; k++) {
        MeanOneSegment.push_back((CumSumData[chp + 1][k] - CumSumData[Changes[j]][k])/(chp - Changes[j] + 1));
      }
      SegmentMeans.push_back(MeanOneSegment);
      chp = Changes[j];
      j++;
    }
    reverse(SegmentMeans.begin(), SegmentMeans.end());
    Changes.pop_back();//remove 0
    reverse(Changes.begin(), Changes.end());
    Changes.pop_back();//remove N

    UnpenalizedCost = VectOfCosts[N] - Penality * (Changes.size());
  }

  List ResAlgoFPOP(){
    List res;
    res["changes"] = GetChanges();
    res["means"] = GetSegmentMeans();
    res["UnpenalizedCost"] = GetUnpenalizedCost();
    res["NumberOfCandidats"] = GetNbOfCandidats();
    return res;
  }

  //Test of  the  UnpenalizedCost  and LastChpt
  bool TestLastChpt(unsigned int* LastChptOfTestX) {
    bool res  = true;
    unsigned int  IndexOfPoint = 0;
    while (res && (IndexOfPoint < N)) {
      if (LastChpt[IndexOfPoint] != LastChptOfTestX[IndexOfPoint] ) {
        res = false;
      }
      IndexOfPoint++;
    }
    return res;
  }

  bool TestUnpenalizedCost(double UnpenalizedCostOfTestX) {
    bool res  = (UnpenalizedCost == UnpenalizedCostOfTestX);
    return res;
  }

  bool TestFPOP(double UnpenalizedCostOfTestX, unsigned int* LastChptOfTestX) {
    bool res  = (UnpenalizedCost == UnpenalizedCostOfTestX);
    if (res) {
      res = TestLastChpt(LastChptOfTestX);
    }
    return res;
  }
};
#endif //FPOP_H

