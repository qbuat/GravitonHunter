//This class is meant to perform the final
// selection of the analysis
// This class inherits from GravitonAnalysis.
#ifndef BkgSampleSelector_h
#define BkgSampleSelector_h

#include <iostream>
#include <utility>
#include <vector>
#include <TTree.h>
#include <TString.h>
#include <TH2.h>
#include <TFile.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include "GravitonSelector.h"

class BkgSampleSelector : public GravitonSelector
{
 private:

  TString m_truthoutput;
  TFile *file_truth;
  TTree *tree_truth;

  //.....................
  double Lead_iso;
  double Lead_weight;
  int    Lead_IsTight;
  int    Lead_IsLoosePrime2;
  int    Lead_IsLoosePrime3;
  int    Lead_IsLoosePrime4;
  int    Lead_IsLoosePrime5;
  //.....................
  double SubLead_iso;
  double SubLead_weight;
  int    SubLead_IsTight;
  int    SubLead_IsLoosePrime2;
  int    SubLead_IsLoosePrime3;
  int    SubLead_IsLoosePrime4;
  int    SubLead_IsLoosePrime5;

  //......................
  double nlo_weight;
  double nlo_weight_xabier;
  double Lead_truth_eta;
  double Lead_truth_phi;
  double Lead_truth_pT;
  int    Lead_truth_status;
  double SubLead_truth_eta;
  double SubLead_truth_phi;
  double SubLead_truth_pT;
  int    SubLead_truth_status;

  double truth_mgg;
  double truth_ptgg;
  double truth_costhetastar;
  double truth_ygg;
  double truth_dphigg;

  //----> Looseprime masks
  unsigned int looseprime2;
  unsigned int looseprime3;
  unsigned int looseprime4;
  unsigned int looseprime5;

  //----------------------
 public :
  BkgSampleSelector();
  BkgSampleSelector(TTree*);
  virtual ~BkgSampleSelector();
  void   Init();
  void   EventLoop(TString outfile,int Nrelaxedcut=0);
  void   TruthEventLoop(TString outfile);

 private :
  void   InitTruthTree();
  void   InitOutTree();
  bool   EventTruthSelector(int entry,bool data);
  void   FillOutputFile();
  void   FillOutputTruthFile();

 public :

  static double GetkFactor_Xabier(double mass);
  static double GetkFactor_40_30(double mass);
  static double GetkFactor_50_50(double mass,TString config="nominal");
};

#endif

