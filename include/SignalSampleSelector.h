// This class is meant to create the various
// RS signal templates
// It inherits from the GravitonAnalysis class 
#ifndef SignalSampleSelector_h
#define SignalSampleSelector_h

#include <iostream>
#include <utility>
#include <math.h>
#include <vector>

#include <TTree.h>
#include <TString.h>
#include <TH2.h>
#include <TFile.h>
#include <TObjArray.h>

#include "GravitonSelector.h"

class SignalSampleSelector : public GravitonSelector
{
 private:


  //.....................
  double Lead_iso;
  double Lead_weight;
  int    Lead_IsTight;
  //.....................
  double SubLead_iso;
  double SubLead_weight;
  int    SubLead_IsTight;

  //......................
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



 public :
  SignalSampleSelector();
  SignalSampleSelector(TTree*);
  virtual ~SignalSampleSelector();
  void     Init();
  void     EventLoop(TString outfile);

 private :

  void    InitOutTree();
  void    FillOutputFile();
  bool    EventTruthSelector(int entry,bool data);
};

#endif

