// This class is used to plot, print results
// It worked independantly from the rest of the code 
// It takes as input a rootfile name and a tree name (2 TString)
// See below which method is used for which rootfile.


#ifndef Results_h
#define Results_h

#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>
#include <TString.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TF1.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TText.h>
#include <TLatex.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TMath.h>
#include <TRandom3.h>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include "ToolsTreeReader.h"
#include "GravitonAnalysis.h"
#include "GoodRunsLists/DQHelperFunctions.h"

using namespace std;
class Results{
 private:
  int         m_nentries;
  TFile*      m_file;
  TTree*      m_tree;
  TreeReader* m_Rd;
 public:
  Results();
  Results(TString,TString);
  virtual ~Results();
  void Init(TString,TString);

  //=== Methods for data.root ===//
  void DrawKinematics();
  //=== Methods for IsoForFitStudyLoosePrime0.root ===//
  void  CreateInputForEvanCode();
  void  CloneInputTree();
  void  CreateListOfEvents();
  void  CompareListOfEvents();
  void  CompareData_rel16_rel17();
  //=== Methods to create input of final results ===//
  void CreateDataInvMassHist();
  void CreateGamGamHist_Fromxsec();
  void CreateInputForBAT_RescaleBK6();
  void SimpleNullHypothesisTest(int NPE);

  //=== Methods for photons.root ===//
  void LCvsJC();
  void DrawHistPt();
  void DrawHisteta();
  void DrawStackedHist();
  void DrawHistMET();
  /* void DrawHistLeadPhPtvsMET(); */

  //== Methods for invmasscleaning.root ===//
  void DrawInvMassStackedHist();
  //== Methods for BinnedBkg and SplotBkg ==//
  void DrawBkgComponents();

  double efficiency(double,double);
  double sigma_efficiency(double,double);


};
#endif // #ifdef Results_h
