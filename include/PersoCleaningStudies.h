// There is also two methods dedicated to the cleaning studies
// This class inherits from GravitonAnalysis.
#ifndef CleaningStudies_h
#define CleaningStudies_h

#include <iostream>
#include <utility>
#include <vector>
#include <TTree.h>
#include <TString.h>
#include <TH2.h>
#include <TFile.h>
#include <TMath.h>

#include "InvMassCalculator.h"
#include "GravitonAnalysis.h"


class CleaningStudies : public GravitonAnalysis
{
 private:

  //--> Boolean for the cut flow 
  bool         m_passTrigger;
  bool         m_inGRL;
  bool         m_passPV;
  bool         m_passPreSel;
  //--------------------------
  TString      m_output;
  TString      m_PU_mcfile;
  TString      m_PU_datafile;
  InvMassCalc* m_IMC;


 public :
  CleaningStudies();
  CleaningStudies(TTree*);
  virtual ~CleaningStudies();
  void     Init(TTree*);
  void     InvMassCleaning(bool data,TString output="DEFAULT");
  void     PhotonStudies(bool data, TString output="DEFAULT");
  //=========Irreducible background ===========//

  //=============================//
  bool     HasJetAssociated(int);

};

#endif

