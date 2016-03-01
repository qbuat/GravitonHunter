//This class is meant to perform the final
// selection of the analysis
// This class inherits from GravitonAnalysis.
#ifndef SinglePhotonSelector_h
#define SinglePhotonSelector_h

#include <iostream>
#include <utility>
#include <vector>
#include <TTree.h>
#include <TString.h>
#include <TH2.h>
#include <TFile.h>
#include <TMath.h>

#include "GravitonAnalysis.h"


class SinglePhotonSelector : public GravitonAnalysis
{
 private:

  int                 m_RunNumber;
  int                 m_EventNumber;
  int                 m_LumiBlock;

  //--> Boolean for the cut flow 
  bool         m_passTrigger;
  bool         m_inGRL;
  bool         m_passPV;
  bool         m_passPreSel;
  //--------------------------
  TString      m_output;
  TString      m_PU_mcfile;
  TString      m_PU_datafile;

  //--> Looseprime mask
  //--> for the isolation root files 
  unsigned int looseprime2;
  unsigned int looseprime3;
  unsigned int looseprime4;
  unsigned int looseprime5;

 public :
  SinglePhotonSelector();
  SinglePhotonSelector(TTree*);
  virtual ~SinglePhotonSelector();
  void     Init(TTree*);
  void     EventLoop( bool data,
		      TString output="DEFAULT",
		      int Nreversedcut=0,
		      TString mctype="" );


};

#endif

