#ifndef PileUpToolHandler_h
#define PileUpToolHandler_h
#include <iostream>
#include <vector>
#include "TCanvas.h"
#include "TTree.h"
#include "TString.h"
#include "ToolsTreeReader.h"
#include "PileupReweighting/TPileupReweighting.h"

class PileUpToolHandler
{
 private:
  TreeReader* m_Rd;
  int m_nentries;
  Root::TPileupReweighting* m_tPileUp;
  TString m_mctype;
  TString m_fileMC;
  TString m_histMC;
  TString m_fileD;
  TString m_histD;
 public:
  PileUpToolHandler();
  PileUpToolHandler(TTree*);
  virtual ~PileUpToolHandler();
  void  Init(TTree*);
  void  SetMCFile(TString);
  void  SetDataFile(TString);
  void  GenerateMCFile(TString);
  TCanvas*  CheckPileUpWeights();

};
#endif

