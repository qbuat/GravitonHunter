// This class inherits from GravitonAnalysis.
// It is meant to perform Photon efficiencies
// measurement for Graviton Analysis

#ifndef PhotonEfficiencies_h
#define PhotonEfficiencies_h

#include <iostream>
#include <utility>
#include <vector>
#include <TString.h>
#include <TH2F.h>
#include "ToolsUtilities.h"
#include "ToolsTruthSelector.h"
#include "GravitonAnalysis.h"

class PhotonEfficiencies : public GravitonAnalysis
{
 private:

  TString m_mcsample;
  std::vector<double>  m_bins_pt;
  std::vector<double>  m_bins_eta;
  std::pair<double,double>  m_mggbin;
  TH2F* m_map_lead_truth;
  TH2F* m_map_lead_reco;
  TH2F* m_map_lead_recotight;
  TH2F* m_map_lead_recotightiso;
  TH2F* m_map_lead_reco_isoapplied;
  TH2F* m_map_lead_recotight_isoapplied;
  TH2F* m_map_sublead_truth;
  TH2F* m_map_sublead_reco;
  TH2F* m_map_sublead_recotight;
  TH2F* m_map_sublead_recotightiso;
  TH2F* m_map_sublead_reco_isoapplied;
  TH2F* m_map_sublead_recotight_isoapplied;

 public :
  PhotonEfficiencies();
  PhotonEfficiencies(TTree*);
  virtual ~PhotonEfficiencies();
  void  Init();
  void  EventLoop(TString MCTYPE = "SMGG");

  void SetMonteCarloSample(TString s="SMGG");
  void SetEtaBins();
  void SetPtBins();
  void CreateOutputFile(TString st="toto_eff.root");

  //== FF method

};

#endif

