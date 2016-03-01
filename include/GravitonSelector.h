//This class is meant to perform the final
// selection of the analysis
// This class inherits from GravitonAnalysis.
#ifndef GravitonSelector_h
#define GravitonSelector_h

#include <iostream>
#include <utility>
#include <vector>
#include <TTree.h>
#include <TString.h>
#include <TH2.h>
#include <TFile.h>
#include <TMath.h>
#include <TLorentzVector.h>

#include "InvMassCalculator.h"
#include "GravitonAnalysis.h"


class GravitonSelector : public GravitonAnalysis
{
 private:

  //--------------------------

  TObjArray*     hist_mgg;
  TObjArray*     hist_mgg_w;

 protected :

  TH1D*          hmgg_final;
  TH1D*          hmgg_final_w;
  TH1F*          hCutFlow;
  TH1F*          hCutFlow_w;

  TString        m_output;
  TFile*         file_out;
  TTree*         tree_out;
  //--> weight var
  double         gen_weight;
  double         m_PUweight;
  double         m_SF_ID_lead;
  double         m_SF_ID_sublead;

  //--> Variables for output tree
  int            m_RunNumber;
  int            m_LumiBlock;
  int            m_EventNumber;
  int            m_NPV;

  double         mgg;
  double         ptgg;
  double         ygg;
  double         deltaphi;
  double         costhetastar;
  double         mgg_smeareddown;
  double         mgg_smearedup;
  double         ptgg_smeareddown;
  double         ptgg_smearedup;
  double         ygg_smeareddown;
  double         ygg_smearedup;
  double         costhetastar_smeareddown;
  double         costhetastar_smearedup;
  double         mymc_m;
  double         weight;
  //.....................
  double         Lead_pT;
  double         Lead_pT_smearedup;
  double         Lead_pT_smeareddown;
  double         Lead_eta;
  double         Lead_eta_PV;
  double         Lead_phi;
  double         Lead_IsConv;
  //.....................
  double         SubLead_pT;
  double         SubLead_pT_smearedup;
  double         SubLead_pT_smeareddown;
  double         SubLead_eta;
  double         SubLead_eta_PV;
  double         SubLead_phi;
  double         SubLead_IsConv;

  //--> TLorentzVector
  TLorentzVector gamgam_lv;
  TLorentzVector Lead_lv;
  TLorentzVector SubLead_lv;

  TLorentzVector gamgam_lv_smearedup;
  TLorentzVector Lead_lv_smearedup;
  TLorentzVector SubLead_lv_smearedup;

  TLorentzVector gamgam_lv_smeareddown;
  TLorentzVector Lead_lv_smeareddown;
  TLorentzVector SubLead_lv_smeareddown;


 public :
  GravitonSelector();
  GravitonSelector(TTree*);
  virtual ~GravitonSelector();
  void   Init();
  void   EventLoop(TString outfile);

 private :
  void   InitOutTree();
  void   FillOutputFile();

 protected : 
  void   ComputeKinematicsOutput(int ilead,int isublead);
  void   ComputeBasicEventInfo();
  void   ComputeIDScaleFactors(int ilead, 
			       int isublead,
			       double ptlead,
			       double ptsublead);
 public :

  void   PrintEventInfo( int RunNumber,
			 int EventNumber ,
			 int LumiBlock ,
			 int ilead ,
			 int isublead,
			 double mggval,
			 double mggcut = 1000. );
  static double CosThetaStar_CS( TLorentzVector v1,
				 TLorentzVector v2 );
				 
};

#endif

