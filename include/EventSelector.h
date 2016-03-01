//This class is meant to perform the final
// selection of the analysis
// A tree is filled with the events passing the preselection
// Global Event variables (mu,NPV,runnumber,LB,...)
// and kinematics variables of the two photons of interests
// and of the diphoton system
// are stored in a TTree.
// This class inherits from GravitonAnalysis.
#ifndef EventSelector_h
#define EventSelector_h

#include <iostream>
#include <utility>
#include <vector>
#include <TTree.h>
#include <TString.h>
#include <TH2.h>
#include <TFile.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include "GravitonAnalysis.h"

class EventSelector : public GravitonAnalysis
{
 private:

  TString m_output;
  TFile*  file_out;
  TTree*  tree_out;


  //--> Cutflow histograms
  TObjArray*     hist_mgg;
  TObjArray*     hist_mgg_w;
  TH1D*          hmgg_final;
  TH1D*          hmgg_final_w;
  TH1F*          hCutFlow;
  TH1F*          hCutFlow_w;

 protected:

  //--> Variables for output tree
  int            m_RunNumber;
  int            m_LumiBlock;
  int            m_EventNumber;
  int            m_NPV;
  double         m_mu;
  double         PV_ID;

  double m_MET_Topo_et;			  
  double m_MET_Topo_sumet;		  
  double m_MET_Topo_sumet_EMB;		  
  double m_MET_Topo_sumet_EME;		  
  double m_MET_Topo_sumet_CentralReg;	  
  double m_MET_Topo_sumet_EndcapRegion;	  
  double m_MET_RefFinal_et;		  
  double m_MET_RefFinal_sumet;		  
  double m_MET_RefFinal_sumet_CentralReg; 
  double m_MET_RefFinal_sumet_EndcapRegion; 



  //--> Variables for the diphoton system 
  double         mgg;
  double         ptgg;
  double         ygg;
  double         deltaphi;
  double         costhetastar;
  double         mgg_smeareddown;
  double         ptgg_smeareddown;
  double         ygg_smeareddown;
  double         costhetastar_smeareddown;
  double         deltaphi_smeareddown;
  double         mgg_smearedup;
  double         ptgg_smearedup;
  double         ygg_smearedup;
  double         costhetastar_smearedup;
  double         deltaphi_smearedup;
  //---> Variables for the leading photon
  double         Lead_pT;
  double         Lead_pT_smearedup;
  double         Lead_pT_smeareddown;
  double         Lead_eta;
  double         Lead_eta_PV;
  double         Lead_phi;
  double         Lead_iso_ptcorr;
  double         Lead_ed_med;
  double         Lead_iso_uncor;
  double         Lead_iso;
  double         Lead_iso_mod;
  int            Lead_IsTight;
  int            Lead_IsConv;
  int            Lead_IsLoosePrime2;
  int            Lead_IsLoosePrime3;
  int            Lead_IsLoosePrime4;
  int            Lead_IsLoosePrime5;
  //---> Variables for the subleading photon
  double         SubLead_pT;
  double         SubLead_pT_smearedup;
  double         SubLead_pT_smeareddown;
  double         SubLead_eta;
  double         SubLead_eta_PV;
  double         SubLead_phi;
  double         SubLead_ed_med;
  double         SubLead_iso_ptcorr;
  double         SubLead_iso_uncor;
  double         SubLead_iso;
  double         SubLead_iso_mod;
  int            SubLead_IsTight;
  int            SubLead_IsConv;
  int            SubLead_IsLoosePrime2;
  int            SubLead_IsLoosePrime3;
  int            SubLead_IsLoosePrime4;
  int            SubLead_IsLoosePrime5;

  //--> Photons cluster variables
  std::vector<float> m_ph_cl_E;
  std::vector<float> m_ph_cl_pt;
  std::vector<float> m_ph_cl_eta;
  std::vector<float> m_ph_cl_phi;
  std::vector<float> m_ph_rawcl_E;
  std::vector<float> m_ph_rawcl_pt;
  std::vector<float> m_ph_rawcl_eta;
  std::vector<float> m_ph_rawcl_phi;
  std::vector<float> m_ph_ED_median;

  //--> Weighting variables
  double         gen_weight;
  double         m_PUweight;
  double         weight;
  double         Lead_weight;
  double         SubLead_weight;



  //--> Variable at truth level
  //--------------------------
  double Lead_truth_eta;
  double Lead_truth_phi;
  double Lead_truth_pT;
  int    Lead_truth_status;
  //--------------------------
  double SubLead_truth_eta;
  double SubLead_truth_phi;
  double SubLead_truth_pT;
  int    SubLead_truth_status;
  //--------------------------
  double truth_mgg;
  double truth_ptgg;
  double truth_costhetastar;
  double truth_ygg;
  double truth_dphigg;
  double mymc_m;
  //-------------------------

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
  EventSelector();
  EventSelector(TTree*);
  virtual ~EventSelector();
  void   SetOutputName(TString st="toto.root") {m_output=st;}

  void   EventLoop();
  bool   InitOutTree();

 private :

  void   FillOutputFile();

 protected :
  void   ComputeBasicEventInfo();
  void   ComputeKinematics(int ilead,int isublead);
  void   ComputeIdDecisions(int ilead,int isublead);
  void   ComputeIsolations(int ilead,int isublead);
  bool   EventTruthSelector(int entry,bool data);
  void   PrintParticleInfo(int index);
  void   PrintEventInfo(int limit = 1000000);
  int    GetOriginalPhotonIndex(int index);
  std::vector<int>   HardScatterPhotonIndices(std::vector<int> v_ph_id);
  bool   IsPhotonFromHardProc(int index);
  bool   IsPhotonFromUE(int index);
  bool   IsPhotonFromRadiation(int index);
  bool   IsPhotonFromFrag(int index);
  float  GetPartonETIso(int index, float cone_size, float parton_minpt_threshold);
  float  ComputeX(int iph);




 public :
  double CosThetaStar_CS( TLorentzVector v1,
			  TLorentzVector v2 );
  void   TestEventRecord(TString filename = "");

};

#endif

