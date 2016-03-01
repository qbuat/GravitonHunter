/// This class provides a framework to perform
/// the bkg fit (invariant mass and isolation)
/// using the RooFit toolkit
/// the input is a TTree with the interesting events
/// and the useful branches of the two photons 
/// (pt,eta,isolation,looseprime and tight decisions) 
/// and the mass of diphoton system

#ifndef BkgFitFramework_h
#define BkgFitFramework_h

#include <iostream>
#include <utility>
#include <TTree.h>
#include <TString.h>
#include <RooRealVar.h>
#include <RooCategory.h>
#include <RooDataSet.h>

class BkgFitFramework
{


 protected:
  
  TTree*              m_tree;
  TString             m_fitparfile;
  TString             m_etacat;
  int                 m_looseprime;
  std::pair<double,double> m_mggbin;
  bool                m_ptdepisocut;

  //------------------------------
  int    m_NPE;
  double m_Isomin;
  double m_Isomax;
  double m_Mggmin;
  double m_Mggmax;

  double m_IsoCut_L;
  double m_IsoCut_SL;
  double m_PtCut_L;
  double m_PtCut_SL;

  TString m_etacut_st;
  TString m_ptcut_st;
  TString m_mggcut_st;
  TString m_isocut_st;

  RooRealVar*  m_pt_L;
  RooRealVar*  m_pt_SL;
  RooRealVar*  m_eta_L;
  RooRealVar*  m_eta_SL;
  RooRealVar*  m_Iso_L;
  RooRealVar*  m_Iso_SL;
  RooRealVar*  m_mgg;
  RooCategory* m_IsTight_L;
  RooCategory* m_IsTight_SL;
  RooCategory* m_IsLoosePrime_L;
  RooCategory* m_IsLoosePrime_SL;
  //---------------------------------------

  RooDataSet* m_DataSet;
  RooDataSet* m_DataSet_looseprimecut;
  RooDataSet* m_DataSet_etacut;
  RooDataSet* m_DataSet_ptcut;
  RooDataSet* m_DataSet_mggcut;
  //----------------------------------------

 public :

  BkgFitFramework();
  BkgFitFramework( TTree*tree,
		   std::pair<double,double> mggbin,
		   TString FitParFile );
  virtual ~BkgFitFramework();
  void   Init( std::pair<double,double> mggbin,
	       TString FitParFile );


  /************* SETTERS METHODS **********************/
  void   SetTree(TTree* tree,int entriesforfit=-1);
  void   SetMggBounds(double min=0.,double max=5000.);
  void   SetIsoBounds(double min=-100.,double max=3000.);
  void   SetEtaCategory(TString cat="NONE"){m_etacat=cat;}
  void   SetPhotonsIsoCut(double l_cut=5.,double sl_cut=5.);
  void   SetPhotonsPtCut(double l_cut=40,double sl_cut=30);
  void   SetParametersFile(TString st="none") {m_fitparfile=st;}
  void   SetNPE(int NPE=1000) {m_NPE = NPE;}
  void   SetLoosePrimeType(int type=0){m_looseprime=type;}
  void   SetPtDependentIsocut(bool b){m_ptdepisocut=b;}
  //--------------------------------------

  static void SetConstantParameters(const RooArgSet& set);

};

#endif

