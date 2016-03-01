#ifndef BkgCategoryMerging_h
#define BkgCategoryMerging_h

#include <iostream>
#include <utility>
#include <vector>
#include <TString.h>
#include <TH1.h>
#include <TGraphAsymmErrors.h>
#include <RooPlot.h>

class BkgCategoryMerging{

 private:

  bool                 m_VERBOSE;
  std::vector<TString> m_filenames;
  //--------------------------------
  TString            m_type;
  //----------------------------
  TH1D*              m_hh_data;
  TGraphAsymmErrors* m_gmgg_data;
  //----------------------------
  TH1D*              m_hh_bkg;
  TGraphAsymmErrors* m_gmgg_bkg;
  //----------------------------
  TH1D*              m_hh_red;
  TGraphAsymmErrors* m_gmgg_red;
  //----------------------------
  TH1D*              m_hh_gj;
  TGraphAsymmErrors* m_gmgg_gj;
  //----------------------------
  TH1D*              m_hh_jg;
  TGraphAsymmErrors* m_gmgg_jg;
  //----------------------------
  TH1D*              m_hh_jj;
  TGraphAsymmErrors* m_gmgg_jj;
  //----------------------------
  TH1D*              m_hh_gg;
  TGraphAsymmErrors* m_gmgg_gg;
  //----------------------------
  TH1D*              m_hh_bkg_syst;

  RooPlot* frame_TiTi_L;
  RooPlot* frame_TiTi_SL;
  RooPlot* frame_PH_L;
  RooPlot* frame_PH_SL;
  RooPlot* frame_Jet_L;
  RooPlot* frame_Jet_SL;


  RooPlot* frame_lead_fake;
  RooPlot* frame_sublead_fake;
  RooPlot* frame_double_fake;

  //-----------------------------------
  TH1D*  m_hh_tot_uncert;
  TH1D*  m_hh_pur_uncert;
  TH1D*  m_hh_red_uncert;
  TH1D*  m_hh_irr_uncert;
  TH1D*  m_hh_pdf_uncert;
  TH1D*  m_hh_iso_uncert;
  TH1D*  m_hh_scale_uncert;





 public:

  BkgCategoryMerging(TString type,TString f0,
		     TString f1="",
		     TString f2="",TString f3="",
		     TString f4="",TString f5="",
		     TString f6="",TString f7="");

  virtual ~BkgCategoryMerging();
  void SetVerbose( bool v) {m_VERBOSE=v;}
  void SetFiles( TString f0,TString f1="",
		 TString f2="",TString f3="",
		 TString f4="",TString f5="",
		 TString f6="",TString f7="");


  void StoreToRootFile(TString st="toto.root");
  void StoreToRootFile_Search(TString st="toto.root");
  void StoreToRootFile_2DFit(TString st="toto.root");
  void StoreToRootFile_RedFit(TString st="toto.root");

 private:

  bool BuildHistograms();
  bool BuildGraphs();
  bool BuildRooPlots_2DFit();
  bool BuildRooPlots_RedFit();
  bool BuildHistograms_Search();


};
#endif // #ifdef BkgCategoryMerging_h
