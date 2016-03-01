#ifndef BkgEstimator_h
#define BkgEstimator_h

#include <iostream>
#include <utility>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TLatex.h>
#include <TLine.h>
#include <TH1.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPad.h>
#include <TGraphAsymmErrors.h>
#include <TObjArray.h>
#include <TMath.h>
#include "ToolsTreeReader.h"
#include "ToolsUtilities.h"


class BkgEstimator{
 private:

  bool               m_VERBOSE;
  /* int                nBins400; */
  /* double*            binning400_log; */
  double             dataNorm;
  //--------------------------------
  TString            name_f_data;
  TString            name_f_red;
  TString            name_f_red_syst;
  TString            name_f_irr;
  TString            name_f_irr_syst;
  TString            name_f_frac_L2;
  TString            name_f_frac_L3;
  TString            name_f_frac_L4;
  TString            name_f_frac_L5;
  TFile*             m_f_data;
  TFile*             m_f_red;
  TFile*             m_f_irr;
  TFile*             m_f_irr_syst;
  TFile*             m_f_frac_L2;
  TFile*             m_f_frac_L3;
  TFile*             m_f_frac_L4;
  TFile*             m_f_frac_L5;
  TFile*             m_f_RS_1250_01;
  TFile*             m_f_RS_1500_01;
  TFile*             m_f_RS_1750_01;
  TFile*             m_f_ADD_2500_GRW;
  TFile*             m_f_ADD_3000_GRW;
  TFile*             m_f_ADD_3500_GRW;
  //----------------------------
  TH1D*              m_hh_data;
  TGraphAsymmErrors* m_gmgg_data;
  //----------------------------
  TGraphAsymmErrors* m_gr_gj;
  TGraphAsymmErrors* m_gr_jg;
  TGraphAsymmErrors* m_gr_jj;
  TH1D*              m_hh_red;
  TH1D*              m_hh_gg;
  TH1D*              m_hh_gj;
  TH1D*              m_hh_jg;
  TH1D*              m_hh_jj;
  //---------------------------------
  TH1D*              m_hh_irr;
  TH1D*              m_hh_irr_syst;
  TH1D*              m_hh_red_syst;
  TH1D*              m_hh_gj_syst;
  TH1D*              m_hh_jg_syst;
  TH1D*              m_hh_jj_syst;
  TH1D*              m_hh_syst;
  TH1D*              m_hh_frac_syst;
  TH1D*              m_hh_frac_stat;
  TH1D*              m_hh_frac_tot;
  //---------------------------------
  TTree*             m_yieldstree_L2;
  TTree*             m_yieldstree_L3;
  TTree*             m_yieldstree_L4;
  TTree*             m_yieldstree_L5;
  //--------------------------------
  TH1D*              m_hh_bkg;
  TH1D*              m_hh_bkg_syst;
  TGraphAsymmErrors* m_gmgg_bkg;
  TGraphAsymmErrors* m_gmgg_red;
  TH1F*              m_hh_signi;
  TH1F*              m_hh_signi_sys;
  //--------------------------------
  TH1D*              m_hh_bkg_400;
  TH1D*              m_hh_data_400;
  TH1D*              m_hh_syst_400;
  TH1D*              m_hh_frac_syst_400;
  TH1D*              m_hh_irr_syst_400;
  TH1D*              m_hh_red_syst_400;
  TH1D*              m_hh_irr_400;
  TH1D*              m_hh_red_400;
  //---------------------------------
  THStack*           m_hs_RS_1250_01;
  THStack*           m_hs_RS_1500_01;
  THStack*           m_hs_RS_1750_01;
  TH1F*              m_hh_RS_1250_01;
  TH1F*              m_hh_RS_1500_01;
  TH1F*              m_hh_RS_1750_01;
  //---------------------------------
  THStack*           m_hs_ADD_2500_GRW;
  THStack*           m_hs_ADD_3000_GRW;
  THStack*           m_hs_ADD_3500_GRW;
  TH1F*              m_hh_ADD_2500_GRW;
  TH1F*              m_hh_ADD_3000_GRW;
  TH1F*              m_hh_ADD_3500_GRW;
  //---------------------------------

 public:
  BkgEstimator();
  virtual ~BkgEstimator();
  void SetDataFile(TString st) {name_f_data=st;}
  void SetIrrFile(TString st) {name_f_irr = st;}
  void SetIrrSystFile(TString st) {name_f_irr_syst = st;}
  void SetRedFile(TString st) {name_f_red = st;}
  void SetRedSystFile(TString st) {name_f_red_syst = st;}
  void SetYieldSystFiles(TString st_L2,
			 TString st_L3,
			 TString st_L5) ;
  void SetYieldFile(TString st) {name_f_frac_L4=st;}
  void SetVerbose(bool verb) {m_VERBOSE=verb;}
  void Init();
  void YieldsSystUncert(bool doplot=false);
  void YieldsStatUncert(bool doplot=false);
  void RedShapeUncert(bool doplot=false);
  void PlotBkgEstimate(bool doRS=false,
		       bool doADD=false,
		       bool doRSADD=false );
  /* void PlotBkgWithSigni(bool dosave); */
  void PlotBkgUncertainty();
  void GetFinalYieldsPerBin();
  void GetFinalYieldsPerBin_NoUncert();
  void CreateBATInput(TString st);
  void CreateBumpHunterInput(TString st);

  std::pair<double,double> GetYield(TString type,
				    TTree* yieldstree,
				    bool doplot = false);

 private:
  void  CreateBinning();
  void  RemoveReducibleStatErrors();
  void  BuildHist400();
  void  BuildRSHists();
  void  BuildADDHists();
  void  BuildTotalBkg();
  void  BuildTotalBkgUncertainty();
  void  BuildTotalBkg400();
  void  BuildBkgErrorsGraph();
  void  BuildDataErrorsGraph();
  void  BuildSignificanceHist();
  TH1D* GetRedShapeSyst( TString bkgtype,
			 bool doplot=false ); 

  TH1D* GetTotalBkgHist( TH1D* hgg,TH1D* hgj,
			 TH1D* hjg,TH1D* hjj,
			 const std::vector<double> yields);
  TH1D* GetRedBkgHist( TH1D* hgj,TH1D* hjg,TH1D* hjj,
		       const std::vector<double> yields);



};
#endif // #ifdef BkgEstimator_h
