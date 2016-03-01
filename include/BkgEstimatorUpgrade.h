#ifndef BkgEstimatorUpgrade_h
#define BkgEstimatorUpgrade_h

#include <iostream>
#include <utility>
#include <vector>
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
#include <RooPlot.h>


class BkgEstimatorUpgrade{
 private:

  bool               m_VERBOSE;
  bool               m_doSys;
  double             dataNorm;
  //--------------------------------
  TString            name_f_data;
  TString            name_f_red;
  TString            name_f_red_L4;
  TString            name_f_red_L5;
  TString            name_f_irr;
  TString            name_f_irr_syst;
  TString            name_f_frac_L2;
  TString            name_f_frac_L3;
  TString            name_f_frac_L4;
  TString            name_f_frac_L5;
  //----------------------------
  TH1D*              m_hh_data;
  TGraphAsymmErrors* m_gmgg_data;
  //----------------------------
  TH1D*              m_hh_bkg;
  TH1D*              m_hh_red;
  TH1D*              m_hh_irr;
  TH1D*              m_hh_gg;
  TH1D*              m_hh_gj;
  TH1D*              m_hh_jg;
  TH1D*              m_hh_jj;
  //---------------------------------

  TH1D*              m_hh_irr_syst;
  TH1D*              m_hh_scale_syst;
  TH1D*              m_hh_iso_syst;
  TH1D*              m_hh_isodatamc_syst;
  TH1D*              m_hh_pdf_syst;
  TH1D*              m_hh_red_syst;
  TH1D*              m_hh_gj_syst;
  TH1D*              m_hh_jg_syst;
  TH1D*              m_hh_jj_syst;
  TH1D*              m_hh_tot_syst;

  TH1D*              m_hh_pur_syst;
  TH1D*              m_hh_fracsyst_syst;
  TH1D*              m_hh_fracstat_syst;
  //--------------------------------


  TH1D*              m_hh_bkg_withsyst;
  TH1D*              m_hh_red_withsyst;
  TGraphAsymmErrors* m_gmgg_bkg;
  TGraphAsymmErrors* m_gmgg_red;

  //--------------------------------
  TH1D*              m_hh_data_search;
  TH1D*              m_hh_bkg_search;
  TH1D*              m_hh_irr_search;
  TH1D*              m_hh_red_search;
  TH1D*              m_hh_tot_syst_search;
  TH1D*              m_hh_pur_syst_search;
  TH1D*              m_hh_irr_syst_search;
  TH1D*              m_hh_red_syst_search;
  //---------------------------------

  TCanvas* canv_YieldsSyst;
  TCanvas* canv_YieldsStat;
  TCanvas* canv_RedShape;  
  TCanvas* canv_IrrShape;  
  RooPlot* frame_Purity;

  //---------------------------------

 public:
  BkgEstimatorUpgrade();
  virtual ~BkgEstimatorUpgrade();

  //--------- Setters---------------------------
  void SetDataFile(TString st) {name_f_data=st;}
  void SetIrrFile(TString st) {name_f_irr = st;}
  void SetIrrSystFile(TString st) {name_f_irr_syst = st;}
  void SetRedFile(TString st) {name_f_red = st;}
  void SetRedSystFiles(TString st_L4,TString st_L5);
  void SetYieldSystFiles(TString st_L2,
			 TString st_L3,
			 TString st_L5) ;
  void SetYieldFile(TString st) {name_f_frac_L4=st;}
  void SetVerbose(bool verb) {m_VERBOSE=verb;}
  void SetUncertFlag(bool b) {m_doSys=b;}

  void Init();

  TCanvas* YieldsSystUncert();
  TCanvas* YieldsStatUncert();
  TCanvas* RedShapeUncert();
  TCanvas* IrrShapeUncert();

  void PrintLoosePrimeResults();
  //--------- Getters---------------------------
  void GetFinalYieldsPerBin();
  void GetFinalYieldsPerBin_NoUncert();
  void CreateBATInput(TString st);
  void CreateBumpHunterInput(TString st);

  std::pair<double,double> GetYield(TString type,
				    TString treefile,
				    bool doplot = false);
  std::pair<double,double> GetPurityInterval(TString treefile,
					     bool doplot=false);



  // private:

  void  RemoveReducibleStatErrors();
  void  BuildHistSearch();
  void  BuildTotalBkg();
  void  BuildTotalBkgUncertainty();
  void  BuildBkgErrorsGraph();
  //  void  BuildDataErrorsGraph();


  TH1D* GetTotalBkgHist( TH1D* hgg,TH1D* hgj,
			 TH1D* hjg,TH1D* hjj,
			 const std::vector<double> yields);

  TH1D* GetTotalBkgHist( TH1D* hgg,TH1D* hgj,
			 TH1D* hjg,TH1D* hjj,
			 const std::vector<double> yields,TString outname);


  TH1D* GetRedBkgHist( TH1D* hgj,TH1D* hjg,TH1D* hjj,
		       const std::vector<double> yields);
  TH1D* GetRedShapeSyst( TString bkgtype ); 
  TH1D* GetScaleShapeSyst();
  TH1D* GetIsoShapeSyst();
  TH1D* GetIsoDataMCShapeSyst();
  TH1D* GetPDFShapeSyst();

};
#endif // #ifdef BkgEstimatorUpgrade_h
