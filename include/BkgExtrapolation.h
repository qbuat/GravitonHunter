// This class perform the mgg fit on the 
// data control samples (reverse-id method)
// It inherits from BkgFitFramework
#ifndef BkgExtrapolation_h
#define BkgExtrapolation_h

#include <TFile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TText.h>
#include <TObjArray.h>
#include <TGraphAsymmErrors.h>
#include <RooRealVar.h>
#include <RooCategory.h>
#include <RooAbsPdf.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooAbsPdf.h>

#include "BkgFitFramework.h"

class BkgExtrapolation : public BkgFitFramework
{
 private:
  
  double              m_Xmin;
  double              m_Xmax;
  bool                m_dofits;
  bool                m_doextrapol;

  ////////////////////////////////////////////// 
  RooDataSet*         m_DataSet_IsoIso; 
  RooDataSet*         m_DataSet_TiTi; 
  RooDataSet*         m_DataSet_LoLo; 
  RooDataSet*         m_DataSet_TiLo; 
  RooDataSet*         m_DataSet_LoTi;
  RooDataSet*         m_DataSet_n_fake;//evan notation

  /////////////////////////////////////////////////
  RooAbsPdf* m_double_fake_pdf;
  RooAbsPdf* m_lead_fake_pdf;
  RooAbsPdf* m_sublead_fake_pdf;
  RooAbsPdf* m_n_fake_pdf;
  /////////////////////////////////////////////////
  RooArgSet* m_double_fake_params;
  RooArgSet* m_lead_fake_params;
  RooArgSet* m_sublead_fake_params;
  RooArgSet* m_n_fake_params;
  ////////////////////////////////////////////////
  RooFitResult* m_FitRes_double_fake;
  RooFitResult* m_FitRes_lead_fake;
  RooFitResult* m_FitRes_sublead_fake;
  RooFitResult* m_FitRes_n_fake;
  ////////////////////////////////////////////////
  RooPlot* frame_double_fake;
  RooPlot* frame_lead_fake;
  RooPlot* frame_sublead_fake;
  RooPlot* frame_n_fake;
  /////////////////////////////////////////////
  TH1D* h_n_fake;
  TH1D* h_lead_fake;
  TH1D* h_sublead_fake;
  TH1D* h_double_fake;
  ///////////////////////////////////////////////
  TGraphAsymmErrors*  graph_n_fake;
  TGraphAsymmErrors*  graph_lead_fake;
  TGraphAsymmErrors*  graph_sublead_fake;
  TGraphAsymmErrors*  graph_double_fake;
  /////////////////////////////////////////////

 public :

  BkgExtrapolation();
  BkgExtrapolation( TTree*tree,
		    std::pair<double,double> mggbin,
		     TString FitParFile );
  virtual ~BkgExtrapolation();

  void  Init( std::pair<double,double> mggbin, 
	      TString FitParFile );
  //--------------- Setters --------------------
  void  SetFitRange( double Xmin=140.,
		     double Xmax=1000. );
  //--------------- Fitters ------------------
  void   Double_fake_fit(bool dofit);
  void   Lead_fake_fit(bool dofit);
  void   Sublead_fake_fit(bool dofit);
  void   n_fake_fit(bool dofit);
  void   Fitter(bool dofits);
  //--------------- Mgg builders ------------
  void   Randomizer(bool VERBOSE=false);
  void   Extrapolate(bool doextrapol);
  //-----------------------------------------

  void   PlotResults(bool dosave=false);
  void   GetChiSquares(int Npar=0);
  void   StoreToRootFile(TString filename);


 private : 

  TH1D* BuildHistogram( const RooAbsPdf& pdf,
			TString name="" );


  TGraphAsymmErrors* FitErrorPropagator( const TH1D& hist,
					 const RooAbsPdf& pdf,
					 const RooAbsData& data,
					 const RooFitResult& fitres,
					 bool VERBOSE = false );


};

#endif
