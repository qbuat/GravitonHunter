// This class perform a 2D template fit on
// the isolation variable
// The input are a tree containing the data 
// and a config file 
// It inherits from BkgFitFramework
// Convention: Ti==Tight;Lo==Loose;L==Leading;SL==SubLeading

#ifndef Bkg2DFit_h
#define Bkg2DFit_h

#include <iostream>
#include <utility>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2.h>
#include <THStack.h>
#include <TObject.h>
#include <TLegend.h>
#include <TLatex.h>
#include <RooRealVar.h>
#include <RooCategory.h>
#include <RooAbsReal.h>
#include <RooDataSet.h>
#include <RooAbsPdf.h>
#include <RooAddPdf.h>
#include <RooProdPdf.h>
#include <RooNDKeysPdf.h>
#include <RooKeysPdf.h>
#include <RooCBShape.h>
#include <RooNovosibirsk.h>
#include <RooGaussian.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooFormulaVar.h>
#include <RooStats/SPlot.h>

#include "BkgFitFramework.h"

class Bkg2DFit : public BkgFitFramework
{
 private:
  
  //-----------------
  bool m_do1Dfits;
  bool m_do2Dfits;
  bool m_doSplot;
  //--------------------
  RooRealVar* Ngamgam;
  RooRealVar* Ngamjet;
  RooRealVar* Njetgam;
  RooRealVar* Njetjet; 
  RooRealVar* iPHPH_TiIso;
  RooRealVar* iPHJET_TiIso;
  RooRealVar* iJETPH_TiIso;
  RooRealVar* iJETJET_TiIso;
  double      m_NgamgamYield;
  double      m_NgamjetYield;
  double      m_NjetgamYield;
  double      m_NjetjetYield;
  double      m_NgamgamYieldError;
  double      m_NgamjetYieldError;
  double      m_NjetgamYieldError;
  double      m_NjetjetYieldError;
  //------------------------------

  std::pair<double,double> m_isonormrange;

  //---------------------------------------
  RooDataSet* m_DataSet_TiTi; 
  RooDataSet* m_DataSet_TiIsoTiIso; 
  RooDataSet* m_DataSet_LoLo; 
  RooDataSet* m_DataSet_Ti_L; 
  RooDataSet* m_DataSet_Ti_SL;
  RooDataSet* m_DataSet_TiLo; 
  RooDataSet* m_DataSet_LoTi;
  //----------------------------------------
  RooDataSet* m_DataSet_TiTi_Ngamgamw; 
  RooDataSet* m_DataSet_TiTi_Ngamjetw; 
  RooDataSet* m_DataSet_TiTi_Njetgamw; 
  RooDataSet* m_DataSet_TiTi_Njetjetw; 
  RooDataSet* m_DataSet_TiIsoTiIso_Ngamgamw; 
  RooDataSet* m_DataSet_TiIsoTiIso_Ngamjetw; 
  RooDataSet* m_DataSet_TiIsoTiIso_Njetgamw; 
  RooDataSet* m_DataSet_TiIsoTiIso_Njetjetw; 

  //---------------------------
  RooAbsPdf* m_LeadJetPdf;
  RooArgSet*      m_JetL;
  //-----------------------------
  RooAbsPdf* m_SubLeadJetPdf;
  RooArgSet*      m_JetSL;
  //-------------------------
  RooNDKeysPdf*   m_JetJetPdf;
  //----------------------------
  RooAbsPdf*      m_LeadPHPdf;
  RooArgSet*      m_PHL;
  //--------------------------------
  RooAbsPdf*      m_SubLeadPHPdf;
  RooArgSet*      m_PHSL;
  //------------------------
  RooProdPdf*     m_PHPHPdf;
  RooProdPdf*     m_PHJetPdf;
  RooProdPdf*     m_JetPHPdf;
  RooAddPdf*      m_TiTiPdf;

  //-----------------------
  RooFitResult* m_FitRes_Jet_L;
  RooFitResult* m_FitRes_Jet_SL;
  RooFitResult* m_FitRes_Ti_L;
  RooFitResult* m_FitRes_Ti_SL;
  RooFitResult* m_FitRes_TiTi;
  //--------------------------
  RooPlot* frame_JetL;
  RooPlot* frame_JetSL;
  RooPlot* frame_PHL;
  RooPlot* frame_PHSL;
  TH2F*    m_hJetJet;
  TH2F*    m_hJetJetFit;

  RooPlot* frame_TiTi_L;
  RooPlot* frame_TiTi_SL;

  //-------------------
  TH1F* hmgg;
  TH1F* hmgg_Ngg;
  TH1F* hmgg_Ngj;
  TH1F* hmgg_Njg;
  TH1F* hmgg_Ngjjg;
  TH1F* hmgg_Njj;
  TH1F* hmgg_red;
  TH1F* hmgg_iso;
  TH1F* hmgg_Ngg_iso;
  TH1F* hmgg_Ngj_iso;
  TH1F* hmgg_Njg_iso;
  TH1F* hmgg_Ngjjg_iso;
  TH1F* hmgg_Njj_iso;
  TH1F* hmgg_red_iso;

 public :

  Bkg2DFit();
  Bkg2DFit( TTree*tree,
	    std::pair<double,double> mggbin,
	    TString FitParFile );
  virtual ~Bkg2DFit();

  void   SetIsoNorm(double Xmin = 10.,double Xmax=25.);

  void   Init( std::pair<double,double> mggbin,
	       TString FitParFile );
  //--------------------------------------
  void   LeadingJetFit(bool doFit=true);
  void   SubLeadingJetFit(bool doFit=true);
  void   JetJetFit(bool doFit=true);
  void   LeadPhotonFit(bool doFit=true);
  void   SubLeadPhotonFit(bool doFit=true);
  void   TightTightFit( bool doFit=true,
			bool doFrame=true );
  //--------------------------------------
  void   Fitter( bool do1Dfits=true,
		 bool do2Dfits=true,
		 bool doSplot=false );
  void   ExtractShapesFrom_sPlot();
  void   PlotResults();
  void   YieldsExtraction();
  void   RandomizeFitResults(int NPE);
  void   RandomizeYieldsResults(int NPE,
				TString out = "none");
  //----------------------------
  void   GetChiSquares(int Npar);
  double GetPurity();
  double GetPurityError();
  double GetSumofYields();
  double GetSumofYields_NoIso();
  double GetSumofYieldsError();
  double GetSumofYieldsError_NoIso();
  double GetNgamgamYield();
  double GetNgamjetYield();
  double GetNjetgamYield();
  double GetNjetjetYield();
  double GetNgamgamYieldError();
  double GetNgamjetYieldError();
  double GetNjetgamYieldError();
  double GetNjetjetYieldError();
  double GetNgamgamYield_NoIso();
  double GetNgamjetYield_NoIso();
  double GetNjetgamYield_NoIso();
  double GetNjetjetYield_NoIso();
  double GetNgamgamYieldError_NoIso();
  double GetNgamjetYieldError_NoIso();
  double GetNjetgamYieldError_NoIso();
  double GetNjetjetYieldError_NoIso();
  TH1F*  GetDataHist();
  TH1F*  GetIsoDataHist();
  RooFitResult* GetTiTiFitResult();
  //---------------------------
  void   Test2DIntegral();
  void   StoreToRootFile(TString out="none");


};

#endif

