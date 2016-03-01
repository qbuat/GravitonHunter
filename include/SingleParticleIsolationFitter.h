
#ifndef SingleParticleIsolationFitter_h
#define SingleParticleIsolationFitter_h

#include <iostream>
#include <utility>
#include <TTree.h>
#include <TString.h>
#include <TH2.h>
#include <TCanvas.h>

#include <RooRealVar.h>
#include <RooCategory.h>
#include <RooDataSet.h>
#include <RooAbsPdf.h>
#include <RooFitResult.h>
#include <RooPlot.h>


class SingleParticleIsolationFitter
{
 private:
  
  TTree*    m_tree;
  int       m_nentries;
  int       m_nentriesforfit;
  TString   m_FitParFile;
  bool      m_dofit;
  bool      m_data;
  //------------------------------
  std::pair<double,double> m_isorange;
  std::pair<double,double> m_isonormrange;
  std::pair<double,double> m_ptrange;
  std::pair<double,double> m_etarange;
  int                 m_conv;
  //-------------------------
  TString m_pt_name;
  TString m_eta_name;
  TString m_iso_name;
  TString m_istight_name;
  TString m_isconv_name;

  //--------------------------
  RooRealVar*  m_pt;
  RooRealVar*  m_eta;
  RooRealVar*  m_iso;
  RooRealVar*  m_istight;
  RooRealVar*  m_isconv;
  //--------------------------
  RooDataSet* m_DataSet;
  RooDataSet* m_DataSet_cut;
  RooDataSet* m_DataSet_Ti;
  RooDataSet* m_DataSet_Lo; 
  //----------------------------------------

  //-----------------------------
  RooArgSet*      m_PHSet;
  RooAbsPdf*      m_TightPdf;
  //----------------------------
  RooAbsPdf*      m_PHPdf;
  RooRealVar*     meanCB;
  RooRealVar*     sigmaCB;
  RooRealVar*     alphaCB;
  RooRealVar*     nCB;
  //-----------------------
  RooFitResult* m_FitRes_Ti;
  //--------------------------
  RooPlot* frame_PH;

  //------------------
  TH1F*  m_h_iso;


 public :

  SingleParticleIsolationFitter();
  SingleParticleIsolationFitter( TTree *tree, TString ParFile );
  virtual ~SingleParticleIsolationFitter();

  //--------- Setters ---------------------------------
  void  SetEntries(TTree* tree,int entries=-1);
  void  SetTree(TTree* tree);
  void  SetIsoBounds(double Xmin=-10.,double Xmax=25.);
  void  SetPtBounds(double Xmin=1.,double Xmax=10000.);
  void  SetEtaBounds(double Xmin=0.,double Xmax=3.);
  void  SetIsConv(int conv=0);
  void  SetParFile(TString ParFile="");
  void  SetPtEtaIsoTightConvNames( TString,
				   TString,
				   TString,
				   TString,
				   TString);

  //---------------------
  void   Init_Vars();
  void   Init_DataSets();

  //--------------------------------------
  void   PhotonFit(bool dofit=true);
  void   Fitter( bool dofit=true);

  //--------------------------------------

  //------ Getters --------------
  bool   GetStreamType() {return m_data;}
  TH1F*  GetPhotonHist() {return m_h_iso;}
  double GetPhotonFitMean() {return meanCB->getVal();}
  double GetPhotonFitMeanError() {return meanCB->getError();}
  double GetPhotonFitWidth() {return sigmaCB->getVal();}
  double GetPhotonFitWidthError() {return sigmaCB->getError();}
  double GetPhotonHistMean() {return m_h_iso->GetMean();}
  double GetPhotonHistRMS() {return m_h_iso->GetRMS();}
  //----------------------------
  void StoreToRootFile(TString st);
  static TH1F*  RooHistToHist(const RooHist& RH);

 private :
  void   SetConstantParameters(const RooArgSet& set);

};

#endif

