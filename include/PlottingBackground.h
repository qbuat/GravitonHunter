#ifndef PlottingBackground_h
#define PlottingBackground_h

#include <iostream>
#include <TCanvas.h>
#include <TString.h>
#include <TFile.h>
#include <TLatex.h>

class PlottingBackground{

 private:
  TFile* m_file2dfit;
  TFile* m_fileredfit;
  TFile* m_filebkg;

  TFile* m_file2dfit_L2;
  TFile* m_file2dfit_L5;

  TFile* m_fileredfit_L3;
  TFile* m_fileredfit_L4;
  TFile* m_fileredfit_L5;

  TString m_cat;

 public:

  PlottingBackground( TString name2dfit="NONE", 
		      TString nameredfit="NONE",
		      TString namebkg="NONE" );
  virtual ~PlottingBackground();
  
  void Init2DFit(TString name2dfit);
  void InitRedFit(TString nameredfit);
  void InitBkg(TString namebkg);

  void Init2DFit_Syst(TString name_L2,TString name_L5);
  void InitRedFit_Syst(TString name_L3,TString name_L4,TString name_L5);


  void SetCategory(TString cat){m_cat=cat;}

  //--------------------------------
  TCanvas* Plot2DFitLeadJet();
  TCanvas* Plot2DFitLeadPh();
  TCanvas* Plot2DFitLeadFinal();
  TCanvas* Plot2DFitSubLeadJet();
  TCanvas* Plot2DFitSubLeadPh();
  TCanvas* Plot2DFitSubLeadFinal();
  TCanvas* Plot2DFitJetJetHist();
  TCanvas* Plot2DFitJetJetCurve();
  //--------------------------------
  TCanvas* Plot2DFitLeadPh_Syst();
  TCanvas* Plot2DFitLeadFinal_Syst();
  TCanvas* Plot2DFitSubLeadPh_Syst();
  TCanvas* Plot2DFitSubLeadFinal_Syst();
  //-----------------------------------
  TCanvas* PlotGJDataFit();
  TCanvas* PlotJGDataFit();
  TCanvas* PlotJJDataFit();
  TCanvas* PlotRedDataFit_syst(TString comp);
  //-----------------------------------
  TCanvas* PlotBkgEstimate();
  TCanvas* PlotBkgEstimate_withratiohist();
  TCanvas* PlotBkgUncertainty();

 private:

};
#endif // #ifdef PlottingBackground_h
