#ifndef BkgKinematics_h
#define BkgKinematics_h

#include <iostream>
#include <map>
#include <utility>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <TH1.h>
#include <TMath.h>
#include <TCanvas.h>

class BkgKinematics{

 private:

  TString m_filedata;
  TString m_filemcgg;
  TString m_etacategory;
  double  m_isocut;
  double m_ptcut_l;
  double m_ptcut_sl;
  double m_purmod;
  //---------------------------------
  TH1D* m_hmgg;
  std::map<TString,TString> m_hname;
  std::map<TString,TString> m_title;
  std::map<TString,int> m_Nbins;
  std::map<TString,double> m_Xmin;
  std::map<TString,double> m_Xmax;
  
  std::map<TString,TH1D*> m_hD;
  std::map<TString,TH1D*> m_hB;
  std::map<TString,TH1D*> m_hI;
  std::map<TString,TH1D*> m_hGJ;
  std::map<TString,TH1D*> m_hJG;
  std::map<TString,TH1D*> m_hJJ;
  //----------------------------------


 public:

  BkgKinematics( TString histsfile);
  BkgKinematics( TString datafile,TString mcfile);
  BkgKinematics( TString datafile,TString mcfile,
		 TString etacat );
  virtual ~BkgKinematics();
  void SetDataFile( TString st){m_filedata=st;}
  void SetMCFile( TString st){m_filemcgg=st;}
  void SetEtaCategory(TString st){m_etacategory=st;}
  void SetIsoCut(double cut) {m_isocut=cut;}
  void SetPtCuts(double c_l,double c_sl){m_ptcut_l=c_l;m_ptcut_sl=c_sl;}
  void SetPurityModifier(double a=0){m_purmod =a;}
  
  void BuildTotalBkgHist(TString var);
  TCanvas* PurityFitter(TString var);
  TCanvas* DrawHists(TString var);
  void SaveToRootFile(TString filename);

 private:

  void InitMaps();
  void InitFromHistsFile(TString histsfile);
  void FillDataHists();
  void FillMCHists();

};
#endif // #ifdef BkgKinematics_h
