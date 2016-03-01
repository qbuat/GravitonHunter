#ifndef BkgShapeBinnedExtractor_h
#define BkgShapeBinnedExtractor_h

#include <iostream>
#include <utility>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH2.h"
#include "THStack.h"
#include "TObject.h"
#include "TLegend.h"
#include "TLatex.h"


class BkgShapeBinnedExtractor
{
 private:
  TFile*      m_file;
  TString     m_filename;
  TTree*      m_tree;
  TFile*      m_file_sysup;
  TString     m_filename_sysup;
  TTree*      m_tree_sysup;
  TFile*      m_file_sysdown;
  TString     m_filename_sysdown;
  TTree*      m_tree_sysdown;
  TString     m_ParamFile;
  bool        m_doSys;
  //-------------------------
  std::vector< std::pair<double,double> > bins;
  TH1F* hmgg;
  TH1F* hmgg_tot;
  TH1F* hmgg_gg;
  TH1F* hmgg_gj;
  TH1F* hmgg_jg;
  TH1F* hmgg_jj;
  TH1F* hmgg_gjjg;
  TH1F* hmgg_tot_sysup;
  TH1F* hmgg_gg_sysup;
  TH1F* hmgg_gj_sysup;
  TH1F* hmgg_jg_sysup;
  TH1F* hmgg_jj_sysup;
  TH1F* hmgg_gjjg_sysup;
  TH1F* hmgg_tot_sysdown;
  TH1F* hmgg_gg_sysdown;
  TH1F* hmgg_gj_sysdown;
  TH1F* hmgg_jg_sysdown;
  TH1F* hmgg_jj_sysdown;
  TH1F* hmgg_gjjg_sysdown;
  //--------------------------
  TH1F* hmgg_iso;
  TH1F* hmgg_tot_iso;
  TH1F* hmgg_gg_iso;
  TH1F* hmgg_gj_iso;
  TH1F* hmgg_jg_iso;
  TH1F* hmgg_jj_iso;
  TH1F* hmgg_gjjg_iso;
  TH1F* hmgg_tot_iso_sysup;
  TH1F* hmgg_gg_iso_sysup;
  TH1F* hmgg_gj_iso_sysup;
  TH1F* hmgg_jg_iso_sysup;
  TH1F* hmgg_jj_iso_sysup;
  TH1F* hmgg_gjjg_iso_sysup;
  TH1F* hmgg_tot_iso_sysdown;
  TH1F* hmgg_gg_iso_sysdown;
  TH1F* hmgg_gj_iso_sysdown;
  TH1F* hmgg_jg_iso_sysdown;
  TH1F* hmgg_jj_iso_sysdown;
  TH1F* hmgg_gjjg_iso_sysdown;
  //---------------------------

 public :

  BkgShapeBinnedExtractor(TString,TString,TString,bool);
  virtual ~BkgShapeBinnedExtractor();
  void Init();
  void FillHists();
  void FillHists_sysup();
  void FillHists_sysdown();
  void StoretoRootFile();
};

#endif

