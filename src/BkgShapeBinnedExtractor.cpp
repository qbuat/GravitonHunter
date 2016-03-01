#define BkgShapeBinnedExtractor_cxx
#include "BkgShapeBinnedExtractor.h"
#include "Bkg2DFit.h"
#include "TError.h"
#include "ToolsTreeReader.h"
///////////////////////////////////////////////////////////////////////////
BkgShapeBinnedExtractor::BkgShapeBinnedExtractor(TString filename,
						 TString filename_sysup,
						 TString filename_sysdown,
						 bool doSys)
//////////////////////////////////////////////////////////////////////////
{
  //Default constructor
  m_filename         = filename;
  m_filename_sysup   = filename_sysup;
  m_filename_sysdown = filename_sysdown;
  m_doSys            = doSys;
  m_ParamFile        = "FitParameters/IsoForFitStudy.config";
  Init();

}
////////////////////////////////////
void BkgShapeBinnedExtractor::Init()
////////////////////////////////////
{
  m_file = new TFile(m_filename,"read");
  m_tree = (TTree*)m_file->Get("tree");
  if(m_doSys){
    m_file_sysup   = new TFile(m_filename_sysup,"read");
    m_tree_sysup   = (TTree*)m_file_sysup->Get("tree");
    m_file_sysdown = new TFile(m_filename_sysdown,"read");
    m_tree_sysdown = (TTree*)m_file_sysdown->Get("tree");
  }

  int    Nbins     = 13;
  double binning[] = {140,160,180,200,220,240,260,
		      280,300,320,340,360,380,400};

  for( int ibin = 0; ibin <Nbins;ibin++){
    std::pair<double,double> mypair;
    mypair.first = binning[ibin];
    mypair.second = binning[ibin+1];
    bins.push_back( mypair );
  }

  hmgg              = new TH1F("hmgg","mgg data",Nbins,binning);
  hmgg_tot          = new TH1F("hmgg_tot","mgg tot",Nbins,binning);
  hmgg_gg           = new TH1F("hmgg_gg","mgg gg",Nbins,binning);
  hmgg_gj           = new TH1F("hmgg_gj","mgg gj",Nbins,binning);
  hmgg_jg           = new TH1F("hmgg_jg","mgg jg",Nbins,binning);
  hmgg_jj           = new TH1F("hmgg_jj","mgg jj",Nbins,binning);
  hmgg_gjjg         = new TH1F("hmgg_gjjg","mgg gj&jg",Nbins,binning);
  hmgg_tot_sysup    = new TH1F("hmgg_tot_sysup","mgg tot sysup",Nbins,binning);
  hmgg_gg_sysup     = new TH1F("hmgg_gg_sysup","mgg gg sysup",Nbins,binning);
  hmgg_gj_sysup     = new TH1F("hmgg_gj_sysup","mgg gj sysup",Nbins,binning);
  hmgg_jg_sysup     = new TH1F("hmgg_jg_sysup","mgg jg sysup",Nbins,binning);
  hmgg_jj_sysup     = new TH1F("hmgg_jj_sysup","mgg jj sysup",Nbins,binning);
  hmgg_gjjg_sysup   = new TH1F("hmgg_gjjg_sysup","mgg gj&jg sysup",Nbins,binning);
  hmgg_tot_sysdown  = new TH1F("hmgg_tot_sysdown","mgg tot sysdown",Nbins,binning);
  hmgg_gg_sysdown   = new TH1F("hmgg_gg_sysdown","mgg gg sysdown",Nbins,binning);
  hmgg_gj_sysdown   = new TH1F("hmgg_gj_sysdown","mgg gj sysdown",Nbins,binning);
  hmgg_jg_sysdown   = new TH1F("hmgg_jg_sysdown","mgg jg sysdown",Nbins,binning);
  hmgg_jj_sysdown   = new TH1F("hmgg_jj_sysdown","mgg jj sysdown",Nbins,binning);
  hmgg_gjjg_sysdown = new TH1F("hmgg_gjjg_sysdown","mgg gj&jg sysdown",Nbins,binning);
  //----------------------------------------------------------------------------------
  hmgg_iso              = new TH1F("hmgg_iso","mgg iso data",Nbins,binning);
  hmgg_tot_iso          = new TH1F("hmgg_tot_iso","mgg iso tot",Nbins,binning);
  hmgg_gg_iso           = new TH1F("hmgg_gg_iso","mgg iso gg",Nbins,binning);
  hmgg_gj_iso           = new TH1F("hmgg_gj_iso","mgg iso gj",Nbins,binning);
  hmgg_jg_iso           = new TH1F("hmgg_jg_iso","mgg iso jg",Nbins,binning);
  hmgg_jj_iso           = new TH1F("hmgg_jj_iso","mgg iso jj",Nbins,binning);
  hmgg_gjjg_iso         = new TH1F("hmgg_gjjg_iso","mgg iso gj&jg",Nbins,binning);
  hmgg_tot_iso_sysup    = new TH1F("hmgg_tot_iso_sysup","mgg tot iso sysup",Nbins,binning);
  hmgg_gg_iso_sysup     = new TH1F("hmgg_gg_iso_sysup","mgg gg iso sysup",Nbins,binning);
  hmgg_gj_iso_sysup     = new TH1F("hmgg_gj_iso_sysup","mgg gj iso sysup",Nbins,binning);
  hmgg_jg_iso_sysup     = new TH1F("hmgg_jg_iso_sysup","mgg jg iso sysup",Nbins,binning);
  hmgg_jj_iso_sysup     = new TH1F("hmgg_jj_iso_sysup","mgg jj isosysup",Nbins,binning);
  hmgg_gjjg_iso_sysup   = new TH1F("hmgg_gjjg_iso_sysup","mgg gj&jg iso_sysup",Nbins,binning);
  hmgg_tot_iso_sysdown  = new TH1F("hmgg_tot_iso_sysdown","mgg tot iso sysdown",Nbins,binning);
  hmgg_gg_iso_sysdown   = new TH1F("hmgg_gg_iso_sysdown","mgg gg iso sysdown",Nbins,binning);
  hmgg_gj_iso_sysdown   = new TH1F("hmgg_gj_iso_sysdown","mgg gj iso sysdown",Nbins,binning);
  hmgg_jg_iso_sysdown   = new TH1F("hmgg_jg_iso_sysdown","mgg jg iso sysdown",Nbins,binning);
  hmgg_jj_iso_sysdown   = new TH1F("hmgg_jj_iso_sysdown","mgg jj iso sysdown",Nbins,binning);
  hmgg_gjjg_iso_sysdown = new TH1F("hmgg_gjjg_iso_sysdown","mgg gj&jg iso sysdown",Nbins,binning);
  //----------------------------------------------------------------------------------------------
}
/////////////////////////////////////////
void BkgShapeBinnedExtractor::FillHists()
/////////////////////////////////////////
{
  int Nbins = bins.size();
  std::cout << "Number of bins = " << Nbins << std::endl;
  Bkg2DFit* RB2DF = new Bkg2DFit(m_tree,bins[0],m_ParamFile);

  //---- Fill data histogram ----//
  int entries = (int)m_tree->GetEntriesFast();
  TreeReader Rd(m_tree);
  for( int entry=0 ; entry< entries ; entry++ ){
    Rd.GetEntry(entry);
    double Iso_L  = Rd.GetVariable("Iso_L");
    double Iso_SL = Rd.GetVariable("Iso_SL");
    int Tight_L   = Rd.GetVariable("IsTight_L");
    int Tight_SL  = Rd.GetVariable("IsTight_SL");
    double mgg    = Rd.GetVariable("mgg");
    hmgg->Fill(mgg);
    if( Tight_L==1 && Tight_SL==1 && 
	Iso_L < 5 && Iso_SL<5 )
      hmgg_iso->Fill(mgg);
  }
  for(int ibin=0 ; ibin < Nbins ; ibin++){
    bool do1Dfits    = true;
    bool do2Dfits    = true;
    bool doSplot     = false;
    RB2DF->Init(bins[ibin],m_ParamFile);
    RB2DF->Fitter(do1Dfits,do2Dfits,doSplot);
    RB2DF->PlotResults();
    hmgg_tot->SetBinContent(ibin+1,RB2DF->GetSumofYields_NoIso()); 
    hmgg_tot->SetBinError(ibin+1,RB2DF->GetSumofYieldsError_NoIso()); 
    hmgg_gg->SetBinContent(ibin+1,RB2DF->GetNgamgamYield_NoIso()); 
    hmgg_gg->SetBinError(ibin+1,RB2DF->GetNgamgamYieldError_NoIso()); 
    hmgg_gj->SetBinContent(ibin+1,RB2DF->GetNgamjetYield_NoIso()); 
    hmgg_gj->SetBinError(ibin+1,RB2DF->GetNgamjetYieldError_NoIso()); 
    hmgg_jg->SetBinContent(ibin+1,RB2DF->GetNjetgamYield_NoIso()); 
    hmgg_jg->SetBinError(ibin+1,RB2DF->GetNjetgamYieldError_NoIso()); 
    hmgg_jj->SetBinContent(ibin+1,RB2DF->GetNjetjetYield_NoIso()); 
    hmgg_jj->SetBinError(ibin+1,RB2DF->GetNjetjetYieldError_NoIso()); 
    //---------------------------------------------------------------
    hmgg_tot_iso->SetBinContent(ibin+1,RB2DF->GetSumofYields()); 
    hmgg_tot_iso->SetBinError(ibin+1,RB2DF->GetSumofYieldsError()); 
    hmgg_gg_iso->SetBinContent(ibin+1,RB2DF->GetNgamgamYield()); 
    hmgg_gg_iso->SetBinError(ibin+1,RB2DF->GetNgamgamYieldError()); 
    hmgg_gj_iso->SetBinContent(ibin+1,RB2DF->GetNgamjetYield()); 
    hmgg_gj_iso->SetBinError(ibin+1,RB2DF->GetNgamjetYieldError()); 
    hmgg_jg_iso->SetBinContent(ibin+1,RB2DF->GetNjetgamYield()); 
    hmgg_jg_iso->SetBinError(ibin+1,RB2DF->GetNjetgamYieldError()); 
    hmgg_jj_iso->SetBinContent(ibin+1,RB2DF->GetNjetjetYield()); 
    hmgg_jj_iso->SetBinError(ibin+1,RB2DF->GetNjetjetYieldError()); 
    //---------------------------------------------------------------
  }
  delete RB2DF;
}
///////////////////////////////////////////////
void BkgShapeBinnedExtractor::FillHists_sysup()
///////////////////////////////////////////////
{
  int Nbins = bins.size();
  std::cout << "Number of bins = " << Nbins << std::endl;
  for(int ibin=0 ; ibin < Nbins ; ibin++){
    Bkg2DFit* RB2DF = new Bkg2DFit(m_tree_sysup,bins[ibin],m_ParamFile);
    bool do1Dfits      = true;
    bool do2Dfits      = true;
    bool doSplot       = false;
    RB2DF->Fitter(do1Dfits,do2Dfits,doSplot);
    RB2DF->PlotResults();
    hmgg_tot_sysup->SetBinContent(ibin+1,RB2DF->GetSumofYields_NoIso()); 
    hmgg_tot_sysup->SetBinError(ibin+1,RB2DF->GetSumofYieldsError_NoIso()); 
    hmgg_gg_sysup->SetBinContent(ibin+1,RB2DF->GetNgamgamYield_NoIso()); 
    hmgg_gg_sysup->SetBinError(ibin+1,RB2DF->GetNgamgamYieldError_NoIso()); 
    hmgg_gj_sysup->SetBinContent(ibin+1,RB2DF->GetNgamjetYield_NoIso()); 
    hmgg_gj_sysup->SetBinError(ibin+1,RB2DF->GetNgamjetYieldError_NoIso()); 
    hmgg_jg_sysup->SetBinContent(ibin+1,RB2DF->GetNjetgamYield_NoIso()); 
    hmgg_jg_sysup->SetBinError(ibin+1,RB2DF->GetNjetgamYieldError_NoIso()); 
    hmgg_jj_sysup->SetBinContent(ibin+1,RB2DF->GetNjetjetYield_NoIso()); 
    hmgg_jj_sysup->SetBinError(ibin+1,RB2DF->GetNjetjetYieldError_NoIso()); 
    //---------------------------------------------------------------
    hmgg_tot_iso_sysup->SetBinContent(ibin+1,RB2DF->GetSumofYields()); 
    hmgg_tot_iso_sysup->SetBinError(ibin+1,RB2DF->GetSumofYieldsError()); 
    hmgg_gg_iso_sysup->SetBinContent(ibin+1,RB2DF->GetNgamgamYield()); 
    hmgg_gg_iso_sysup->SetBinError(ibin+1,RB2DF->GetNgamgamYieldError()); 
    hmgg_gj_iso_sysup->SetBinContent(ibin+1,RB2DF->GetNgamjetYield()); 
    hmgg_gj_iso_sysup->SetBinError(ibin+1,RB2DF->GetNgamjetYieldError()); 
    hmgg_jg_iso_sysup->SetBinContent(ibin+1,RB2DF->GetNjetgamYield()); 
    hmgg_jg_iso_sysup->SetBinError(ibin+1,RB2DF->GetNjetgamYieldError()); 
    hmgg_jj_iso_sysup->SetBinContent(ibin+1,RB2DF->GetNjetjetYield()); 
    hmgg_jj_iso_sysup->SetBinError(ibin+1,RB2DF->GetNjetjetYieldError()); 
    //---------------------------------------------------------------
  }
}
/////////////////////////////////////////////////
void BkgShapeBinnedExtractor::FillHists_sysdown()
/////////////////////////////////////////////////
{
  int Nbins = bins.size();
  std::cout << "Number of bins = " << Nbins << std::endl;
  for(int ibin=0 ; ibin < Nbins ; ibin++){
    Bkg2DFit* RB2DF = new Bkg2DFit(m_tree_sysdown,bins[ibin],m_ParamFile);
    bool do1Dfits    = true;
    bool do2Dfits    = true;
    bool doSplot     = false;
    RB2DF->Fitter(do1Dfits,do2Dfits,doSplot);
    RB2DF->PlotResults();
    hmgg_tot_sysdown->SetBinContent(ibin+1,RB2DF->GetSumofYields_NoIso()); 
    hmgg_tot_sysdown->SetBinError(ibin+1,RB2DF->GetSumofYieldsError_NoIso()); 
    hmgg_gg_sysdown->SetBinContent(ibin+1,RB2DF->GetNgamgamYield_NoIso()); 
    hmgg_gg_sysdown->SetBinError(ibin+1,RB2DF->GetNgamgamYieldError_NoIso()); 
    hmgg_gj_sysdown->SetBinContent(ibin+1,RB2DF->GetNgamjetYield_NoIso()); 
    hmgg_gj_sysdown->SetBinError(ibin+1,RB2DF->GetNgamjetYieldError_NoIso()); 
    hmgg_jg_sysdown->SetBinContent(ibin+1,RB2DF->GetNjetgamYield_NoIso()); 
    hmgg_jg_sysdown->SetBinError(ibin+1,RB2DF->GetNjetgamYieldError_NoIso()); 
    hmgg_jj_sysdown->SetBinContent(ibin+1,RB2DF->GetNjetjetYield_NoIso()); 
    hmgg_jj_sysdown->SetBinError(ibin+1,RB2DF->GetNjetjetYieldError_NoIso()); 
    //---------------------------------------------------------------
    hmgg_tot_iso_sysdown->SetBinContent(ibin+1,RB2DF->GetSumofYields()); 
    hmgg_tot_iso_sysdown->SetBinError(ibin+1,RB2DF->GetSumofYieldsError()); 
    hmgg_gg_iso_sysdown->SetBinContent(ibin+1,RB2DF->GetNgamgamYield()); 
    hmgg_gg_iso_sysdown->SetBinError(ibin+1,RB2DF->GetNgamgamYieldError()); 
    hmgg_gj_iso_sysdown->SetBinContent(ibin+1,RB2DF->GetNgamjetYield()); 
    hmgg_gj_iso_sysdown->SetBinError(ibin+1,RB2DF->GetNgamjetYieldError()); 
    hmgg_jg_iso_sysdown->SetBinContent(ibin+1,RB2DF->GetNjetgamYield()); 
    hmgg_jg_iso_sysdown->SetBinError(ibin+1,RB2DF->GetNjetgamYieldError()); 
    hmgg_jj_iso_sysdown->SetBinContent(ibin+1,RB2DF->GetNjetjetYield()); 
    hmgg_jj_iso_sysdown->SetBinError(ibin+1,RB2DF->GetNjetjetYieldError()); 
    //---------------------------------------------------------------
  }
}
///////////////////////////////////////////////
void BkgShapeBinnedExtractor::StoretoRootFile()
///////////////////////////////////////////////
{
  int mggmin  = bins[0].first;
  int mggmax  = bins[bins.size()-1].second;
  TString ext = Form("_mgg%d%d",mggmin,mggmax); 
  TFile fout("rootfiles/BinnedBkg"+ext+".root","RECREATE");
  fout.Add( hmgg );
  fout.Add( hmgg_tot );
  fout.Add( hmgg_gg );
  fout.Add( hmgg_gj );
  fout.Add( hmgg_jg );
  fout.Add( hmgg_jj );
  fout.Add( hmgg_gjjg );
  if( m_doSys ){
    fout.Add( hmgg_tot_sysup );
    fout.Add( hmgg_gg_sysup );
    fout.Add( hmgg_gj_sysup );
    fout.Add( hmgg_jg_sysup );
    fout.Add( hmgg_jj_sysup );
    fout.Add( hmgg_gjjg_sysup );
    fout.Add( hmgg_tot_sysdown );
    fout.Add( hmgg_gg_sysdown );
    fout.Add( hmgg_gj_sysdown );
    fout.Add( hmgg_jg_sysdown );
    fout.Add( hmgg_jj_sysdown );
    fout.Add( hmgg_gjjg_sysdown );
  }
  fout.Add( hmgg_iso );
  fout.Add( hmgg_tot_iso );
  fout.Add( hmgg_gg_iso );
  fout.Add( hmgg_gj_iso );
  fout.Add( hmgg_jg_iso );
  fout.Add( hmgg_jj_iso );
  fout.Add( hmgg_gjjg_iso );
  if( m_doSys ){
    fout.Add( hmgg_tot_iso_sysup );
    fout.Add( hmgg_gg_iso_sysup );
    fout.Add( hmgg_gj_iso_sysup );
    fout.Add( hmgg_jg_iso_sysup );
    fout.Add( hmgg_jj_iso_sysup );
    fout.Add( hmgg_gjjg_iso_sysup );
    fout.Add( hmgg_tot_iso_sysdown );
    fout.Add( hmgg_gg_iso_sysdown );
    fout.Add( hmgg_gj_iso_sysdown );
    fout.Add( hmgg_jg_iso_sysdown );
    fout.Add( hmgg_jj_iso_sysdown );
    fout.Add( hmgg_gjjg_iso_sysdown );
  }
  fout.Write();
  fout.Close();

}
///////////////////////////////////////////////////
BkgShapeBinnedExtractor::~BkgShapeBinnedExtractor()
//////////////////////////////////////////////////
{
  //destructor
  delete m_tree;
  delete m_tree_sysup;
  delete m_tree_sysdown;
  delete m_file;
  delete m_file_sysup;
  delete m_file_sysdown;
  delete hmgg;
  delete hmgg_tot;
  delete hmgg_gg;
  delete hmgg_gj;
  delete hmgg_jg;
  delete hmgg_jj;
  delete hmgg_gjjg;
  delete hmgg_tot_sysup;
  delete hmgg_gg_sysup;
  delete hmgg_gj_sysup;
  delete hmgg_jg_sysup;
  delete hmgg_jj_sysup;
  delete hmgg_gjjg_sysup;
  delete hmgg_tot_sysdown;
  delete hmgg_gg_sysdown;
  delete hmgg_gj_sysdown;
  delete hmgg_jg_sysdown;
  delete hmgg_jj_sysdown;
  delete hmgg_gjjg_sysdown;
  delete hmgg_iso;
  delete hmgg_tot_iso;
  delete hmgg_gg_iso;
  delete hmgg_gj_iso;
  delete hmgg_jg_iso;
  delete hmgg_jj_iso;
  delete hmgg_gjjg_iso;
  delete hmgg_tot_iso_sysup;
  delete hmgg_gg_iso_sysup;
  delete hmgg_gj_iso_sysup;
  delete hmgg_jg_iso_sysup;
  delete hmgg_jj_iso_sysup;
  delete hmgg_gjjg_iso_sysup;
  delete hmgg_tot_iso_sysdown;
  delete hmgg_gg_iso_sysdown;
  delete hmgg_gj_iso_sysdown;
  delete hmgg_jg_iso_sysdown;
  delete hmgg_jj_iso_sysdown;
  delete hmgg_gjjg_iso_sysdown;
}
