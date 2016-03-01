#include "BkgEstimatorUpgrade.h"
#include <TFile.h>
#include <TTree.h>
#include <TError.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooMsgService.h>

#include "BkgUncertaintyBuilder.h"
#include "ToolsCommons.h"
#include "ToolsSignificanceHist.h"
#include "ToolsExtendedCanvas.h"
#include "ToolsTreeReader.h"
#include "ToolsUtilities.h"

//////////////////////////////////////////
BkgEstimatorUpgrade::BkgEstimatorUpgrade()
/////////////////////////////////////////
{
  //Default constructor
  SetVerbose(false);
  SetUncertFlag(true);
}
///////////////////////////////////////////
BkgEstimatorUpgrade::~BkgEstimatorUpgrade()
///////////////////////////////////////////
{
  delete frame_Purity;

}

//////////////////////////////////////////////////////////////////////////////////////
void BkgEstimatorUpgrade::SetYieldSystFiles(TString st_L2,TString st_L3,TString st_L5)
//////////////////////////////////////////////////////////////////////////////////////
{
  name_f_frac_L2 = st_L2;
  name_f_frac_L3 = st_L3;
  name_f_frac_L5 = st_L5;
}

///////////////////////////////////////////////////////////////////////
void BkgEstimatorUpgrade::SetRedSystFiles(TString st_L4,TString st_L5)
//////////////////////////////////////////////////////////////////////
{
  name_f_red_L4 = st_L4;
  name_f_red_L5 = st_L5;

}

////////////////////////////////
void BkgEstimatorUpgrade::Init()
////////////////////////////////
{
  //--> First you need to call the file setters
  if(m_VERBOSE) std::cout << "Read the files";
  //----------------------------------------------------------------
  TFile* f_data         = new TFile(name_f_data,"read");
  TFile* f_red          = new TFile(name_f_red,"read");
  TFile* f_irr          = new TFile(name_f_irr,"read");
  //----------------------------------------------------------------
  if(m_VERBOSE) std::cout << "... done !"<< std::endl;

  if(m_VERBOSE) std::cout << "Get the objects";
  //---------------------------------------------------------------
  m_hh_data       = (TH1D*)f_data->Get("mgg_data");
  m_gmgg_data     = AnalysisTools::DataGraph(*m_hh_data);
  m_gmgg_data->SetName("gmgg_data");
  //---------------------------------------------------------------
  m_hh_gj = (TH1D*)f_red->Get("h_subleading_fake");
  m_hh_jg = (TH1D*)f_red->Get("h_leading_fake");
  m_hh_jj = (TH1D*)f_red->Get("h_double_fake");
  m_hh_gj->SetName("bkg_gammajet_gg_full"); 
  m_hh_jg->SetName("bkg_jetgamma_gg_full");
  m_hh_jj->SetName("bkg_jetjet_gg_full"); 
  m_hh_jj->SetTitle("bkg_jetjet_gg_full");
  m_hh_gj->SetTitle("bkg_gammajet_gg_full");
  m_hh_jg->SetTitle("bkg_jetgamma_gg_full");
  //--------------------------------------------------------------------
  // m_hh_irr      = (TH1D*)m_f_irr->Get("irreducible_shape_nokfac");
  m_hh_irr      = (TH1D*)f_irr->Get("irreducible_shape");
  m_hh_gg       = (TH1D*)m_hh_irr->Clone("bkg_irreducible_gg_full");
  m_hh_gg->SetTitle(m_hh_gg->GetName());
  //--------------------------------------------------------------------
  //---------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------------
  if(m_VERBOSE) std::cout << "... done !"<< std::endl;


  //-----------------------------------------------------------------------------------------------------------
  dataNorm = m_hh_data->Integral(Commons::norm_bin.first,Commons::norm_bin.second);
  //---------------------------------------------------------
  if(m_VERBOSE) std::cout << "Set yields";
  std::pair<double,double> y_gg = GetYield("Ngg",name_f_frac_L4);
  std::pair<double,double> y_gj = GetYield("Ngj",name_f_frac_L4);
  std::pair<double,double> y_jg = GetYield("Njg",name_f_frac_L4);
  std::pair<double,double> y_jj = GetYield("Njj",name_f_frac_L4);
  Commons::SetYields(y_gg.first,y_gj.first,y_jg.first,y_jj.first,true);
  Commons::SetYieldsStatError(y_gg.second,y_gj.second,y_jg.second,y_jj.second);

  if(m_VERBOSE) std::cout << "Build total bkg \n";
  BuildTotalBkg();//Build the total background 
  if(m_VERBOSE) std::cout << "... done !"<< std::endl;
  if(m_VERBOSE) std::cout << "Remove stat errors bars on the red hist \n";
  RemoveReducibleStatErrors();//Set stat errors of the red hists to 0
  if(m_VERBOSE) std::cout << "... done !"<< std::endl;
  if(m_VERBOSE) std::cout << "Build total bkg uncert ... " << std::endl;
  if(m_doSys)
    BuildTotalBkgUncertainty();//Build the total bkg syst
  if(m_VERBOSE) std::cout << "Build Search Histograms: \n";
  BuildHistSearch();//Create the histograms for BAT
  if(m_VERBOSE) std::cout << "... done !"<< std::endl;
  if(m_VERBOSE) std::cout << "Build Background graphs \n";
  if(m_doSys)
    BuildBkgErrorsGraph();
  if(m_VERBOSE) std::cout << "... done !"<< std::endl;


}

/////////////////////////////////////////////////////
void BkgEstimatorUpgrade::RemoveReducibleStatErrors()
/////////////////////////////////////////////////////
{
  int Nbins = m_hh_red->GetNbinsX();
  for( int ibin=1; ibin<Nbins+1;ibin++){
    m_hh_red->SetBinError(ibin,0);
    m_hh_gj->SetBinError(ibin,0);
    m_hh_jg->SetBinError(ibin,0);
    m_hh_gj->SetBinError(ibin,0);
  }
}
///////////////////////////////////////////
void BkgEstimatorUpgrade::BuildHistSearch()
///////////////////////////////////////////
{
  std::cout << "data, " ;
  m_hh_data_search = AnalysisTools::GetTruncatedHist(*m_hh_data,Commons::nBins_search,Commons::binning_search);
  m_hh_data_search->SetName("hh_data");
  m_hh_data_search->SetTitle(m_hh_data_search->GetName());

  std::cout << "irreducible, " ;
  m_hh_irr_search = AnalysisTools::GetTruncatedHist(*m_hh_irr,Commons::nBins_search,Commons::binning_search);
  m_hh_irr_search->SetName("bkg_irreductible_gg");
  m_hh_irr_search->SetTitle(m_hh_irr_search->GetName());

  std::cout << "reducible, " ;
  m_hh_red_search = AnalysisTools::GetTruncatedHist(*m_hh_red,Commons::nBins_search,Commons::binning_search);
  m_hh_red_search->SetName("bkg_reductible_gg");
  m_hh_red_search->SetTitle(m_hh_red_search->GetName());

  std::cout << "total bkg, " ;
  m_hh_bkg_search = AnalysisTools::GetTruncatedHist(*m_hh_bkg,Commons::nBins_search,Commons::binning_search);
  m_hh_bkg_search->SetName("bkg_total_gg");
  m_hh_bkg_search->SetTitle(m_hh_bkg_search->GetName());

  if(m_doSys){
    std::cout << "total uncert, " ;
    m_hh_tot_syst_search = AnalysisTools::GetTruncatedHist(*m_hh_tot_syst,Commons::nBins_search,Commons::binning_search);
    m_hh_tot_syst_search->SetName("bkg_total_syst_gg");
    m_hh_tot_syst_search->SetTitle(m_hh_tot_syst_search->GetName());
    
    std::cout << "purity uncert, " ;
    m_hh_pur_syst_search = AnalysisTools::GetTruncatedHist(*m_hh_pur_syst,Commons::nBins_search,Commons::binning_search);
    m_hh_pur_syst_search->SetName("bkg_purity_syst_gg");
    m_hh_pur_syst_search->SetTitle(m_hh_pur_syst_search->GetName());
    
    std::cout << "reducible uncert, " ;
    m_hh_red_syst_search = AnalysisTools::GetTruncatedHist(*m_hh_red_syst,Commons::nBins_search,Commons::binning_search);
    m_hh_red_syst_search->SetName("bkg_reducible_syst_gg");
    m_hh_red_syst_search->SetTitle(m_hh_red_syst_search->GetName());
    
    std::cout << "irreducible uncert." ;
    m_hh_irr_syst_search = AnalysisTools::GetTruncatedHist(*m_hh_irr_syst,Commons::nBins_search,Commons::binning_search);
    m_hh_irr_syst_search->SetName("bkg_irreducible_syst_gg");
    m_hh_irr_syst_search->SetTitle(m_hh_irr_syst_search->GetName());
    std::cout << std::endl;
  }
}

/////////////////////////////////////////
void BkgEstimatorUpgrade::BuildTotalBkg()
////////////////////////////////////////
{

  //  m_hh_bkg= GetTotalBkgHist(m_hh_gg,m_hh_gj,m_hh_jg,m_hh_jj,Commons::yields);
  m_hh_bkg= GetTotalBkgHist(m_hh_gg,m_hh_gj,m_hh_jg,m_hh_jj,Commons::yields,"nomhists.root");
  m_hh_bkg->SetName("bkg_total_gg_full");
  m_hh_bkg->SetTitle("bkg_total_gg_full");

  std::vector<double> yields_red;
  yields_red.push_back(Commons::yields[1]);
  yields_red.push_back(Commons::yields[2]);
  yields_red.push_back(Commons::yields[3]);
  double sumyields = 0;
  for(int iy=0;iy<(int)Commons::yields.size();iy++) sumyields += Commons::yields[iy];
  double purity = Commons::yields[0]/sumyields;

  m_hh_red = GetRedBkgHist(m_hh_gj,m_hh_jg,m_hh_jj,yields_red);
  m_hh_red->SetName("bkg_reducible_gg_full");
  m_hh_red->SetTitle("bkg_reducible_gg_full");
  m_hh_red->Scale(1-purity);

}

////////////////////////////////////////////////////
void BkgEstimatorUpgrade::BuildTotalBkgUncertainty()
////////////////////////////////////////////////////
{
  //--------------------------------------------------------------
  // PrintLoosePrimeResults();
  if(m_VERBOSE) std::cout << "Compute bkg uncertainty components";
  canv_YieldsSyst = YieldsSystUncert();//Fill m_hh_fracsyst_syst
  canv_YieldsStat = YieldsStatUncert();//Fill m_hh_fracstat_syst
  canv_RedShape = RedShapeUncert();//Fill m_hh_red_syst
  canv_IrrShape = IrrShapeUncert();// Fill m_hh_irr_syst_***
  if(m_VERBOSE) std::cout << "... done !" << std::endl;
  //--------------------------------------------------------------

  //--------------------------------------------------------------
  if(m_VERBOSE) std::cout << "Compute bkg uncertainty ";
  m_hh_pur_syst = AnalysisTools::HistAddedUncert(*m_hh_fracsyst_syst,*m_hh_fracstat_syst);
  m_hh_pur_syst->SetName("bkg_purity_syst_gg_full");
  m_hh_pur_syst->SetTitle(m_hh_pur_syst->GetName());
  m_hh_tot_syst = AnalysisTools::HistAddedUncert(*m_hh_irr_syst,*m_hh_pur_syst);
  m_hh_tot_syst = AnalysisTools::HistAddedUncert(*m_hh_tot_syst,*m_hh_red_syst);
  m_hh_tot_syst->SetName("bkg_total_syst_gg_full");
  m_hh_tot_syst->SetTitle(m_hh_tot_syst->GetName());
  if(m_VERBOSE) std::cout << "... done !" << std::endl;
  //--------------------------------------------------------------


  //--------------------------------------------------------------
  if(m_VERBOSE) std::cout << "Compute bkg histogram with systematics included ";
  //---- Build total background with systematics included for BumpHunter
  m_hh_bkg_withsyst = new TH1D( "bkg_total_gg_full_withsyst","bkg_total_gg_full_withsyst",
				Commons::nBins,Commons::binning );
  m_hh_red_withsyst = new TH1D( "bkg_reducible_gg_full_withsyst","bkg_reducible_gg_full_withsyst",
				Commons::nBins,Commons::binning );
  for( int ibin=1; ibin <= Commons::nBins ;ibin++){
    //--> Total Background
    double val_tot    = m_hh_bkg->GetBinContent(ibin);
    double e_syst_tot = m_hh_tot_syst->GetBinContent(ibin)*val_tot;
    double e_stat_tot = m_hh_bkg->GetBinError(ibin);
    double e_sum_tot  = sqrt(e_syst_tot*e_syst_tot+e_stat_tot*e_stat_tot);
    m_hh_bkg_withsyst->SetBinContent(ibin,val_tot);
    m_hh_bkg_withsyst->SetBinError(ibin,e_sum_tot);
    //--> Reducible Background
    double val_red    = m_hh_red->GetBinContent(ibin);
    double e_syst_red = m_hh_red_syst->GetBinContent(ibin)*val_tot;
    double e_stat_red = m_hh_red->GetBinError(ibin);
    double e_sum_red  = sqrt(e_stat_red*e_stat_red+e_syst_red*e_syst_red);
    m_hh_red_withsyst->SetBinContent(ibin,val_red);
    m_hh_red_withsyst->SetBinError(ibin,e_sum_red);
  }
  if(m_VERBOSE) std::cout << "... done !" << std::endl;
  //--------------------------------------------------------------

}
///////////////////////////////////////////////
void BkgEstimatorUpgrade::BuildBkgErrorsGraph()
//////////////////////////////////////////////
{
  m_gmgg_bkg = AnalysisTools::BkgGraph(*m_hh_bkg_withsyst);
  m_gmgg_red = AnalysisTools::BkgGraph(*m_hh_red_withsyst);
  m_gmgg_bkg->SetName("gmgg_bkg");
  m_gmgg_red->SetName("gmgg_red");
  m_gmgg_bkg->SetTitle("gmgg_bkg");
  m_gmgg_red->SetTitle("gmgg_red");
}



/////////////////////////////////////////////////
void BkgEstimatorUpgrade::GetFinalYieldsPerBin()
////////////////////////////////////////////////
{
  int Nbins = m_hh_data->GetNbinsX();
  std::cout << "Bin"     << " | "
	    << "low edge"<< " | " 
	    << "Bakground expectation" << " | " 
	    << "Data" << std::endl;

  for( int ibin=1; ibin< Nbins ; ibin++){
    std::cout << ibin  << " | "
	      << m_hh_data->GetBinLowEdge(ibin)     << " | " 
	      << m_hh_bkg_withsyst->GetBinContent(ibin) << "+/-"
	      << m_hh_bkg_withsyst->GetBinError(ibin)   << " | "
	      << m_hh_data->GetBinContent(ibin)     << std::endl;
  }

  int binfirst [] = {20,48,54,58,63,66,69,72,75,77,79,81,83};
  int binlast  [] = {48,54,58,63,66,69,72,75,77,79,81,83,100};
  int Nbins_Int   = 13;
  for(int ibin = 0; ibin <Nbins_Int;ibin++){
    double irrInt  = m_hh_irr->Integral(binfirst[ibin],binlast[ibin]-1);
    double redInt  = m_hh_red->Integral(binfirst[ibin],binlast[ibin]-1);
    double bkgInt  = m_hh_bkg->Integral(binfirst[ibin],binlast[ibin]-1);
    double dataInt = m_hh_data->Integral(binfirst[ibin],binlast[ibin]-1);
    double bkgerrInt = 0;
    for( int i =binfirst[ibin];i<binlast[ibin]; i++){
      double BC = m_hh_bkg->GetBinContent(i);
      double rel_err = m_hh_tot_syst->GetBinContent(i);
      bkgerrInt += BC*rel_err;
    }

    std::cout << "$[" 
	      << (int)m_hh_data->GetBinLowEdge(binfirst[ibin]) <<","
	      << (int)m_hh_data->GetBinLowEdge(binlast[ibin]) << "]" 
	      << "$ & $" << irrInt << "pm" << "" //need to compute the errors
	      << "$ & $" << redInt << "pm" << "" //need to compute the errors
	      << "$ & $" << bkgInt << "pm" << bkgerrInt
	      << "$ & $" << dataInt<< "$ " 
	      << std::endl;
  }
}
/////////////////////////////////////////////////////////
void BkgEstimatorUpgrade::GetFinalYieldsPerBin_NoUncert()
/////////////////////////////////////////////////////////
{

  std::cout << " Bin |\t reducible |\t irreducible |\t total |\t data |\t data-total" << std::endl;
  for( int ibin = 1; ibin<=m_hh_bkg->GetNbinsX();ibin++){
    std::cout << Form("[%1.1f,%1.1f]",m_hh_bkg->GetBinLowEdge(ibin),m_hh_bkg->GetBinLowEdge(ibin+1))
	      << " |\t "
	      << m_hh_red->GetBinContent(ibin) 
	      << "|\t "
	      << m_hh_gg->GetBinContent(ibin) 
	      << "|\t "
	      << m_hh_bkg->GetBinContent(ibin) 
	      << "|\t "
	      << m_hh_data->GetBinContent(ibin) 
	      << "|\t "
	      << m_hh_data->GetBinContent(ibin)-m_hh_bkg->GetBinContent(ibin)  
	      << std::endl;
  }


}

////////////////////////////////////////////////
TCanvas* BkgEstimatorUpgrade::YieldsSystUncert()
////////////////////////////////////////////////
{
  
  // double Yield_L2[] = {9646,3086,1222,1707};
  // double Yield_L3[] = {10120,2986,1207,1379};
  // double Yield_L4[] = {10812,2640,1057,1167};
  // double Yield_L5[] = {11647,2265,886,890};

  //-------> Construction of the yields uncertainty from a relative uncertainty on the purity
  double alpha = Commons::purity_relative_uncert;//--> Arbitrary 10 %
  double beta = alpha*Commons::yields[0]/(Commons::yields[1]+Commons::yields[2]+Commons::yields[3]);
  std::vector<double> Yields_up;
  std::vector<double> Yields_do;
  Yields_up.push_back(Commons::yields[0]*(1+alpha));
  Yields_up.push_back(Commons::yields[1]*(1-beta));
  Yields_up.push_back(Commons::yields[2]*(1-beta));
  Yields_up.push_back(Commons::yields[3]*(1-beta));
  Yields_do.push_back(Commons::yields[0]*(1-alpha));
  Yields_do.push_back(Commons::yields[1]*(1+beta));
  Yields_do.push_back(Commons::yields[2]*(1+beta));
  Yields_do.push_back(Commons::yields[3]*(1+beta));

  TH1D* hh_bkg_up = GetTotalBkgHist(m_hh_gg,m_hh_gj,m_hh_jg,m_hh_jj,Yields_up);
  TH1D* hh_bkg_do = GetTotalBkgHist(m_hh_gg,m_hh_gj,m_hh_jg,m_hh_jj,Yields_do); 

  BkgUncertaintyBuilder *BUB;
  BUB = new BkgUncertaintyBuilder(*m_hh_bkg,*hh_bkg_up);
  TH1D* huncert_up = BUB->GetRelUncertHist();
  TH1D* huncert_frac = BUB->GetBkgFracHist();
  BUB = new BkgUncertaintyBuilder(*m_hh_bkg,*hh_bkg_do);
  TH1D* huncert_do = BUB->GetRelUncertHist();


  TH1D* huncert_envelop = (TH1D*)huncert_frac->Clone("uncert_envelop");
  for(int ibin=0;ibin<huncert_envelop->GetNbinsX()+1;ibin++){
    if( ibin < Commons::norm_bin.second) continue;
    if( huncert_envelop->GetBinContent(ibin)<huncert_envelop->GetBinContent(ibin-1))
      huncert_envelop->SetBinContent(ibin,huncert_envelop->GetBinContent(ibin-1));
  }
  huncert_envelop->SetLineStyle(kDashed);

  TCanvas* c = new TCanvas("cyieldsyst","Yields systematic uncertainties",800,600);
  c->cd();
  huncert_up->GetXaxis()->SetRange(Commons::norm_bin.first,Commons::binning[Commons::nBins]);
  huncert_up->GetYaxis()->SetRangeUser(-0.2,0.2);
  huncert_up->Draw("HIST");
  huncert_do->Draw("SAMEHIST");
  huncert_envelop->Draw("SAMEHIST");

  m_hh_fracsyst_syst = huncert_frac;//huncert_envelop;

  return c;
}

////////////////////////////////////////////////
TCanvas* BkgEstimatorUpgrade::YieldsStatUncert()
///////////////////////////////////////////////
{

  //---->  Determine the 1sigma band on the purity
  std::pair<double,double> PurityInterval = GetPurityInterval(name_f_frac_L4);


  //----> Iterate over the yield tree to get the 
  //----> set of yields giving the closest value of 
  //----> the -1sigma and +1sigma purity values 
  TFile* file0 = TFile::Open(name_f_frac_L4,"READ");
  TTree* yieldstree = (TTree*)file0->Get("yieldstree");
 
  int i_do = -999;
  int i_up = -999;
  double dist_do = 1e18;
  double dist_up = 1e18;

  TreeReader Rd(yieldstree);
  int nentries = yieldstree->GetEntriesFast();
  for( int ipe = 0 ; ipe<nentries;ipe++){
    Rd.GetEntry(ipe);
    double purity = Rd.GetVariable("Ngg/(Ngg+Ngj+Njg+Njj)");
    if( fabs(purity-PurityInterval.first) < dist_do){
      dist_do = fabs(purity-PurityInterval.first);
      i_do = ipe;
    }
    if( fabs(purity-PurityInterval.second) < dist_up){
      dist_up = fabs(purity-PurityInterval.second);
      i_up = ipe;
    }
  }

  Rd.GetEntry(i_do);
  std::vector<double> Yields_do;
  Yields_do.push_back(Rd.GetVariable("Ngg"));
  Yields_do.push_back(Rd.GetVariable("Ngj"));
  Yields_do.push_back(Rd.GetVariable("Njg"));
  Yields_do.push_back(Rd.GetVariable("Njj"));
  Rd.GetEntry(i_up);
  std::vector<double> Yields_up;
  Yields_up.push_back(Rd.GetVariable("Ngg"));
  Yields_up.push_back(Rd.GetVariable("Ngj"));
  Yields_up.push_back(Rd.GetVariable("Njg"));
  Yields_up.push_back(Rd.GetVariable("Njj"));

  TH1D* hh_bkg_up = GetTotalBkgHist(m_hh_gg,m_hh_gj,m_hh_jg,m_hh_jj,Yields_up);
  TH1D* hh_bkg_do = GetTotalBkgHist(m_hh_gg,m_hh_gj,m_hh_jg,m_hh_jj,Yields_do); 

  BkgUncertaintyBuilder *BUB;
  BUB = new BkgUncertaintyBuilder(*m_hh_bkg,*hh_bkg_up);
  TH1D* huncert_up = BUB->GetRelUncertHist();
  m_hh_fracstat_syst = BUB->GetBkgFracHist();
  BUB = new BkgUncertaintyBuilder(*m_hh_bkg,*hh_bkg_do);
  TH1D* huncert_do = BUB->GetRelUncertHist();

  TCanvas* c = new TCanvas("cyieldsstat","Yields statistical uncertainties",800,600);
  c->cd();
  huncert_up->GetXaxis()->SetRange(Commons::norm_bin.first,Commons::binning[Commons::nBins]);
  huncert_up->GetYaxis()->SetRangeUser(-0.01,0.01);
  huncert_up->Draw("HIST");
  huncert_do->Draw("SAMEHIST");
  return c;


}
//////////////////////////////////////////////
TCanvas* BkgEstimatorUpgrade::RedShapeUncert()
/////////////////////////////////////////////
{
  m_hh_gj_syst = GetRedShapeSyst("gj");
  m_hh_jg_syst = GetRedShapeSyst("jg");
  m_hh_jj_syst = GetRedShapeSyst("jj");
  //----------------------------------------------
  m_hh_red_syst = (TH1D*)m_hh_gj_syst->Clone("bkg_reducible_syst_gg_full");

  m_hh_red_syst = AnalysisTools::HistAddedUncert(*m_hh_gj_syst,*m_hh_jg_syst);
  m_hh_red_syst = AnalysisTools::HistAddedUncert(*m_hh_red_syst,*m_hh_jj_syst);
  m_hh_red_syst->SetName("bkg_reducible_syst_gg_full");
  m_hh_red_syst->SetTitle(m_hh_red_syst->GetName());


  TCanvas* c = new TCanvas("credshape","Red Shape Uncert",800,600);
  c->cd();
  m_hh_red_syst->GetXaxis()->SetRange(Commons::norm_bin.first,Commons::binning[Commons::nBins]);
  m_hh_red_syst->Draw("HIST");
  m_hh_gj_syst->Draw("sameHIST");
  m_hh_jg_syst->Draw("sameHIST");
  m_hh_jj_syst->Draw("sameHIST");
  c->BuildLegend();
  return c;
}
//////////////////////////////////////////////
TCanvas* BkgEstimatorUpgrade::IrrShapeUncert()
//////////////////////////////////////////////
{
  m_hh_scale_syst = GetScaleShapeSyst();
  m_hh_iso_syst   = GetIsoShapeSyst();
  m_hh_isodatamc_syst   = GetIsoDataMCShapeSyst();
  m_hh_pdf_syst   = GetPDFShapeSyst();
  //----------------------------------------------




  m_hh_irr_syst = AnalysisTools::HistAddedUncert(*m_hh_scale_syst,*m_hh_iso_syst);
  m_hh_irr_syst = AnalysisTools::HistAddedUncert(*m_hh_irr_syst,*m_hh_pdf_syst);
  m_hh_irr_syst = AnalysisTools::HistAddedUncert(*m_hh_irr_syst,*m_hh_isodatamc_syst);
  m_hh_irr_syst->SetName("bkg_irreducible_syst_gg_full"); 
  m_hh_irr_syst->SetTitle(m_hh_irr_syst->GetName());

  TCanvas* c = new TCanvas("cirrshape","Irreducible shape uncert");
  c->cd();
  m_hh_irr_syst->GetXaxis()->SetRange(Commons::norm_bin.first,Commons::binning[Commons::nBins]);
  m_hh_irr_syst->Draw("HIST");
  m_hh_scale_syst->Draw("SAMEHIST");
  m_hh_iso_syst->Draw("sameHIST");
  m_hh_isodatamc_syst->Draw("sameHIST");
  m_hh_pdf_syst->Draw("sameHIST");
  c->BuildLegend();
  return c;
}
///////////////////////////////////////////////////////////
TH1D* BkgEstimatorUpgrade::GetRedShapeSyst(TString bkgtype)
///////////////////////////////////////////////////////////
{
  TFile *f4 = new TFile(name_f_red_L4,"read");
  TFile *f5 = new TFile(name_f_red_L5,"read");
  TString histname;
  if(bkgtype=="gj") histname = "h_subleading_fake";
  else if(bkgtype=="jg") histname = "h_leading_fake";
  else if(bkgtype=="jj") histname = "h_double_fake";
  else std::cout << "WRONG BKG TYPE !!!!" << std::endl;


  TH1D* h4  = (TH1D*)f4->Get(histname);
  TH1D* h5 = (TH1D*)f5->Get(histname);

  TH1D* hsys;
  if(bkgtype=="gj") hsys = GetTotalBkgHist(m_hh_gg,h4,m_hh_jg,m_hh_jj,Commons::yields);
  else if(bkgtype=="jg") hsys = GetTotalBkgHist(m_hh_gg,m_hh_gj,h4,m_hh_jj,Commons::yields);
  else if(bkgtype=="jj") hsys = GetTotalBkgHist(m_hh_gg,m_hh_gj,m_hh_jg,h5,Commons::yields);

  BkgUncertaintyBuilder *BUB;
  BUB = new BkgUncertaintyBuilder(*m_hh_bkg,*hsys);

  TH1D* htemp = BUB->GetBkgFracHist();
  htemp->SetName("bkg_"+bkgtype+"_syst_gg_full");
  htemp->SetTitle("bkg_"+bkgtype+"_syst_gg_full");
  return htemp;

}

/////////////////////////////////////////////
TH1D* BkgEstimatorUpgrade::GetScaleShapeSyst()
/////////////////////////////////////////////
{
  TFile* f_irr  = new TFile(name_f_irr,"read");
  TH1D* h_scale = (TH1D*)f_irr->Get("irreducible_shape_scale");

  //  TH1D* hbkg_scale = GetTotalBkgHist(h_scale,m_hh_gj,m_hh_jg,m_hh_jj,Commons::yields);
  TH1D* hbkg_scale = GetTotalBkgHist(h_scale,m_hh_gj,m_hh_jg,m_hh_jj,Commons::yields,"scalehists.root");


  BkgUncertaintyBuilder *BUB;
  BUB = new BkgUncertaintyBuilder(*m_hh_bkg,*hbkg_scale);
  TH1D* htemp = BUB->GetBkgFracHist();
  htemp->SetName("bkg_scale_syst_gg_full");
  htemp->SetTitle("bkg_scale_syst_gg_full");
  return htemp;
}
////////////////////////////////////////////
TH1D* BkgEstimatorUpgrade::GetPDFShapeSyst()
////////////////////////////////////////////
{
  TFile* f_irr   = new TFile(name_f_irr,"read");
  TH1D* h_pdfset = (TH1D*)f_irr->Get("irreducible_shape_pdfset");
  TH1D* h_pdfeig = (TH1D*)f_irr->Get("irreducible_shape_pdfeig");


  std::cout<<"+++++++++++++++++++ PDFSHAPE ++++++++++++++++++++++"<<std::endl;

//  TH1D* hbkg_pdfset = GetTotalBkgHist(h_pdfset,m_hh_gj,m_hh_jg,m_hh_jj,Commons::yields);
//  TH1D* hbkg_pdfeig = GetTotalBkgHist(h_pdfeig,m_hh_gj,m_hh_jg,m_hh_jj,Commons::yields);
  TH1D* hbkg_pdfset = GetTotalBkgHist(h_pdfset,m_hh_gj,m_hh_jg,m_hh_jj,Commons::yields,"pdfsethists.root");
  TH1D* hbkg_pdfeig = GetTotalBkgHist(h_pdfeig,m_hh_gj,m_hh_jg,m_hh_jj,Commons::yields,"pdfeighists.root");


  BkgUncertaintyBuilder *BUB;
  BUB = new BkgUncertaintyBuilder(*m_hh_bkg,*hbkg_pdfset);
  TH1D* htemp_pdfset = BUB->GetBkgFracHist();
  BUB = new BkgUncertaintyBuilder(*m_hh_bkg,*hbkg_pdfeig);
  TH1D* htemp_pdfeig = BUB->GetBkgFracHist();

  // TFile * fpdf = new TFile("filepdf"+name_f_irr,"RECREATE");
  // m_hh_gg->Write();
  // m_hh_bkg->Write();
  // //  hbkg_pdfset->Write();
  // hbkg_pdfeig->Write();
  // h_pdfeig->Write();

  // //  htemp_pdfset->Write();
  // htemp_pdfeig->Write();
  // fpdf->Close();
  


  TH1D* htemp = AnalysisTools::HistAddedUncert(*htemp_pdfset,*htemp_pdfeig);
  htemp->SetName("bkg_pdf_syst_gg_full") ; 
  htemp->SetTitle("bkg_pdf_syst_gg_full") ; 
  return htemp;
}
////////////////////////////////////////////
TH1D* BkgEstimatorUpgrade::GetIsoShapeSyst()
////////////////////////////////////////////
{
  TFile* f_irr = new TFile(name_f_irr,"read");
  TH1D* h_iso  = (TH1D*)f_irr->Get("irreducible_shape_iso");

  //  TH1D* hbkg_iso = GetTotalBkgHist(h_iso,m_hh_gj,m_hh_jg,m_hh_jj,Commons::yields);
  TH1D* hbkg_iso = GetTotalBkgHist(h_iso,m_hh_gj,m_hh_jg,m_hh_jj,Commons::yields,"isohists.root");

  BkgUncertaintyBuilder *BUB = new BkgUncertaintyBuilder(*m_hh_bkg,*hbkg_iso);
  TH1D* htemp = BUB->GetBkgFracHist();
  htemp->SetName("bkg_iso_syst_gg_full");
  htemp->SetTitle("bkg_iso_syst_gg_full");
  return htemp;
}


////////////////////////////////////////////
TH1D* BkgEstimatorUpgrade::GetIsoDataMCShapeSyst()
////////////////////////////////////////////
{
  TFile* f_irr = new TFile(name_f_irr,"read");
  TH1D* h_iso  = (TH1D*)f_irr->Get("irreducible_shape_isodatamc");

  //  TH1D* hbkg_iso = GetTotalBkgHist(h_iso,m_hh_gj,m_hh_jg,m_hh_jj,Commons::yields);
  TH1D* hbkg_iso = GetTotalBkgHist(h_iso,m_hh_gj,m_hh_jg,m_hh_jj,Commons::yields,"isodatamchists.root");

  BkgUncertaintyBuilder *BUB = new BkgUncertaintyBuilder(*m_hh_bkg,*hbkg_iso);
  TH1D* htemp = BUB->GetBkgFracHist();
  htemp->SetName("bkg_isodatamc_syst_gg_full");
  htemp->SetTitle("bkg_isodatamc_syst_gg_full");
  return htemp;
}


////////////////////////////////////////////////////
void BkgEstimatorUpgrade::CreateBATInput(TString st)
///////////////////////////////////////////////////
{
  m_hh_bkg_search->SetName("bkg_total_gg");
  m_hh_bkg_search->SetTitle("bkg_total_gg");

  TObjArray Hlist(0);
  Hlist.Add(m_hh_data_search);
  Hlist.Add(m_hh_bkg_search);
  Hlist.Add(m_hh_red_search);
  Hlist.Add(m_hh_irr_search);
  if(m_doSys){
    Hlist.Add(m_hh_tot_syst_search);
    Hlist.Add(m_hh_red_syst_search);
    Hlist.Add(m_hh_irr_syst_search);
    Hlist.Add(m_hh_pur_syst_search);
  }
  TFile fout(st,"RECREATE");
  Hlist.Write();
  fout.Close();
}
///////////////////////////////////////////////////////////
void BkgEstimatorUpgrade::CreateBumpHunterInput(TString st)
///////////////////////////////////////////////////////////
{
  TObjArray Hlist1(0);

  Hlist1.Add(m_hh_data);
  Hlist1.Add(m_hh_bkg);
  Hlist1.Add(m_hh_red);
  Hlist1.Add(m_hh_gg);
  Hlist1.Add(m_hh_gj);
  Hlist1.Add(m_hh_jg);
  Hlist1.Add(m_hh_jj);
  if(m_doSys){
    Hlist1.Add(m_hh_tot_syst);
    Hlist1.Add(m_hh_irr_syst);
    Hlist1.Add(m_hh_red_syst);
    Hlist1.Add(m_hh_pur_syst);
    Hlist1.Add(m_hh_scale_syst);
    Hlist1.Add(m_hh_iso_syst);
    Hlist1.Add(m_hh_isodatamc_syst);
    Hlist1.Add(m_hh_pdf_syst);
    Hlist1.Add(m_hh_gj_syst);
    Hlist1.Add(m_hh_jg_syst);
    Hlist1.Add(m_hh_jj_syst);
    Hlist1.Add(m_hh_fracsyst_syst);
    Hlist1.Add(m_hh_fracstat_syst);
    Hlist1.Add(m_hh_bkg_withsyst);
    Hlist1.Add(m_hh_red_withsyst);
  }
  Hlist1.Add(m_gmgg_data);
  if(m_doSys){
    Hlist1.Add(m_gmgg_bkg);
    Hlist1.Add(m_gmgg_red);
  }
  if(m_doSys){
    Hlist1.Add(canv_YieldsSyst);
    Hlist1.Add(canv_YieldsStat);
    Hlist1.Add(canv_RedShape);  
    Hlist1.Add(canv_IrrShape);    
  }
  Hlist1.Add(frame_Purity);
  TFile fout1(st,"RECREATE");
  Hlist1.Write();
  fout1.Close();
  


}
/////////////////////////////////////////////////////////////////////////////
TH1D* BkgEstimatorUpgrade::GetTotalBkgHist( TH1D* hgg,TH1D* hgj,TH1D* hjg,TH1D* hjj,
					    const std::vector<double> yields)
//////////////////////////////////////////////////////////////////////
{
  if( yields.size() !=4 )
    Fatal( "BkgEstimatorUpgrade::GetTotalBkgHist",
	   "Invalid size of the yields vector !!" );

  double sumyields = 0; std::vector<double> frac;
  for(int iy=0;iy<(int)yields.size();iy++) sumyields += yields[iy];
  for(int iy=0;iy<(int)yields.size();iy++) frac.push_back(yields[iy]/sumyields);


  std::cout<<"=========== BkgEstimatorUpgrade::GetTotalBkgHist for "<<hgg->GetName()<<" ============="<<std::endl;

  std::cout<<"yields[0] "<<yields[0]<<" yields[1] "<<yields[1]<<" yields[2] "<<yields[2]<<" yields[3] "<<yields[3]<<std::endl;
  std::cout<<"sumyields "<<sumyields <<std::endl;
  std::cout<<"frac[0] "<<frac[0] <<" frac[1] "<<frac[1] <<" frac[2] "<<frac[2] <<" frac[3] "<<frac[3] <<std::endl;

  std::cout<<"norm bin 1 "<<Commons::norm_bin.first<<" normbin 2 " <<Commons::norm_bin.second<<std::endl;
  std::cout<<"gg integral in norm region = "<<hgg->Integral(Commons::norm_bin.first,Commons::norm_bin.second)<<std::endl;
  std::cout<<"gg histogram scaled by "<<frac[0]*dataNorm/hgg->Integral(Commons::norm_bin.first,Commons::norm_bin.second)<<std::endl;



  hgg->Scale(frac[0]*dataNorm/hgg->Integral(Commons::norm_bin.first,Commons::norm_bin.second));
  hgj->Scale(frac[1]*dataNorm/hgj->Integral(Commons::norm_bin.first,Commons::norm_bin.second));
  hjg->Scale(frac[2]*dataNorm/hjg->Integral(Commons::norm_bin.first,Commons::norm_bin.second));
  hjj->Scale(frac[3]*dataNorm/hjj->Integral(Commons::norm_bin.first,Commons::norm_bin.second));

  TH1D* htemp = (TH1D*)hgg->Clone("hh_totbkg");
  htemp->Add(hgj);htemp->Add(hjg);htemp->Add(hjj);
  return htemp;
}




/////////////////////////////////////////////////////////////////////////////
TH1D* BkgEstimatorUpgrade::GetTotalBkgHist( TH1D* hgg,TH1D* hgj,TH1D* hjg,TH1D* hjj,
					    const std::vector<double> yields, TString outname)
//////////////////////////////////////////////////////////////////////
{
  if( yields.size() !=4 )
    Fatal( "BkgEstimatorUpgrade::GetTotalBkgHist",
	   "Invalid size of the yields vector !!" );

  double sumyields = 0; std::vector<double> frac;
  for(int iy=0;iy<(int)yields.size();iy++) sumyields += yields[iy];
  for(int iy=0;iy<(int)yields.size();iy++) frac.push_back(yields[iy]/sumyields);


  std::cout<<"=========== BkgEstimatorUpgrade::GetTotalBkgHist for "<<hgg->GetName()<<" ============="<<std::endl;
//
//  std::cout<<"yields[0] "<<yields[0]<<" yields[1] "<<yields[1]<<" yields[2] "<<yields[2]<<" yields[3] "<<yields[3]<<std::endl;
//  std::cout<<"sumyields "<<sumyields <<std::endl;
//  std::cout<<"frac[0] "<<frac[0] <<" frac[1] "<<frac[1] <<" frac[2] "<<frac[2] <<" frac[3] "<<frac[3] <<std::endl;
//
//  std::cout<<"norm bin 1 "<<Commons::norm_bin.first<<" normbin 2 " <<Commons::norm_bin.second<<std::endl;
//  std::cout<<"gg integral in norm region = "<<hgg->Integral(Commons::norm_bin.first,Commons::norm_bin.second)<<std::endl;
//  std::cout<<"gg histogram scaled by "<<frac[0]*dataNorm/hgg->Integral(Commons::norm_bin.first,Commons::norm_bin.second)<<std::endl;



  hgg->Scale(frac[0]*dataNorm/hgg->Integral(Commons::norm_bin.first,Commons::norm_bin.second));
  hgj->Scale(frac[1]*dataNorm/hgj->Integral(Commons::norm_bin.first,Commons::norm_bin.second));
  hjg->Scale(frac[2]*dataNorm/hjg->Integral(Commons::norm_bin.first,Commons::norm_bin.second));
  hjj->Scale(frac[3]*dataNorm/hjj->Integral(Commons::norm_bin.first,Commons::norm_bin.second));


  TH1D* htemp = (TH1D*)hgg->Clone("hh_totbkg");
  htemp->Add(hgj);htemp->Add(hjg);htemp->Add(hjj);
//  TFile * f = new TFile(outname,"RECREATE");
//  hgg->Write();
//  hgj->Write();
//  hjg->Write();
//  hjj->Write();
//  f->Close();


  return htemp;
}






///////////////////////////////////////////////////////////////////////
TH1D* BkgEstimatorUpgrade::GetRedBkgHist( TH1D* hgj,TH1D* hjg,TH1D* hjj,
					  const std::vector<double> yields)
//////////////////////////////////////////////////////////////////////
{
  if( yields.size() !=3 )
    Fatal( "BkgEstimatorUpgrade::GetRedBkgHist",
	   "Invalid size of the yields vector !!");

  double sumyields = 0; std::vector<double> frac;
  for(int iy=0;iy<(int)yields.size();iy++) sumyields += yields[iy];
  for(int iy=0;iy<(int)yields.size();iy++) frac.push_back(yields[iy]/sumyields);

  hgj->Scale(frac[0]*dataNorm/hgj->Integral(Commons::norm_bin.first,Commons::norm_bin.second));
  hjg->Scale(frac[1]*dataNorm/hjg->Integral(Commons::norm_bin.first,Commons::norm_bin.second));
  hjj->Scale(frac[2]*dataNorm/hjj->Integral(Commons::norm_bin.first,Commons::norm_bin.second));

  TH1D* htemp = (TH1D*)hgj->Clone("hh_redbkg");
  htemp->Add(hjg);htemp->Add(hjj);
  return htemp;
}
//////////////////////////////////////////////////////////////////////////////////////////////////
std::pair<double,double> BkgEstimatorUpgrade::GetYield(TString type,TString treefile,bool doplot)
//////////////////////////////////////////////////////////////////////////////////////////////////
{
  TFile* file0 = TFile::Open(treefile,"READ");
  TTree* yieldstree = (TTree*)file0->Get("yieldstree");
  double Xmin = 1e18;
  double Xmax = -1e18;
  int nentries = (int)yieldstree->GetEntriesFast();
  TreeReader Rd(yieldstree);
  for( int ipe = 0 ; ipe<nentries;ipe++){
    Rd.GetEntry(ipe);
    double N = Rd.GetVariable(type);
    if( N<Xmin) Xmin=N;
    if( N>Xmax) Xmax=N;
  }
  RooRealVar N(type,type,Xmin,Xmax);
  double half= (Xmax-Xmin)/2.;
  RooDataSet data("data","data",RooArgSet(N),RooFit::Import(*yieldstree) );
  RooRealVar mean("mean","mean",Xmin+half,Xmin,Xmax) ;
  RooRealVar width("width","width",half/4.,1,half) ;
  RooGaussian gauss("gauss","gaussian PDF",N,mean,width) ;
  RooMsgService::instance().setSilentMode(true);
  gauss.fitTo(data);

  if(doplot){
    RooPlot *frame = N.frame();
    data.plotOn(frame,RooFit::Name("data"));
    gauss.plotOn(frame,RooFit::Name("gauss"),RooFit::LineColor(2));
    ExtendedCanvas *c = new ExtendedCanvas("c"+type,type,800,600,2);
    TPad* p1 = (TPad*)c->GetPad(1);
    TPad* p2 = (TPad*)c->GetPad(2);
    p1->cd();
    frame->Draw();
    p2->cd();
    SignificanceHist SH(*frame->getHist("data"),*frame->getCurve("gauss"));
    TH1F* hsigni = SH.GetChiHist(3);
    hsigni->GetXaxis()->SetTitle(frame->GetXaxis()->GetTitle());
    hsigni->Draw("HIST");
    std::cout << "Chi2/Ndf = " 
	      << frame->chiSquare("gauss","data",2)
	      << std::endl;
  }

  std::pair<double,double> temp(mean.getVal(),width.getVal());
  if(doplot) std::cout << type << " = " << temp.first << "+/-" << temp.second << std::endl;

  return temp;
}

/////////////////////////////////////////////////////////////////////////////////////////////
std::pair<double,double> BkgEstimatorUpgrade::GetPurityInterval(TString treefile,bool doplot)
/////////////////////////////////////////////////////////////////////////////////////////////
{
  //---> Perform a fit of the purity and determine the 1sigma band interval 


  TFile* file0 = TFile::Open(treefile,"READ");
  TTree* yieldstree = (TTree*)file0->Get("yieldstree");
  double Xmin = 1e18;
  double Xmax = -1e18;
  int nentries = (int)yieldstree->GetEntriesFast();
  TreeReader Rd(yieldstree);

  double purity = 0;
  TFile* file = new TFile("rootfiles/toto.root","RECREATE");
  file->cd();
  TTree * tree_fit = new TTree("tree_fit","tree_fit"); 
  tree_fit->Branch("Purity",&purity,"Purity/D");
  for( int ipe = 0 ; ipe<nentries;ipe++){
    Rd.GetEntry(ipe);
    purity = Rd.GetVariable("Ngg/(Ngg+Ngj+Njg+Njj)");
    tree_fit->Fill();
    if( purity<Xmin) Xmin=purity;
    if( purity>Xmax) Xmax=purity;
  }

  RooRealVar Purity("Purity","Purity",Xmin,Xmax );
  double half= (Xmax-Xmin)/2.;
  RooDataSet data("data","data",RooArgSet(Purity),RooFit::Import(*tree_fit) );
  RooRealVar mean("mean","mean",Xmin+half,Xmin,Xmax) ;
  RooRealVar width("width","width",half/4.,0.001,half) ;
  RooGaussian gauss("gauss","gaussian PDF",Purity,mean,width) ;
  RooMsgService::instance().setSilentMode(true);
  gauss.fitTo(data);

  frame_Purity = Purity.frame();
  frame_Purity->SetName("Purity");
  data.plotOn(frame_Purity,RooFit::Name("data"));
  gauss.plotOn(frame_Purity,RooFit::Name("gauss"),RooFit::LineColor(2));

  std::pair<double,double> temp(mean.getVal()-width.getVal(),mean.getVal()+width.getVal());
  return temp;
}

///////////////////////////////////////////////////
void BkgEstimatorUpgrade::PrintLoosePrimeResults()
///////////////////////////////////////////////////
{
  std::cout << "Type    |       Ngg       |       Ngj       |       Njg       |       Njj       |" << std::endl;
  std::pair<double,double> yield = GetYield("Ngg",name_f_frac_L2);
  std::cout << "Loose'2 | " << yield.first << "+/-" << yield.second << " | ";
  yield = GetYield("Ngj",name_f_frac_L2);
  std::cout << yield.first << "+/-" << yield.second << " | ";
  yield = GetYield("Njg",name_f_frac_L2);
  std::cout << yield.first << "+/-" << yield.second << " | ";
  yield = GetYield("Njj",name_f_frac_L2);
  std::cout << yield.first << "+/-" << yield.second << " | " << std::endl;
  yield = GetYield("Ngg",name_f_frac_L3);
  std::cout << "Loose'3 | " << yield.first << "+/-" << yield.second << " | ";
  yield = GetYield("Ngj",name_f_frac_L3);
  std::cout << yield.first << "+/-" << yield.second << " | ";
  yield = GetYield("Njg",name_f_frac_L3);
  std::cout << yield.first << "+/-" << yield.second << " | ";
  yield = GetYield("Njj",name_f_frac_L3);
  std::cout << yield.first << "+/-" << yield.second << " | " << std::endl;
  yield = GetYield("Ngg",name_f_frac_L4);
  std::cout << "Loose'4 | " << yield.first << "+/-" << yield.second << " | ";
  yield = GetYield("Ngj",name_f_frac_L4);
  std::cout << yield.first << "+/-" << yield.second << " | ";
  yield = GetYield("Njg",name_f_frac_L4);
  std::cout << yield.first << "+/-" << yield.second << " | ";
  yield = GetYield("Njj",name_f_frac_L4);
  std::cout << yield.first << "+/-" << yield.second << " | " << std::endl;
  yield = GetYield("Ngg",name_f_frac_L5);
  std::cout << "Loose'5 | " << yield.first << "+/-" << yield.second << " | ";
  yield = GetYield("Ngj",name_f_frac_L5);
  std::cout << yield.first << "+/-" << yield.second << " | ";
  yield = GetYield("Njg",name_f_frac_L5);
  std::cout << yield.first << "+/-" << yield.second << " | ";
  yield = GetYield("Njj",name_f_frac_L5);
  std::cout << yield.first << "+/-" << yield.second << " |" << std::endl;
}
