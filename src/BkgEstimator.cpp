#include "BkgEstimator.h"
#include <TError.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooPlot.h>
#include <RooMsgService.h>

#include "BkgUncertaintyBuilder.h"
#include "ToolsCommons.h"
#include "ToolsSignificanceHist.h"
#include "ToolsExtendedCanvas.h"

///////////////////////////////////
BkgEstimator::BkgEstimator()
///////////////////////////////////
{
  //Default constructor
  CreateBinning();
  SetVerbose(false);
}
/////////////////////////////////////
BkgEstimator::~BkgEstimator()
////////////////////////////////////
{
}

///////////////////////////////////////////////////////////////////////////////
void BkgEstimator::SetYieldSystFiles(TString st_L2,TString st_L3,TString st_L5)
///////////////////////////////////////////////////////////////////////////////
{
  name_f_frac_L2 = st_L2;
  name_f_frac_L3 = st_L3;
  name_f_frac_L5 = st_L5;
}


/////////////////////////////////////
void BkgEstimator::Init()
////////////////////////////////////
{
  //--> First you need to call the file setters
  if(m_VERBOSE) std::cout << "Read the files";
  //----------------------------------------------------------------
  m_f_data         = new TFile(name_f_data,"read");
  m_f_red          = new TFile(name_f_red,"read");
  m_f_irr          = new TFile(name_f_irr,"read");
  m_f_irr_syst     = new TFile(name_f_irr_syst,"read");
  m_f_frac_L2      = new TFile(name_f_frac_L2,"read");
  m_f_frac_L3      = new TFile(name_f_frac_L3,"read");
  m_f_frac_L4      = new TFile(name_f_frac_L4,"read");
  m_f_frac_L5      = new TFile(name_f_frac_L5,"read");
  // m_f_RS_1250_01   = new TFile(name_f_RS_1250_01,"read");
  // m_f_RS_1500_01   = new TFile(name_f_RS_1500_01,"read");
  // m_f_RS_1750_01   = new TFile(name_f_RS_1750_01,"read");
  // m_f_ADD_2500_GRW = new TFile(name_f_ADD_2500_GRW,"read");
  // m_f_ADD_3000_GRW = new TFile(name_f_ADD_3000_GRW,"read");
  // m_f_ADD_3500_GRW = new TFile(name_f_ADD_3500_GRW,"read");
  if(m_VERBOSE) std::cout << "... done !"<< std::endl;

  if(m_VERBOSE) std::cout << "Get the object";
  //---------------------------------------------------------------
  m_hh_data  = (TH1D*)m_f_data->Get("mgg_data");
  //---------------------------------------------------------------
  m_gr_gj         = (TGraphAsymmErrors*)m_f_red->Get("graph_h_subleading_fake");
  m_gr_jg         = (TGraphAsymmErrors*)m_f_red->Get("graph_h_leading_fake");
  m_gr_jj         = (TGraphAsymmErrors*)m_f_red->Get("graph_h_double_fake");
  m_hh_gj         = (TH1D*)m_f_red->Get("h_subleading_fake");
  m_hh_jg         = (TH1D*)m_f_red->Get("h_leading_fake");
  m_hh_jj         = (TH1D*)m_f_red->Get("h_double_fake");
  //--------------------------------------------------------------------
  // m_hh_irr      = (TH1D*)m_f_irr->Get("irreducible_shape_nokfac");
  m_hh_irr      = (TH1D*)m_f_irr->Get("irreducible_shape");
  m_hh_gg       = (TH1D*)m_hh_irr->Clone("m_hh_gg");
  //--------------------------------------------------------------------
  m_hh_irr_syst = (TH1D*)m_f_irr_syst->Get("hh_irreducible_uncert");
  m_hh_syst     = (TH1D*)m_f_irr_syst->Get("hh_total_uncert");//MUST BE RECOMPUTED
  //---------------------------------------------------------------------
  m_yieldstree_L2  = (TTree*)m_f_frac_L2 ->Get("yieldstree");
  m_yieldstree_L3  = (TTree*)m_f_frac_L3 ->Get("yieldstree");
  m_yieldstree_L4  = (TTree*)m_f_frac_L4 ->Get("yieldstree");
  m_yieldstree_L5  = (TTree*)m_f_frac_L5 ->Get("yieldstree");
  //----------------------------------------------------------------------------------------------------------
  if(m_VERBOSE) std::cout << "... done !"<< std::endl;

  if(m_VERBOSE) std::cout << "Define the search histograms";
  m_hh_data_400      = new TH1D("hh_data","hh_data",Commons::nBins_search,Commons::binning_search);
  m_hh_irr_400       = new TH1D("hh_irr_400","hh_irr_400",Commons::nBins_search,Commons::binning_search);
  m_hh_red_400       = new TH1D("hh_red_400","hh_red_400",Commons::nBins_search,Commons::binning_search);
  m_hh_bkg_400       = new TH1D("hh_bkg_400","hh_bkg_400",Commons::nBins_search,Commons::binning_search);
  m_hh_syst_400      = new TH1D("bkg_total_syst_gg","bkg_total_syst_gg",Commons::nBins_search,Commons::binning_search);
  m_hh_frac_syst_400 = new TH1D("bkg_purity_syst_gg","bkg_purity_syst_gg",Commons::nBins_search,Commons::binning_search);
  m_hh_red_syst_400  = new TH1D("bkg_reducible_syst_gg","bkg_reducible_syst_gg",Commons::nBins_search,Commons::binning_search);
  m_hh_irr_syst_400  = new TH1D("bkg_irreducible_syst_gg","bkg_irreducible_syst_gg",Commons::nBins_search,Commons::binning_search);
  //-----------------------------------------------------------------------------------------------------------
  dataNorm = m_hh_data->Integral(Commons::norm_bin.first,Commons::norm_bin.second);
  //---------------------------------------------------------
  // m_hh_RS_1250_01  = (TH1F*)m_f_RS_1250_01->Get("hmgg_log2");
  // m_hh_RS_1500_01  = (TH1F*)m_f_RS_1500_01->Get("hmgg_log2");
  // m_hh_RS_1750_01  = (TH1F*)m_f_RS_1750_01->Get("hmgg_log2");
  //-------------------------------------------------------------
  // m_hh_ADD_2500_GRW  = (TH1F*)m_f_ADD_2500_GRW->Get("hmgg_log2");
  // m_hh_ADD_3000_GRW  = (TH1F*)m_f_ADD_3000_GRW->Get("hmgg_log2");
  // m_hh_ADD_3500_GRW  = (TH1F*)m_f_ADD_3500_GRW->Get("hmgg_log2");
  //-------------------------------------------------------------
  if(m_VERBOSE) std::cout << "... done !"<< std::endl;

  if(m_VERBOSE) std::cout << "Set yields";
  std::pair<double,double> y_gg = GetYield("Ngg",m_yieldstree_L4);
  std::pair<double,double> y_gj = GetYield("Ngj",m_yieldstree_L4);
  std::pair<double,double> y_jg = GetYield("Njg",m_yieldstree_L4);
  std::pair<double,double> y_jj = GetYield("Njj",m_yieldstree_L4);
  Commons::SetYields(y_gg.first,y_gj.first,y_jg.first,y_jj.first,true);
  Commons::SetYieldsStatError(y_gg.second,y_gj.second,y_jg.second,y_jj.second);
  if(m_VERBOSE) std::cout << "Build total bkg";
  BuildTotalBkg();//Build the total background 
  if(m_VERBOSE) std::cout << "... done !"<< std::endl;
  if(m_VERBOSE) std::cout << "Remove stat errors bars on the red hist";
  RemoveReducibleStatErrors();//Set stat errors of the red hists to 0
  if(m_VERBOSE) std::cout << "... done !"<< std::endl;
  if(m_VERBOSE) std::cout << "Build total bkg uncert ... " << std::endl;
  BuildTotalBkgUncertainty();//Build the total bkg syst
  if(m_VERBOSE) std::cout << "Build search hist";
  BuildHist400();//Create the histograms for BAT
  if(m_VERBOSE) std::cout << "... done !"<< std::endl;


  // BuildRSHists();//Build RS signal examples
  // BuildADDHists();//Build ADD signal example

  BuildBkgErrorsGraph();
  BuildDataErrorsGraph();
  BuildSignificanceHist();

}
//////////////////////////////////
void BkgEstimator::CreateBinning()
//////////////////////////////////
{
  // std::vector<double> vec_bin400;
  // for (int i=1 ;i<=Commons::nBins ;i++){
  //   if(Commons::binning[i]>400.){
  //     vec_bin400.push_back(Commons::binning[i]);
  //   }
  // }
  // nBins400 = vec_bin400.size()-1;
  // binning400_log = new double[nBins400+1];
  // for (int i=0 ;i<=nBins400 ;i++){
  //   binning400_log[i] = vec_bin400[i];
  // }
}
///////////////////////////////////////////////
void BkgEstimator::RemoveReducibleStatErrors()
///////////////////////////////////////////////
{
  int Nbins = m_hh_red->GetNbinsX();
  for( int ibin=1; ibin<Nbins+1;ibin++){
    m_hh_red->SetBinError(ibin,0);
    m_hh_gg->SetBinError(ibin,0);
    m_hh_gj->SetBinError(ibin,0);
    m_hh_jg->SetBinError(ibin,0);
    m_hh_gj->SetBinError(ibin,0);
  }
}
/////////////////////////////////
void BkgEstimator::BuildHist400()
/////////////////////////////////
{
  int Nbins_data = (int)m_hh_data->GetNbinsX();
  int ibin_400 =1;
  for( int ibin=1 ; ibin <= Nbins_data;ibin++){
    if( m_hh_data->GetBinLowEdge(ibin) < 400) continue;
    m_hh_data_400->SetBinContent(ibin_400,m_hh_data->GetBinContent(ibin) );
    m_hh_red_400->SetBinContent(ibin_400,m_hh_red->GetBinContent(ibin) );
    m_hh_irr_400->SetBinContent(ibin_400,m_hh_irr->GetBinContent(ibin) );
    m_hh_syst_400->SetBinContent(ibin_400,m_hh_syst->GetBinContent(ibin) );
    m_hh_frac_syst_400->SetBinContent(ibin_400,m_hh_frac_syst->GetBinContent(ibin) );
    m_hh_irr_syst_400->SetBinContent(ibin_400,m_hh_irr_syst->GetBinContent(ibin) );
    m_hh_red_syst_400->SetBinContent(ibin_400,m_hh_red_syst->GetBinContent(ibin) );
    m_hh_red_400->SetBinError(ibin_400,0. );
    m_hh_irr_400->SetBinError(ibin_400,m_hh_irr->GetBinError(ibin) );
    ibin_400++;
  }
}
//////////////////////////////////
void BkgEstimator::BuildTotalBkg()
//////////////////////////////////
{



  m_hh_bkg= GetTotalBkgHist(m_hh_gg,m_hh_gj,m_hh_jg,m_hh_jj,Commons::yields);
  m_hh_bkg->SetName("bkg_total_gg_140");
  m_hh_bkg->SetTitle("bkg_total_gg");

  std::vector<double> yields_red;
  yields_red.push_back(Commons::yields[1]);
  yields_red.push_back(Commons::yields[2]);
  yields_red.push_back(Commons::yields[3]);
  double sumyields = 0;
  for(int iy=0;iy<(int)Commons::yields.size();iy++) sumyields += Commons::yields[iy];
  double purity = Commons::yields[0]/sumyields;

  m_hh_red = GetRedBkgHist(m_hh_gj,m_hh_jg,m_hh_jj,yields_red);
  m_hh_red->SetName("bkg_red_gg_140");
  m_hh_red->SetTitle("bkg_red_gg");
  m_hh_red->Scale(1-purity);

}
/////////////////////////////////
void BkgEstimator::BuildRSHists()
/////////////////////////////////
{
  double lumi = Commons::int_lumi_fb ;
  double kfac = 1.75;
  double xsec_1250_01 = 36.6;
  double xsec_1500_01 = 10.2;
  double xsec_1750_01 = 3.31;
  m_hh_RS_1250_01->Scale(kfac*xsec_1250_01*lumi/10000);
  m_hh_RS_1500_01->Scale(kfac*xsec_1500_01*lumi/10000);
  m_hh_RS_1750_01->Scale(kfac*xsec_1750_01*lumi/10000);
  m_hs_RS_1250_01 = new THStack("hs_RS_1250_01","hs");
  m_hs_RS_1250_01->Add(m_hh_bkg);
  m_hs_RS_1250_01->Add(m_hh_RS_1250_01);
  m_hh_RS_1250_01->SetLineColor(2);
  m_hs_RS_1500_01 = new THStack("hs_RS_1500_01","hs");
  m_hs_RS_1500_01->Add(m_hh_bkg);
  m_hs_RS_1500_01->Add(m_hh_RS_1500_01);
  m_hh_RS_1500_01->SetLineColor(3);
  m_hs_RS_1750_01 = new THStack("hs_RS_1750_01","hs");
  m_hs_RS_1750_01->Add(m_hh_bkg);
  m_hs_RS_1750_01->Add(m_hh_RS_1750_01);
  m_hh_RS_1750_01->SetLineColor(4);
}

/////////////////////////////////
void BkgEstimator::BuildADDHists()
/////////////////////////////////
{
  int binfirst = 53;
  int binlast  = 54;
  double bkgnorm = m_hh_bkg->Integral(binfirst,binlast);
  if(m_VERBOSE){
    std::cout << std::endl;
    std::cout << "ADD normalization : [ "
	      << Commons::binning[binfirst] << ", "
	      << Commons::binning[binlast] << " ]" 
	      << std::endl;
  }
  double int_2500_GRW = m_hh_ADD_2500_GRW->Integral();
  double int_3000_GRW = m_hh_ADD_3000_GRW->Integral();
  double int_3500_GRW = m_hh_ADD_3500_GRW->Integral();
  m_hh_ADD_2500_GRW->Scale(1./int_2500_GRW);
  m_hh_ADD_3000_GRW->Scale(1./int_3000_GRW);
  m_hh_ADD_3500_GRW->Scale(1./int_3500_GRW);
  double int_2500_GRW_n = m_hh_ADD_2500_GRW->Integral(binfirst,binlast);
  double int_3000_GRW_n = m_hh_ADD_3000_GRW->Integral(binfirst,binlast);
  double int_3500_GRW_n = m_hh_ADD_3500_GRW->Integral(binfirst,binlast);
  m_hh_ADD_2500_GRW->Scale(bkgnorm/int_2500_GRW_n);
  m_hh_ADD_3000_GRW->Scale(bkgnorm/int_3000_GRW_n);
  m_hh_ADD_3500_GRW->Scale(bkgnorm/int_3500_GRW_n);

  //------------------------------------------------------
  m_hs_ADD_2500_GRW = new THStack("hs_ADD_2500_GRW","hs");
  m_hs_ADD_2500_GRW->Add(m_hh_ADD_2500_GRW);
  m_hh_ADD_2500_GRW->SetLineColor(2);
  //------------------------------------------------------
  m_hs_ADD_3000_GRW = new THStack("hs_ADD_3000_GRW","hs");
  m_hs_ADD_3000_GRW->Add(m_hh_ADD_3000_GRW);
  m_hh_ADD_3000_GRW->SetLineColor(3);
  //------------------------------------------------------
  m_hs_ADD_3500_GRW = new THStack("hs_ADD_3500_GRW","hs");
  m_hs_ADD_3500_GRW->Add(m_hh_ADD_3500_GRW);
  m_hh_ADD_3500_GRW->SetLineColor(4);
  //------------------------------------------------------
}
/////////////////////////////////////////////
void BkgEstimator::BuildTotalBkgUncertainty()
/////////////////////////////////////////////
{
  if(m_VERBOSE) std::cout << "Compute bkg uncertainty components";
  YieldsSystUncert();//Fill m_hh_frac_syst
  YieldsStatUncert();//Fill m_hh_frac_stat
  RedShapeUncert();//Fill m_hh_red_syst
  if(m_VERBOSE) std::cout << "... done !" << std::endl;

  if(m_VERBOSE) std::cout << "Compute bkg uncertainty ";
  for(int ibin = 1 ;ibin<=Commons::nBins;ibin++){
    double frac_stat = m_hh_frac_stat->GetBinContent(ibin);
    double frac_syst = m_hh_frac_syst->GetBinContent(ibin);
    double irr_shape = m_hh_irr_syst->GetBinContent(ibin);
    double red_shape = m_hh_red_syst->GetBinContent(ibin);

    double data_stat = 0.0001;
    double error = sqrt(frac_stat*frac_stat+
			frac_syst*frac_syst+
			irr_shape*irr_shape+
			red_shape*red_shape+
			data_stat*data_stat
			);
    m_hh_syst->SetBinContent(ibin,error);
    m_hh_syst->SetBinError(ibin,0);
  }
  if(m_VERBOSE) std::cout << "... done !" << std::endl;
  if(m_VERBOSE) std::cout << "Compute bkg histogram with systematics included ";
  //---- Build total background with systematics included for BumpHunter
  m_hh_bkg_syst = new TH1D("hh_bkg_syst_gg","hh_bkg_syst_gg",Commons::nBins,Commons::binning);
  if(m_VERBOSE) std::cout << "... done !" << std::endl;
  if(m_VERBOSE) std::cout << "ibin | value    | syst error | stat error | tot error | relative error|" << std::endl;
  for( int ibin=1; ibin <= Commons::nBins ;ibin++){
    double val      = m_hh_bkg->GetBinContent(ibin);
    double syst_err = m_hh_syst->GetBinContent(ibin)*val;
    double stat_err = m_hh_bkg->GetBinError(ibin);
    //--> Add errors in quadrature
    double tot_err  = sqrt(stat_err*stat_err+syst_err*syst_err);
    if(m_VERBOSE && ibin> 19) {
      std::cout << "-----------------------------------"
		<< "-----------------------------------"
		<< std::endl;
      std::cout << ibin     << "       " 
		<< val      << "    "
		<< syst_err << "    " 
		<< stat_err << "    " 
		<< tot_err  << "    "
		<< tot_err/val 
		<< std::endl;
      m_hh_bkg_syst->SetBinContent(ibin,val);
      m_hh_bkg_syst->SetBinError(ibin,tot_err);
    }
  }
  
}
////////////////////////////////////////
void BkgEstimator::BuildBkgErrorsGraph()
///////////////////////////////////////
{
  m_gmgg_bkg = new TGraphAsymmErrors();
  m_gmgg_red = new TGraphAsymmErrors();
  m_gmgg_bkg->SetName("gmgg_bkg");
  m_gmgg_red->SetName("gmgg_red");
  m_gmgg_bkg->SetTitle("gmgg_bkg");
  m_gmgg_red->SetTitle("gmgg_red");

  for (int ibin = 1 ; ibin <= m_hh_bkg_syst->GetNbinsX()+1;ibin++){
    double value  = m_hh_bkg_syst->GetBinContent(ibin);
    double erry   = m_hh_bkg_syst->GetBinError(ibin);
    double errx   = m_hh_bkg_syst->GetBinWidth(ibin)/2;
    if(value !=0){
      m_gmgg_bkg->SetPoint(ibin-1, m_hh_bkg_syst->GetBinCenter(ibin), value);
      m_gmgg_bkg->SetPointError(ibin-1, errx, errx, erry, erry);
    }      
    double value_red  = m_hh_red->GetBinContent(ibin);
    double erry_red   = m_hh_red->GetBinContent(ibin)*m_hh_red_syst->GetBinContent(ibin);
    double errx_red   = m_hh_red->GetBinWidth(ibin)/2;
    if(value_red !=0){
      m_gmgg_red->SetPoint(ibin-1, m_hh_red->GetBinCenter(ibin), value_red);
      m_gmgg_red->SetPointError(ibin-1, errx_red, errx_red, erry_red, erry_red);
    }      
  }
}
//////////////////////////////////////////
void BkgEstimator::BuildDataErrorsGraph()
/////////////////////////////////////////
{
  m_gmgg_data = new TGraphAsymmErrors();
  m_gmgg_data->SetName("gmgg_data");
  m_gmgg_data->SetTitle("gmgg_data");

  for (int ibin = 1 ; ibin <= m_hh_data->GetNbinsX()+1;ibin++){
    double value = m_hh_data->GetBinContent(ibin);
    if( value!=0){
      double y1 = value + 1.0;
      double d1 = 1.0 - 1.0/(9.0*y1) + 1.0/(3*sqrt(y1));
      double error_poisson_up = y1*d1*d1*d1 - value;
      double y2 = value;
      double d2 = 1.0 - 1.0/(9.0*y2) - 1.0/(3.0*sqrt(y2));
      double error_poisson_down = value - y2*d2*d2*d2;
      m_gmgg_data->SetPoint(ibin-1, m_hh_data->GetBinCenter(ibin), value);
      m_gmgg_data->SetPointError(ibin-1, 0, 0, error_poisson_down, error_poisson_up);
    }      
  }
}
///////////////////////////////////////////
void BkgEstimator::BuildSignificanceHist()
/////////////////////////////////////////
{
  //---> from arXiv:1111.2062v2 
  //---------------------------------------------------//
  //---------- NO SYSTEMATIC --------------------------//
  //---------------------------------------------------//
  SignificanceHist SH(*(TH1F*)m_hh_data,*(TH1F*)m_hh_bkg);
  m_hh_signi = SH.GetSignificanceHist(4);
  m_hh_signi->SetName("hh_signi");
  m_hh_signi->SetTitle("hh_signi");
  //---------------------------------------------------------------------//
  //--------------------- WITH SYSTEMATIC -------------------------------//
  //---------------------------------------------------------------------//
  SignificanceHist SH_withsyst(*(TH1F*)m_hh_data,*(TH1F*)m_hh_bkg_syst);
  m_hh_signi_sys = SH_withsyst.GetSignificanceHist_withsyst();
  m_hh_signi_sys->SetName("hh_signi_sys");
  m_hh_signi_sys->SetTitle("hh_signi_sys");

}
///////////////////////////////////////////////////////////////////////
void BkgEstimator::PlotBkgEstimate(bool doRS, bool doADD, bool doRSADD)
///////////////////////////////////////////////////////////////////////
{
  //-------------------------------------------------------------------------
  m_hh_bkg->GetXaxis()->SetRange(Commons::norm_bin.first,Commons::nBins);
  m_hh_bkg->GetYaxis()->SetRangeUser(0.0001,5000);
  m_hh_bkg->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  m_hh_bkg->GetYaxis()->SetTitle("Events/bin");
  m_hh_bkg->GetXaxis()->SetMoreLogLabels();
  m_hh_bkg->GetXaxis()->SetNoExponent();
  m_hh_bkg->SetLineColor(4);
  m_gmgg_bkg->SetFillColor(kOrange);
  m_gmgg_bkg->SetLineColor(kOrange);
  m_gmgg_bkg->SetMarkerSize(0);
  m_hh_red->SetLineColor(4);
  m_hh_red->SetLineStyle(2);
  m_gmgg_red->SetFillColor(kYellow);
  m_gmgg_red->SetLineColor(kYellow);
  //---------------------------------------------------------------------------
  TLegend *leg = new TLegend(0.2,0.2,0.4,0.4);
  leg->SetFillColor(0);
  leg->AddEntry(m_gmgg_data,"2012 data","lp");
  leg->AddEntry(m_hh_bkg,"Total Background","f");
  leg->AddEntry(m_hh_red,"Reducible Background","f");
  leg->AddEntry(m_gmgg_bkg,"syst #oplus stat (total)","f");
  if(doRS){
    leg->AddEntry(m_hh_RS_1250_01,"RS k/M_{Pl} = 0.1, m_{G} = 1.25 TeV","f");
    leg->AddEntry(m_hh_RS_1500_01,"RS k/M_{Pl} = 0.1, m_{G} = 1.5 TeV","f");
    leg->AddEntry(m_hh_RS_1750_01,"RS k/M_{Pl} = 0.1, m_{G} = 1.75 TeV","f");
  }
  if(doADD){
    leg->AddEntry(m_hh_ADD_2500_GRW,"ADD, GRW, m_{S}=2.5 TeV","f");
    leg->AddEntry(m_hh_ADD_3000_GRW,"ADD, GRW, m_{S}=3.0 TeV","f");
    leg->AddEntry(m_hh_ADD_3500_GRW,"ADD, GRW, m_{S}=3.5 TeV","f");
  }
  if(doRSADD){
    leg->AddEntry(m_hh_RS_1500_01,"RS k/M_{Pl} = 0.1, m_{G} = 1.5 TeV","f");
    leg->AddEntry(m_hh_ADD_2500_GRW,"ADD, GRW, m_{S} = 2.5 TeV","f");
  }
  //--------------------------------------------------------------------------
  TLatex *t = new TLatex();
  TString writetext1 = "#font[72]{ATLAS} Internal";
  TString writetext2 = Form("#int L dt = %1.1f fb^{-1}",Commons::int_lumi_fb);
  TString writetext3 = "#sqrt{s} = 8 TeV";
  t->SetTextAlign(13);
  t->SetTextColor(kBlack);
  t->SetTextSize(0.035);


  //---------------------------------------------------------------------//
  TCanvas * cbkg = new TCanvas("cbkgestimate","Bkg estimation",800,600);
  cbkg->SetLogx();cbkg->SetLogy();cbkg->cd();
  m_hh_bkg->Draw("HIST][");
  if(doRS){
    m_hs_RS_1250_01->Draw("sameHIST");
    m_hs_RS_1500_01->Draw("sameHIST");
    m_hs_RS_1750_01->Draw("sameHIST");
  }
  if(doADD){
    m_hs_ADD_2500_GRW->Draw("sameHIST");
    m_hs_ADD_3000_GRW->Draw("sameHIST");
    m_hs_ADD_3500_GRW->Draw("sameHIST");
  }
  if(doRSADD){
    m_hs_RS_1500_01->Draw("sameHIST");
    m_hs_ADD_2500_GRW->Draw("sameHIST");
  }
  TH1D* hh_bkg2 = (TH1D*)m_hh_bkg->Clone("hh_bkg2");
  hh_bkg2->SetFillColor(10);
  hh_bkg2->Draw("sameHIST");
  m_hh_bkg->Draw("sameaxis");
  // m_gmgg_bkg->Draw("sameE2");
  m_hh_bkg->Draw("sameHIST");
  m_gmgg_data->Draw("sameP");
  m_hh_red->Draw("sameHIST][");
  leg->Draw("same");
  t->DrawLatex(1500,2000,writetext1);
  t->DrawLatex(1500,500,writetext2);
  t->DrawLatex(1500,100,writetext3);

  //-----------------------------------------------------------------------------
  TCanvas *cbkg_signi= AnalysisTools::Get2PadPlot("cbkg_signi","bkg with significance");
  TPad* p1 = (TPad*)cbkg_signi->GetPad(1);
  p1->SetLogy();p1->SetLogx();
  TPad* p2 = (TPad*)cbkg_signi->GetPad(2);
  p2->SetLogx();
  p1->cd();
  m_hh_bkg->Draw("HIST][");
  // m_gmgg_bkg->Draw("sameE2");
  // m_hh_bkg->Draw("sameHIST");
  m_gmgg_data->Draw("sameP");
  m_hh_red->Draw("sameHIST][");
  leg->Draw("same");
  t->DrawLatex(1500,2000,writetext1);
  t->DrawLatex(1500,500,writetext2);
  t->DrawLatex(1500,100,writetext3);
  p1->Update();
  p2->cd();
  m_hh_signi->GetXaxis()->SetRange(Commons::norm_bin.first,Commons::nBins);
  m_hh_signi->Draw("HIST");
  m_hh_signi->GetXaxis()->SetMoreLogLabels();
  m_hh_signi->GetXaxis()->SetNoExponent();
  p2->Update();
}

///////////////////////////////////////
void BkgEstimator::PlotBkgUncertainty()
///////////////////////////////////////
{
  TCanvas *c1 = new TCanvas("cbkguncert","Bkg uncertainty",800,600);
  c1->cd();
  m_hh_frac_stat->SetLineColor(2);
  m_hh_frac_syst->SetLineColor(3);
  m_hh_red_syst->SetLineColor(4);
  m_hh_irr_syst->SetLineColor(6);
  m_hh_frac_stat->GetXaxis()->SetRange(Commons::norm_bin.first,Commons::nBins);
  m_hh_frac_stat->GetYaxis()->SetRange(-0.01,0.40);
  m_hh_syst->Draw("HIST][");
  m_hh_frac_stat->Draw("sameHIST][");
  m_hh_frac_syst->Draw("sameHIST][");
  m_hh_irr_syst->Draw("sameHIST][");
  m_hh_red_syst->Draw("sameHIST][");
  TLegend *leg = new TLegend(0.2,0.7,0.5,0.9);
  leg->SetFillColor(0);
  leg->AddEntry(m_hh_syst,"Total uncertainty","l");
  leg->AddEntry(m_hh_irr_syst,"Irreducible background","l");
  leg->AddEntry(m_hh_frac_syst,"Yields estimate","l");
  leg->AddEntry(m_hh_red_syst,"Reducible background","l");
  leg->Draw("same");
  TLatex *t = new TLatex();
  TString writetext1 = "#font[72]{ATLAS} For Approval";
  t->SetTextAlign(13);
  t->SetTextColor(kBlack);
  t->SetTextSize(0.035);
  t->DrawLatex(2000,0.3,writetext1);

}
/////////////////////////////////////////
void BkgEstimator::GetFinalYieldsPerBin()
/////////////////////////////////////////
{
  int Nbins = m_hh_data->GetNbinsX();
  std::cout << "Bin"     << " | "
	    << "low edge"<< " | " 
	    << "Bakground expectation" << " | " 
	    << "Data" << std::endl;

  for( int ibin=1; ibin< Nbins ; ibin++){
    std::cout << ibin  << " | "
	      << m_hh_data->GetBinLowEdge(ibin)     << " | " 
	      << m_hh_bkg_syst->GetBinContent(ibin) << "+/-"
	      << m_hh_bkg_syst->GetBinError(ibin)   << " | "
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
      double rel_err = m_hh_syst->GetBinContent(i);
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
///////////////////////////////////////////////////
void BkgEstimator::GetFinalYieldsPerBin_NoUncert()
//////////////////////////////////////////////////
{

  std::cout << " Bin | reducible | irreducible | total | data | data-total" << std::endl;
  for( int ibin = 1; ibin<=m_hh_bkg->GetNbinsX();ibin++){
    std::cout << Form("[%1.1f,%1.1f]",m_hh_bkg->GetBinLowEdge(ibin),m_hh_bkg->GetBinLowEdge(ibin+1))
	      << " | "
	      << m_hh_red->GetBinContent(ibin) 
	      << "| "
	      << m_hh_gg->GetBinContent(ibin) 
	      << "| "
	      << m_hh_bkg->GetBinContent(ibin) 
	      << "| "
	      << m_hh_data->GetBinContent(ibin) 
	      << "| "
	      << m_hh_data->GetBinContent(ibin)-m_hh_bkg->GetBinContent(ibin)  
	      << std::endl;
  }


}






/////////////////////////////////////////////////
void BkgEstimator::YieldsSystUncert(bool doplot)
////////////////////////////////////////////////
{
  
  std::pair<double,double> y_gg_L2 = GetYield("Ngg",m_yieldstree_L2);
  std::pair<double,double> y_gj_L2 = GetYield("Ngj",m_yieldstree_L2);
  std::pair<double,double> y_jg_L2 = GetYield("Njg",m_yieldstree_L2);
  std::pair<double,double> y_jj_L2 = GetYield("Njj",m_yieldstree_L2);

  std::pair<double,double> y_gg_L3 = GetYield("Ngg",m_yieldstree_L3);
  std::pair<double,double> y_gj_L3 = GetYield("Ngj",m_yieldstree_L3);
  std::pair<double,double> y_jg_L3 = GetYield("Njg",m_yieldstree_L3);
  std::pair<double,double> y_jj_L3 = GetYield("Njj",m_yieldstree_L3);


  std::pair<double,double> y_gg_L5 = GetYield("Ngg",m_yieldstree_L5);
  std::pair<double,double> y_gj_L5 = GetYield("Ngj",m_yieldstree_L5);
  std::pair<double,double> y_jg_L5 = GetYield("Njg",m_yieldstree_L5);
  std::pair<double,double> y_jj_L5 = GetYield("Njj",m_yieldstree_L5);

  double Yield_L2[] = {y_gg_L2.first,y_gj_L2.first,y_jg_L2.first,y_jj_L2.first};
  double Yield_L3[] = {y_gg_L3.first,y_gj_L3.first,y_jg_L3.first,y_jj_L3.first};
  double Yield_L4[] = {Commons::yields[0],Commons::yields[1],Commons::yields[2],Commons::yields[3]};
  double Yield_L5[] = {y_gg_L5.first,y_gj_L5.first,y_jg_L5.first,y_jj_L5.first};

  // double Yield_L2[] = {9646,3086,1222,1707};
  // double Yield_L3[] = {10120,2986,1207,1379};
  // double Yield_L4[] = {10812,2640,1057,1167};
  // double Yield_L5[] = {11647,2265,886,890};

  std::vector<double> Yield_L2_vec(Yield_L2,Yield_L2+4);
  std::vector<double> Yield_L3_vec(Yield_L3,Yield_L3+4);
  std::vector<double> Yield_L4_vec(Yield_L4,Yield_L4+4);
  std::vector<double> Yield_L5_vec(Yield_L5,Yield_L5+4);

 
  TH1D* hh_bkg_L2 = GetTotalBkgHist(m_hh_gg,m_hh_gj,m_hh_jg,m_hh_jj,Yield_L2_vec);
  TH1D* hh_bkg_L3 = GetTotalBkgHist(m_hh_gg,m_hh_gj,m_hh_jg,m_hh_jj,Yield_L3_vec); 
  TH1D* hh_bkg_L4 = GetTotalBkgHist(m_hh_gg,m_hh_gj,m_hh_jg,m_hh_jj,Yield_L4_vec); 
  TH1D* hh_bkg_L5 = GetTotalBkgHist(m_hh_gg,m_hh_gj,m_hh_jg,m_hh_jj,Yield_L5_vec); 

  hh_bkg_L2->SetName("hh_bkg_L2"); hh_bkg_L2->SetLineColor(2);
  hh_bkg_L3->SetName("hh_bkg_L3"); hh_bkg_L3->SetLineColor(3);
  hh_bkg_L4->SetName("hh_bkg_L4"); hh_bkg_L4->SetLineColor(4);
  hh_bkg_L5->SetName("hh_bkg_L5"); hh_bkg_L5->SetLineColor(5);

  
  TH1D* hh_ratio_L2_nom = (TH1D*)hh_bkg_L2->Clone("ratio_L2_nom");
  TH1D* hh_ratio_L3_nom = (TH1D*)hh_bkg_L3->Clone("ratio_L3_nom");
  TH1D* hh_ratio_L4_nom = (TH1D*)hh_bkg_L4->Clone("ratio_L4_nom");
  TH1D* hh_ratio_L5_nom = (TH1D*)hh_bkg_L5->Clone("ratio_L5_nom");
  TH1D* hh_ratio_nom_nom = (TH1D*)m_hh_bkg->Clone("ratio_nom_nom");
  hh_ratio_L2_nom->Divide(m_hh_bkg);
  hh_ratio_L3_nom->Divide(m_hh_bkg);
  hh_ratio_L4_nom->Divide(m_hh_bkg);
  hh_ratio_L5_nom->Divide(m_hh_bkg);
  hh_ratio_nom_nom->Divide(m_hh_bkg);
  hh_ratio_L2_nom->SetLineColor(2);
  hh_ratio_L3_nom->SetLineColor(3);
  hh_ratio_L4_nom->SetLineColor(4);
  hh_ratio_L5_nom->SetLineColor(5);

  m_hh_frac_syst = (TH1D*)hh_ratio_L2_nom->Clone("frac_envelop_syst");
  for(int ibin=1;ibin<=m_hh_frac_syst->GetNbinsX();ibin++){
    double BC = m_hh_frac_syst->GetBinContent(ibin);
    m_hh_frac_syst->SetBinError(ibin,0);
    m_hh_frac_syst->SetBinContent(ibin,fabs(BC-1));
    if( fabs(m_hh_frac_syst->GetBinContent(ibin)-1)< 
	fabs(hh_ratio_L3_nom->GetBinContent(ibin)-1))
      m_hh_frac_syst->SetBinContent(ibin,fabs(hh_ratio_L3_nom->GetBinContent(ibin)-1) ) ;
    if( fabs(m_hh_frac_syst->GetBinContent(ibin)-1)< 
	fabs(hh_ratio_L5_nom->GetBinContent(ibin)-1))
      m_hh_frac_syst->SetBinContent(ibin,fabs(hh_ratio_L5_nom->GetBinContent(ibin)-1) );
  }
  if(doplot){
    TLegend * leg = new TLegend(0.6,0.7,0.8,0.9);
    leg->SetFillColor(0);
    leg->AddEntry(hh_ratio_L2_nom,"Loose' 2","l");
    leg->AddEntry(hh_ratio_L3_nom,"Loose' 3","l");
    leg->AddEntry(hh_ratio_nom_nom,"Nominal","l");
    leg->AddEntry(hh_ratio_L5_nom,"Loose' 5","l");
    leg->Draw("same");
    TCanvas * c = new TCanvas("c","c",800,600);
    c->cd();
    hh_bkg_L2->GetXaxis()->SetRange(Commons::norm_bin.first,Commons::nBins);
    hh_bkg_L2->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    hh_bkg_L2->Draw("HIST");
    hh_bkg_L3->Draw("sameHIST");
    // hh_bkg_L4->Draw("sameHIST");
    hh_bkg_L5->Draw("sameHIST");
    m_hh_bkg->Draw("sameHIST");
    leg->Draw("same");
    TCanvas* c1 = new TCanvas("c1","c1",800,600);
    c1->cd();
    hh_ratio_L2_nom->GetXaxis()->SetRange(Commons::norm_bin.first,Commons::nBins);
    hh_ratio_L2_nom->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    hh_ratio_L2_nom->GetYaxis()->SetRangeUser(0.7,1.3);
    hh_ratio_L2_nom->Draw("HIST][");
    hh_ratio_L3_nom->Draw("sameHIST][");
    hh_ratio_L4_nom->Draw("sameHIST][");
    hh_ratio_L5_nom->Draw("sameHIST][");
    hh_ratio_nom_nom->Draw("sameHIST][");
    leg->Draw("same");
  }

}
////////////////////////////////////////////////
void BkgEstimator::YieldsStatUncert(bool doplot)
///////////////////////////////////////////////
{
  double Yield_nom[] = { Commons::yields[0],
			 Commons::yields[1],
			 Commons::yields[2],
			 Commons::yields[3] };

  double Yield_err[] = { Commons::yields_er[0],
			 Commons::yields_er[1],
			 Commons::yields_er[2],
			 Commons::yields_er[3] };


  int nentries = (int)m_yieldstree_L4->GetEntriesFast();
  std::vector<TH1D*> h_bkg_err;
  std::vector<TH1D*> h_ratio_err_nom;
  TreeReader Rd(m_yieldstree_L4);

  for( int ipe = 0 ; ipe<nentries;ipe++){
    Rd.GetEntry(ipe);
    std::vector<double> Yield;
    Yield.push_back((int)Rd.GetVariable("Ngg"));
    Yield.push_back((int)Rd.GetVariable("Ngj"));
    Yield.push_back((int)Rd.GetVariable("Njg"));
    Yield.push_back((int)Rd.GetVariable("Njj"));
    bool IsIn1sigma = true;
    for( int i= 0;i<4;i++){
      if( Yield[i]>Yield_nom[i]+Yield_err[i] ||
	  Yield[i]<Yield_nom[i]-Yield_err[i] ){
	IsIn1sigma = false;
	break;
      }
    }
    
    if(!IsIn1sigma) continue;
 
    TH1D* h_bkg_temp = GetTotalBkgHist(m_hh_gg,m_hh_gj,m_hh_jg,m_hh_jj,Yield);
    h_bkg_temp->SetName( Form("h_bkg_err%d",ipe) );
    h_bkg_temp->SetTitle( Form("h_bkg_err%d",ipe) );
    h_bkg_temp->SetLineColor(2);


    TH1D* h_ratio_temp_nom = (TH1D*)h_bkg_temp->Clone( Form("h_ratio_temp_nom%d",ipe) );
    h_ratio_temp_nom->Divide(m_hh_bkg);
    h_ratio_temp_nom->SetLineColor(2);
    h_bkg_err.push_back(h_bkg_temp);
    h_ratio_err_nom.push_back(h_ratio_temp_nom);
  }

  m_hh_frac_stat = (TH1D*)h_ratio_err_nom[0]->Clone("frac_envelop_stat");
  for(int ibin=1; ibin <=m_hh_frac_stat->GetNbinsX();ibin++){
    double BC = m_hh_frac_stat->GetBinContent(ibin);
    m_hh_frac_stat->SetBinError(ibin,0);
    m_hh_frac_stat->SetBinContent(ibin,fabs(BC-1));
      for(int ipe = 1; ipe< (int)h_ratio_err_nom.size();ipe++){
	double BC_new = h_ratio_err_nom[ipe]->GetBinContent(ibin);
	if(fabs(BC_new-1)> fabs(m_hh_frac_stat->GetBinContent(ibin)-1))
	  m_hh_frac_stat->SetBinContent(ibin,fabs(BC_new-1));	 
      }
  }
  
  if(doplot){
    TH1D* hh_ratio_nom_nom = (TH1D*)m_hh_bkg->Clone("ratio_nom_nom");
    hh_ratio_nom_nom->Divide(m_hh_bkg);
    TCanvas* c1 = new TCanvas("cratio_stat","ratio stat",800,600);
    c1->cd();
    hh_ratio_nom_nom->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    hh_ratio_nom_nom->GetXaxis()->SetRange(Commons::norm_bin.first,Commons::nBins);
    hh_ratio_nom_nom->GetYaxis()->SetRangeUser(0.98,1.02);
    hh_ratio_nom_nom->Draw("HIST][");
    for( int ipe = 0 ; ipe<(int)h_ratio_err_nom.size();ipe++){
      h_ratio_err_nom[ipe]->Draw("sameHIST][");
    }
    hh_ratio_nom_nom->Draw("sameHIST][");
    TLegend* leg = new TLegend(0.2,0.8,0.6,0.9);
    leg->SetFillColor(0);
    leg->AddEntry(h_ratio_err_nom[0],"Yields randomized","l");
    leg->Draw("same");
  }
}
//////////////////////////////////////////////
void BkgEstimator::RedShapeUncert(bool doplot)
//////////////////////////////////////////////
{
  m_hh_gj_syst = GetRedShapeSyst("gj",doplot);
  m_hh_jg_syst = GetRedShapeSyst("jg",doplot);
  m_hh_jj_syst = GetRedShapeSyst("jj",doplot);
  //----------------------------------------------
  m_hh_red_syst = (TH1D*)m_hh_gj_syst->Clone("hh_red_syst");
  m_hh_red_syst->SetTitle("hh_red_syst");

  for( int ibin=1;ibin<=m_hh_gj_syst->GetNbinsX();ibin++ ){
    double gj_BC = m_hh_gj_syst->GetBinContent(ibin);
    double jg_BC = m_hh_jg_syst->GetBinContent(ibin);
    double jj_BC = m_hh_jj_syst->GetBinContent(ibin);
    m_hh_red_syst->SetBinContent(ibin, sqrt(gj_BC*gj_BC+jg_BC*jg_BC+jj_BC*jj_BC)); 
  }

  if(doplot){
    TCanvas *c = new TCanvas();
    c->cd();
    m_hh_red_syst->Draw("HIST");
    m_hh_gj_syst->Draw("sameHIST");
    m_hh_jg_syst->Draw("sameHIST");
    m_hh_jj_syst->Draw("sameHIST");
  }

}
////////////////////////////////////////////////////////////////
TH1D* BkgEstimator::GetRedShapeSyst(TString bkgtype,bool doplot)
////////////////////////////////////////////////////////////////
{
  TFile *f4 = new TFile(name_f_red_syst,"read");

  TString histname;
  if(bkgtype=="gj") histname = "h_subleading_fake";
  else if(bkgtype=="jg") histname = "h_leading_fake";
  else if(bkgtype=="jj") histname = "h_double_fake";
  else std::cout << "WRONG BKG TYPE !!!!" << std::endl;


  TH1D* h4 = (TH1D*)f4->Get(histname);
  // TH1D* h5 = (TH1D*)f5->Get(histname);
  TH1D* hsys;
  if(bkgtype=="gj") hsys = GetTotalBkgHist(m_hh_gg,h4,m_hh_jg,m_hh_jj,Commons::yields);
  else if(bkgtype=="jg") hsys = GetTotalBkgHist(m_hh_gg,m_hh_gj,h4,m_hh_jj,Commons::yields);
  else if(bkgtype=="jj") hsys = GetTotalBkgHist(m_hh_gg,m_hh_gj,m_hh_jg,h4,Commons::yields);


  BkgUncertaintyBuilder *BUB;
  BUB = new BkgUncertaintyBuilder(*m_hh_bkg,*hsys);

  TH1D* htemp = BUB->GetBkgFracHist();
  htemp->SetName("frachist_"+bkgtype);
  htemp->SetTitle("frachist_"+bkgtype);

  if( doplot){
    TCanvas *c = new TCanvas("c_"+bkgtype,"c_"+bkgtype,800,600);
    c->cd();
    htemp->Draw("HIST");
  }
  return htemp;



}

/////////////////////////////////////////////
void BkgEstimator::CreateBATInput(TString st)
/////////////////////////////////////////////
{
  m_hh_bkg_400->SetName("bkg_total_gg");
  m_hh_bkg_400->SetTitle("bkg_total_gg");

  TObjArray Hlist(0);
  Hlist.Add(m_hh_data_400);
  Hlist.Add(m_hh_bkg_400);
  Hlist.Add(m_hh_red_400);
  Hlist.Add(m_hh_irr_400);
  Hlist.Add(m_hh_syst_400);
  Hlist.Add(m_hh_red_syst_400);
  Hlist.Add(m_hh_irr_syst_400);
  Hlist.Add(m_hh_frac_syst_400);
  TFile fout(st,"RECREATE");
  Hlist.Write();
  fout.Close();
}
////////////////////////////////////////////////////
void BkgEstimator::CreateBumpHunterInput(TString st)
////////////////////////////////////////////////////
{
  TObjArray Hlist1(0);

  Hlist1.Add(m_hh_data);
  Hlist1.Add(m_hh_bkg);
  Hlist1.Add(m_hh_red);
  Hlist1.Add(m_hh_gg);
  Hlist1.Add(m_hh_gj);
  Hlist1.Add(m_hh_jg);
  Hlist1.Add(m_hh_jj);
  Hlist1.Add(m_hh_irr_syst);
  Hlist1.Add(m_hh_red_syst);
  Hlist1.Add(m_hh_gj_syst);
  Hlist1.Add(m_hh_jg_syst);
  Hlist1.Add(m_hh_jj_syst);
  Hlist1.Add(m_hh_frac_syst);
  Hlist1.Add(m_hh_frac_stat);
  Hlist1.Add(m_hh_frac_tot);
  Hlist1.Add(m_hh_bkg_syst);
  Hlist1.Add(m_gr_gj);
  Hlist1.Add(m_gr_jg);
  Hlist1.Add(m_gr_jj);
  Hlist1.Add(m_gmgg_data);
  Hlist1.Add(m_gmgg_bkg);
  Hlist1.Add(m_gmgg_red);

  TFile fout1(st,"RECREATE");
  Hlist1.Write();
  fout1.Close();
}
/////////////////////////////////////////////////////////////////////////////
TH1D* BkgEstimator::GetTotalBkgHist( TH1D* hgg,TH1D* hgj,TH1D* hjg,TH1D* hjj,
				     const std::vector<double> yields)
//////////////////////////////////////////////////////////////////////
{
  if( yields.size() !=4 )
    Fatal( "BkgEstimator::GetTotalBkgHist",
	   "Invalid size of the yields vector !!" );

  double sumyields = 0; std::vector<double> frac;
  for(int iy=0;iy<(int)yields.size();iy++) sumyields += yields[iy];
  for(int iy=0;iy<(int)yields.size();iy++) frac.push_back(yields[iy]/sumyields);

  hgg->Scale(frac[0]*dataNorm/hgg->Integral(Commons::norm_bin.first,Commons::norm_bin.second));
  hgj->Scale(frac[1]*dataNorm/hgj->Integral(Commons::norm_bin.first,Commons::norm_bin.second));
  hjg->Scale(frac[2]*dataNorm/hjg->Integral(Commons::norm_bin.first,Commons::norm_bin.second));
  hjj->Scale(frac[3]*dataNorm/hjj->Integral(Commons::norm_bin.first,Commons::norm_bin.second));

  TH1D* htemp = (TH1D*)hgg->Clone("hh_totbkg");
  htemp->Add(hgj);htemp->Add(hjg);htemp->Add(hjj);
  return htemp;
}
///////////////////////////////////////////////////////////////////////
TH1D* BkgEstimator::GetRedBkgHist( TH1D* hgj,TH1D* hjg,TH1D* hjj,
				   const std::vector<double> yields)
//////////////////////////////////////////////////////////////////////
{
  if( yields.size() !=3 )
    Fatal( "BkgEstimator::GetRedBkgHist",
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
///////////////////////////////////////////////////////////////////////////////////////////
std::pair<double,double> BkgEstimator::GetYield(TString type,TTree* yieldstree,bool doplot)
//////////////////////////////////////////////////////////////////////////////////////////
{
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
  RooMsgService::instance().setSilentMode(!m_VERBOSE);
  gauss.fitTo(data);

  if(doplot){
    RooPlot *frame = N.frame();
    data.plotOn(frame,RooFit::Name("data"));
    gauss.plotOn(frame,RooFit::Name("gauss"),RooFit::LineColor(2));
    // TCanvas *c = new TCanvas();
    // c->cd();
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
