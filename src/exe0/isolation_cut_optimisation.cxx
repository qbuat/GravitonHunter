#include <iostream>
#include <TSystem.h>
#include <TError.h>
#include <map>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TMath.h>
#include <TAxis.h>
#include <TH1.h>

#include "ToolsCommons.h"
#include "ToolsTreeReader.h"
#include "ToolsChainMaker.h"

using std::cerr;

std::pair<double,double> GetIsoCut(double efficiency,double ptmin,double ptmax,TString L_or_SL="L", double isocut = 5.);
int main(int argc, char ** argv){
  

  //------------- PROTECTION AGAINST WRONG NUMBER OF PARAMETERS ------------
  if( argc < 2 ){
    std::cerr << "Wrong usage ! "
	      << argv[0]
	      << " outfile_label isocut\n" ;
    exit(1);
  }
  //------------------------------------------------------------------------

  //------------- INPUT PARAMETERS VALUES ------------------------------------
  TString outfile_label = argv[1]; //label used to name the file
  double isocut = atof(argv[2]);
  //--------------------------------------------------------------------------
  std::cerr << "Processing : " << argv[0] << " " 
	    << outfile_label << "\n";

  Commons::Setup();

  //--> Initialize the file used for this study (flat rs gg file)
  TFile* file = TFile::Open("/sps/atlas/q/qbuat/OUTDS/GravitonAnalysis8TeV/sig_mc12_Ggg_flat_TOPO.root");
  TTree* tree = (TTree*)file->Get("tree");
  int nentries = tree->GetEntries();
  TreeReader Rd(tree);

  const int Nptbins = 17;
  double ptbins[Nptbins+1] = {50,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1700,2000};
  double Ntight_L[Nptbins];
  double Nisotight_L[Nptbins];
  double Ntight_SL[Nptbins];
  double Nisotight_SL[Nptbins];
  for( int i=0;i<Nptbins;i++){
    Ntight_L[i]=0;    
    Nisotight_L[i]=0; 
    Ntight_SL[i]=0;   
    Nisotight_SL[i]=0;
  }

  for( int entry=0;entry<nentries;entry++){//entries loop
    Rd.GetEntry(entry);
    if( (int)Rd.GetVariable("IsTight_L") == 0 ) continue;
    if( (int)Rd.GetVariable("IsTight_SL") == 0 ) continue;
    double weight = Rd.GetVariable("weight");

    for( int i=0;i<Nptbins;i++){
      if( Rd.GetVariable("pT_L") >= ptbins[i] && Rd.GetVariable("pT_L")<ptbins[i+1])
	Ntight_L[i] += weight;
      if( Rd.GetVariable("pT_SL") >= ptbins[i] && Rd.GetVariable("pT_SL")<ptbins[i+1])
	Ntight_SL[i] += weight;
    }

    //    if( Rd.GetVariable("Iso_L") - Commons::GetTopoIsoPtcorr(Rd.GetVariable("pT_L") ) < isocut ){
    if( Rd.GetVariable("Iso_L") < isocut ){
      for( int i=0;i<Nptbins;i++){
	if( Rd.GetVariable("pT_L") >= ptbins[i] && Rd.GetVariable("pT_L")<ptbins[i+1])
	  Nisotight_L[i] += weight;
      }
    } 

    //    if( Rd.GetVariable("Iso_SL") - Commons::GetTopoIsoPtcorr(Rd.GetVariable("pT_SL") ) < isocut ){
    if( Rd.GetVariable("Iso_SL") < isocut ){
      for( int i=0;i<Nptbins;i++){
	if( Rd.GetVariable("pT_SL") >= ptbins[i] && Rd.GetVariable("pT_SL")<ptbins[i+1])
	  Nisotight_SL[i] += weight;
      }
    } 

  }//entries loop


  TGraphErrors *gr_effl = new TGraphErrors();
  TGraphErrors *gr_effsl = new TGraphErrors();
  for(int i=0 ;i < Nptbins ; i++){
    double eff_L = Nisotight_L[i]/Ntight_L[i];
    double err_eff_L = sqrt(eff_L*(1-eff_L)/Ntight_L[i]);
    gr_effl->SetPoint(i,0.5*(ptbins[i+1]+ptbins[i]),eff_L);
    gr_effl->SetPointError(i,0.5*(ptbins[i+1]-ptbins[i]),err_eff_L);

    double eff_SL = Nisotight_SL[i]/Ntight_SL[i];
    double err_eff_SL = sqrt(eff_SL*(1-eff_SL)/Ntight_SL[i]);
    gr_effsl->SetPoint(i,0.5*(ptbins[i+1]+ptbins[i]),eff_SL);
    gr_effsl->SetPointError(i,0.5*(ptbins[i+1]-ptbins[i]),err_eff_SL);
  }

  double eff_targeted_value_l = Nisotight_L[0]/Ntight_L[0];//0.96;
  double eff_targeted_value_sl = Nisotight_SL[0]/Ntight_SL[0];//0.91;

  std::cout<<"Targeted values: "<<std::endl;
  std::cout<<"Leading photon: "<<eff_targeted_value_l<<std::endl;
  std::cout<<"Subleading photon: "<<eff_targeted_value_sl<<std::endl;


  TGraphErrors* gr_isocut_l  = new TGraphErrors();
  TGraphErrors* gr_isocut_sl = new TGraphErrors();
  for( int i=0; i<Nptbins; i++){
    std::pair<double,double> a = GetIsoCut(eff_targeted_value_l,ptbins[i],ptbins[i+1],"L", isocut);
    std::cout << "IsoCut = " << a.first << "+/-" << a.second << std::endl;
    gr_isocut_l->SetPoint(i,0.5*(ptbins[i]+ptbins[i+1]),a.first);
    double precision_uncert_l = a.second;
    double x_l,y_l;
    gr_effl->GetPoint(i,x_l,y_l);
    double stateff_uncert_l = gr_effl->GetErrorY(i)/y_l*a.first;
    double error_l = sqrt( pow(precision_uncert_l,2) + pow(stateff_uncert_l,2));
    gr_isocut_l->SetPointError(i,0.5*(ptbins[i]-ptbins[i+1]),error_l);

    std::pair<double,double> b = GetIsoCut(eff_targeted_value_sl,ptbins[i],ptbins[i+1],"SL", isocut);
    std::cout << "IsoCut = " << b.first << "+/-" << b.second << std::endl;
    gr_isocut_sl->SetPoint(i,0.5*(ptbins[i]+ptbins[i+1]),b.first);
    double precision_uncert_sl = b.second;
    double x_sl,y_sl;
    gr_effl->GetPoint(i,x_sl,y_sl);
    double stateff_uncert_sl = gr_effsl->GetErrorY(i)/y_sl*b.first;
    double error_sl = sqrt( pow(precision_uncert_sl,2) + pow(stateff_uncert_sl,2));
    gr_isocut_sl->SetPointError(i,0.5*(ptbins[i]-ptbins[i+1]),error_sl);
  }

  gr_effl->SetName("efficiency_leading_photon");
  gr_effl->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
  gr_effl->GetYaxis()->SetTitle("Isolation Efficiency (#Iso+Tight/#Tight)");

  gr_effsl->SetName("efficiency_subleading_photon");
  gr_effsl->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
  gr_effsl->GetYaxis()->SetTitle("Isolation Efficiency (#Iso+Tight/#Tight)");

  gr_isocut_l->SetName("isolation_cut_leading_photon");
  gr_isocut_l->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
  gr_isocut_l->GetYaxis()->SetTitle("Isolation Cut [GeV]");

  gr_isocut_sl->SetName("isolation_cut_subleading_photon");
  gr_isocut_sl->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
  gr_isocut_sl->GetYaxis()->SetTitle("Isolation Cut [GeV]");


  TH1F* htargeted_effs = new TH1F("htargeted_effs","htargeted_effs",2,0,2);
  htargeted_effs->SetBinContent(1,eff_targeted_value_l);
  htargeted_effs->SetBinContent(2,eff_targeted_value_sl);

  TString filest = "./rootfiles/IsolationCutOptimisation"+outfile_label+".root";
  TFile fout(filest,"RECREATE");
  fout.Add(gr_effl);
  fout.Add(gr_effsl);
  fout.Add(gr_isocut_l);
  fout.Add(gr_isocut_sl);
  fout.Add(htargeted_effs);
  fout.Write();
  fout.Close();


}


std::pair<double,double> GetIsoCut(double efficiency,double ptmin,double ptmax,TString L_or_SL, double isocut)
{

  TFile* file = TFile::Open("/sps/atlas/q/qbuat/OUTDS/GravitonAnalysis8TeV/sig_mc12_Ggg_flat_TOPO.root");

  TTree* tree = (TTree*)file->Get("tree");
  int nentries = tree->GetEntries();
  TreeReader Rd(tree);

  TString st_pt;
  TString st_iso;
  if( L_or_SL == "L" ){
    st_pt = "pT_L" ;
    st_iso = "Iso_L" ;
  } else if( L_or_SL = "SL" ){
    st_pt = "pT_SL";
    st_iso = "Iso_SL";
  }else Fatal("GetIsoCut","Wrong L_or_SL choice");
      
  double isocut_test = isocut;
  while ( true ){
    std::cout << "test the value " << isocut_test ;
    double Ntight    = 0;
    double Nisotight = 0;
    for( int entry=0;entry<nentries;entry++){//entries loop
      Rd.GetEntry(entry);
      if( (int)Rd.GetVariable("IsTight_L") == 0 ) continue;
      if( (int)Rd.GetVariable("IsTight_SL") == 0 ) continue;
      double weight = Rd.GetVariable("weight");
      if( Rd.GetVariable(st_pt) < ptmin ) continue;
      if( Rd.GetVariable(st_pt)>= ptmax ) continue;
      Ntight += weight;
      if( Rd.GetVariable(st_iso) <isocut_test ) Nisotight += weight;
    }					
    double eff = Nisotight/Ntight;
    std::cout << ", the efficiency is " << eff << std::endl;
    if( eff < efficiency ) isocut_test += 0.1;
    else{
      // break;
      double err_iso = isocut_test*fabs(eff-efficiency)/efficiency;
      std::pair<double,double> isocut_witherror(isocut_test,err_iso);
      return isocut_witherror;
    }
  }
}




