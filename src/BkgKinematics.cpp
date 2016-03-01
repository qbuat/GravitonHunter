#include "BkgKinematics.h"
#include <TError.h>
#include <TLegend.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h> 
#include <RooAddPdf.h> 
#include <RooFitResult.h>
#include <RooPlot.h>  

#include "ToolsTreeReader.h"
#include "ToolsUtilities.h"
#include "ToolsCommons.h"
#include "ToolsExtendedCanvas.h"
#include "ToolsSignificanceHist.h"

/////////////////////////////////////////////////////////////
BkgKinematics::BkgKinematics(TString datafile,TString mcfile)
/////////////////////////////////////////////////////////////
{
  SetMCFile(mcfile);
  SetDataFile(datafile);
  SetEtaCategory("NONE");
  SetIsoCut(5);
  SetPtCuts(50,50);
  InitMaps();
  std::cout << "Enter FillDataHists()" << std::endl;
  FillDataHists();
  std::cout << "Enter FillMCHists()" << std::endl;
  FillMCHists();
}
////////////////////////////////////////////////////////////////////////////
BkgKinematics::BkgKinematics(TString datafile,TString mcfile,TString etacat)
///////////////////////////////////////////////////////////////////////////
{
  SetMCFile(mcfile);
  SetDataFile(datafile);
  SetEtaCategory(etacat);
  SetPurityModifier(0.);
  InitMaps();
  SetPtCuts(50,50);
  SetIsoCut(5);
  std::cout << "Enter FillDataHists()" << std::endl;
  FillDataHists();
  std::cout << "Enter FillMCHists()" << std::endl;
  FillMCHists();
  std::map<TString,TString>::iterator it;
  for ( it=m_hname.begin() ; it != m_hname.end(); it++ )
    BuildTotalBkgHist( (*it).second );
}
//////////////////////////////////////////////
BkgKinematics::BkgKinematics(TString histsfile)
//////////////////////////////////////////////
{ InitFromHistsFile(histsfile); }


/////////////////////////////////////
BkgKinematics::~BkgKinematics()
////////////////////////////////////
{  delete m_hmgg; }
//////////////////////////////
void BkgKinematics::InitMaps()
//////////////////////////////
{
  m_hname["pT"] = "pT";
  m_hname["pT_L"] = "pT_L";
  m_hname["pT_SL"] = "pT_SL";
  m_hname["eta_L"] = "eta_L";
  m_hname["eta_SL"] = "eta_SL";
  m_hname["mgg"] = "mgg";
  m_hname["ptgg"] = "ptgg";
  m_hname["costhetastar"] = "costhetastar";
  m_hname["deltaphi"] = "deltaphi";

  m_title["pT"] = "p_{T}^{#gamma} [GeV]";
  m_title["pT_L"] = "p_{T}^{#gamma,leading} [GeV]";
  m_title["pT_SL"] = "p_{T}^{#gamma,subleading} [GeV]";
  m_title["eta_L"] = "#eta^{#gamma,leading}";
  m_title["eta_SL"] = "#eta^{#gamma,subleading}";
  m_title["mgg"] = "m_{#gamma#gamma} [GeV]";
  m_title["ptgg"] = "p_{T,#gamma#gamma} [GeV]";
  m_title["costhetastar"] = "cos(#theta*)";
  m_title["deltaphi"] = "#Delta #phi_{#gamma#gamma}";

  m_Nbins["pT"] = 110;
  m_Nbins["pT_L"] = 110;
  m_Nbins["pT_SL"] = 110;
  m_Nbins["eta_L"] = 60;
  m_Nbins["eta_SL"] = 60;
  m_Nbins["mgg"] = 130;
  m_Nbins["ptgg"] = 50;
  m_Nbins["costhetastar"] = 200;

  m_Xmin["pT"] = 30;  m_Xmax["pT"] = 260;
  m_Xmin["pT_L"] = 40;  m_Xmax["pT_L"] = 260;
  m_Xmin["pT_SL"] = 30;  m_Xmax["pT_SL"] = 250;
  m_Xmin["eta_L"] = -3;  m_Xmax["eta_L"] = 3;
  m_Xmin["eta_SL"] = -3;  m_Xmax["eta_SL"] = 3;
  // m_Xmin["mgg"] = 140;  m_Xmax["mgg"] = 400;
  m_Xmin["mgg"] = Commons::binning[Commons::norm_bin.first];  
  m_Xmax["mgg"] = Commons::binning[Commons::norm_bin.second];
  m_Xmin["ptgg"] = 0; m_Xmax["ptgg"] = 600;
  m_Xmin["costhetastar"] = -1;  m_Xmax["costhetastar"] = 1;

  m_hD["pT"] = new TH1D(m_hname["pT"]+"D","",m_Nbins["pT"],m_Xmin["pT"],m_Xmax["pT"]);
  m_hD["pT_L"] = new TH1D(m_hname["pT_L"]+"D","",m_Nbins["pT_L"],m_Xmin["pT_L"],m_Xmax["pT_L"]);
  m_hD["pT_SL"] = new TH1D(m_hname["pT_SL"]+"D","",m_Nbins["pT_SL"],m_Xmin["pT_SL"],m_Xmax["pT_SL"]);
  m_hD["eta_L"] = new TH1D(m_hname["eta_L"]+"D","",m_Nbins["eta_L"],m_Xmin["eta_L"],m_Xmax["eta_L"]);
  m_hD["eta_SL"] = new TH1D(m_hname["eta_SL"]+"D","",m_Nbins["eta_SL"],m_Xmin["eta_SL"],m_Xmax["eta_SL"]);
  m_hD["mgg"] = new TH1D(m_hname["mgg"]+"D","",m_Nbins["mgg"],m_Xmin["mgg"],m_Xmax["mgg"]);
  m_hD["ptgg"] = new TH1D(m_hname["ptgg"]+"D","",m_Nbins["ptgg"],m_Xmin["ptgg"],m_Xmax["ptgg"]);
  m_hD["costhetastar"] = new TH1D(m_hname["costhetastar"]+"D","",m_Nbins["costhetastar"],m_Xmin["costhetastar"],m_Xmax["costhetastar"]);
  //--------------------------------------------------------------------------------------------------------------------------------------
  m_hB["pT"] = new TH1D(m_hname["pT"]+"B","",m_Nbins["pT"],m_Xmin["pT"],m_Xmax["pT"]);
  m_hB["pT_L"] = new TH1D(m_hname["pT_L"]+"B","",m_Nbins["pT_L"],m_Xmin["pT_L"],m_Xmax["pT_L"]);
  m_hB["pT_SL"] = new TH1D(m_hname["pT_SL"]+"B","",m_Nbins["pT_SL"],m_Xmin["pT_SL"],m_Xmax["pT_SL"]);
  m_hB["eta_L"] = new TH1D(m_hname["eta_L"]+"B","",m_Nbins["eta_L"],m_Xmin["eta_L"],m_Xmax["eta_L"]);
  m_hB["eta_SL"] = new TH1D(m_hname["eta_SL"]+"B","",m_Nbins["eta_SL"],m_Xmin["eta_SL"],m_Xmax["eta_SL"]);
  m_hB["mgg"] = new TH1D(m_hname["mgg"]+"B","",m_Nbins["mgg"],m_Xmin["mgg"],m_Xmax["mgg"]);
  m_hB["ptgg"] = new TH1D(m_hname["ptgg"]+"B","",m_Nbins["ptgg"],m_Xmin["ptgg"],m_Xmax["ptgg"]);
  m_hB["costhetastar"] = new TH1D(m_hname["costhetastar"]+"B","",m_Nbins["costhetastar"],m_Xmin["costhetastar"],m_Xmax["costhetastar"]);
  //--------------------------------------------------------------------------------------------------------------------------------------
  m_hI["pT"] = new TH1D(m_hname["pT"]+"I","",m_Nbins["pT"],m_Xmin["pT"],m_Xmax["pT"]);
  m_hI["pT_L"] = new TH1D(m_hname["pT_L"]+"I","",m_Nbins["pT_L"],m_Xmin["pT_L"],m_Xmax["pT_L"]);
  m_hI["pT_SL"] = new TH1D(m_hname["pT_SL"]+"I","",m_Nbins["pT_SL"],m_Xmin["pT_SL"],m_Xmax["pT_SL"]);
  m_hI["eta_L"] = new TH1D(m_hname["eta_L"]+"I","",m_Nbins["eta_L"],m_Xmin["eta_L"],m_Xmax["eta_L"]);
  m_hI["eta_SL"] = new TH1D(m_hname["eta_SL"]+"I","",m_Nbins["eta_SL"],m_Xmin["eta_SL"],m_Xmax["eta_SL"]);
  m_hI["mgg"] = new TH1D(m_hname["mgg"]+"I","",m_Nbins["mgg"],m_Xmin["mgg"],m_Xmax["mgg"]);
  m_hI["ptgg"] = new TH1D(m_hname["ptgg"]+"I","",m_Nbins["ptgg"],m_Xmin["ptgg"],m_Xmax["ptgg"]);
  m_hI["costhetastar"] = new TH1D(m_hname["costhetastar"]+"I","",m_Nbins["costhetastar"],m_Xmin["costhetastar"],m_Xmax["costhetastar"]);
  //---------------------------------------------------------------------------------------------------------------------------------------
  m_hGJ["pT"] = new TH1D(m_hname["pT"]+"GJ","",m_Nbins["pT"],m_Xmin["pT"],m_Xmax["pT"]);
  m_hGJ["pT_L"] = new TH1D(m_hname["pT_L"]+"GJ","",m_Nbins["pT_L"],m_Xmin["pT_L"],m_Xmax["pT_L"]);
  m_hGJ["pT_SL"] = new TH1D(m_hname["pT_SL"]+"GJ","",m_Nbins["pT_SL"],m_Xmin["pT_SL"],m_Xmax["pT_SL"]);
  m_hGJ["eta_L"] = new TH1D(m_hname["eta_L"]+"GJ","",m_Nbins["eta_L"],m_Xmin["eta_L"],m_Xmax["eta_L"]);
  m_hGJ["eta_SL"] = new TH1D(m_hname["eta_SL"]+"GJ","",m_Nbins["eta_SL"],m_Xmin["eta_SL"],m_Xmax["eta_SL"]);
  m_hGJ["mgg"] = new TH1D(m_hname["mgg"]+"GJ","",m_Nbins["mgg"],m_Xmin["mgg"],m_Xmax["mgg"]);
  m_hGJ["ptgg"] = new TH1D(m_hname["ptgg"]+"GJ","",m_Nbins["ptgg"],m_Xmin["ptgg"],m_Xmax["ptgg"]);
  m_hGJ["costhetastar"] = new TH1D(m_hname["costhetastar"]+"GJ","",m_Nbins["costhetastar"],m_Xmin["costhetastar"],m_Xmax["costhetastar"]);
  //---------------------------------------------------------------------------------------------------------------------------------------
  m_hJG["pT"] = new TH1D(m_hname["pT"]+"JG","",m_Nbins["pT"],m_Xmin["pT"],m_Xmax["pT"]);
  m_hJG["pT_L"] = new TH1D(m_hname["pT_L"]+"JG","",m_Nbins["pT_L"],m_Xmin["pT_L"],m_Xmax["pT_L"]);
  m_hJG["pT_SL"] = new TH1D(m_hname["pT_SL"]+"JG","",m_Nbins["pT_SL"],m_Xmin["pT_SL"],m_Xmax["pT_SL"]);
  m_hJG["eta_L"] = new TH1D(m_hname["eta_L"]+"JG","",m_Nbins["eta_L"],m_Xmin["eta_L"],m_Xmax["eta_L"]);
  m_hJG["eta_SL"] = new TH1D(m_hname["eta_SL"]+"JG","",m_Nbins["eta_SL"],m_Xmin["eta_SL"],m_Xmax["eta_SL"]);
  m_hJG["mgg"] = new TH1D(m_hname["mgg"]+"JG","",m_Nbins["mgg"],m_Xmin["mgg"],m_Xmax["mgg"]);
  m_hJG["ptgg"] = new TH1D(m_hname["ptgg"]+"JG","",m_Nbins["ptgg"],m_Xmin["ptgg"],m_Xmax["ptgg"]);
  m_hJG["costhetastar"] = new TH1D(m_hname["costhetastar"]+"JG","",m_Nbins["costhetastar"],m_Xmin["costhetastar"],m_Xmax["costhetastar"]);
  //---------------------------------------------------------------------------------------------------------------------------------------
  m_hJJ["pT"] = new TH1D(m_hname["pT"]+"JJ","",m_Nbins["pT"],m_Xmin["pT"],m_Xmax["pT"]);
  m_hJJ["pT_L"] = new TH1D(m_hname["pT_L"]+"JJ","",m_Nbins["pT_L"],m_Xmin["pT_L"],m_Xmax["pT_L"]);
  m_hJJ["pT_SL"] = new TH1D(m_hname["pT_SL"]+"JJ","",m_Nbins["pT_SL"],m_Xmin["pT_SL"],m_Xmax["pT_SL"]);
  m_hJJ["eta_L"] = new TH1D(m_hname["eta_L"]+"JJ","",m_Nbins["eta_L"],m_Xmin["eta_L"],m_Xmax["eta_L"]);
  m_hJJ["eta_SL"] = new TH1D(m_hname["eta_SL"]+"JJ","",m_Nbins["eta_SL"],m_Xmin["eta_SL"],m_Xmax["eta_SL"]);
  m_hJJ["mgg"] = new TH1D(m_hname["mgg"]+"JJ","",m_Nbins["mgg"],m_Xmin["mgg"],m_Xmax["mgg"]);
  m_hJJ["ptgg"] = new TH1D(m_hname["ptgg"]+"JJ","",m_Nbins["ptgg"],m_Xmin["ptgg"],m_Xmax["ptgg"]);
  m_hJJ["costhetastar"] = new TH1D(m_hname["costhetastar"]+"JJ","",m_Nbins["costhetastar"],m_Xmin["costhetastar"],m_Xmax["costhetastar"]);
  m_hmgg = new TH1D("hmgg","hmgg",Commons::nBins,Commons::binning);

}


/////////////////////////////////////////////////////////
void BkgKinematics::InitFromHistsFile(TString histsfile)
////////////////////////////////////////////////////////
{
  TFile *f = new TFile(histsfile,"read");
  m_hname["pT"] = "pT";
  m_hname["pT_L"] = "pT_L";
  m_hname["pT_SL"] = "pT_SL";
  m_hname["eta_L"] = "eta_L";
  m_hname["eta_SL"] = "eta_SL";
  m_hname["mgg"] = "mgg";
  m_hname["ptgg"] = "ptgg";
  m_hname["costhetastar"] = "costhetastar";

  std::map<TString,TString>::iterator it;
  for ( it=m_hname.begin() ; it != m_hname.end(); it++ ){
    m_hD[(*it).first] = (TH1D*)f->Get( (*it).second+"D" );
    m_hB[(*it).first] = (TH1D*)f->Get( "hB"+(*it).second);
    m_hI[(*it).first] = (TH1D*)f->Get( (*it).second+"I" );
    m_hGJ[(*it).first] = (TH1D*)f->Get( (*it).second+"GJ" );
    m_hJG[(*it).first] = (TH1D*)f->Get( (*it).second+"JG" );
    m_hJJ[(*it).first] = (TH1D*)f->Get( (*it).second+"JJ" );
  }

}
////////////////////////////////////
void BkgKinematics::FillDataHists()
////////////////////////////////////
{
  std::cout << m_ptcut_l << "/" << m_ptcut_sl << std::endl;
  TFile * f = new TFile(m_filedata,"read");
  TChain * trD = (TChain*)f->Get("tree");
  TreeReader Rd_D(trD);
  int nentriesD = (int)trD->GetEntriesFast();
  //---------------------------------------------------
  for(int entry=0;entry<nentriesD;entry++){
    AnalysisTools::Processing(entry,nentriesD,(int)nentriesD/100);
    Rd_D.GetEntry(entry);

    double mgg = Rd_D.GetVariable("mgg");
    double eta_L = Rd_D.GetVariable("eta_L");
    double eta_SL = Rd_D.GetVariable("eta_SL");

    if( !Commons::EtaCategory(eta_L,eta_SL,m_etacategory) ) continue;
    if(mgg<Commons::binning[Commons::norm_bin.first]) continue;
    if(mgg>Commons::binning[Commons::norm_bin.second]) continue;
    if( Rd_D.GetVariable("pT_L")<m_ptcut_l) continue;
    if( Rd_D.GetVariable("pT_SL")<m_ptcut_sl) continue;
    if( Rd_D.GetVariable("Iso_L")>m_isocut) continue;
    if( Rd_D.GetVariable("Iso_SL")>m_isocut) continue;
    std::map<TString,TH1D*>::iterator it;

    if( (int)Rd_D.GetVariable("IsTight_L") ==1 && 
	(int)Rd_D.GetVariable("IsTight_SL")==1){
      for ( it=m_hD.begin() ; it != m_hD.end(); it++ )
	if( (*it).first != "pT" )
	  (*it).second->Fill( Rd_D.GetVariable((*it).first) );
	else{
	  (*it).second->Fill( Rd_D.GetVariable("pT_L") );
	  (*it).second->Fill( Rd_D.GetVariable("pT_SL") );
	}	  
      m_hmgg->Fill(mgg);
    }else if( (int)Rd_D.GetVariable("IsTight_L") ==0 && 
	      (int)Rd_D.GetVariable("IsTight_SL")==1)
      for ( it=m_hJG.begin() ; it != m_hJG.end(); it++ )
	if( (*it).first != "pT" )
	  (*it).second->Fill( Rd_D.GetVariable((*it).first) );
	else{
	  (*it).second->Fill( Rd_D.GetVariable("pT_L") );
	  (*it).second->Fill( Rd_D.GetVariable("pT_SL") );
	}	  
    else if( (int)Rd_D.GetVariable("IsTight_L") ==1 && 
	     (int)Rd_D.GetVariable("IsTight_SL")==0)
      for ( it=m_hGJ.begin() ; it != m_hGJ.end(); it++ )
	if( (*it).first != "pT" )
	  (*it).second->Fill( Rd_D.GetVariable((*it).first) );
	else{
	  (*it).second->Fill( Rd_D.GetVariable("pT_L") );
	  (*it).second->Fill( Rd_D.GetVariable("pT_SL") );
	}	  
    else if( (int)Rd_D.GetVariable("IsTight_L") ==0 && 
	     (int)Rd_D.GetVariable("IsTight_SL")==0)
      for ( it=m_hJJ.begin() ; it != m_hJJ.end(); it++ )
	if( (*it).first != "pT" )
	  (*it).second->Fill( Rd_D.GetVariable((*it).first) );
	else{
	  (*it).second->Fill( Rd_D.GetVariable("pT_L") );
	  (*it).second->Fill( Rd_D.GetVariable("pT_SL") );
	}	  

  }
  //------------------------------------------------------------

}
////////////////////////////////////
void BkgKinematics::FillMCHists()
////////////////////////////////////
{
  TFile* fbkg = new TFile(m_filemcgg,"read");
  TTree * trB = (TTree*)fbkg->Get("tree");
  TreeReader Rd_B(trB);
  int nentriesB = (int)trB->GetEntriesFast();
  //----------------------------------------------------
  for(int entry=0;entry<nentriesB;entry++){
    AnalysisTools::Processing(entry,nentriesB,10000);
    Rd_B.GetEntry(entry);

    double w_1 = Rd_B.GetVariable("weight");
    double w_g = Rd_B.GetVariable("gen_weight");
    double w_nlo = Commons::GetkFactor_40_30(Rd_B.GetVariable("truth_mgg"));
    double mgg = Rd_B.GetVariable("mgg");
    double eta_L = Rd_B.GetVariable("eta_L");
    double eta_SL = Rd_B.GetVariable("eta_SL");

    if( !Commons::EtaCategory(eta_L,eta_SL,m_etacategory) ) continue;
    if(mgg<Commons::binning[Commons::norm_bin.first]) continue;
    if(mgg>Commons::binning[Commons::norm_bin.second]) continue;
    if( (int)Rd_B.GetVariable("IsTight_L")!=1) continue;
    if( (int)Rd_B.GetVariable("IsTight_SL")!=1) continue;
    if( Rd_B.GetVariable("pT_L") < m_ptcut_l) continue;
    if( Rd_B.GetVariable("pT_SL")< m_ptcut_sl) continue;
    if( Rd_B.GetVariable("Iso_L")>m_isocut) continue;
    if( Rd_B.GetVariable("Iso_SL")>m_isocut) continue;

    std::map<TString,TH1D*>::iterator it;
    for ( it=m_hI.begin() ; it != m_hI.end(); it++ )
      if( (*it).first != "pT" )
	(*it).second->Fill( Rd_B.GetVariable((*it).first),w_1*w_g*w_nlo );
      else{
	(*it).second->Fill( Rd_B.GetVariable("pT_L"),w_1*w_g*w_nlo );
	(*it).second->Fill( Rd_B.GetVariable("pT_SL"),w_1*w_g*w_nlo );
	}	  

  }
  //--------------------------------------------------------------------------

}

///////////////////////////////////////////////////
void BkgKinematics::BuildTotalBkgHist(TString var)
//////////////////////////////////////////////////
{
  double sumyields = 0;
  for(int i=0;i<(int)Commons::yields.size();i++) sumyields+=Commons::yields[i]; 
  double dataNorm = m_hmgg->Integral(Commons::norm_bin.first,Commons::norm_bin.second);

  double scalepurgg = (1+m_purmod);
  double scalepurother = (sumyields-(1+m_purmod)*Commons::yields[0])/(sumyields-Commons::yields[0]);
  if( var != "pT"){
    m_hI[var]->Scale((scalepurgg)*Commons::yields[0]/sumyields*dataNorm/m_hI[var]->Integral());
    m_hGJ[var]->Scale((scalepurother)*Commons::yields[1]/sumyields*dataNorm/m_hGJ[var]->Integral());
    m_hJG[var]->Scale((scalepurother)*Commons::yields[2]/sumyields*dataNorm/m_hJG[var]->Integral());
    m_hJJ[var]->Scale((scalepurother)*Commons::yields[3]/sumyields*dataNorm/m_hJJ[var]->Integral());
  }else {
    m_hI[var]->Scale((scalepurgg)*2*Commons::yields[0]/sumyields*dataNorm/m_hI[var]->Integral());
    m_hGJ[var]->Scale((scalepurother)*2*Commons::yields[1]/sumyields*dataNorm/m_hGJ[var]->Integral());
    m_hJG[var]->Scale((scalepurother)*2*Commons::yields[2]/sumyields*dataNorm/m_hJG[var]->Integral());
    m_hJJ[var]->Scale((scalepurother)*2*Commons::yields[3]/sumyields*dataNorm/m_hJJ[var]->Integral());
  }

  m_hB[var] = (TH1D*)m_hI[var]->Clone("hB"+var);
  m_hB[var]->GetXaxis()->SetTitle(m_title[var]);
  m_hB[var]->Add(m_hGJ[var]);m_hB[var]->Add(m_hJG[var]);m_hB[var]->Add(m_hJJ[var]);
  m_hD[var]->GetXaxis()->SetTitle(m_title[var]);
  // m_hB[var]->GetXaxis()->SetRange(Commons::norm_bin.first,Commons::nBins);
  // m_hD[var]->GetXaxis()->SetRange(Commons::norm_bin.first,Commons::nBins);
  std::cout << "---- " << m_hname[var] << " ----" << std::endl;
  std::cout << "dataNorm = " << dataNorm  << std::endl;
  std::cout << "m_hD->Integral() =" << m_hD[var]->Integral() << std::endl;

}
///////////////////////////////////////////////
TCanvas* BkgKinematics::DrawHists(TString var)
//////////////////////////////////////////////
{


  m_hB[var]->SetLineColor(4);
  m_hI[var]->SetLineColor(4);m_hI[var]->SetLineStyle(kDashed);
  m_hGJ[var]->SetLineColor(3);m_hGJ[var]->SetLineStyle(kDashed);
  m_hJG[var]->SetLineColor(2);m_hJG[var]->SetLineStyle(kDashed);
  m_hJJ[var]->SetLineColor(6);m_hJJ[var]->SetLineStyle(kDashed);

  ExtendedCanvas *cbkg_signi = new ExtendedCanvas("cbkg_signi"+var,"bkg with significance "+var,800,600,2);
  TPad* p1 = (TPad*)cbkg_signi->GetPad(1);
  TPad* p2 = (TPad*)cbkg_signi->GetPad(2);

  double maxhist = m_hB[var]->GetBinContent(m_hB[var]->GetMaximumBin());
  if(  m_hD[var]->GetBinContent(m_hD[var]->GetMaximumBin())> maxhist )
    maxhist = m_hD[var]->GetBinContent(m_hD[var]->GetMaximumBin());

  
  m_hD[var]->GetYaxis()->SetRangeUser(1e-8,maxhist+0.1*maxhist);


  p1->cd();m_hD[var]->Draw("PE");m_hB[var]->Draw("sameHIST");
  m_hI[var]->Draw("sameHIST");m_hGJ[var]->Draw("sameHIST");
  m_hJG[var]->Draw("sameHIST");m_hJJ[var]->Draw("sameHIST");




  TLegend *leg = new TLegend(0.5,0.82,0.78,0.95);
  leg->SetFillColor(0);
  leg->SetNColumns(3);
  leg->AddEntry(m_hD[var],"Data","lp");
  leg->AddEntry(m_hB[var],"total bkg","l");
  leg->AddEntry(m_hI[var],"#gamma#gamma bkg","l");
  leg->AddEntry(m_hGJ[var],"#gammaj bkg","l");
  leg->AddEntry(m_hJG[var],"j#gamma bkg","l");
  leg->AddEntry(m_hJJ[var],"jj bkg","l");
  leg->Draw("same");
  p1->RedrawAxis();
  if(m_etacategory != "allcat")
    cbkg_signi->SetEtaCategoryLabel(0.25,0.88,m_etacategory);
  p1->Update();
  p2->cd();
  SignificanceHist SH(*(TH1F*)m_hD[var],*(TH1F*)m_hB[var]);
  TH1F* hh_signi = SH.GetChiHist(5);
  hh_signi->SetName("hh_signi");hh_signi->SetTitle("hh_signi");
  hh_signi->Draw("HIST");
  p2->Update();
  return cbkg_signi;

}

/////////////////////////////////////////////////////
void BkgKinematics::SaveToRootFile(TString filename)
/////////////////////////////////////////////////////
{
  TFile f(filename,"RECREATE");
  std::map<TString,TH1D*>::iterator it;
  for ( it=m_hI.begin() ; it != m_hI.end(); it++ )
    f.Add( (*it).second );
  for ( it=m_hD.begin() ; it != m_hD.end(); it++ )
    f.Add( (*it).second );
  for ( it=m_hB.begin() ; it != m_hB.end(); it++ )
    f.Add( (*it).second );
  for ( it=m_hGJ.begin() ; it != m_hGJ.end(); it++ )
    f.Add( (*it).second );
  for ( it=m_hJG.begin() ; it != m_hJG.end(); it++ )
    f.Add( (*it).second );
  for ( it=m_hJJ.begin() ; it != m_hJJ.end(); it++ )
    f.Add( (*it).second );
  f.Write(); f.Close();

}
/////////////////////////////////////////////
TCanvas* BkgKinematics::PurityFitter(TString var)
/////////////////////////////////////////////
{
  TH1D* hred = (TH1D*)m_hGJ[var]->Clone("hred");
  hred->Add(m_hJG[var]);
  hred->Add(m_hJJ[var]);

  RooRealVar x("x",m_title[var],m_hD[var]->GetBinLowEdge(1),m_hD[var]->GetBinLowEdge(m_hD[var]->GetNbinsX()));
  RooDataHist data_H("data","data",x,m_hD[var]);
  RooDataHist irr_H("irr","irr",x,m_hI[var]);
  RooDataHist red_H("red","red",x,hred);

  RooHistPdf pdf_irr("pdf_irr","pdf_irr",x,irr_H);
  RooHistPdf pdf_red("pdf_red","pdf_red",x,red_H);

  RooRealVar purity("purity","purity",0.5,1);
  RooAddPdf  pdf_tot("pdf_tot","pdf_tot", pdf_irr,pdf_red,purity);

  RooFitResult* fitres = pdf_tot.fitTo(data_H,RooFit::Save(kTRUE));
  fitres->Print("v"); 
  RooPlot * frame = x.frame();
  data_H.plotOn(frame,RooFit::Name("data"));
  pdf_tot.plotOn(frame,RooFit::LineColor(4),RooFit::Name("bkg"));
  pdf_tot.plotOn(frame,RooFit::Components(pdf_irr),RooFit::LineColor(4),RooFit::LineStyle(kDashed));
  pdf_tot.plotOn(frame,RooFit::Components(pdf_red),RooFit::LineColor(2),RooFit::LineStyle(kDashed));
  data_H.plotOn(frame,RooFit::Name("data2"));
  pdf_tot.paramOn(frame);


  SignificanceHist SID( *(RooHist*)frame->getHist("data"),
			*(RooCurve*)frame->getCurve("bkg") );
  TH1F* hsigni = SID.GetSignificanceHist(5);
  hsigni->GetXaxis()->SetTitle(x.GetTitle());
  ExtendedCanvas *c = new ExtendedCanvas("c"+var,"c"+var,800,600,2);
  TPad* p1 = (TPad*)c->GetPad(1);
  TPad* p2 = (TPad*)c->GetPad(2);
  p1->cd();
  frame->Draw();
  if(m_etacategory != "allcat")
    c->SetEtaCategoryLabel(0.25,0.88,m_etacategory);
  p1->Update();
  p2->cd();
  hsigni->Draw("HIST");
  p2->Update();
  return c;
}
