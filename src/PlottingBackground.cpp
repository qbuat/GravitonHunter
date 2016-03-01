#include "PlottingBackground.h"
#include <TH1.h>
#include <TPad.h>
#include <RooPlot.h>
#include <RooCurve.h>
#include <RooHist.h>
#include <TLegend.h>

#include "ToolsExtendedCanvas.h"
#include "ToolsUtilities.h"
#include "ToolsSignificanceHist.h"
#include "ToolsCommons.h"

////////////////////////////////////////////////////////////////////////////////////////////
PlottingBackground::PlottingBackground(TString name2dfit,TString nameredfit,TString namebkg)
////////////////////////////////////////////////////////////////////////////////////////////
{
  Init2DFit(name2dfit);
  InitRedFit(nameredfit);
  InitBkg(namebkg);
  Init2DFit_Syst("NONE","NONE");
  InitRedFit_Syst("NONE","NONE","NONE");

  m_cat = "NONE";

}
/////////////////////////////////////////
PlottingBackground::~PlottingBackground()
////////////////////////////////////////
{
  delete m_file2dfit;
  delete m_fileredfit;
  delete m_filebkg;
  delete m_fileredfit_L3;
  delete m_fileredfit_L4;
  delete m_fileredfit_L5;
  delete m_file2dfit_L2;
  delete m_file2dfit_L5;

}


/////////////////////////////////////////////////////
void PlottingBackground::Init2DFit(TString name2dfit)
/////////////////////////////////////////////////////
{
  if(name2dfit=="NONE") m_file2dfit = 0;
  else  m_file2dfit = new TFile(name2dfit,"read");
}
////////////////////////////////////////////////////////////////////////
void PlottingBackground::Init2DFit_Syst(TString name_L2,TString name_L5)
////////////////////////////////////////////////////////////////////////
{
  if(name_L2=="NONE") m_file2dfit_L2 = 0;
  else  m_file2dfit_L2 = new TFile(name_L2,"read");
  if(name_L5=="NONE") m_file2dfit_L5 = 0;
  else  m_file2dfit_L5 = new TFile(name_L5,"read");
}
//////////////////////////////////////////////////////
void PlottingBackground::InitRedFit(TString nameredfit)
//////////////////////////////////////////////////////
{
  if(nameredfit=="NONE") m_fileredfit = 0;
  else  m_fileredfit = new TFile(nameredfit,"read");
}
/////////////////////////////////////////////////////////////////////////////////////////
void PlottingBackground::InitRedFit_Syst(TString name_L3,TString name_L4,TString name_L5)
////////////////////////////////////////////////////////////////////////////////////////
{
  if(name_L3=="NONE") m_fileredfit_L3 = 0;
  else  m_fileredfit_L3 = new TFile(name_L3,"read");

  if(name_L4=="NONE") m_fileredfit_L4 = 0;
  else  m_fileredfit_L4 = new TFile(name_L4,"read");

  if(name_L5=="NONE") m_fileredfit_L5 = 0;
  else  m_fileredfit_L5 = new TFile(name_L5,"read");
}
/////////////////////////////////////////////////
void PlottingBackground::InitBkg(TString namebkg)
////////////////////////////////////////////////
{
  if(namebkg=="NONE")    m_filebkg = 0;
  else     m_filebkg = new TFile(namebkg,"read");
}

///////////////////////////////////////////////
TCanvas* PlottingBackground::Plot2DFitLeadJet()
//////////////////////////////////////////////
{
  RooPlot* frame_JetL   = (RooPlot*)m_file2dfit->Get("JetL");
  SignificanceHist SJL( *(RooHist*)frame_JetL->getHist("data_JL"),
                        *(RooCurve*)frame_JetL->getCurve("pdf_JL"));
  TH1F* signi_JL = SJL.GetSignificanceHist(4);
  signi_JL->GetXaxis()->SetTitle(frame_JetL->GetXaxis()->GetTitle());

  ExtendedCanvas *cJL = new ExtendedCanvas("cJL","Lead Jet",800,600,2);
  TPad* pJL1 = (TPad*)cJL->GetPad(1);
  TPad* pJL2 = (TPad*)cJL->GetPad(2);
  pJL1->cd(); frame_JetL->Draw();
  cJL->SetLumiLabel(0.72,0.83,Commons::int_lumi_fb);
  cJL->SetCdmLabel(0.75,0.73);
  cJL->SetEtaCategoryLabel(0.18,0.90,m_cat);
  pJL1->Update();
  pJL2->cd(); signi_JL->Draw("HIST");
  pJL2->Update();
  return cJL;

}
//////////////////////////////////////////////
TCanvas* PlottingBackground::Plot2DFitLeadPh()
/////////////////////////////////////////////
{
  RooPlot* frame_PHL    = (RooPlot*)m_file2dfit->Get("PHL");
  SignificanceHist SPHL( *(RooHist*)frame_PHL->getHist("data_TiL"),
                        *(RooCurve*)frame_PHL->getCurve("pdf_TiL"));
  TH1F* signi_PHL = SPHL.GetSignificanceHist(4);
  signi_PHL->GetXaxis()->SetTitle(frame_PHL->GetXaxis()->GetTitle());

  ExtendedCanvas *cPHL = new ExtendedCanvas("cPHL","Lead Photon",800,600,2);
  TPad* pPHL1 = (TPad*)cPHL->GetPad(1);
  TPad* pPHL2 = (TPad*)cPHL->GetPad(2);
  pPHL1->cd(); frame_PHL->Draw();
  cPHL->SetLumiLabel(0.72,0.83,Commons::int_lumi_fb);
  cPHL->SetCdmLabel(0.75,0.73);
  cPHL->SetEtaCategoryLabel(0.18,0.90,m_cat);
  pPHL1->Update();
  pPHL2->cd(); signi_PHL->Draw("HIST");
  pPHL2->Update();
  return cPHL;

}
////////////////////////////////////////////////
TCanvas* PlottingBackground::Plot2DFitLeadFinal()
////////////////////////////////////////////////
{
  RooPlot* frame_TiTiL  = (RooPlot*)m_file2dfit->Get("TiTiL");

  TString names_TiTi_L[5] = { "pdf_tot","pdf_phph",
			      "pdf_phjet","pdf_jetph",
                              "pdf_jetjet" };
  TString signs_TiTi_L[5] = { "#gamma#gamma+#gammaj+j#gamma+jj",
                              "#gamma#gamma","#gammaj","j#gamma","jj" };

  TLegend* leg_TiTi_L = new TLegend(0.7,0.25,0.9,0.65);
  leg_TiTi_L->SetFillColor(0);
  for( int ileg=0;ileg<5;ileg++)
    leg_TiTi_L->AddEntry( frame_TiTiL->findObject(names_TiTi_L[ileg].Data()),
                          signs_TiTi_L[ileg],"l" );


  ExtendedCanvas *cTiTiL = new ExtendedCanvas("cTiTiL","Lead Tight final",800,600,2);
  TPad* pTiTiL1 = (TPad*)cTiTiL->GetPad(1);
  TPad* pTiTiL2 = (TPad*)cTiTiL->GetPad(2);
  pTiTiL1->cd();
  frame_TiTiL->Draw();
  cTiTiL->SetLumiLabel(0.72,0.83,Commons::int_lumi_fb);
  cTiTiL->SetCdmLabel(0.75,0.73);
  cTiTiL->SetEtaCategoryLabel(0.18,0.90,m_cat);
  pTiTiL1->Update();
  leg_TiTi_L->Draw("same");
  pTiTiL2->cd();
  SignificanceHist STiTiL( *(RooHist*)frame_TiTiL->getHist("data_L"),
                           *(RooCurve*)frame_TiTiL->getCurve("pdf_tot"));
  TH1F* signi_TiTiL = STiTiL.GetSignificanceHist(4);
  signi_TiTiL->GetXaxis()->SetTitle(frame_TiTiL->GetXaxis()->GetTitle());
  signi_TiTiL->Draw("HIST");
  pTiTiL2->Update();
  return cTiTiL;

} 
/////////////////////////////////////////////////
TCanvas* PlottingBackground::Plot2DFitSubLeadJet()
////////////////////////////////////////////////
{
  RooPlot* frame_JetSL  = (RooPlot*)m_file2dfit->Get("JetSL");
  SignificanceHist SJSL( *(RooHist*)frame_JetSL->getHist("data_JSL"),
                         *(RooCurve*)frame_JetSL->getCurve("pdf_JSL"));
  TH1F* signi_JSL = SJSL.GetSignificanceHist(4);
  signi_JSL->GetXaxis()->SetTitle(frame_JetSL->GetXaxis()->GetTitle());

  ExtendedCanvas *cJSL = new ExtendedCanvas("cSJL", "SubLead Jet",800,600,2);
  TPad* pJSL1 = (TPad*)cJSL->GetPad(1);
  TPad* pJSL2 = (TPad*)cJSL->GetPad(2);
  pJSL1->cd(); frame_JetSL->Draw();
  cJSL->SetLumiLabel(0.72,0.83,Commons::int_lumi_fb);
  cJSL->SetCdmLabel(0.75,0.73);
  cJSL->SetEtaCategoryLabel(0.18,0.90,m_cat);
  pJSL1->Update();
  pJSL2->cd(); signi_JSL->Draw("HIST");
  pJSL2->Update();
  return cJSL;

}
////////////////////////////////////////////////
TCanvas* PlottingBackground::Plot2DFitSubLeadPh()
////////////////////////////////////////////////
{
  RooPlot* frame_PHSL   = (RooPlot*)m_file2dfit->Get("PHSL");
  SignificanceHist SPHSL( *(RooHist*)frame_PHSL->getHist("data_TiSL"),
                        *(RooCurve*)frame_PHSL->getCurve("pdf_TiSL"));
  TH1F* signi_PHSL = SPHSL.GetSignificanceHist(4);
  signi_PHSL->GetXaxis()->SetTitle(frame_PHSL->GetXaxis()->GetTitle());

  ExtendedCanvas *cPHSL = new ExtendedCanvas("cPHSL","Sublead Photon",800,600,2);
  TPad* pPHSL1 = (TPad*)cPHSL->GetPad(1);
  TPad* pPHSL2 = (TPad*)cPHSL->GetPad(2);
  pPHSL1->cd(); frame_PHSL->Draw();
  cPHSL->SetLumiLabel(0.72,0.83,Commons::int_lumi_fb);
  cPHSL->SetCdmLabel(0.75,0.73);
  cPHSL->SetEtaCategoryLabel(0.18,0.90,m_cat);
  pPHSL1->Update();
  pPHSL2->cd(); signi_PHSL->Draw("HIST");
  pPHSL2->Update();
  return cPHSL;

}
///////////////////////////////////////////////////
TCanvas* PlottingBackground::Plot2DFitSubLeadFinal()
///////////////////////////////////////////////////
{
  RooPlot* frame_TiTiSL = (RooPlot*)m_file2dfit->Get("TiTiSL");

  TLegend* leg_TiTi_SL = new TLegend(0.7,0.3,0.9,0.7);
  leg_TiTi_SL->SetFillColor(0);
  TString names_TiTi_SL[5] = { "pdf_tot","pdf_phph",
                               "pdf_jetph","pdf_phjet",
                               "pdf_jetjet" };
  TString signs_TiTi_SL[5] = { "#gamma#gamma+#gammaj+j#gamma+jj",
                               "#gamma#gamma","j#gamma","#gammaj","jj" };
  for( int ileg=0;ileg<5;ileg++)
    leg_TiTi_SL->AddEntry( frame_TiTiSL->findObject(names_TiTi_SL[ileg].Data()),
                           signs_TiTi_SL[ileg],"l" );

  ExtendedCanvas *cTiTiSL = new ExtendedCanvas("cTiTiSL","Sublead Tight Final",800,600,2);
  TPad* pTiTiSL1 = (TPad*)cTiTiSL->GetPad(1);
  TPad* pTiTiSL2 = (TPad*)cTiTiSL->GetPad(2);
  pTiTiSL1->cd();
  frame_TiTiSL->Draw();
  leg_TiTi_SL->Draw("same");
  cTiTiSL->SetLumiLabel(0.72,0.83,Commons::int_lumi_fb);
  cTiTiSL->SetCdmLabel(0.75,0.73);
  cTiTiSL->SetEtaCategoryLabel(0.18,0.90,m_cat);
  pTiTiSL1->Update();
  pTiTiSL2->cd();
  SignificanceHist STiTiSL( *(RooHist*)frame_TiTiSL->getHist("data_SL"),
                            *(RooCurve*)frame_TiTiSL->getCurve("pdf_tot"));
  TH1F* signi_TiTiSL = STiTiSL.GetSignificanceHist(4);
  signi_TiTiSL->GetXaxis()->SetTitle(frame_TiTiSL->GetXaxis()->GetTitle());
  signi_TiTiSL->Draw("HIST");
  pTiTiSL2->Update();
  return cTiTiSL;

}
//////////////////////////////////////////////////
TCanvas* PlottingBackground::Plot2DFitJetJetHist()
/////////////////////////////////////////////////
{
  TH2F* hJetJet = (TH2F*)m_file2dfit->Get("JetJetHist");
  TCanvas * cJJ = new TCanvas("cJJ","Jet-Jet Histogram",800,600);
  cJJ->cd(); hJetJet->Draw("lego");
  return cJJ;

}
//////////////////////////////////////////////////
TCanvas* PlottingBackground::Plot2DFitJetJetCurve()
//////////////////////////////////////////////////
{
  TH2F* hJetJetFit = (TH2F*)m_file2dfit->Get("JetJetPdf");
  TCanvas* cJJF = new TCanvas("cJJF","Jet-Jet Smoothed hist",800,600);
  cJJF->cd(); hJetJetFit->Draw("lego");
  return cJJF;

}
///////////////////////////////////////////
TCanvas* PlottingBackground::PlotGJDataFit()
///////////////////////////////////////////
{


  std::cout<<m_fileredfit->GetName()<<std::endl;
  RooPlot* frame_sublead_fake = (RooPlot*)m_fileredfit->Get("sublead_fake");
  SignificanceHist SISL( *(RooHist*)frame_sublead_fake->getHist("TiLo"),
			 *(RooCurve*)frame_sublead_fake->getCurve("sublead_fake_pdf") );
  TH1F* signi_sublead = SISL.GetSignificanceHist(4);

  ExtendedCanvas * csublead = new ExtendedCanvas("csublead","sublead fake",800,600,2);
  TPad* psublead1 = (TPad*)csublead->GetPad(1);
  TPad* psublead2 = (TPad*)csublead->GetPad(2);
  psublead1->cd(); psublead1->SetLogy();
  //  frame_sublead_fake->GetYaxis()->SetRangeUser(1e-5,20000);
  frame_sublead_fake->GetYaxis()->SetRangeUser(0.01,10000);
  frame_sublead_fake->Draw();
  csublead->SetGenericLabel(0.4,0.88,"Tight-AntiTight");
  csublead->SetEtaCategoryLabel(0.4,0.80,m_cat);
  csublead->SetLumiLabel(0.72,0.83,Commons::int_lumi_fb);
  csublead->SetCdmLabel(0.75,0.73);
  psublead1->Update();
  psublead2->cd();
  signi_sublead->GetXaxis()->SetTitle(frame_sublead_fake->GetXaxis()->GetTitle());
  signi_sublead->Draw("HIST");
  psublead2->Update();
  return csublead;

}
/////////////////////////////////////////////
TCanvas*  PlottingBackground::PlotJGDataFit()
////////////////////////////////////////////
{
  RooPlot* frame_lead_fake =(RooPlot*)m_fileredfit->Get("lead_fake");
  SignificanceHist SIL( *(RooHist*)frame_lead_fake->getHist("LoTi"),
			*(RooCurve*)frame_lead_fake->getCurve("lead_fake_pdf") );
  TH1F* signi_lead = SIL.GetSignificanceHist(4);
  signi_lead->GetXaxis()->SetTitle(frame_lead_fake->GetXaxis()->GetTitle());
 
  ExtendedCanvas * clead = new ExtendedCanvas("clead","lead fake",800,600,2);
  TPad* plead1 = (TPad*)clead->GetPad(1);
  TPad* plead2 = (TPad*)clead->GetPad(2);
  plead1->cd(); plead1->SetLogy();
  //  frame_lead_fake->GetYaxis()->SetRangeUser(1e-5,20000);
  frame_lead_fake->GetYaxis()->SetRangeUser(0.01,10000);
  frame_lead_fake->Draw();
  clead->SetGenericLabel(0.4,0.88,"AntiTight-Tight");
  clead->SetEtaCategoryLabel(0.4,0.80,m_cat);
  clead->SetLumiLabel(0.72,0.83,Commons::int_lumi_fb);
  clead->SetCdmLabel(0.75,0.73);
  plead1->Update(); plead2->cd();
  signi_lead->Draw("HIST");
  plead2->Update();
  return clead;

}
////////////////////////////////////////////
TCanvas*  PlottingBackground::PlotJJDataFit()
////////////////////////////////////////////
{
  RooPlot* frame_double_fake =(RooPlot*)m_fileredfit->Get("double_fake");
  SignificanceHist SID( *(RooHist*)frame_double_fake->getHist("LoLo"),
			*(RooCurve*)frame_double_fake->getCurve("double_fake_pdf") );
  TH1F* signi_double = SID.GetSignificanceHist(4);

  //  frame_double_fake->GetYaxis()->SetRangeUser(1e-5,20000);
  frame_double_fake->GetYaxis()->SetRangeUser(0.01,10000);
  signi_double->GetXaxis()->SetTitle(frame_double_fake->GetXaxis()->GetTitle());

  ExtendedCanvas *cdouble = new ExtendedCanvas("cdouble","double fake",800,600,2);
  TPad* pdouble1 = (TPad*)cdouble->GetPad(1);
  TPad* pdouble2 = (TPad*)cdouble->GetPad(2);
  pdouble1->cd(); pdouble1->SetLogy();
  frame_double_fake->Draw();
  cdouble->SetGenericLabel(0.4,0.88,"AntiTight-AntiTight");
  cdouble->SetEtaCategoryLabel(0.4,0.80,m_cat);
  cdouble->SetLumiLabel(0.72,0.83,Commons::int_lumi_fb);
  cdouble->SetCdmLabel(0.75,0.73);
  pdouble1->Update();
  pdouble2->cd();
  signi_double->Draw("HIST");
  pdouble2->Update();
  return cdouble;

}
/////////////////////////////////////////////////////////////
TCanvas* PlottingBackground::PlotRedDataFit_syst(TString comp)
////////////////////////////////////////////////////////////
{
  TH1D* hnom = (TH1D*)m_fileredfit->Get(comp);
  TH1D* hsys_L3 = (TH1D*)m_fileredfit_L3->Get(comp);
  TH1D* hsys_L4 = (TH1D*)m_fileredfit_L4->Get(comp);
  TH1D* hsys_L5 = (TH1D*)m_fileredfit_L5->Get(comp);

  hnom->Scale(1./hnom->Integral(Commons::norm_bin.first,Commons::norm_bin.second));
  hsys_L3->Scale(1./hsys_L3->Integral(Commons::norm_bin.first,Commons::norm_bin.second));
  hsys_L4->Scale(1./hsys_L4->Integral(Commons::norm_bin.first,Commons::norm_bin.second));
  hsys_L5->Scale(1./hsys_L5->Integral(Commons::norm_bin.first,Commons::norm_bin.second));

  hsys_L3->Divide(hnom);
  hsys_L4->Divide(hnom);
  hsys_L5->Divide(hnom);
  hnom->Divide(hnom);

  hsys_L3->SetLineColor(kBlue);
  hsys_L4->SetLineColor(kRed);
  hsys_L5->SetLineColor(kGreen);

  hnom->GetYaxis()->SetRangeUser(-1,3);
  TCanvas* c1 = new TCanvas();
  c1->cd();
  hnom->Draw("HIST");
  hsys_L3->Draw("SAMEHIST");
  hsys_L4->Draw("SAMEHIST");
  hsys_L5->Draw("SAMEHIST");
  c1->BuildLegend();
  return c1;
}






/////////////////////////////////////////////
TCanvas* PlottingBackground::PlotBkgEstimate()
////////////////////////////////////////////
{
  TH1D* hh_dat = (TH1D*)m_filebkg->Get("mgg_data");
  TH1D* hh_bkg = (TH1D*)m_filebkg->Get("bkg_total_gg_full");
  TH1D* hh_red = (TH1D*)m_filebkg->Get("bkg_reducible_gg_full");

  // TGraphAsymmErrors* gmgg_bkg = (TGraphAsymmErrors*)m_filebkg->Get("gmgg_bkg");
  // TGraphAsymmErrors* gmgg_red = (TGraphAsymmErrors*)m_filebkg->Get("gmgg_red");
  TGraphAsymmErrors* gmgg_dat = (TGraphAsymmErrors*)m_filebkg->Get("gmgg_data");
  std::cout << "Histograms and graphs are loaded ... "<< std::endl;
  


  //--> Recompute the bkg graph (to draw the error bars) from the histograms
  TH1D* hh_bkguncert = (TH1D*)m_filebkg->Get("bkg_total_syst_gg_full");
  TH1D* hh_bkg_withsyst = (TH1D*)hh_bkg->Clone("bkg_total_gg_full_withsyst");
  for(int ibin= 0; ibin < hh_bkg_withsyst->GetNbinsX()+1;ibin++){
    hh_bkg_withsyst->SetBinError(ibin, 
				 sqrt( pow(hh_bkg_withsyst->GetBinError(ibin),2)+
				       pow(hh_bkg->GetBinContent(ibin)*hh_bkguncert->GetBinContent(ibin),2) )  );
  }
  TGraphAsymmErrors* gmgg_bkg = AnalysisTools::BkgGraph(*hh_bkg_withsyst);


  TH1D* huncert_red = (TH1D*)m_filebkg->Get("bkg_reducible_syst_gg_full");
  TH1D* hh_redbkg_withsyst = (TH1D*)hh_red->Clone("bkg_red_gg_full_withsyst");
  //--> Recompute purity
  double purity = (hh_bkg->Integral(Commons::norm_bin.first,Commons::norm_bin.second)-hh_red->Integral(Commons::norm_bin.first,Commons::norm_bin.second));
  purity *= 1./hh_bkg->Integral(Commons::norm_bin.first,Commons::norm_bin.second);
  double antipurity = hh_red->Integral(Commons::norm_bin.first,Commons::norm_bin.second)/ hh_bkg->Integral(Commons::norm_bin.first,Commons::norm_bin.second);
  std::cout << purity << " / " << 1-antipurity << std::endl;

  double antipurity_relative_uncert = Commons::purity_relative_uncert*purity/antipurity; //antipurity relative uncert on the red bkg

  for( int ibin = 0; ibin<=hh_redbkg_withsyst->GetNbinsX()+1;ibin++){

    double red_shape_uncert = huncert_red->GetBinContent(ibin)*hh_bkg->GetBinContent(ibin);
    double red_uncert = sqrt( pow(red_shape_uncert,2) +//shape
                              pow(antipurity_relative_uncert*hh_red->GetBinContent(ibin),2) );//antipurity
    hh_redbkg_withsyst->SetBinError(ibin, red_uncert);
  }
  TGraphAsymmErrors* gmgg_red = AnalysisTools::BkgGraph(*hh_redbkg_withsyst);


  SignificanceHist SH(*(TH1F*)hh_dat,*(TH1F*)hh_bkg);
  TH1F* hh_signi = SH.GetSignificanceHist(4);
  hh_signi->SetName("hh_signi_bkg");
  hh_signi->SetTitle("hh_signi_bkg");

  std::cout << "The significance is computed ... "<< std::endl;

  hh_bkg->GetXaxis()->SetRange(Commons::norm_bin.first,Commons::nBins);
  hh_bkg->GetYaxis()->SetRangeUser(0.004,10000);
  hh_bkg->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  hh_bkg->GetYaxis()->SetTitle("Events/bin");
  hh_bkg->GetXaxis()->SetMoreLogLabels();
  hh_bkg->GetXaxis()->SetNoExponent();
  hh_bkg->SetLineColor(4);
  std::cout << "The background histogram style is set ... "<< std::endl;
  hh_signi->GetXaxis()->SetRange(Commons::norm_bin.first,Commons::nBins);  
  hh_signi->GetXaxis()->SetMoreLogLabels();
  hh_signi->GetXaxis()->SetNoExponent();

  gmgg_bkg->SetFillColor(kOrange);
  gmgg_bkg->SetLineColor(kOrange);
  gmgg_bkg->SetMarkerSize(0);

  hh_red->SetLineColor(4);
  hh_red->SetLineStyle(2);
  gmgg_red->SetFillColor(kYellow);
  gmgg_red->SetLineColor(kYellow);
  std::cout << "The other histograms/graphs style is set ... "<< std::endl;

  //---------------------------------------------------------------------------
  TLegend *leg = new TLegend(0.2,0.2,0.4,0.4);
  leg->SetFillColor(0);leg->SetTextSize(0.045);
  leg->AddEntry(gmgg_dat,"2012 data","lp");
  leg->AddEntry(hh_bkg,"Total Background","f");
  leg->AddEntry(hh_red,"Reducible Background","f");
  leg->AddEntry(gmgg_red,"syst #oplus stat (reducible)","f");
  leg->AddEntry(gmgg_bkg,"syst #oplus stat (total)","f");
  //--------------------------------------------------------------------------
  std::cout << "The legend is built ... "<< std::endl;

  //------------------------------------------------------------------------------------
  ExtendedCanvas *cbkg_signi = new ExtendedCanvas("cbkg_signi","bkg with significance",800,600,2);
  TPad* p1 = (TPad*)cbkg_signi->GetPad(1);
  p1->SetLogy();p1->SetLogx();
  TPad* p2 = (TPad*)cbkg_signi->GetPad(2);
  p2->SetLogx();
  p1->cd();
  hh_bkg->Draw("HIST][");
  gmgg_bkg->Draw("sameE2");
  hh_bkg->Draw("sameHIST][");
  gmgg_red->Draw("sameE2");
  hh_red->Draw("sameHIST][");
  gmgg_dat->Draw("sameP");


  leg->Draw("same");
  cbkg_signi->SetAtlasLabel(0.52,0.88);
  cbkg_signi->SetLumiLabel(0.72,0.73,Commons::int_lumi_fb);
  cbkg_signi->SetCdmLabel(0.75,0.63);
  cbkg_signi->SetEtaCategoryLabel(0.18,0.05,m_cat);
  p1->RedrawAxis();
  p1->Update();
  p2->cd();
  hh_signi->Draw("HIST");
  p2->RedrawAxis();
  p2->Update();






  std::cout << " Bin \t\t reducible \t irreducible \t\t total \t\t\t data \t\t data-total \t\t Significance" << std::endl;
  for( int ibin = 1; ibin<=hh_bkg->GetNbinsX();ibin++){
    if( ibin< Commons::norm_bin.first ) continue;
    if( ibin == Commons::norm_bin.first ) std::cout << "---- \t\t Start of the control region \t\t ---" << std::endl;
    std::cout << Form("[%1.1f,%1.1f]",hh_bkg->GetBinLowEdge(ibin),hh_bkg->GetBinLowEdge(ibin+1))
	      << "\t"
	      << Form("%1.3f",hh_red->GetBinContent(ibin) )
	      << "\t\t"
	      << Form("%1.3f",hh_bkg->GetBinContent(ibin) - hh_red->GetBinContent(ibin) )
	      << "\t\t"
	      << Form("%1.3f",hh_bkg->GetBinContent(ibin) )
	      << "\t\t"
	      << Form("%1.3f",hh_dat->GetBinContent(ibin) )
	      << "\t\t"
	      << Form("%1.3f",hh_dat->GetBinContent(ibin)-hh_bkg->GetBinContent(ibin) )
	      << "\t\t"
	      << Form("%1.3f",hh_signi->GetBinContent(ibin) )
	      << std::endl;
    if( ibin == Commons::norm_bin.second ) std::cout << "---- \t\t End of the control region \t\t ---" << std::endl;

  }

  return cbkg_signi;

}

/////////////////////////////////////////////////////////////
TCanvas* PlottingBackground::PlotBkgEstimate_withratiohist()
////////////////////////////////////////////////////////////
{
  TH1D* hh_dat = (TH1D*)m_filebkg->Get("mgg_data");
  TH1D* hh_bkg = (TH1D*)m_filebkg->Get("bkg_total_gg_full");
  TH1D* hh_red = (TH1D*)m_filebkg->Get("bkg_reducible_gg_full");
  // TGraphAsymmErrors* gmgg_bkg = (TGraphAsymmErrors*)m_filebkg->Get("gmgg_bkg");
  // TGraphAsymmErrors* gmgg_red = (TGraphAsymmErrors*)m_filebkg->Get("gmgg_red");
  TGraphAsymmErrors* gmgg_dat = (TGraphAsymmErrors*)m_filebkg->Get("gmgg_data");
  std::cout << "Histograms and graphs are loaded ... "<< std::endl;
  
  SignificanceHist SH(*(TH1F*)hh_dat,*(TH1F*)hh_bkg);
  TH1F* hh_ratio = SH.GetRatioHist(0.5);
  hh_ratio->SetName("hh_ratio_bkg");
  hh_ratio->SetTitle("hh_ratio_bkg");
  std::cout << "The ratio is computed ... "<< std::endl;


  TH1D* hh_bkguncert = (TH1D*)m_filebkg->Get("bkg_total_syst_gg_full");
  TH1D* hh_bkg_withsyst = (TH1D*)hh_bkg->Clone("bkg_total_gg_full_withsyst");
  for(int ibin= 0; ibin < hh_bkg_withsyst->GetNbinsX()+1;ibin++){
    hh_bkg_withsyst->SetBinError(ibin, 
				 sqrt( pow(hh_bkg_withsyst->GetBinError(ibin),2)+
				       pow(hh_bkg->GetBinContent(ibin)*hh_bkguncert->GetBinContent(ibin),2) )  );
  }
  TGraphAsymmErrors* gmgg_bkg = AnalysisTools::BkgGraph(*hh_bkg_withsyst);


  TH1D* huncert_red = (TH1D*)m_filebkg->Get("bkg_reducible_syst_gg_full");
  TH1D* hh_redbkg_withsyst = (TH1D*)hh_red->Clone("bkg_red_gg_full_withsyst");
  //--> Recompute purity
  double purity = (hh_bkg->Integral(Commons::norm_bin.first,Commons::norm_bin.second)-hh_red->Integral(Commons::norm_bin.first,Commons::norm_bin.second));
  purity *= 1./hh_bkg->Integral(Commons::norm_bin.first,Commons::norm_bin.second);
  double antipurity = hh_red->Integral(Commons::norm_bin.first,Commons::norm_bin.second)/ hh_bkg->Integral(Commons::norm_bin.first,Commons::norm_bin.second);
  std::cout << purity << " / " << 1-antipurity << std::endl;

  double antipurity_relative_uncert = Commons::purity_relative_uncert*purity/antipurity; //antipurity relative uncert on the red bkg

  for( int ibin = 0; ibin<=hh_redbkg_withsyst->GetNbinsX()+1;ibin++){

    double red_shape_uncert = huncert_red->GetBinContent(ibin)*hh_bkg->GetBinContent(ibin);
    double red_uncert = sqrt( pow(red_shape_uncert,2) +//shape
                              pow(antipurity_relative_uncert*hh_red->GetBinContent(ibin),2) );//antipurity
    hh_redbkg_withsyst->SetBinError(ibin, red_uncert);
  }
  TGraphAsymmErrors* gmgg_red = AnalysisTools::BkgGraph(*hh_redbkg_withsyst);


  TGraphAsymmErrors* gr_uncertband = new TGraphAsymmErrors();
  for(int ibin=1;ibin<(hh_bkguncert->GetNbinsX()+1);ibin++){
    gr_uncertband->SetPoint(ibin-1,hh_bkguncert->GetBinCenter(ibin),1);
    gr_uncertband->SetPointError(ibin-1,hh_bkguncert->GetBinWidth(ibin)/2,hh_bkguncert->GetBinWidth(ibin)/2,
  				 hh_bkguncert->GetBinContent(ibin),hh_bkguncert->GetBinContent(ibin));
  }
  gr_uncertband->SetFillColor(kOrange);
  hh_bkg->GetXaxis()->SetRange(Commons::norm_bin.first,Commons::nBins);
  hh_bkg->GetYaxis()->SetRangeUser(0.004,10000);
  hh_bkg->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  hh_bkg->GetYaxis()->SetTitle("Events/bin");
  hh_bkg->GetXaxis()->SetMoreLogLabels();
  hh_bkg->GetXaxis()->SetNoExponent();
  hh_bkg->SetLineColor(4);
  std::cout << "The background histogram style is set ... "<< std::endl;

  hh_ratio->GetXaxis()->SetRange(Commons::norm_bin.first,Commons::nBins);  
  hh_ratio->GetXaxis()->SetMoreLogLabels();
  hh_ratio->GetXaxis()->SetNoExponent();

  gmgg_bkg->SetFillColor(kOrange);
  gmgg_bkg->SetLineColor(kOrange);
  gmgg_bkg->SetMarkerSize(0);



  hh_red->SetLineColor(4);
  hh_red->SetLineStyle(2);
  gmgg_red->SetFillColor(kYellow);
  gmgg_red->SetLineColor(kYellow);
  gmgg_red->SetMarkerSize(0);



  std::cout << "The other histograms/graphs style is set ... "<< std::endl;

  //---------------------------------------------------------------------------
  TLegend *leg = new TLegend(0.2,0.2,0.4,0.4);
  leg->SetFillColor(0);leg->SetTextSize(0.045);
  leg->AddEntry(gmgg_dat,"2012 data","lp");
  leg->AddEntry(hh_bkg,"Total Background","f");
  leg->AddEntry(hh_red,"Reducible Background","f");
  leg->AddEntry(gmgg_red,"syst #oplus stat (reducible)","f");
  leg->AddEntry(gmgg_bkg,"syst #oplus stat (total)","f");
  //--------------------------------------------------------------------------
  std::cout << "The legend is built ... "<< std::endl;

  //------------------------------------------------------------------------------------
  ExtendedCanvas *cbkg_ratio = new ExtendedCanvas("cbkg_ratio","bkg with ratio",800,600,2);
  TPad* p1 = (TPad*)cbkg_ratio->GetPad(1);
  p1->SetLogy();p1->SetLogx();
  TPad* p2 = (TPad*)cbkg_ratio->GetPad(2);
  p2->SetLogx();
  p1->cd();
  hh_bkg->Draw("HIST][");
  gmgg_red->Draw("sameE2");
  hh_red->Draw("sameHIST][");
  gmgg_bkg->Draw("sameE2");
  hh_bkg->Draw("sameHIST][");
  gmgg_dat->Draw("sameP");
  p1->RedrawAxis();

  leg->Draw("same");
  cbkg_ratio->SetAtlasLabel(0.52,0.88);
  cbkg_ratio->SetLumiLabel(0.72,0.73,Commons::int_lumi_fb);
  cbkg_ratio->SetCdmLabel(0.75,0.63);
  cbkg_ratio->SetEtaCategoryLabel(0.18,0.05,m_cat);

  p1->Update();
  p2->cd();
  hh_ratio->SetMarkerSize(0.5);
  hh_ratio->Draw("PE");
  gr_uncertband->Draw("sameE2");
  hh_ratio->Draw("samePE");
  p2->RedrawAxis();
  p2->Update();

  return cbkg_ratio;

}


//////////////////////////////////////////////////
TCanvas* PlottingBackground::PlotBkgUncertainty()
//////////////////////////////////////////////////
{

  TH1D* hh_pur_syst = (TH1D*)m_filebkg->Get("bkg_purity_syst_gg_full");
  TH1D* hh_irr_syst = (TH1D*)m_filebkg->Get("bkg_irreducible_syst_gg_full");
  TH1D* hh_red_syst = (TH1D*)m_filebkg->Get("bkg_reducible_syst_gg_full");
  TH1D* hh_tot_syst = (TH1D*)m_filebkg->Get("bkg_total_syst_gg_full");


  std::cout<<hh_irr_syst->GetMaximum()<<std::endl;



  ExtendedCanvas *c1 = new ExtendedCanvas("cbkguncert","Bkg uncertainty",800,600);
  c1->cd();
  hh_pur_syst->SetLineColor(3);
  hh_red_syst->SetLineColor(4);
  hh_irr_syst->SetLineColor(6);
  hh_tot_syst->GetXaxis()->SetRange(Commons::norm_bin.first,Commons::nBins);
  hh_tot_syst->GetYaxis()->SetRange(-0.01,0.40);
  hh_tot_syst->Draw("HIST][");
  hh_pur_syst->Draw("sameHIST][");
  hh_irr_syst->Draw("sameHIST][");
  hh_red_syst->Draw("sameHIST][");
  TLegend *leg = new TLegend(0.2,0.7,0.5,0.9);
  leg->SetFillColor(0);
  leg->AddEntry(hh_tot_syst,"Total uncertainty","l");
  leg->AddEntry(hh_irr_syst,"Irreducible background","l");
  leg->AddEntry(hh_pur_syst,"Yields estimate","l");
  leg->AddEntry(hh_red_syst,"Reducible background","l");
  leg->Draw("same");
  c1->SetAtlasLabel(0.52,0.88);
  c1->SetCdmLabel(0.75,0.63);
  c1->SetEtaCategoryLabel(0.18,0.05,m_cat);

  return c1;
}


///////////////////////////////////////////////////
TCanvas* PlottingBackground::Plot2DFitLeadPh_Syst()
///////////////////////////////////////////////////
{
  ExtendedCanvas *c1 = new ExtendedCanvas("cPHL_syst","Lead Photon with syst",800,600);
  RooPlot* frame_PHL    = (RooPlot*)m_file2dfit->Get("PHL");
  RooPlot* frame_PHL_L2 = (RooPlot*)m_file2dfit_L2->Get("PHL");
  RooPlot* frame_PHL_L5 = (RooPlot*)m_file2dfit_L5->Get("PHL");
  RooCurve* Curv_PH_L2  = (RooCurve*)frame_PHL_L2->getCurve("pdf_PHL");
  RooCurve* Curv_PH_L5  = (RooCurve*)frame_PHL_L5->getCurve("pdf_PHL");
  RooCurve* Curv_J_L2  = (RooCurve*)frame_PHL_L2->getCurve("pdf_JL");
  RooCurve* Curv_J_L5  = (RooCurve*)frame_PHL_L5->getCurve("pdf_JL");

  frame_PHL->Draw();
  Curv_PH_L2->Draw("same");
  Curv_PH_L5->Draw("same");
  Curv_J_L2->Draw("same");
  Curv_J_L5->Draw("same");
  c1->SetLumiLabel(0.72,0.83,Commons::int_lumi_fb);
  c1->SetCdmLabel(0.75,0.73);
  c1->SetEtaCategoryLabel(0.18,0.90,m_cat);
  return c1;

}
//////////////////////////////////////////////////////
TCanvas* PlottingBackground::Plot2DFitLeadFinal_Syst()
//////////////////////////////////////////////////////
{
  ExtendedCanvas *c1 = new ExtendedCanvas("cTiTiL_syst","Lead Final with syst",800,600);
  return c1;
}
//////////////////////////////////////////////////////
TCanvas* PlottingBackground::Plot2DFitSubLeadPh_Syst()
//////////////////////////////////////////////////////
{

  ExtendedCanvas *c1 = new ExtendedCanvas("cPHSL_syst","Sublead Photon with syst",800,600);
  RooPlot* frame_PHSL    = (RooPlot*)m_file2dfit->Get("PHSL");
  RooPlot* frame_PHSL_L2 = (RooPlot*)m_file2dfit_L2->Get("PHSL");
  RooPlot* frame_PHSL_L5 = (RooPlot*)m_file2dfit_L5->Get("PHSL");
  RooCurve* Curv_PH_L2  = (RooCurve*)frame_PHSL_L2->getCurve("pdf_PHSL");
  RooCurve* Curv_PH_L5  = (RooCurve*)frame_PHSL_L5->getCurve("pdf_PHSL");
  RooCurve* Curv_J_L2  = (RooCurve*)frame_PHSL_L2->getCurve("pdf_JSL");
  RooCurve* Curv_J_L5  = (RooCurve*)frame_PHSL_L5->getCurve("pdf_JSL");

  frame_PHSL->Draw();
  Curv_PH_L2->Draw("same");
  Curv_PH_L5->Draw("same");
  Curv_J_L2->Draw("same");
  Curv_J_L5->Draw("same");
  c1->SetLumiLabel(0.72,0.83,Commons::int_lumi_fb);
  c1->SetCdmLabel(0.75,0.73);
  c1->SetEtaCategoryLabel(0.18,0.90,m_cat);
  return c1;

}
/////////////////////////////////////////////////////////
TCanvas* PlottingBackground::Plot2DFitSubLeadFinal_Syst()
/////////////////////////////////////////////////////////
{
  ExtendedCanvas *c1 = new ExtendedCanvas("cTiTiSL_syst","Sublead Final with syst",800,600);
  return c1;

}
