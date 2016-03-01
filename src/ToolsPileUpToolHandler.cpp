#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMath.h"
#include "ToolsPileUpToolHandler.h"
#include "ToolsUtilities.h"

///////////////////////////////////////
PileUpToolHandler::PileUpToolHandler()
//////////////////////////////////////
{
  //Default constructor
  m_nentries=0;
  m_mctype="";
  m_fileMC="";
  m_histMC="";
  m_fileD="";
  m_histD="";
}
/////////////////////////////////////////////////
PileUpToolHandler::PileUpToolHandler(TTree* tree)
////////////////////////////////////////////////
{
  Init(tree);
}
///////////////////////////////////////
PileUpToolHandler::~PileUpToolHandler()
//////////////////////////////////////
{
  //Destructor
  delete m_Rd;
  delete m_tPileUp;
}
//////////////////////////////////////////
void PileUpToolHandler::Init(TTree* tree)
/////////////////////////////////////////
{
  if (tree == 0) {
    std::cout << "You need to specify an input tree !!!"<<std::endl;
  }
  m_Rd       = new TreeReader(tree);
  m_nentries = (int)tree->GetEntriesFast();
  m_tPileUp  = 0;

  m_fileD  = "";
  m_histD  = "avgintperbx";
  m_fileMC="";
  m_histMC="";


}
//////////////////////////////////////////////////////////////////////////////
void PileUpToolHandler::SetMCFile(TString filename) { m_fileMC  = filename; }
/////////////////////////////////////////////////////////////////////////////
void PileUpToolHandler::SetDataFile(TString filename) {m_fileD = filename;}
///////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////
void PileUpToolHandler::GenerateMCFile(TString toolname)
////////////////////////////////////////////////////////
{
  m_tPileUp = new Root::TPileupReweighting(toolname);
  m_tPileUp->UsePeriodConfig("MC12a");
  m_tPileUp->initialize();
  for(int jentry=0;jentry<m_nentries;jentry++){
    m_Rd->GetEntry(jentry);

    //================== Trigger Requirement ===================//
    bool m_passTrigger    = (bool)m_Rd->GetVariable("EF_g35_loose_g25_loose");
    //__________________________________________________________//
    if( !m_passTrigger ) continue;

    double mu             = m_Rd->GetVariable("averageIntPerXing");
    int mcRunNumber       = (int)m_Rd->GetVariable("RunNumber");
    int mc_channel_number = (int)m_Rd->GetVariable("mc_channel_number");
    int mcevt_w00         = (int)m_Rd->GetVariable("mcevt_weight[0][0]");
    m_tPileUp->Fill(mcRunNumber,mc_channel_number,mcevt_w00,mu);
    AnalysisTools myAT;
    myAT.ShowNumberOfProcessedEvents(jentry,m_nentries,1000);
  }
  m_tPileUp->WriteToFile();
}
////////////////////////////////////////////
TCanvas* PileUpToolHandler::CheckPileUpWeights()
///////////////////////////////////////////
{

  m_tPileUp = new Root::TPileupReweighting("MyPUTool");
  m_tPileUp->AddConfigFile(m_fileMC);
  // m_tPileUp->SetDataScaleFactors(1./1.11);
  m_tPileUp->AddLumiCalcFile(m_fileD); 
  m_tPileUp->SetUnrepresentedDataAction(2);
  m_tPileUp->initialize();

  TFile* f = new TFile(m_fileD,"read");
  TH1F* hD = (TH1F*)f->Get("avgintperbx");
  hD->Rebin(10);

  TH1F* hMC1 = new TH1F("hMC1","hMC1",50,0,50);
  TH1F* hMC2 = new TH1F("hMC2","hMC2",50,0,50);
  float sumofweights=0;

  hD->Sumw2();
  hMC1->Sumw2();
  hMC2->Sumw2();

  for(int i=0; i<m_nentries;i++){
    m_Rd->GetEntry(i);
    double mu             = m_Rd->GetVariable("averageIntPerXing");
    int mcRunNumber       = (int)m_Rd->GetVariable("RunNumber");
    int mc_channel_number = (int)m_Rd->GetVariable("mc_channel_number");
    float PUWeight        = m_tPileUp->GetCombinedWeight(mcRunNumber,mc_channel_number,mu);
    hMC1->Fill( mu,PUWeight );
    hMC2->Fill( mu );
    sumofweights=sumofweights+PUWeight;
  }
  std::cout << "nentries : " << m_nentries << std ::endl;
  std::cout << "Comparison of the number of events" << std::endl;
  std::cout << "Before Reweighting = " << m_nentries <<std::endl;
  std::cout << "After Reweighting  = " << (int)sumofweights << std::endl;
  std::cout << "Difference in %    = " <<fabs(sumofweights-m_nentries)/m_nentries*100 << std::endl;

  hMC1->SetLineColor(2);
  hMC1->SetMarkerColor(2);
  hMC2->SetLineColor(3);
  hMC2->SetMarkerColor(3);
  double int_D = hD->Integral();
  hD->Scale(1./int_D);
  double int_MC1 = hMC1->Integral();
  hMC1->Scale(1./int_MC1);
  double int_MC2 = hMC2->Integral();
  hMC2->Scale(1./int_MC2);
  hD->GetYaxis()->SetTitle("Entries/1.00");
  hD->GetYaxis()->SetRangeUser(0,0.2);

  TLegend * leg = new TLegend(0.2,0.5,0.4,0.8);
  leg->SetTextSize(0.042);
  leg->SetFillColor(0);
  leg->AddEntry(hD,"Data 2012","lp");
  leg->AddEntry(hMC1,"MC12a with reweighting","lp");
  leg->AddEntry(hMC2,"MC12a without reweighting","lp");

  TCanvas * c1 = new TCanvas("c1","c1",800,600);
  c1->cd();
  hD->Draw("PE");
  hMC1->Draw("samePE");
  hMC2->Draw("samePE");
  leg->Draw("same");
  return c1;

}
