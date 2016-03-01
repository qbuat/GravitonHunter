#include "PersoResults.h"
#include "ToolsCommons.h"
///////////////////////////////////////////////////
Results::Results(TString filename,TString treename)
///////////////////////////////////////////////////
{
  if (filename=="") {
    std::cout << "No input file specified "<<std::endl;
  }
  if (treename=="") {
    std::cout << "No input tree specified"<<std::endl;
  }
  
  Init(filename,treename);
}
///////////////////////////////////
Results::Results()
///////////////////////////////////
{
  //Default constructor
  m_nentries=0;
  m_file=new TFile();
  m_tree=new TTree();

}
/////////////////////////////////////
Results::~Results()
////////////////////////////////////
{
  delete m_Rd;
  delete m_tree;
  delete m_file;
}


////////////////////////////////////////////////////
void Results::Init(TString filename,TString treename)
////////////////////////////////////////////////////
{
  m_file=new TFile(filename,"read");
  m_tree=(TTree*)m_file->Get(treename);
  m_nentries=m_tree->GetEntriesFast();
  m_Rd=new TreeReader(m_tree);
}
//////////////////////////////
void Results::DrawKinematics()
//////////////////////////////
{
  TH1F* hmgg_log   = new TH1F("hmgg_log","mgg",Commons::nBins,Commons::binning);
  TH1F* hleadpt    = new TH1F("hleadpt","leadpt",100,15,900);
  TH1F* hsubleadpt = new TH1F("hsubleadpt","subleadpt",100,15,900);
  TH1F* hdphi      = new TH1F("hdphi","dphi",60,-0.5,3.5);
  TH1F* hdeta      = new TH1F("hdeta","deta",50,-5,5);
  std::cout << "RunNumber" << " | " << "EventNumber" << " | " << "LumiBlock" << " | " << "mgg" << std::endl;

  for(int entry = 0 ; entry < m_nentries ; entry++){
    m_Rd->GetEntry(entry);
    int RunNumber      = (int)m_Rd->GetVariable("RunNumber");
    int EventNumber    = (int)m_Rd->GetVariable("EventNumber");
    int LumiBlock      = (int)m_Rd->GetVariable("LumiBlock");
    double mgg         = 1000*m_Rd->GetVariable("mgg_old");
    double Lead_pT     = m_Rd->GetVariable("Lead_pT"); 
    double Lead_eta    = m_Rd->GetVariable("Lead_eta"); 
    double SubLead_pT  = m_Rd->GetVariable("SubLead_pT"); 
    double SubLead_eta = m_Rd->GetVariable("SubLead_eta"); 
    double Dphi       = m_Rd->GetVariable("Deltaphi");
    std::cout << RunNumber << " | " << EventNumber << " | " << LumiBlock << " | "<< mgg << std::endl;
    hmgg_log->Fill(mgg);
    if(mgg>140){
      hleadpt->Fill(Lead_pT);
      hsubleadpt->Fill(SubLead_pT);
      hdeta->Fill(Lead_eta-SubLead_eta);
      hdphi->Fill(Dphi);
    }
  }

  TCanvas *c1 = new TCanvas("c1","mgg",800,800);
  c1->SetLogx();
  c1->SetLogy();
  c1->cd();
  hmgg_log->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  hmgg_log->GetXaxis()->SetRangeUser(140,2100);
  hmgg_log->GetYaxis()->SetTitle("Events/ bin");
  hmgg_log->GetYaxis()->SetRangeUser(0.0001,3000);
  hmgg_log->Draw("PE");

  TCanvas *c1_1 = new TCanvas("c1_1","mgg_1",800,800);
  c1_1->SetLogx();
  c1_1->SetLogy();
  c1_1->cd();
  TGraphAsymmErrors* gmgg = new TGraphAsymmErrors();
  for (int ibin = 1 ; ibin < hmgg_log->GetNbinsX()+1;ibin++){
    double value = hmgg_log->GetBinContent(ibin);
    if( value!=0){
      double y1 = value + 1.0;
      double d1 = 1.0 - 1.0/(9.0*y1) + 1.0/(3*sqrt(y1));
      double error_poisson_up = y1*d1*d1*d1 - value;
      double y2 = value;
      double d2 = 1.0 - 1.0/(9.0*y2) - 1.0/(3.0*sqrt(y2));
      double error_poisson_down = value - y2*d2*d2*d2;
      gmgg->SetPoint(ibin-1, hmgg_log->GetBinCenter(ibin), value);
      gmgg->SetPointError(ibin-1, 0, 0, error_poisson_down, error_poisson_up);
    }      
    std::cout << "bin n " << ibin << " : mgg = " << hmgg_log->GetBinCenter(ibin) << std::endl;
  }
  gmgg->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  gmgg->GetXaxis()->SetRangeUser(140,2100);
  gmgg->GetYaxis()->SetTitle("Events/ bin");
  gmgg->GetYaxis()->SetRangeUser(0.0001,3000);
  gmgg->Draw("AP");
  TCanvas *c2 = new TCanvas("c2","leadpt",800,800);
  c2->SetLogy();
  c2->cd();
  hleadpt->GetXaxis()->SetTitle("p_{T}^{lead #gamma} [GeV]");
  hleadpt->GetYaxis()->SetTitle("Events/ bin");
  hleadpt->Draw("PE");
  TCanvas *c3 = new TCanvas("c3","subleadpt",800,800);
  c3->SetLogy();
  c3->cd();
  hsubleadpt->GetXaxis()->SetTitle("p_{T}^{sublead #gamma} [GeV]");
  hsubleadpt->GetYaxis()->SetTitle("Events/ bin");
  hsubleadpt->Draw("PE");
  TCanvas *c4 = new TCanvas("c4","deta",800,800);
  c4->SetLogy();
  c4->cd();
  hdeta->GetXaxis()->SetTitle("#Delta #eta");
  hdeta->GetYaxis()->SetTitle("Events/ bin");
  hdeta->Draw("PE");
  TCanvas *c5 = new TCanvas("c5","dphi",800,800);
  c5->SetLogy();
  c5->cd();
  hdphi->GetXaxis()->SetTitle("#Delta #phi");
  hdphi->GetYaxis()->SetTitle("Events/ bin");
  hdphi->Draw("PE");


}
////////////////////////////////////////
void Results::SimpleNullHypothesisTest(int NPE)
////////////////////////////////////////
{
  double mean = 3.22502857408520427;
  //double mean = 3.89373670246682924;
  double Ndata = 6;
  TRandom3* Randomizer = new TRandom3();
  int Nbelow = 0;
  int Nup = 0;

  for(int ipe=0;ipe<NPE;ipe++){
    double temp = Randomizer->PoissonD(mean);
    if(temp > Ndata)
      Nup++;
    if(temp <= Ndata)
      Nbelow++;
  }
  std::cout << "NPE    : " << NPE    << std::endl;
  std::cout << "up     : " << Nup    << std::endl;
  std::cout << "below  : " << Nbelow << std::endl;
  std::cout << "Pup :  : " << 100*(double)Nup/(double)NPE << std::endl;
}
///////////////////////////////////////////
void Results::CreateInputForBAT_RescaleBK6()
///////////////////////////////////////////
{
  TFile *fD = new TFile("/sps/atlas/q/qbuat/OUTDS/GravitonAnalysis/data.root","read");
  TH1F *hh_data = (TH1F*)fD->Get("hmgg_log");
  hh_data->GetXaxis()->SetRangeUser(450,3000);
  hh_data->SetName("hh_data"); 
  hh_data->SetTitle("hh_data"); 
  double scalefactor = 2.36019;
  TFile *fBkg = new TFile("/afs/in2p3.fr/home/q/qbuat/BAT/inputfiles/Bkg_gg_Templates_IsoTemplates_B_K6_newCaloC.root","read");
  TH1D* bkg_total_gg_old        = (TH1D*)fBkg->Get("bkg_total_gg");
  TH1D* bkg_reducible_gg_old    = (TH1D*)fBkg->Get("bkg_reducible_gg");
  TH1D* bkg_irreducible_gg_old  = (TH1D*)fBkg->Get("bkg_irreducible_gg");
  TH1D* bkg_total_syst_gg       = (TH1D*)fBkg->Get("bkg_total_syst_gg");
  bkg_total_syst_gg->SetName("bkg_total_syst_gg");
  TH1D* bkg_reducible_syst_gg   = (TH1D*)fBkg->Get("bkg_reducible_syst_gg");
  bkg_reducible_syst_gg->SetName("bkg_reducible_syst_gg");
  TH1D* bkg_irreducible_syst_gg = (TH1D*)fBkg->Get("bkg_irreducible_syst_gg");
  bkg_irreducible_syst_gg->SetName("bkg_irreducible_syst_gg");
  TH1D* bkg_purity_syst_gg      = (TH1D*)fBkg->Get("bkg_purity_syst_gg");
  bkg_purity_syst_gg->SetName("bkg_purity_syst_gg");

  TH1D*  bkg_total_gg       = (TH1D*)bkg_total_gg_old->Clone();
  bkg_total_gg->SetName("bkg_total_gg");
  TH1D*  bkg_reducible_gg   = (TH1D*)bkg_reducible_gg_old->Clone();
  bkg_reducible_gg->SetName("bkg_reducible_gg");
  TH1D*  bkg_irreducible_gg = (TH1D*)bkg_irreducible_gg_old->Clone();
  bkg_irreducible_gg->SetName("bkg_irreducible_gg");
  bkg_total_gg->Sumw2();
  bkg_total_gg->Scale(scalefactor);
  bkg_reducible_gg->Sumw2();
  bkg_reducible_gg->Scale(scalefactor);
  bkg_irreducible_gg->Sumw2();
  bkg_irreducible_gg->Scale(scalefactor);

  TObjArray Hlist(0);
  Hlist.Add(bkg_total_gg);
  Hlist.Add(bkg_reducible_gg);
  Hlist.Add(bkg_irreducible_gg);
  Hlist.Add(bkg_total_syst_gg);
  Hlist.Add(bkg_reducible_syst_gg);
  Hlist.Add(bkg_irreducible_syst_gg);
  Hlist.Add(bkg_purity_syst_gg);
  Hlist.Add(hh_data);

  TFile fout("Bkg_gg_Template_B-M.root","RECREATE");
  Hlist.Write();
  fout.Close();

}
//////////////////////////////////////
void Results::CreateInputForEvanCode()
//////////////////////////////////////
{
  int Nbins    = 300;
  double Xmin  = 0.;
  double Xmax  = 3000.;
  TH1F * h_datafinal_mass = new TH1F("h_datafinal_mass","h_datafinal_mass",Nbins,Xmin,Xmax);
  TH1F * h_n_fake_mass = new TH1F("h_n_fake_mass","h_n_fake_mass",Nbins,Xmin,Xmax);
  TH1F * h_double_fake_mass = new TH1F("h_double_fake_mass","h_double_fake_mass",Nbins,Xmin,Xmax);
  TH1F * h_leading_fake_mass = new TH1F("h_leading_fake_mass","h_leading_fake_mass",Nbins,Xmin,Xmax);
  TH1F * h_subleading_fake_mass = new TH1F("h_subleading_fake_mass","h_subleading_fake_mass",Nbins,Xmin,Xmax);
  const int nBins_log = 100;
  double xmin_GeV     = 70;
  double xmax_GeV     = 3000;
  double logxmin_GeV  = log10(xmin_GeV);
  double logxmax_GeV  = log10(xmax_GeV);
  double binwidth_GeV = (logxmax_GeV-logxmin_GeV)/double(nBins_log);
  double binning_log[nBins_log+1];
  binning_log[0] = xmin_GeV;
  for (int i=1 ;i<=nBins_log ;i++){
    binning_log[i] = pow(10,logxmin_GeV+i*binwidth_GeV);
  }
  TH1F* h_datafinal_log_mass = new TH1F("h_datafinal_log_mass","h_datafinal_log_mass",nBins_log,binning_log);
  TH1F* h_n_fake_log_mass = new TH1F("h_n_fake_log_mass","h_n_fake_log_mass",nBins_log,binning_log);
  TH1F* h_double_fake_log_mass = new TH1F("h_double_fake_log_mass","h_double_fake_log_mass",nBins_log,binning_log);
  TH1F* h_leading_fake_log_mass = new TH1F("h_leading_fake_log_mass","h_leading_fake_log_mass",nBins_log,binning_log);
  TH1F* h_subleading_fake_log_mass = new TH1F("h_subleading_fake_log_mass","h_subleading_fake_log_mass",nBins_log,binning_log);

  for( int jentry=0 ; jentry<m_nentries ; jentry++){
    m_Rd->GetEntry(jentry);
    int RunNumber  = (int)m_Rd->GetVariable("RunNumber"); 
    double Iso_L   = m_Rd->GetVariable("Iso_L");
    double Iso_SL  = m_Rd->GetVariable("Iso_SL");
    int IsTight_L  = (int)m_Rd->GetVariable("IsTight_L"); 
    int IsTight_SL = (int)m_Rd->GetVariable("IsTight_SL"); 
    double mgg     = m_Rd->GetVariable("mgg");
    std::cout << jentry<< " / " << RunNumber << std::endl;
    // if( RunNumber <  188902) continue;
    if (Iso_L <5 && Iso_SL < 5){
      if( IsTight_L!=1 && IsTight_SL!=1 ){
	h_double_fake_mass->Fill(mgg);
	h_double_fake_log_mass->Fill(mgg);
      }
      else if( IsTight_L!=1 && IsTight_SL==1 ){
	h_leading_fake_mass->Fill(mgg);
	h_leading_fake_log_mass->Fill(mgg);
      }
      else if( IsTight_L==1 && IsTight_SL!=1 ){
	h_subleading_fake_mass->Fill(mgg);
	h_subleading_fake_log_mass->Fill(mgg);
      }
      else if( IsTight_L==1 && IsTight_SL==1 ){
	h_datafinal_mass->Fill(mgg);
	h_datafinal_log_mass->Fill(mgg);
      }
      if( IsTight_L!=1 || IsTight_SL!=1){
	h_n_fake_mass->Fill(mgg);
	h_n_fake_log_mass->Fill(mgg);
      }
    }
  }

  std::cout << "TOTO" << std::endl;
  TObjArray Hlist(0);
  Hlist.Add(h_n_fake_mass);
  Hlist.Add(h_double_fake_mass);
  Hlist.Add(h_leading_fake_mass);
  Hlist.Add(h_subleading_fake_mass);
  Hlist.Add(h_n_fake_log_mass);
  Hlist.Add(h_double_fake_log_mass);
  Hlist.Add(h_leading_fake_log_mass);
  Hlist.Add(h_subleading_fake_log_mass);
  Hlist.Add(h_datafinal_mass);
  Hlist.Add(h_datafinal_log_mass);
  std::cout << "TOTO2" << std::endl;
  TFile *fout = new TFile("hists_reducible.root","RECREATE");
  Hlist.Write();
  std::cout << "TOTO3" << std::endl;
  fout->Close();
  std::cout << "TOTO4" << std::endl;

}
//////////////////////////////////////
void Results::CloneInputTree()
//////////////////////////////////////
{
  // DQ::SetXMLFile("/afs/in2p3.fr/home/q/qbuat/GravToGamGam/GravInvMass/run/data11_7TeV.periodAllYear_DetStatus-v32-pro09_CoolRunQuery-00-04-00_Eg_standard.xml");
  std::set<unsigned long long> B_K6_events;
  GravitonAnalysis GA;
  B_K6_events = GA.build_eeset("list_rel16.txt");
  TFile * f = new TFile("/sps/atlas/q/qbuat/OUTDS/GravitonAnalysis/IsoForFitStudyLooseprime0_BK6_b.root","RECREATE");
  TTree * clonetree =(TTree*)m_tree->CloneTree(0);

  for( int jentry=0 ; jentry<m_nentries ; jentry++){
    m_Rd->GetEntry(jentry);
    m_tree->GetEntry(jentry);
    if ( (int)m_Rd->GetVariable("RunNumber")>187815 ) continue;
    // int RunNumber    = (int)m_Rd->GetVariable("RunNumber"); 
    // int EventNumber  = (int)m_Rd->GetVariable("EventNumber"); 
    // int LumiBlock    = (int)m_Rd->GetVariable("LumiBlock"); 
    // if( GA.EventInZprimme(B_K6_events,RunNumber,EventNumber) ) continue;
    // if( !DQ::PassRunLB(RunNumber,LumiBlock) ) continue;
    // if( RunNumber >  187815) continue;
    clonetree->Fill();
  }
  clonetree->Write();
  f->Close();
}

//////////////////////////////////////
void Results::CreateListOfEvents()
//////////////////////////////////////
{
  std::ofstream file_list;
  file_list.open("list_rel16.txt");
  for( int jentry=0 ; jentry<m_nentries ; jentry++){
    m_Rd->GetEntry(jentry);
    m_tree->GetEntry(jentry);
    int RunNumber  = (int)m_Rd->GetVariable("RunNumber"); 
    int EventNumber  = (int)m_Rd->GetVariable("EventNumber"); 
    int LumiBlock  = (int)m_Rd->GetVariable("LumiBlock"); 
    file_list << RunNumber << " " << EventNumber << " "<< LumiBlock << std::endl;
  }
  file_list.close();
}
///////////////////////////////////////
void Results::CompareData_rel16_rel17()
//////////////////////////////////////
{
  TFile *f16       = new TFile("/afs/in2p3.fr/home/q/qbuat/GravToGamGam/GravInvMass/run/IsoForFitStudyLooseprime0.root","read");
  TFile *f17       = new TFile("/sps/atlas/q/qbuat/OUTDS/GravitonAnalysis/IsoForFitStudyLooseprime0_BK6rel16events.root","read");
  TTree *t16       = (TTree*)f16->Get("tree");
  TTree *t17       = (TTree*)f17->Get("tree");
  TreeReader *Rd16 = new TreeReader(t16);
  TreeReader *Rd17 = new TreeReader(t17);
  int nentries16   = t16->GetEntries();
  int nentries17   = t17->GetEntries();
  TH1F* hmgg_16_n   = new TH1F("hmgg_16_n","mgg",Commons::nBins,Commons::binning);
  TH1F* hmgg_17_n   = new TH1F("hmgg_17_n","mgg",Commons::nBins,Commons::binning);
  TH1F* hmgg_16_l   = new TH1F("hmgg_16_l","mgg",Commons::nBins,Commons::binning);
  TH1F* hmgg_17_l   = new TH1F("hmgg_17_l","mgg",Commons::nBins,Commons::binning);
  TH1F* hmgg_16_s   = new TH1F("hmgg_16_s","mgg",Commons::nBins,Commons::binning);
  TH1F* hmgg_17_s   = new TH1F("hmgg_17_s","mgg",Commons::nBins,Commons::binning);
  TH1F* hmgg_16_d   = new TH1F("hmgg_16_d","mgg",Commons::nBins,Commons::binning);
  TH1F* hmgg_17_d   = new TH1F("hmgg_17_d","mgg",Commons::nBins,Commons::binning);
  for (int jentry=0 ;jentry< nentries16;jentry++){
    Rd16->GetEntry(jentry);
    int IsTight_L  = (int)Rd16->GetVariable("IsTight_L");
    int IsTight_SL = (int)Rd16->GetVariable("IsTight_SL");
    double Iso_L   = Rd16->GetVariable("Iso_L");
    double Iso_SL  = Rd16->GetVariable("Iso_SL");
    double mgg     = Rd16->GetVariable("mgg");
    if( Iso_L<5 && Iso_SL<5){
      if(IsTight_L==0 || IsTight_SL==0){
	hmgg_16_n->Fill(mgg);
      }
      if(IsTight_L==1 && IsTight_SL==0){
	hmgg_16_s->Fill(mgg);
      }
      if(IsTight_L==0 && IsTight_SL==1){
	hmgg_16_l->Fill(mgg);
      }
      if(IsTight_L==0 && IsTight_SL==0){
	hmgg_16_d->Fill(mgg);
      }
    }
  }
  for (int jentry=0 ;jentry< nentries17;jentry++){
    Rd17->GetEntry(jentry);
    if ( (int)Rd17->GetVariable("RunNumber")>187815 ) continue;
    int IsTight_L  = (int)Rd17->GetVariable("IsTight_L");
    int IsTight_SL = (int)Rd17->GetVariable("IsTight_SL");
    double Iso_L   = Rd17->GetVariable("Iso_L");
    double Iso_SL  = Rd17->GetVariable("Iso_SL");
    double mgg     = Rd17->GetVariable("mgg");
    if( Iso_L<5 && Iso_SL<5){
      if(IsTight_L==0 || IsTight_SL==0){
	hmgg_17_n->Fill(mgg);
      }
      if(IsTight_L==1 && IsTight_SL==0){
	hmgg_17_s->Fill(mgg);
      }
      if(IsTight_L==0 && IsTight_SL==1){
	hmgg_17_l->Fill(mgg);
      }
      if(IsTight_L==0 && IsTight_SL==0){
	hmgg_17_d->Fill(mgg);
      }
    }
  }
  hmgg_16_n->Sumw2();
  hmgg_17_n->Sumw2();
  hmgg_16_l->Sumw2();
  hmgg_17_l->Sumw2();
  hmgg_16_s->Sumw2();
  hmgg_17_s->Sumw2();
  hmgg_16_d->Sumw2();
  hmgg_17_d->Sumw2();

  // hmgg_16_n->SetNormFactor(1);
  // hmgg_17_n->SetNormFactor(1);
  // hmgg_16_l->SetNormFactor(1);
  // hmgg_17_l->SetNormFactor(1);
  // hmgg_16_s->SetNormFactor(1);
  // hmgg_17_s->SetNormFactor(1);
  // hmgg_16_d->SetNormFactor(1);
  // hmgg_17_d->SetNormFactor(1);

  TH1F* hr_n =(TH1F*)hmgg_17_n->Clone();
  hr_n->Divide(hmgg_16_n);
  hr_n->Sumw2();
  TH1F* hr_l =(TH1F*)hmgg_17_l->Clone();
  hr_l->Divide(hmgg_16_l);
  hr_l->Sumw2();
  TH1F* hr_s =(TH1F*)hmgg_17_s->Clone();
  hr_s->Divide(hmgg_16_s);
  hr_s->Sumw2();
  TH1F* hr_d =(TH1F*)hmgg_17_d->Clone();
  hr_d->Divide(hmgg_16_d);
  hr_d->Sumw2();

  int rel17color = 2;
  // int rel16color = 1;
  TCanvas *c_n = new TCanvas("c_n","c_n",800,800);
  c_n->Divide(1,2);
  c_n->cd(1);

  hmgg_17_n->SetLineColor(rel17color);
  hmgg_17_n->SetMarkerColor(rel17color);
  hmgg_17_n->Draw("PE");
  hmgg_16_n->Draw("samePE");
  c_n->cd(2);
  hr_n->Draw("PE");
  TCanvas *c_l = new TCanvas("c_l","c_l",800,800);
  c_l->Divide(1,2);
  c_l->cd(1);
  hmgg_17_l->SetLineColor(rel17color);
  hmgg_17_l->SetMarkerColor(rel17color);
  hmgg_17_l->Draw("PE");
  hmgg_16_l->Draw("samePE");
  c_l->cd(2);
  hr_l->Draw("PE");
  TCanvas *c_s = new TCanvas("c_s","c_s",800,800);
  c_s->Divide(1,2);
  c_s->cd(1);
  hmgg_16_s->Draw("PE");
  hmgg_17_s->SetLineColor(rel17color);
  hmgg_17_s->SetMarkerColor(rel17color);
  hmgg_17_s->Draw("PE");
  hmgg_16_s->Draw("samePE");
  c_s->cd(2);
  hr_s->Draw("PE");
  TCanvas *c_d = new TCanvas("c_d","c_d",800,800);
  c_d->Divide(1,2);
  c_d->cd(1);
  hmgg_17_d->SetLineColor(rel17color);
  hmgg_17_d->SetMarkerColor(rel17color);
  hmgg_17_d->Draw("PE");
  hmgg_16_d->Draw("samePE");
  c_d->cd(2);
  hr_d->Draw("PE");

}
//////////////////////////////////////
void Results::CompareListOfEvents()
//////////////////////////////////////
{

  std::ifstream rel16;
  rel16.open("list_rel16.txt");
  std::string in16;
  int Nevt = 0;
  if (rel16.is_open()) {
    while (!rel16.eof()) {
      getline(rel16,in16);
      bool islineinrel17 = false;
      
      std::ifstream rel17;
      rel17.open("list_rel17.txt");
      std::string in17;
      if (rel17.is_open()) {
	while (!rel17.eof()) {
	  getline(rel17,in17);
	  // std::cout << in17 << std::endl;
	  if(in17==in16){
	    islineinrel17 = true;
	  }
	}
      }
      if (islineinrel17==false){
	Nevt++;
	std::cout << in16 << std::endl;
      }
    }
  }
  std::cout << "Number of events = " << Nevt << std::endl;
}
/////////////////////////////////////
void Results::CreateDataInvMassHist()
/////////////////////////////////////
{
  TH1D* hD = new TH1D("hmgg_data","hmgg_data",
		      Commons::nBins,Commons::binning);

  TString filename = "/sps/atlas/q/qbuat/OUTDS/GravitonAnalysis/data.root";
  TFile *f0 = new TFile(filename,"read");
  TTree *tree = (TTree*)f0->Get("tree");
  int nentries = (int)tree->GetEntriesFast();
  TreeReader *Rd = new TreeReader(tree);
  for(int entry=0;entry<nentries;entry++){
    Rd->GetEntry(entry);
    double mgg = Rd->GetVariable("mgg_old");
    hD->Fill(mgg);
  }

  TFile fout("datahist_B_M.root","RECREATE");
  fout.Add(hD);
  fout.Write();
  fout.Close();
}
/////////////////////////////////////////
void Results::CreateGamGamHist_Fromxsec()
////////////////////////////////////////
{

  TString path = "/afs/in2p3.fr/home/q/qbuat/GravToGamGam/GravInvMass-rel17-MC11/run/rootfiles/";
  // TFile *f0 = new TFile(path+"input_Pythiagamgam15_irr.root","read");
  // TFile *f0_fix = new TFile(path+"input_Pythiagamgam15_irr_fix.root","read");
  // TFile *f1 = new TFile(path+"input_Pythiagamgam15_highmass_irr.root","read");
  // TFile *f2 = new TFile(path+"input_Pythiagamgam15_M800_irr.root","read");
  // TFile *f3 = new TFile(path+"input_Pythiagamgam15_M1500_irr.root","read");

  TFile *f0 = new TFile(path+"input_Pythiagamgam15_irr_linearbinning.root","read");
  TFile *f0_fix = new TFile(path+"input_Pythiagamgam15_irr_linearbinning_fix.root","read");
  TFile *f1 = new TFile(path+"input_Pythiagamgam15_highmass_irr_linearbinning.root","read");
  TFile *f2 = new TFile(path+"input_Pythiagamgam15_M800_irr_linearbinning.root","read");
  TFile *f3 = new TFile(path+"input_Pythiagamgam15_M1500_irr_linearbinning.root","read");

  TH1D* h0_r  = (TH1D*)f0->Get("hmgg_log2_truthcut_nlo");
  TH1D* h1_r  = (TH1D*)f1->Get("hmgg_log2_truthcut_nlo");
  TH1D* h2_r  = (TH1D*)f2->Get("hmgg_log2_truthcut_nlo");
  TH1D* h3_r  = (TH1D*)f3->Get("hmgg_log2_truthcut_nlo");

  // TH1D* h0_r  = (TH1D*)f0->Get("hmgg_log2_truthcut");
  // TH1D* h1_r  = (TH1D*)f1->Get("hmgg_log2_truthcut");
  // TH1D* h2_r  = (TH1D*)f2->Get("hmgg_log2_truthcut");
  // TH1D* h3_r  = (TH1D*)f3->Get("hmgg_log2_truthcut");

  TH1D* h0_fix  = (TH1D*)f0_fix->Get("htruthmgg_log_truthcut");
  TH1D* h0  = (TH1D*)f0->Get("htruthmgg_log_truthcut");
  TH1D* h1  = (TH1D*)f1->Get("htruthmgg_log_truthcut");
  TH1D* h2  = (TH1D*)f2->Get("htruthmgg_log_truthcut");
  TH1D* h3  = (TH1D*)f3->Get("htruthmgg_log_truthcut");

  h1->SetLineColor(2);
  h2->SetLineColor(3);
  h3->SetLineColor(4);

  h1_r->SetLineColor(2);
  h2_r->SetLineColor(3);
  h3_r->SetLineColor(4);

  std::cout << "h0 entries : " << (int)h0->GetEntries() << std::endl;
  std::cout << "h1 entries : " << (int)h1->GetEntries() << std::endl;
  std::cout << "h2 entries : " << (int)h2->GetEntries() << std::endl;
  std::cout << "h3 entries : " << (int)h3->GetEntries() << std::endl;

  // double Nevent_massbin[] = {52020 ,71676 ,33424,9039};//log binning
  double Nevent_massbin[] = {80460,71676,33424,9039};//linear binning

  // double Nsimevents[]     = {500000,200000,80000,20000};
  double Nsimevents[]     = {499950,199948,79999,20000};
  double xsec[]           = {8.6811e-01,3.3321e-03,2.8570e-05,1.5445e-06};
  double filtereff[]      = {0.13240,0.41824,0.30355,0.21329};
  double scale[4];
  for( int i=0; i<4;i++){
    scale[i] = Nevent_massbin[i]/Nsimevents[i]*xsec[i]*filtereff[i]*10000/Nsimevents[i];
  }

  
  
  // h0->Scale(scale[0]);
  h1->Scale(scale[1]);
  h2->Scale(scale[2]);
  h3->Scale(scale[3]);
  // double I0 = h0_fix->Integral(16,26);//log binning
  // double I1 = h1->Integral(16,26);//log binning
  double I0 = h0_fix->Integral(22,28);
  double I1 = h1->Integral(22,28);
  double scale_fix = I1/I0;
  std::cout << "HERE6 : " << scale_fix << std::endl;
  h0->Scale(scale_fix);


  h0_r->Sumw2();
  h1_r->Sumw2();
  h2_r->Sumw2();
  h3_r->Sumw2();

  h0_r->Scale(scale_fix);
  // h0_r->Scale(scale[0]);
  h1_r->Scale(scale[1]);
  h2_r->Scale(scale[2]);
  h3_r->Scale(scale[3]);



  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  c1->cd();
  c1->SetLogy();
  h0->GetYaxis()->SetRangeUser(0.0000001,2);
  h0->Draw("HIST");
  h1->Draw("sameHIST");
  h2->Draw("sameHIST");
  h3->Draw("sameHIST");

  TCanvas *c2 = new TCanvas("c2","c2",800,800);
  c2->cd();
  c2->SetLogy();
  h0_r->GetYaxis()->SetRangeUser(0.0000001,2);
  h0_r->Draw("HIST");
  h1_r->Draw("sameHIST");
  h2_r->Draw("sameHIST");
  h3_r->Draw("sameHIST");

  TH1D* h_r = (TH1D*)h0_r->Clone("h_r");
  h_r->Add(h1_r);
  h_r->Add(h2_r);
  h_r->Add(h3_r);
  h_r->Sumw2();
  TCanvas *c3 = new TCanvas("c3","c3",800,800);
  c3->cd();
  h_r->Draw("HIST");

  TFile *f = new TFile("irreducible_shape.root","RECREATE");
  h_r->Write("irreducible_shape");
  f->Close();
}



/////////////////////////
void Results::LCvsJC()
/////////////////////////
{

  int NkilledbyLC=0;
  int NkilledbyLCnew=0;
  int NkilledbyJC=0;
  int NpassLC=0;
  int NpassJC=0;
  int NkilledbyLCtight=0;
  int NkilledbyLCnewtight=0;
  int NkilledbyJCtight=0;
  int NpassLCtight=0;
  int NpassJCtight=0;
  int NpassBaseline=0;
  int NpassBaselineTight=0;
  TH1F *hpt_goodLC=new TH1F("hpt_goodLC","hpt_goodLC",10,0,1000);
  TH1F *hpt_badLC=new TH1F("hpt_badLC","hpt_badLC",10,0,1000);

  TH1F *hpt_goodLC2=new TH1F("hpt_goodLC2","hpt_goodLC2",10,0,1000);
  TH1F *hpt_badLC2=new TH1F("hpt_badLC2","hpt_badLC2",10,0,1000);

  TH1F *hpt_goodLC_T=new TH1F("hpt_goodLC_T","hpt_goodLC_T",10,0,1000);
  TH1F *hpt_badLC_T=new TH1F("hpt_badLC_T","hpt_badLC_T",10,0,1000);

  TH1F *heta_goodLC=new TH1F("heta_goodLC","heta_goodLC",14,-2.8,2.8);
  TH1F *heta_badLC=new TH1F("heta_badLC","heta_badLC",14,-2.8,2.8);

  TH1F *hpt_goodJC=new TH1F("hpt_goodJC","hpt_goodJC",10,0,1000);
  TH1F *hpt_badJC=new TH1F("hpt_badJC","hpt_badJC",10,0,1000);

  TH1F *heta_goodJC=new TH1F("heta_goodJC","heta_goodJC",14,-2.8,2.8);
  TH1F *heta_badJC=new TH1F("heta_badJC","heta_badJC",14,-2.8,2.8);


  //== LOOP OVER PHOTONS ===// 
  for(int entry=0;entry<m_nentries;entry++){
    m_Rd->GetEntry(entry);
    if((int)m_Rd->GetVariable("larError")>1) continue;
    NpassBaseline++;
    if((int)m_Rd->GetVariable("phistight")==1){
      NpassBaselineTight++;
    }
    if((int)m_Rd->GetVariable("phPassLC")==0){
      NkilledbyLC++;
      // if( m_Rd->GetVariable("phreta")>0.98 ||
      // 	  m_Rd->GetVariable("phrphi")>1.){
      if( (int)m_Rd->GetVariable("phPassLC2")==0){
	NkilledbyLCnew++;
      }
      hpt_badLC->Fill( m_Rd->GetVariable("phclpt")/1000. );
      heta_badLC->Fill( m_Rd->GetVariable("phetas2") );
      if((int)m_Rd->GetVariable("phistight")==1){
	hpt_badLC_T->Fill( m_Rd->GetVariable("phclpt")/1000. );
	NkilledbyLCtight++;
	if( m_Rd->GetVariable("phreta")>0.98 ||
	    m_Rd->GetVariable("phrphi")>1.){
	  NkilledbyLCnewtight++;
	}
      }
    }else{
      NpassLC++;
      hpt_goodLC->Fill( m_Rd->GetVariable("phclpt")/1000. );
      heta_goodLC->Fill( m_Rd->GetVariable("phetas2") );
      if( (int)m_Rd->GetVariable("phistight") ==1) {
	NpassLCtight++;
	hpt_goodLC_T->Fill( m_Rd->GetVariable("phclpt")/1000. );
      }
    }
    if((int)m_Rd->GetVariable("phPassLC2")==0)
      hpt_badLC2->Fill( m_Rd->GetVariable("phclpt")/1000. );
    else
      hpt_goodLC2->Fill( m_Rd->GetVariable("phclpt")/1000. );
      
    if((int)m_Rd->GetVariable("phPassJC")==0){
      NkilledbyJC++;
      hpt_badJC->Fill( m_Rd->GetVariable("phclpt")/1000. );
      heta_badJC->Fill( m_Rd->GetVariable("phetas2") );
      if((int)m_Rd->GetVariable("phistight")==1){
	NkilledbyJCtight++;
      }
    }else{
      NpassJC++;
      hpt_goodJC->Fill( m_Rd->GetVariable("phclpt")/1000. );
      heta_goodJC->Fill( m_Rd->GetVariable("phetas2") );
     if((int)m_Rd->GetVariable("phistight")==1){
	NpassJCtight++;
      }
    }
  }
  

  ////End of LOOP over photons //////////
  std::cout << "======= Baseline Selection ====" << std::endl;
  std::cout << "Ntot = "<< NpassBaseline << std::endl;
  std::cout << "LAr Cleaning new"<< std::endl;
  std::cout << "Killed by LArCleaning new: "<< NkilledbyLCnew << " --- "<< 100*efficiency(NkilledbyLCnew,NpassBaseline)<< std::endl;
  std::cout << "LAr Cleaning"<< std::endl;  
  std::cout << "Passing LArCleaning   : "<< NpassLC << std::endl;
  std::cout << "Killed by LArCleaning : "<< NkilledbyLC <<" --- "<< 100*efficiency(NkilledbyLC,NpassBaseline)<< std::endl;
  std::cout << "Jet Cleaning"<< std::endl;  
  std::cout << "Passing JetCleaning   : "<< NpassJC << std::endl;
  std::cout << "Killed by JetCleaning : "<< NkilledbyJC <<" --- "<< 100*efficiency(NkilledbyJC,NpassBaseline)<<std::endl;
  std::cout << "======= Baseline Selection +Tight ====" << std::endl;
  std::cout << "Ntot = "<< NpassBaselineTight << std::endl;
  std::cout << "LAr Cleaning new"<< std::endl;
  std::cout << "Killed by LArCleaning new: "<< NkilledbyLCnewtight <<" --- "<< 100*efficiency(NkilledbyLCnewtight,NpassBaselineTight) << std::endl;
  std::cout << "LAr Cleaning"<< std::endl;  
  std::cout << "Passing LArCleaning   : "<< NpassLCtight << std::endl;
  std::cout << "Killed by LArCleaning : "<< NkilledbyLCtight <<" --- " << 100*efficiency(NkilledbyLCtight,NpassBaselineTight)<< std::endl;
  std::cout << "Jet Cleaning"<< std::endl;  
  std::cout << "Passing JetCleaning   : "<< NpassJCtight << std::endl;
  std::cout << "Killed by JetCleaning : "<< NkilledbyJCtight <<" --- " << 100*efficiency(NkilledbyJCtight,NpassBaselineTight)<< std::endl;

  double xpt[50];
  double e_xpt[50];
  double ypt[50];
  double e_ypt[50];
  int npt=0;
  int jpt=0;
  double xpta[50];
  double e_xpta[50];
  double ypta[50];
  double e_ypta[50];
  int npta=0;
  int jpta=0;
  double xptb[50];
  double e_xptb[50];
  double yptb[50];
  double e_yptb[50];
  int nptb=0;
  int jptb=0;
  double xpt_T[50];
  double e_xpt_T[50];
  double ypt_T[50];
  double e_ypt_T[50];
  int npt_T=0;
  int jpt_T=0;

  for(int i=0;i<50;i++){
    xpt[i]=0;
    e_xpt[i]=0;
    ypt[i]=0;
    e_ypt[i]=0;
    xpta[i]=0;
    e_xpta[i]=0;
    ypta[i]=0;
    e_ypta[i]=0;
    xptb[i]=0;
    e_xptb[i]=0;
    yptb[i]=0;
    e_yptb[i]=0;
    xpt_T[i]=0;
    e_xpt_T[i]=0;
    ypt_T[i]=0;
    e_ypt_T[i]=0;
  }


  for(int i=0; i<= 25 ; i++){
    if(hpt_badLC->GetBinContent(i)  !=0 &&
       hpt_goodLC->GetBinContent(i) !=0 ){
      double eff=100*efficiency( hpt_badLC->GetBinContent(i),hpt_goodLC->GetBinContent(i)+hpt_badLC->GetBinContent(i) );
      double sig_eff=100*sigma_efficiency( hpt_badLC->GetBinContent(i),hpt_goodLC->GetBinContent(i)+hpt_badLC->GetBinContent(i) );
      jpt++;
      xpt[jpt]=hpt_goodLC->GetBinCenter(i);
      e_xpt[jpt]=hpt_goodLC->GetBinWidth(i)/2.;
      ypt[jpt]=eff;
      e_ypt[jpt]=sig_eff;
      npt++;
    }
    if(hpt_badLC2->GetBinContent(i)  !=0 &&
       hpt_goodLC2->GetBinContent(i) !=0 ){
      double eff=100*efficiency( hpt_badLC2->GetBinContent(i),hpt_goodLC2->GetBinContent(i)+hpt_badLC2->GetBinContent(i) );
      double sig_eff=100*sigma_efficiency( hpt_badLC2->GetBinContent(i),hpt_goodLC2->GetBinContent(i)+hpt_badLC2->GetBinContent(i) );
      jptb++;
      xptb[jpt]=hpt_goodLC2->GetBinCenter(i);
      e_xptb[jpt]=hpt_goodLC2->GetBinWidth(i)/2.;
      yptb[jpt]=eff;
      e_yptb[jpt]=sig_eff;
      nptb++;
    }
    if(hpt_badLC_T->GetBinContent(i)  !=0 &&
       hpt_goodLC_T->GetBinContent(i) !=0 ){
       double eff=100*efficiency( hpt_badLC_T->GetBinContent(i),hpt_goodLC_T->GetBinContent(i)+hpt_badLC_T->GetBinContent(i) );
      double sig_eff=100*sigma_efficiency( hpt_badLC_T->GetBinContent(i),hpt_goodLC_T->GetBinContent(i)+hpt_badLC_T->GetBinContent(i) );
      jpt_T++;
      xpt_T[jpt]=hpt_goodLC_T->GetBinCenter(i);
      e_xpt_T[jpt]=hpt_goodLC_T->GetBinWidth(i)/2.;
      ypt_T[jpt]=eff;
      e_ypt_T[jpt]=sig_eff;
      npt_T++;
    }
    if(hpt_badJC->GetBinContent(i)  !=0 &&
       hpt_goodJC->GetBinContent(i) !=0 ){
      double eff=100*efficiency( hpt_badJC->GetBinContent(i),hpt_goodJC->GetBinContent(i)+hpt_badJC->GetBinContent(i) );
      double sig_eff=100*sigma_efficiency( hpt_badJC->GetBinContent(i),hpt_goodJC->GetBinContent(i)+hpt_badJC->GetBinContent(i) );
      jpta++;
      xpta[jpta]=hpt_goodJC->GetBinCenter(i);
      e_xpta[jpta]=hpt_goodJC->GetBinWidth(i)/2.;
      ypta[jpta]=eff;
      e_ypta[jpta]=sig_eff;
      npta++;
    }
  }

  for( int i =0 ; i<npt ; i++){
    std::cout << "xpt["<<i<<"] = "<< xpt[i] << std::endl;
    std::cout << "ypt["<<i<<"] = "<< ypt[i] << std::endl;
  }
  for( int i =0 ; i<nptb ; i++){
    std::cout << "xptb["<<i<<"] = "<< xptb[i] << std::endl;
    std::cout << "yptb["<<i<<"] = "<< yptb[i] << std::endl;
  }
  TGraphErrors *grpt   = new TGraphErrors(npt,xpt,ypt,e_xpt,e_ypt);
  TGraphErrors *grpt_T = new TGraphErrors(npt_T,xpt_T,ypt_T,e_xpt_T,e_ypt_T);
  TGraphErrors *grpta  = new TGraphErrors(npta,xpta,ypta,e_xpta,e_ypta);
  TGraphErrors *grptb  = new TGraphErrors(nptb,xptb,yptb,e_xptb,e_yptb);

  grpt->GetXaxis()->SetTitle("p_{T} [GeV]");
  grpt->GetYaxis()->SetTitle("Killing Rate (%)");
  grpt->GetXaxis()->SetRangeUser(10,700);
  grpt->GetYaxis()->SetRangeUser(0.0001,4);
  grpt->SetLineColor(3);
  grpt->SetMarkerColor(3);
  grpt->SetMarkerStyle(22);

  grpta->GetXaxis()->SetTitle("p_{T} [GeV]");
  grpta->GetYaxis()->SetTitle("Killing Rate (%)");
  grpta->GetXaxis()->SetRangeUser(10,700);
  grpta->GetYaxis()->SetRangeUser(0.0001,4);
  grpta->SetLineColor(2);
  grpta->SetMarkerColor(2);
  grpta->SetMarkerStyle(24);

  grptb->GetXaxis()->SetTitle("p_{T} [GeV]");
  grptb->GetYaxis()->SetTitle("Killing Rate (%)");
  grptb->GetXaxis()->SetRangeUser(10,700);
  grptb->GetYaxis()->SetRangeUser(0.0001,4);
  grptb->SetLineColor(4);
  grptb->SetMarkerColor(4);
  grptb->SetMarkerStyle(25);

  grpt_T->GetXaxis()->SetTitle("p_{T} [GeV]");
  grpt_T->GetYaxis()->SetTitle("Killing Rate (%)");
  grpt_T->GetXaxis()->SetRangeUser(10,700);
  grpt_T->GetYaxis()->SetRangeUser(0.0001,4);
  grpt_T->SetLineColor(4);
  grpt_T->SetMarkerColor(4);

  TLegend *leg0 = new TLegend(0.2,0.8,0.5,0.9);
  leg0->SetFillColor(0);
  leg0->AddEntry(grpt,"killed by LArCleaning","lp");
  leg0->AddEntry(grptb,"killed by the PhotonCleaning","lp");
  // leg0->AddEntry(grpta,"killed by JetCleaning","lp");

  TCanvas *c1=new TCanvas("c1","c1",800,800);
  c1->cd();
  c1->SetLogy();
  grpt->Draw("AP");
  // grpta->Draw("sameP");
  grptb->Draw("sameP");
  leg0->Draw("same");
  // grpt_T->Draw("sameP");

  TCanvas *c1_bis=new TCanvas("c1_bis","c1_bis",800,800);
  c1_bis->cd();
  c1_bis->SetLogy();
  grpt->Draw("AP");
  grpt_T->Draw("sameP");

  double xeta[50];
  double e_xeta[50];
  double yeta[50];
  double e_yeta[50];
  int neta=0;
  int jeta=0;

  double xeta_a[50];
  double e_xeta_a[50];
  double yeta_a[50];
  double e_yeta_a[50];
  int neta_a=0;
  int jeta_a=0;

  for(int i=0;i<50;i++){
    xeta[i]=0;
    e_xeta[i]=0;
    yeta[i]=0;
    e_yeta[i]=0;
    xeta_a[i]=0;
    e_xeta_a[i]=0;
    yeta_a[i]=0;
    e_yeta_a[i]=0;
  }
  for(int i=0; i<= 14 ; i++){
    if(heta_badLC->GetBinContent(i)  != 0 &&
       heta_goodLC->GetBinContent(i) != 0 ){
      double eff=100*efficiency( heta_badLC->GetBinContent(i),heta_goodLC->GetBinContent(i)+heta_badLC->GetBinContent(i) );
      double sig_eff=100*sigma_efficiency( heta_badLC->GetBinContent(i),heta_goodLC->GetBinContent(i)+heta_badLC->GetBinContent(i) );
      jeta++;
      xeta[jeta]=heta_goodLC->GetBinCenter(i);
      e_xeta[jeta]=heta_goodLC->GetBinWidth(i)/2.;
      yeta[jeta]=eff;
      e_yeta[jeta]=sig_eff;
      neta++;
   }
    if(heta_badJC->GetBinContent(i)  != 0 &&
       heta_goodJC->GetBinContent(i) != 0 ){
      double eff=100*efficiency( heta_badJC->GetBinContent(i),heta_goodJC->GetBinContent(i)+heta_badJC->GetBinContent(i) );
      double sig_eff=100*sigma_efficiency( heta_badJC->GetBinContent(i),heta_goodJC->GetBinContent(i)+heta_badJC->GetBinContent(i) );
      jeta_a++;
      xeta_a[jeta_a]=heta_goodJC->GetBinCenter(i);
      e_xeta_a[jeta_a]=heta_goodJC->GetBinWidth(i)/2.;
      yeta_a[jeta_a]=eff;
      e_yeta_a[jeta_a]=sig_eff;
      neta_a++;
      }

  }

  TGraphErrors *greta=new TGraphErrors(neta,xeta,yeta,e_xeta,e_yeta);
  TGraphErrors *greta_a=new TGraphErrors(neta_a,xeta_a,yeta_a,e_xeta_a,e_yeta_a);

  TCanvas *c2=new TCanvas("c2","c2",700,700);
  c2->cd();
  greta->GetXaxis()->SetTitle("#eta");
  greta->GetYaxis()->SetTitle("Loss of photons (%)");
  greta->GetYaxis()->SetRangeUser(0.01,0.26);
  greta->SetLineColor(3);
  greta->SetMarkerColor(3);
  greta_a->GetXaxis()->SetTitle("#eta");
  greta_a->GetYaxis()->SetTitle("Loss of photons (%)");
  greta_a->GetYaxis()->SetRangeUser(0.01,0.26);
  greta_a->SetLineColor(2);
  greta_a->SetMarkerColor(2);
  greta->Draw("AP");
  greta_a->Draw("sameP");


}

///////////////////////////////
void Results::DrawHistPt()
///////////////////////////////
{

  TH1F *hpt = new TH1F("hpt","hpt",20,0,500);
  hpt->GetXaxis()->SetTitle("p_{T} [GeV]");
  TH1F *hpt_badLC=new TH1F("hpt_badLC","hpt_badLC",20,0,500);
  hpt_badLC->GetXaxis()->SetTitle("p_{T} [GeV]");
  TH1F *hpt_badJC=new TH1F("hpt_badJC","hpt_badJC",20,0,500);
  hpt_badJC->GetXaxis()->SetTitle("p_{T} [GeV]");
  TH1F *hpt_ratio=new TH1F("hpt_ratio","hpt_ratio",20,0,500);
  hpt_ratio->GetXaxis()->SetTitle("p_{T} [GeV]");
  TH1F *hpt_ratio1=new TH1F("hpt_ratio1","hpt_ratio1",20,0,500);
  hpt_ratio1->GetXaxis()->SetTitle("p_{T} [GeV]");

  for(int iph=0;iph<m_nentries;iph++){
    m_Rd->GetEntry(iph);
    hpt->Fill((float)m_Rd->GetVariable("phclpt")/1000.);
    if( (int)m_Rd->GetVariable("phPassLC2")==0){
      hpt_badLC->Fill( (float)m_Rd->GetVariable("phclpt")/1000. );
    }
    if( (int)m_Rd->GetVariable("phPassJC")==0){
      hpt_badJC->Fill( (float)m_Rd->GetVariable("phclpt")/1000. );
    }

  }


  hpt_badJC->SetLineColor(2);
  hpt_badJC->SetMarkerColor(2);

  hpt_badLC->SetLineColor(3);
  hpt_badLC->SetMarkerColor(3);

  TCanvas *c2=new TCanvas("c2","c2",700,700);
  c2->cd();
  hpt_badLC->Sumw2();
  hpt_badJC->Sumw2();
  hpt->Sumw2();
  hpt_ratio->Divide(hpt_badLC,hpt);
  hpt_ratio1->Divide(hpt_badJC,hpt);
  hpt_ratio->SetLineColor(3);
  hpt_ratio->SetMarkerColor(3);
  hpt_ratio1->SetLineColor(2);
  hpt_ratio1->SetMarkerColor(2);
  hpt_ratio->GetXaxis()->SetRangeUser(0,225);
  hpt_ratio1->GetXaxis()->SetRangeUser(0,225);
  hpt_ratio->Draw("PE");
  hpt_ratio1->Draw("samePE");

  hpt_badLC->SetNormFactor(1);
  hpt_badJC->SetNormFactor(1);
  hpt->SetNormFactor(1);
  TCanvas *c1=new TCanvas("c1","c1",700,700);
  c1->cd();
  hpt->Draw("PE");
  hpt_badLC->Draw("samePE");
  hpt_badJC->Draw("samePE");

 }
///////////////////////////////
void Results::DrawHisteta()
///////////////////////////////
{

  TH1F *heta = new TH1F("heta","heta",14,-2.8,2.8);
  heta->GetXaxis()->SetTitle("#eta");
  TH1F *heta_badLC=new TH1F("heta_badLC","heta_badLC",14,-2.8,2.8);
  heta_badLC->GetXaxis()->SetTitle("#eta");
  TH1F *heta_badJC=new TH1F("heta_badJC","heta_badJC",14,-2.8,2.8);
  heta_badJC->GetXaxis()->SetTitle("#eta");
  TH1F *heta_ratio=new TH1F("heta_ratio","heta_ratio",14,-2.8,2.8);
  heta_ratio->GetXaxis()->SetTitle("#eta");
  TH1F *heta_ratio1=new TH1F("heta_ratio1","heta_ratio1",14,-2.8,2.8);
  heta_ratio1->GetXaxis()->SetTitle("#eta");

  for(int iph=0;iph<m_nentries;iph++){
    m_Rd->GetEntry(iph);
    heta->Fill((float)m_Rd->GetVariable("phetas2") );
    if( (int)m_Rd->GetVariable("phPassLC2")==0){
      heta_badLC->Fill( (float)m_Rd->GetVariable("phetas2") );
    }
    if( (int)m_Rd->GetVariable("phPassJC")==0){
      heta_badJC->Fill( (float)m_Rd->GetVariable("phetas2") );
    }

  }


  heta_badJC->SetLineColor(2);
  heta_badJC->SetMarkerColor(2);

  heta_badLC->SetLineColor(3);
  heta_badLC->SetMarkerColor(3);

 TCanvas *c2=new TCanvas("c2","c2",700,700);
  c2->cd();
  heta_badLC->Sumw2();
  heta_badJC->Sumw2();
  heta->Sumw2();
  heta_ratio->Divide(heta_badLC,heta);
  heta_ratio1->Divide(heta_badJC,heta);
  heta_ratio->SetLineColor(3);
  heta_ratio->SetMarkerColor(3);
  heta_ratio1->SetLineColor(2);
  heta_ratio1->SetMarkerColor(2);
//   heta_ratio->GetXaxis()->SetRangeUser(0,225);
//   heta_ratio1->GetXaxis()->SetRangeUser(0,225);
  heta_ratio->Draw("PE");
  heta_ratio1->Draw("samePE");

  heta_badLC->SetNormFactor(1);
  heta_badJC->SetNormFactor(1);
  heta->SetNormFactor(1);
  TCanvas *c1=new TCanvas("c1","c1",700,700);
  c1->cd();
  heta_badLC->Draw("PE");
  heta_badJC->Draw("samePE");
  heta->Draw("samePE");
 
 }

///////////////////////////////
void Results::DrawStackedHist()
///////////////////////////////
{
  
  TH1F *hpt_goodLC=new TH1F("hpt_goodLC","hpt_goodLC",40,0,500);
  TH1F *hpt_badLC=new TH1F("hpt_badLC","hpt_badLC",40,0,500);

  TH1F *heta_goodLC=new TH1F("heta_goodLC","heta_goodLC",56,-2.8,2.8);
  TH1F *heta_badLC=new TH1F("heta_badLC","heta_badLC",56,-2.8,2.8);

  TH1F *hpt_goodJC=new TH1F("hpt_goodJC","hpt_goodJC",40,0,500);
  TH1F *hpt_badJC=new TH1F("hpt_badJC","hpt_badJC",40,0,500);

  TH1F *heta_goodJC=new TH1F("heta_goodJC","heta_goodJC",56,-2.8,2.8);
  TH1F *heta_badJC=new TH1F("heta_badJC","heta_badJC",56,-2.8,2.8);


  //== LOOP OVER PHOTONS ===// 
  for(int entry=0;entry<m_nentries;entry++){
    m_Rd->GetEntry(entry);
    if((int)m_Rd->GetVariable("larError")!=0) continue;
    if((int)m_Rd->GetVariable("phistight")==0) continue;

    if((int)m_Rd->GetVariable("phPassLC2")==0){
      hpt_badLC->Fill( m_Rd->GetVariable("phclpt")/1000. );
      heta_badLC->Fill( m_Rd->GetVariable("phetas2") );
    }else{
      hpt_goodLC->Fill( m_Rd->GetVariable("phclpt")/1000. );
      heta_goodLC->Fill( m_Rd->GetVariable("phetas2") );
    }
    
    if((int)m_Rd->GetVariable("phPassJC")==0){
      hpt_badJC->Fill( m_Rd->GetVariable("phclpt")/1000. );
      heta_badJC->Fill( m_Rd->GetVariable("phetas2") );
    }else{
      hpt_goodJC->Fill( m_Rd->GetVariable("phclpt")/1000. );
      heta_goodJC->Fill( m_Rd->GetVariable("phetas2") );
    }
  }
  
  hpt_badLC->SetLineColor(3);
  hpt_badLC->SetFillColor(3);
  hpt_badJC->SetLineColor(2);
  hpt_badJC->SetFillColor(2);
  hpt_goodLC->SetFillColor(15);
  hpt_goodJC->SetFillColor(15);

  heta_badLC->SetLineColor(3);
  heta_badLC->SetFillColor(3);
  heta_badJC->SetLineColor(2);
  heta_badJC->SetFillColor(2);
  heta_goodLC->SetFillColor(15);
  heta_goodJC->SetFillColor(15);

  THStack *hspt_LC=new THStack("hspt_LC","hspt");
  hspt_LC->SetMinimum(1);
  hspt_LC->Add(hpt_badLC);
  hspt_LC->Add(hpt_goodLC);

  THStack *hspt_JC=new THStack("hspt_JC","hspt");
  hspt_JC->SetMinimum(1);
  hspt_JC->Add(hpt_badJC);
  hspt_JC->Add(hpt_goodJC);

  TCanvas *c1=new TCanvas("c1","c1",700,700);
  c1->cd();
  c1->SetLogy();
  hspt_LC->Draw();
  hspt_LC->GetXaxis()->SetTitle("p_{T} [GeV]");
  hspt_LC->Draw();

  TCanvas *c2=new TCanvas("c2","c2",700,700);
  c2->cd();
  c2->SetLogy();
  hspt_JC->Draw();
  hspt_JC->GetXaxis()->SetTitle("p_{T} [GeV]");
  hspt_JC->Draw();

  THStack *hseta_LC=new THStack("hseta_LC","hseta");
  hseta_LC->SetMinimum(1);
  hseta_LC->Add(heta_badLC);
  hseta_LC->Add(heta_goodLC);

  THStack *hseta_JC=new THStack("hseta_JC","hseta");
  hseta_JC->SetMinimum(1);
  hseta_JC->Add(heta_badJC);
  hseta_JC->Add(heta_goodJC);

  TCanvas *c3=new TCanvas("c3","c3",700,700);
  c3->cd();
  c3->SetLogy();
  hseta_LC->Draw();
  hseta_LC->GetXaxis()->SetTitle("#eta");
  hseta_LC->Draw();

  TCanvas *c4=new TCanvas("c4","c4",700,700);
  c4->cd();
  c4->SetLogy();
  hseta_JC->Draw();
  hseta_JC->GetXaxis()->SetTitle("#eta");
  hseta_JC->Draw();



}
///////////////////////////////
void Results::DrawHistMET()
///////////////////////////////
{

  TH1F *hMET_goodLC = new TH1F("hMET_goodLC","hMET_goodLC",31,-10,300);
  hMET_goodLC->GetXaxis()->SetTitle("MET [GeV]");
  TH1F *hMET_badLC=new TH1F("hMET_badLC","hMET_badLC",31,-10,300);
  hMET_badLC->GetXaxis()->SetTitle("MET [GeV]");
  TH1F *hMET_goodJC = new TH1F("hMET_goodJC","hMET_goodJC",31,-10,300);
  hMET_goodJC->GetXaxis()->SetTitle("MET [GeV]");
  TH1F *hMET_badJC=new TH1F("hMET_badJC","hMET_badJC",31,-10,300);
  hMET_badJC->GetXaxis()->SetTitle("MET [GeV]");

  int myEvent=0; 
  bool newEvent=false;
  for(int iph=0;iph<m_nentries;iph++){
    m_Rd->GetEntry(iph);
    newEvent=false;
    if(myEvent!=m_Rd->GetVariable("EventNumber")){
      myEvent=m_Rd->GetVariable("EventNumber");
      newEvent=true;
    }
    if((int)m_Rd->GetVariable("larError")!=0) continue;
    if(newEvent){
      if( (int)m_Rd->GetVariable("phPassLC2")==0){
	hMET_badLC->Fill( (float)m_Rd->GetVariable("MET")/1000. );
      }else{
	hMET_goodLC->Fill((float)m_Rd->GetVariable("MET")/1000.);
      }
      if( (int)m_Rd->GetVariable("phPassJC")==0){
	hMET_badJC->Fill( (float)m_Rd->GetVariable("MET")/1000. );
      }else{
	hMET_goodJC->Fill((float)m_Rd->GetVariable("MET")/1000.);
      }
    }

  }

  hMET_badLC->SetLineColor(3);
  hMET_badLC->SetFillColor(3);
  hMET_badJC->SetLineColor(2);
  hMET_badJC->SetFillColor(2);
  hMET_goodLC->SetFillColor(15);
  hMET_goodJC->SetFillColor(15);

  THStack *hsMET_LC=new THStack("hsMET_LC","hsMET");
  hsMET_LC->SetMinimum(1);
  hsMET_LC->Add(hMET_badLC);
  hsMET_LC->Add(hMET_goodLC);
  THStack *hsMET_JC=new THStack("hsMET_JC","hsMET");
  hsMET_JC->SetMinimum(1);
  hsMET_JC->Add(hMET_badJC);
  hsMET_JC->Add(hMET_goodJC);

  TCanvas *c1=new TCanvas("c1","c1",700,700);
  c1->cd();
  c1->SetLogy();
  hsMET_LC->Draw();
  hsMET_LC->GetXaxis()->SetTitle("MET [GeV]");
  hsMET_LC->Draw();

  TCanvas *c2=new TCanvas("c2","c2",700,700);
  c2->cd();
  c2->SetLogy();
  hsMET_JC->Draw();
  hsMET_JC->GetXaxis()->SetTitle("MET [GeV]");
  hsMET_JC->Draw();

 }

//////////////////////////////////////
void Results::DrawInvMassStackedHist()
//////////////////////////////////////
{
  TH1F *hgood=new TH1F("hgood","hgood",92,80,1000);
  TH1F *hbad=new TH1F("hbad","hbad",92,80,1000);
  hbad->SetFillColor(3);
  hbad->SetLineColor(3);
  hgood->SetFillColor(15);
  for(int entry=0;entry<m_nentries;entry++){
    m_Rd->GetEntry(entry);
    if((int)m_Rd->GetVariable("larError")!=0) continue;
    if( (int)m_Rd->GetVariable("Lead_PassLC2")==0 
	||(int)m_Rd->GetVariable("SubLead_PassLC2")==0){
      hbad->Fill(m_Rd->GetVariable("mgg"));
    }else{
      hgood->Fill(m_Rd->GetVariable("mgg"));
    }
  }
  THStack *hs=new THStack("hs","hs");
  hs->Add(hbad);
  hs->Add(hgood);
  TCanvas *c1=new TCanvas("c1","c1",700,700);
  c1->cd();
  c1->SetLogy();
  hs->Draw();
  hs->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
  hs->Draw();

}
/////////////////////////////////
void Results::DrawBkgComponents()
/////////////////////////////////
{
  TFile* fBin   = new TFile("rootfiles/BinnedBkg_mgg140400.root","read");
  TH1F* hmgg_Bin = (TH1F*)fBin->Get("hmgg_gg");
  TH1F* hmgj_Bin = (TH1F*)fBin->Get("hmgg_gj");
  TH1F* hmjg_Bin = (TH1F*)fBin->Get("hmgg_jg");
  TH1F* hmjj_Bin = (TH1F*)fBin->Get("hmgg_jj");
  TH1F* hmgg_iso_Bin = (TH1F*)fBin->Get("hmgg_gg_iso");
  TH1F* hmgj_iso_Bin = (TH1F*)fBin->Get("hmgg_gj_iso");
  TH1F* hmjg_iso_Bin = (TH1F*)fBin->Get("hmgg_jg_iso");
  TH1F* hmjj_iso_Bin = (TH1F*)fBin->Get("hmgg_jj_iso");
  TFile* fSplot = new TFile("rootfiles/SplotBkg_mgg140400.root","read");
  TH1F* hmgg_Splot = (TH1F*)fSplot->Get("hmgg_gg");
  TH1F* hmgj_Splot = (TH1F*)fSplot->Get("hmgg_gj");
  TH1F* hmjg_Splot = (TH1F*)fSplot->Get("hmgg_jg");
  TH1F* hmjj_Splot = (TH1F*)fSplot->Get("hmgg_jj");
  TH1F* hmgg_iso_Splot = (TH1F*)fSplot->Get("hmgg_gg_iso");
  TH1F* hmgj_iso_Splot = (TH1F*)fSplot->Get("hmgg_gj_iso");
  TH1F* hmjg_iso_Splot = (TH1F*)fSplot->Get("hmgg_jg_iso");
  TH1F* hmjj_iso_Splot = (TH1F*)fSplot->Get("hmgg_jj_iso");

  TCanvas *cTi = new TCanvas("cTi","Components TiTi",800,800);
  cTi->Divide(2,2);
  cTi->cd(1);
  hmgg_Bin->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  hmgg_Bin->Draw("HIST");
  hmgg_Splot->Draw("samePE");
  cTi->cd(2);
  hmgj_Bin->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  hmgj_Bin->Draw("HIST");
  hmgj_Splot->Draw("samePE");
  cTi->cd(3);
  hmjg_Bin->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  hmjg_Bin->Draw("HIST");
  hmjg_Splot->Draw("samePE");
  cTi->cd(4);
  hmjj_Bin->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  hmjj_Bin->Draw("HIST");
  hmjj_Splot->Draw("samePE");
  cTi->SaveAs("./plots/ComponentsTiTi.eps");
  TCanvas *cTiIso = new TCanvas("cTiIso","Components TiIsoTiIso",800,800);
  cTiIso->Divide(2,2);
  cTiIso->cd(1);
  hmgg_iso_Bin->Draw("HIST");
  hmgg_iso_Splot->Draw("samePE");
  cTiIso->cd(2);
  hmgj_iso_Bin->Draw("HIST");
  hmgj_iso_Splot->Draw("samePE");
  cTiIso->cd(3);
  hmjg_iso_Bin->Draw("HIST");
  hmjg_iso_Splot->Draw("samePE");
  cTiIso->cd(4);
  hmjj_iso_Bin->Draw("HIST");
  hmjj_iso_Splot->Draw("samePE");
  cTiIso->SaveAs("./plots/ComponentsTiIsoTiIso.eps");

}
//////////////////////////////////////////////////////
double Results::efficiency(double k,double n)
//////////////////////////////////////////////////////
{
  double temp=(k+1)/(n+2);
  return temp;
}
////////////////////////////////////////////////////////////
double Results::sigma_efficiency(double k,double n)
////////////////////////////////////////////////////////////
{
  double r1=(k+1)/(n+2);
  double r2=(k+2)/(n+3);
  double r3=(k+1)/(n+2);
  double temp=sqrt(r1*r2-r3*r3);
  return temp;
}
