#define IsolationStudies_cxx
#include "TFile.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH2.h"
#include "PersoIsolationStudies.h"
#include "InvMassCalculator.h"
#include "ToolsUtilities.h"
#include "GoodRunsLists/DQHelperFunctions.h"


////////////////////////////////////////////////////////////////////////
IsolationStudies::IsolationStudies(TTree* tree) : GravitonAnalysis(tree)
////////////////////////////////////////////////////////////////////////
{
   if (tree == 0) {
     std::cout << "You need to specify an input tree !!!"<<std::endl;
   }
   Init(tree);
}
/////////////////////////////////////////////////////////
IsolationStudies::IsolationStudies() : GravitonAnalysis()
/////////////////////////////////////////////////////////
{
  //Default constructor
  m_Npt = 0;

  m_passTrigger = false;
  m_inGRL       = false;
  m_passPV      = false;
  m_passPreSel  = false;
}
/////////////////////////////////////
IsolationStudies::~IsolationStudies()
/////////////////////////////////////
{
  //Default destructor

}
////////////////////////////////////////
void IsolationStudies::Init(TTree* tree)
////////////////////////////////////////
{
  m_passTrigger = false;
  m_inGRL       = false;
  m_passPV      = false;
  m_passPreSel  = false;

  double ptbins[]     = {27.5,32.5,37.5,45.,55.,65.,75.,85.,95.,150.};
  double err_ptbins[] = { 2.5, 2.5, 2.5, 5., 5., 5., 5., 5., 5., 50.};
  m_bin_pt.assign(ptbins,ptbins+10);
  m_err_bin_pt.assign(err_ptbins,err_ptbins+10);

  m_Npt=m_bin_pt.size();
  Nph_isotight_pt.assign(m_Npt,0);
  Nph_tight_pt.assign(m_Npt,0);

}
/////////////////////////////////////////////////
void IsolationStudies::CutOptimisation(bool data)
////////////////////////////////////////////////
{
  double Nevent=0;
  double Ntot_TL=0;
  double Ntot_notTL=0;
  double Ntot_TSL=0;
  double Ntot_notTSL=0;
//   int Ncut=10;
//   double x[]  ={1,2.5,3.5,4.5,5.5,6.5,8.5,12.5,17.5,22.5};
//   double e_x[]={1,0.5,0.5,0.5,0.5,0.5,1.5,2.5,2.5,2.5};
//   double cut[]={2,3,4,5,6,7,10,15,20,25};
  int Ncut=22;
  double x[]  ={0.25,0.75,1.25,1.75,2.1,2.3,2.5,2.7,2.9,3.1,
		3.3,3.5,3.7,3.9,4.25,4.75,5.5,6.5,8.5,12.5,17.5,22.5};
  double e_x[]={0.25,0.25,0.25,0.25,0.1,0.1,0.1,0.1,0.1,0.1,
		0.1,0.1,0.1,0.1,0.25,0.25,0.5,0.5,1.5,2.50,2.50,2.50};
  double cut[]={0.50,1.00,1.50,2.00,2.2,2.4,2.6,2.8,3.0,3.2,
		3.4,3.6,3.8,4.0,4.50,5.00,6.0,7.0,10.,15.0,20.0,25.0};

  int Nbins=45;
  double Xmin=-10;
  double Xmax=35;
  TH1F *hL=new TH1F("hL","hL",Nbins,Xmin,Xmax);
  hL->GetXaxis()->SetTitle("Etcone40_corrected [GeV]");
  TH1F *hT_L=new TH1F("hT_L","hT_L",Nbins,Xmin,Xmax);
  hT_L->GetXaxis()->SetTitle("Etcone40_corrected [GeV]");
  hT_L->SetLineColor(3);
  TH1F *hnotT_L=new TH1F("hnotT_L","hnotT_L",Nbins,Xmin,Xmax);
  hnotT_L->GetXaxis()->SetTitle("Etcone40_corrected [GeV]");
  hnotT_L->SetLineColor(4);
  TH1F *hSL=new TH1F("hSL","hSL",Nbins,Xmin,Xmax);
  hSL->GetXaxis()->SetTitle("Etcone40_corrected [GeV]");
  TH1F *hT_SL=new TH1F("hT_SL","hT_SL",Nbins,Xmin,Xmax);
  hT_SL->GetXaxis()->SetTitle("Etcone40_corrected [GeV]");
  hT_SL->SetLineColor(3);
  TH1F *hnotT_SL=new TH1F("hnotT_SL","hnotT_SL",Nbins,Xmin,Xmax);
  hnotT_SL->GetXaxis()->SetTitle("Etcone40_corrected [GeV]");
  hnotT_SL->SetLineColor(4);

  double N_T_L[Ncut];
  double N_notT_L[Ncut];
  double rej_L[Ncut];
  double e_rej_L[Ncut];
  double eff_L[Ncut];
  double e_eff_L[Ncut];
  double SoverB_L[Ncut];
  double e_SoverB_L[Ncut];

  double N_T_SL[Ncut];
  double N_notT_SL[Ncut];
  double rej_SL[Ncut];
  double e_rej_SL[Ncut];
  double eff_SL[Ncut];
  double e_eff_SL[Ncut];
  double SoverB_SL[Ncut];
  double e_SoverB_SL[Ncut];

  double N_eventS[Ncut][Ncut];
  double N_eventB[Ncut][Ncut];


  for(int icut=0;icut<Ncut;icut++){
  N_T_L[icut]=0;
  N_notT_L[icut]=0;
  rej_L[icut]=0;
  e_rej_L[icut]=0;
  eff_L[icut]=0;
  e_eff_L[icut]=0;
  SoverB_L[icut]=0;
  e_SoverB_L[icut]=0;

  N_T_SL[icut]=0;
  N_notT_SL[icut]=0;
  rej_SL[icut]=0;
  e_rej_SL[icut]=0;
  eff_SL[icut]=0;
  e_eff_SL[icut]=0;
  SoverB_SL[icut]=0;
  e_SoverB_SL[icut]=0;
  }
  m_data        = data;//Choose between data and mc
  double weight = 1;
  //== PileUp weight ==//
  double PUweight = 1;
  //==Declare the GRL object==//
  DQ::SetXMLFile(m_GRL);
  //__________________________//
  
  
  //======================================================================//  
  //================ Start the Loop over all entries =====================//
  //======================================================================//  
  for(int jentry=0;jentry<m_nentries;jentry++){
    m_Rd->GetEntry(jentry);
    
    //==================== PU Weight for MC ONLY ========================//
    if(!m_data){
      PUweight = m_pileupTool->getPileupWeight(m_Rd->GetVariable("averageIntPerXing"),
					     (int)m_Rd->GetVariable("RunNumber"),
					     (int)m_Rd->GetVariable("mc_channel_number") 
					     );
    }
    //___________________________________________________________________//
 
    //===== clean all flags ======//     
    m_passTrigger       = false ;
    m_inGRL             = false ;
    m_passPV            = false ;
    m_passPreSel        = false ;
    //____________________________//
 
   //================= GRL Selection ==================//
    int myRunNumber=(int)m_Rd->GetVariable("RunNumber");
    int myLumiBlock=(int)m_Rd->GetVariable("lbn");
    if(m_data){
      m_inGRL = DQ::PassRunLB(myRunNumber ,myLumiBlock);
    }else{
      m_inGRL=true;
    }
    //_________________________________________________//
    
    //=================== Trigger Requirement ====================//
    m_passTrigger    = (bool)m_Rd->GetVariable("EF_2g20_loose");
    //____________________________________________________________//
 
    
    //==================== Primary Vertex Requirement =========================//
    unsigned int n_PV=(unsigned int)m_Rd->GetVariable("@PV_nTracks.size()");
    for(unsigned int i=0;i<n_PV;i++){
      if( (int)m_Rd->GetVariable(Form("PV_nTracks[%d]",i))>2){
	  m_passPV=true;
	  break;
      }
    }
    //_________________________________________________________________________//
    

    //============================= Diphotons selection ============================//
    if( !m_passTrigger ) continue;
    if( !m_inGRL )       continue; 
    if( !m_passPV )      continue; 

    double pt_lead=-99999; double pt_sublead=-99999;
    int ilead=-99999; int isublead=-99999;
    m_passPreSel=EventPreSelectionOK(&ilead,&isublead,&pt_lead,&pt_sublead);
    weight=PUweight;
    
    if( !m_passPreSel )  continue;
    if(pt_lead<40 || pt_sublead<40) continue;
    double mggGEV_old =0;
    double E1      = m_Rd->GetVariable( Form("ph_E[%d]",ilead) );
    double eta1    = m_Rd->GetVariable( Form("ph_etas1[%d]",ilead) );
    double phi1    = m_Rd->GetVariable( Form("ph_phi[%d]",ilead) );
    double E2      = m_Rd->GetVariable( Form("ph_E[%d]",isublead) );
    double eta2    = m_Rd->GetVariable( Form("ph_etas1[%d]",isublead) );
    double phi2    = m_Rd->GetVariable( Form("ph_phi[%d]",isublead) );
    double PV_ID   = m_Rd->GetVariable("PV_z[0]");
    mggGEV_old     = (1./1000.)*InvMassCalc::GetCorrectedInvMass( E1,eta1,phi1,
								  E2,eta2,phi2,
								  PV_ID );
    if( mggGEV_old<120) continue;
 
    ++Nevent;
    hL->Fill( m_Rd->GetVariable(Form("ph_Etcone40_corrrected[%d]",ilead))/1000. );
    hSL->Fill( m_Rd->GetVariable(Form("ph_Etcone40_corrected[%d]",isublead))/1000. ); 
    if( PhotonIsTightOK(ilead) ){
      ++Ntot_TL;
      hT_L->Fill( m_Rd->GetVariable(Form("ph_Etcone40_corrected[%d]",ilead))/1000. );
    }else{
      ++Ntot_notTL;
      hnotT_L->Fill( m_Rd->GetVariable(Form("ph_Etcone40_corrected[%d]",ilead))/1000. );
    }
    if( PhotonIsTightOK(isublead) ){
      ++Ntot_TSL;
      hT_SL->Fill( m_Rd->GetVariable(Form("ph_Etcone40_corrected[%d]",isublead))/1000. );
    }else{
      ++Ntot_notTSL;
      hnotT_SL->Fill( m_Rd->GetVariable(Form("ph_Etcone40_corrected[%d]",isublead))/1000. );
    }

    double Iso_L=  m_Rd->GetVariable(Form("ph_Etcone40_corrected[%d]",ilead))/1000;   
    double Iso_SL= m_Rd->GetVariable(Form("ph_Etcone40_corrected[%d]",isublead))/1000;   

    for(int icut=0;icut<Ncut;icut++){
      if( PhotonIsTightOK(ilead) ){
	if(Iso_L<cut[icut]){
	  N_T_L[icut]=N_T_L[icut]+1;
	}
      }else{
	if(Iso_L<cut[icut]){
	  N_notT_L[icut]=N_notT_L[icut]+1;
	}
      }
      if( PhotonIsTightOK(isublead) ){
	if(Iso_SL<cut[icut]){
	  N_T_SL[icut]=N_T_SL[icut]+1;
	}
      }else{
	if(Iso_SL<cut[icut]){
	  N_notT_SL[icut]=N_notT_SL[icut]+1;
	}
      }
    }
    for(int icutL=0;icutL<Ncut;icutL++){
      for(int icutSL=0;icutSL<Ncut;icutSL++){
	if (Iso_L< cut[icutL] && Iso_SL< cut[icutSL]){
	  if( PhotonIsTightOK(ilead) &&
	      PhotonIsTightOK(isublead) ){
	    N_eventS[icutL][icutSL]= N_eventS[icutL][icutSL]+1;
	  }else{
	    N_eventB[icutL][icutSL]= N_eventB[icutL][icutSL]+1;
	  }
	}
      }
    }
    //_________________________________________________________________________________//
    AnalysisTools myAT;
    myAT.ShowNumberOfProcessedEvents(jentry,m_nentries);
  }
  //========================================================================//
  //=============== End of the loop over all entries =======================//
  //========================================================================//
  TCanvas *c1=new TCanvas("c1","c1",1024,800);
  c1->Divide(3,2);
  c1->cd(1);
  hT_L->Draw();
  c1->cd(2);
  hnotT_L->Draw();
  c1->cd(3);
  hL->Draw("PE");
  hT_L->Draw("same");
  hnotT_L->Draw("same");
  c1->cd(4);
  hT_SL->Draw();
  c1->cd(5);
  hnotT_SL->Draw();
  c1->cd(6);
  hSL->Draw("PE");
  hT_SL->Draw("same");
  hnotT_SL->Draw("same");

  for(int icut=0;icut<Ncut;icut++){
    std::cout << x[icut] << std::endl;
    std::cout << N_T_L[icut]  << " " << N_notT_L[icut]<< " " 
	      << N_T_SL[icut] << " " << N_notT_SL[icut] 
	      << std::endl;
    AnalysisTools myAT;
    eff_L[icut]       = myAT.efficiency(N_T_L[icut],Ntot_TL);
    e_eff_L[icut]     = myAT.sigma_efficiency(N_T_L[icut],Ntot_TL);
    rej_L[icut]       = myAT.efficiency(Ntot_notTL-N_notT_L[icut],Ntot_notTL);
    e_rej_L[icut]     = myAT.sigma_efficiency(Ntot_notTL-N_notT_L[icut],Ntot_notTL);
    SoverB_L[icut]    = N_T_L[icut]/sqrt(N_notT_L[icut]);
    e_SoverB_L[icut]  = 0.5;
    eff_SL[icut]      = myAT.efficiency(N_T_SL[icut],Ntot_TSL);
    e_eff_SL[icut]    = myAT.sigma_efficiency(N_T_SL[icut],Ntot_TSL);
    rej_SL[icut]      = myAT.efficiency(Ntot_notTSL-N_notT_SL[icut],Ntot_notTSL);
    e_rej_SL[icut]    = myAT.sigma_efficiency(Ntot_notTSL-N_notT_SL[icut],Ntot_notTSL);
    SoverB_SL[icut]   = N_T_SL[icut]/sqrt(N_notT_SL[icut]);
    e_SoverB_SL[icut] = 0.001;
    std::cout << eff_L[icut] << " " 
	      << e_eff_L[icut] << " " 
	      << rej_L[icut] << " "
	      << e_rej_L[icut]
	      << std::endl;
  }
  TH2F* h2D=new TH2F("h2D","h2D",21,x,21,x);
  h2D->GetXaxis()->SetTitle("Leading Cut"); 
  h2D->GetYaxis()->SetTitle("SubLeading Cut"); 
  for(int icutL=0;icutL<Ncut;icutL++){
    for(int icutSL=0;icutSL<Ncut;icutSL++){
      double SoverB=N_eventS[icutL][icutSL]/sqrt(N_eventB[icutL][icutSL]) ;
      h2D->SetBinContent(icutL+1,icutSL+1,SoverB);
    }
  } 
  TGraphErrors *gR_effL=new TGraphErrors(Ncut,x,eff_L,e_x,e_eff_L);
  gR_effL->GetXaxis()->SetTitle("Cut Value [GeV]");
  gR_effL->SetLineColor(3);
  gR_effL->SetMarkerColor(3);
  TGraphErrors *gR_rejL=new TGraphErrors(Ncut,x,rej_L,e_x,e_rej_L);
  gR_rejL->GetXaxis()->SetTitle("Cut Value [GeV]");
  gR_rejL->SetLineColor(4);
  gR_rejL->SetMarkerColor(4);
  TGraphErrors *gR_SoverBL=new TGraphErrors(Ncut,x,SoverB_L,e_x,e_SoverB_L);
  gR_SoverBL->GetXaxis()->SetTitle("Cut Value [GeV]");
  gR_SoverBL->SetLineColor(2);
  gR_SoverBL->SetMarkerColor(2);
  gR_SoverBL->GetYaxis()->SetTitle("#frac{S}{#sqrt{B}}");


  TGraphErrors *gR_effSL=new TGraphErrors(Ncut,x,eff_SL,e_x,e_eff_SL);
  gR_effSL->GetXaxis()->SetTitle("Cut Value [GeV]");
  gR_effSL->SetLineColor(3);
  gR_effSL->SetMarkerColor(3);
  TGraphErrors *gR_rejSL=new TGraphErrors(Ncut,x,rej_SL,e_x,e_rej_SL);
  gR_rejSL->GetXaxis()->SetTitle("Cut Value [GeV]");
  gR_rejSL->SetLineColor(4);
  gR_rejSL->SetMarkerColor(4);
  TGraphErrors *gR_SoverBSL=new TGraphErrors(Ncut,x,SoverB_SL,e_x,e_SoverB_SL);
  gR_SoverBSL->GetXaxis()->SetTitle("Cut Value [GeV]");
  gR_SoverBSL->GetYaxis()->SetTitle("#frac{S}{#sqrt{B}}");

  gR_SoverBSL->SetLineColor(2);
  gR_SoverBSL->SetMarkerColor(2);


  TCanvas *c2=new TCanvas("c2","L",800,800);
  c2->Divide(1,2);
  c2->cd(1);
  gR_effL->GetYaxis()->SetRangeUser(0,1);
  gR_effL->Draw("AP");
  gR_rejL->Draw("sameP");
  c2->cd(2);
  gR_SoverBL->Draw("AP");
  TCanvas *c3=new TCanvas("c3","SL",800,800);
  c3->Divide(1,2);
  c3->cd(1);
  gR_effSL->GetYaxis()->SetRangeUser(0,1);
  gR_effSL->Draw("AP");
  gR_rejSL->Draw("sameP");
  c3->cd(2);
  gR_SoverBSL->Draw("AP");
 
 TCanvas *c4=new TCanvas("c4","2D",800,800);
 c4->cd();
 h2D->Draw("colz");

}
//////////////////////////////////////////////////////////////////////////////
void IsolationStudies::IsoPerNPVbins(bool data,TString output,TString mctype)
//////////////////////////////////////////////////////////////////////////////
{

  m_data  = data;//Choose between data and mc
  m_mctype= mctype;// Choose mc type "mc10a" or "mc10b"

  int NIsobins=45;
  double XIsomin=-10;
  double XIsomax=35;

  std::vector< std::pair<double,double> > bins;
  std::pair<double,double> bin0(0,2);
  std::pair<double,double> bin1(2,4);
  std::pair<double,double> bin2(4,6);
  std::pair<double,double> bin3(6,8);
  std::pair<double,double> bin4(8,10);
  bins.push_back(bin0);
  bins.push_back(bin1);
  bins.push_back(bin2);
  bins.push_back(bin3);
  bins.push_back(bin4);
//   int Nbins=bins.size();
  std::vector<TH1F*>  HistVec;
  std::vector<TTree*> TreeVec;

  double Iso_L_bin[]     = {0,0,0,0,0};
  double weight_L_bin[]  = {0,0,0,0,0};
  int IsTight_L_bin[]    = {0,0,0,0,0};
  double Iso_SL_bin[]    = {0,0,0,0,0};
  double weight_SL_bin[] = {0,0,0,0,0};
  int IsTight_SL_bin[]   = {0,0,0,0,0};

  for(unsigned int ibin=0;ibin<bins.size();ibin++){
    TString name=Form("h%d%d",(int)bins[ibin].first,(int)bins[ibin].second);
    std::cout << name << std::endl;
    TH1F* htemp=new TH1F(name,name,NIsobins,XIsomin,XIsomax);
    htemp->GetXaxis()->SetTitle("Etcone40_corrected [GeV]");
    HistVec.push_back(htemp);
    TString treename=Form("tree%d",(int)ibin);
    TString Branchname_L          = Form("Iso_L_bin%d",(int)ibin);
    TString Branchname_IsTight_L  = Form("IsTight_L_bin%d",(int)ibin);
    TString Branchname_W_L        = Form("weight_L_bin%d",(int)ibin);
    TString Branchname_SL         = Form("Iso_SL_bin%d",(int)ibin);
    TString Branchname_IsTight_SL = Form("IsTight_SL_bin%d",(int)ibin);
    TString Branchname_W_SL       = Form("weight_SL_bin%d",(int)ibin);
    TTree *tree=new TTree(treename,"Iso vs NPV");
    tree->Branch(Branchname_L,&Iso_L_bin[ibin],Branchname_L+"/D");
    tree->Branch(Branchname_IsTight_L,&IsTight_L_bin[ibin],Branchname_IsTight_L+"/I");
    tree->Branch(Branchname_W_L,&weight_L_bin[ibin],Branchname_W_L+"/D");
    tree->Branch(Branchname_SL,&Iso_SL_bin[ibin],Branchname_SL+"/D");
    tree->Branch(Branchname_IsTight_SL,&IsTight_SL_bin[ibin],Branchname_IsTight_SL+"/I");
    tree->Branch(Branchname_W_SL,&weight_SL_bin[ibin],Branchname_W_SL+"/D");
    TreeVec.push_back(tree);
//     delete htemp;
//     delete tree;
  }

  
 //======weight ======//
  double weight   = 1;
  double PUweight = 1;
  //==================//

  //============== Declare output variables ===============//

  //== Declare the GRL object ==//
  DQ::SetXMLFile(m_GRL);
  //____________________________//
  
  
  //======================================================================//  
  //================ Start the Loop over all entries =====================//
  //======================================================================//  
  for(int jentry=0;jentry<m_nentries;jentry++){
    m_Rd->GetEntry(jentry);
    
    //======== PU Weight && energy smearing for MC ONLY ==================//
    if(!m_data){
      m_Rescale.SetRandomSeed((int)m_Rd->GetVariable("EventNumber"));
      PUweight = m_pileupTool->getPileupWeight(m_Rd->GetVariable("averageIntPerXing"),
					     (int)m_Rd->GetVariable("RunNumber"),
					     (int)m_Rd->GetVariable("mc_channel_number") 
					     );
    }
    //__________________________________________________________________//
 
   // ==== clean all flags =====//     
    m_passTrigger       = false ;
    m_inGRL             = false ;
    m_passPV            = false ;
    m_passPreSel        = false ;
    //__________________________//
 
   //================= GRL Selection ====================//
    int myRunNumber=(int)m_Rd->GetVariable("RunNumber");
 //    int myEventNumber=(int)m_Rd->GetVariable("EventNumber");
    int myLumiBlock=(int)m_Rd->GetVariable("lbn");
    if(m_data){
      m_inGRL = DQ::PassRunLB(myRunNumber ,myLumiBlock);
    }else{
      m_inGRL=true;
    }
    //__________________________________________________//
    
 

    //================== Trigger Requirement ===================//
    //For Photon D3PDs
    m_passTrigger    = (bool)m_Rd->GetVariable("EF_2g20_loose");
    //__________________________________________________________//
 
    
    //======================== Primary Vertex Requirement ======================//
    unsigned int n_PV=(unsigned int)m_Rd->GetVariable("@PV_nTracks.size()");
    int NPV=(int)n_PV;
    if( !data )
      NPV = (int)m_Rd->GetVariable("@mc_vx_x.size()");
    for(unsigned int i=0;i<n_PV;i++){
      if( (int)m_Rd->GetVariable(Form("PV_nTracks[%d]",i))>2){
	  m_passPV=true;
	  break;
      }
    }
    //_________________________________________________________________________//
    

    //================================ Diphotons selection ====================================//
    if(!m_passTrigger) continue;
    if (!m_inGRL)      continue; 
    if( !m_passPV )    continue; 
    if( (int)m_Rd->GetVariable("larError")==2 ) continue;
    if( (int)m_Rd->GetVariable("tileError")==2 ) continue;

    double pt_lead=-99999; double pt_sublead=-99999;
    int ilead=-99999; int isublead=-99999;
    m_passPreSel=EventPreSelectionOK(&ilead,&isublead,&pt_lead,&pt_sublead);
    weight=PUweight;

    if( !m_passPreSel ) continue;
 //    if( !PhotonIsTightOK(ilead) )    continue;
//     if( !PhotonIsTightOK(isublead) ) continue;


 //    if( !PhotonIsolation_toolOK(ilead) )    continue;
//     if( !PhotonIsolation_toolOK(isublead) ) continue;
    
    float Iso_L;float Iso_SL;
    if( m_isotype == "CONE"){
      Iso_L  = PhotonIsolation_tool(ilead)/1000.;
      Iso_SL = PhotonIsolation_tool(isublead)/1000.;
    }else if( m_isotype == "TOPO"){
      Iso_L  = PhotonTopoIsolation_tool(ilead)/1000.;
      Iso_SL = PhotonTopoIsolation_tool(isublead)/1000.;
    }else Fatal("IsolationStudies::IsoPerNPVbins","Wrong isolation type !!!");      

    for(unsigned int i=0;i<bins.size();i++){
      if(NPV>bins[i].first &&
	 NPV<=bins[i].second){
	HistVec[i]->Fill(Iso_L,PUweight);
	HistVec[i]->Fill(Iso_SL,PUweight);
	Iso_L_bin[i]  = Iso_L;
	Iso_SL_bin[i] = Iso_SL;
	if( PhotonIsTightOK(ilead) )
	  IsTight_L_bin[i]=1;
	else
	  IsTight_L_bin[i]=0;

	if( PhotonIsTightOK(isublead) )
	  IsTight_SL_bin[i]=1;
	else
	  IsTight_SL_bin[i]=0;

	weight_L_bin[i]  = PUweight;
	weight_SL_bin[i] = PUweight;
	TreeVec[i]->Fill();

      }
    }
    AnalysisTools myAT;
    myAT.ShowNumberOfProcessedEvents(jentry,m_nentries);
  }//End of the loop over all entries
  //========================================================================//
  //========== End of the loop over all entries ============================//
  //========================================================================//



  TObjArray Hlist(0);
  for(unsigned int ibin=0;ibin<bins.size();ibin++){
    Hlist.Add(HistVec[ibin]);
    Hlist.Add(TreeVec[ibin]);
  }
  TFile *f1=new TFile(output,"RECREATE");  
  Hlist.Write();
  f1->Close();

  delete f1;

}
