#include "GravitonSelector.h"
#include "ToolsUtilities.h"
#include "ToolsCommons.h"
#include <TError.h>

/////////////////////////////////////////////////////////////////////////
GravitonSelector::GravitonSelector(TTree* tree) : GravitonAnalysis(tree)
/////////////////////////////////////////////////////////////////////////
{
   if (tree == 0) {
     std::cout << "You need to specify an input tree !!!"<<std::endl;
   }
   Init();
}
/////////////////////////////////////////////////////////
GravitonSelector::GravitonSelector() : GravitonAnalysis()
/////////////////////////////////////////////////////////
{
  //Default constructor
}
//////////////////////////////////////
GravitonSelector::~GravitonSelector()
/////////////////////////////////////
{
  //Default destructor

}
/////////////////////////////
void GravitonSelector::Init()
/////////////////////////////
{
  //========  weights ======//
  m_PUweight      = 1;
  m_SF_ID_lead    = 1;
  m_SF_ID_sublead = 1;

  m_output = "toto.root";

}
/////////////////////////////////////
void GravitonSelector::InitOutTree()
/////////////////////////////////////
{
  //-----> Output variables  
  file_out   = new TFile(m_output,"RECREATE");
  tree_out   = new TTree("diphoton","diphoton");
  hist_mgg   = new TObjArray();
  hist_mgg_w = new TObjArray();

  //------> Event informations
  tree_out->Branch("RunNumber",&m_RunNumber,"RunNumber/I");
  tree_out->Branch("EventNumber",&m_EventNumber,"EventNumber/I");
  tree_out->Branch("LumiBlock",&m_LumiBlock,"LumiBlock/I");
  tree_out->Branch("NPV",&m_NPV,"NPV/I");
  //------> Leader photon informations
  tree_out->Branch("Lead_pT",&Lead_pT,"Lead_pT/D");
  tree_out->Branch("Lead_pT_smearedup",&Lead_pT_smearedup,"Lead_pT_smearedup/D");
  tree_out->Branch("Lead_pT_smeareddown",&Lead_pT_smeareddown,"Lead_pT_smeareddown/D");
  tree_out->Branch("Lead_phi",&Lead_phi,"Lead_phi/D");
  tree_out->Branch("Lead_eta",&Lead_eta,"Lead_eta/D");
  tree_out->Branch("Lead_eta_PV",&Lead_eta_PV,"Lead_eta_PV/D");
  tree_out->Branch("Lead_IsConv",&Lead_IsConv,"Lead_IsConv/I");
  //------> Subleader photon informations
  tree_out->Branch("SubLead_pT",&SubLead_pT,"SubLead_pT/D");
  tree_out->Branch("SubLead_pT_smearedup",&SubLead_pT_smearedup,"SubLead_pT_smearedup/D");
  tree_out->Branch("SubLead_pT_smeareddown",&SubLead_pT_smeareddown,"SubLead_pT_smeareddown/D");
  tree_out->Branch("SubLead_phi",&SubLead_phi,"SubLead_phi/D");
  tree_out->Branch("SubLead_eta",&SubLead_eta,"SubLead_eta/D");
  tree_out->Branch("SubLead_eta_PV",&SubLead_eta_PV,"SubLead_eta_PV/D");
  tree_out->Branch("SubLead_IsConv",&SubLead_IsConv,"SubLead_IsConv/I");
  //------> Invariant mass calculation
  tree_out->Branch("mgg",&mgg,"mgg/D");
  tree_out->Branch("mgg_smearedup",&mgg_smearedup,"mgg_smearedup/D");
  tree_out->Branch("mgg_smeareddown",&mgg_smeareddown,"mgg_smeareddown/D");
  //------> ptgg calculation
  tree_out->Branch("ptgg",&ptgg,"ptgg/D");
  tree_out->Branch("ptgg_smearedup",&ptgg_smearedup,"ptgg_smearedup/D");
  tree_out->Branch("ptgg_smeareddown",&ptgg_smeareddown,"ptgg_smeareddown/D");
  //------> ygg calculation
  tree_out->Branch("ygg",&ygg,"ygg/D");
  tree_out->Branch("ygg_smearedup",&ygg_smearedup,"ygg_smearedup/D");
  tree_out->Branch("ygg_smeareddown",&ygg_smeareddown,"ygg_smeareddown/D");
  //------> Cos(theta*) calculation
  tree_out->Branch("costhetastar",&costhetastar,"costhetastar/D");
  tree_out->Branch("costhetastar_smearedup",&costhetastar_smearedup,"costhetastar_smearedup/D");
  tree_out->Branch("costhetastar_smeareddown",&costhetastar_smeareddown,"costhetastar_smeareddown/D");
  //------> Deltaphi calculation
  tree_out->Branch("deltaphi",&deltaphi,"deltaphi/D");
  //----> True mass (MC only)
  tree_out->Branch("mc_m",&mymc_m,"mc_m/D");
  //----> Event weight
  tree_out->Branch("weight",&weight,"weight/D");

}
/////////////////////////////////////////////////
void GravitonSelector::EventLoop(TString outfile)
////////////////////////////////////////////////
{
  m_output = outfile;
  InitOutTree();

  //============== Declare output variables ===============//
  int nchecks = 11;
  double n[nchecks];
  double n_w[nchecks];
  hCutFlow   = new TH1F("CutFlow","CutFlow",nchecks,0,nchecks);
  hCutFlow_w = new TH1F("CutFlow_w","CutFlow_w",nchecks,0,nchecks);
  //-----------------------------------------------------------------
  hmgg_final   = new TH1D("hmgg","mgg",Commons::nBins,Commons::binning);
  hmgg_final_w = new TH1D("hmgg_w","weighted mgg",Commons::nBins,Commons::binning);
  TH1F* hmgg[nchecks];
  TH1F* hmgg_w[nchecks];
  //--------------------------------------------------------------------------------------//
  
  
  //=========================== Cut Flow checks =========================//
  for(int ic=0 ;ic<nchecks;ic++){
    n[ic]=0;n_w[ic]=0;
    hmgg[ic]   = new TH1F( Form("hmgg%d",ic),Form("mgg%d",ic),
			   Commons::nBins,Commons::binning );
    hmgg_w[ic] = new TH1F( Form("hmgg_w%d",ic),Form("weighted mgg%d",ic),
			   Commons::nBins,Commons::binning );
  }
  //____________________________________________________________________//

  
  //==============================================================================//  
  //==================== Start the Loop over all entries =========================//
  //==============================================================================//  
  for(int jentry=0;jentry<m_nentries;jentry++){
    m_Rd->GetEntry(jentry);
    //--> Print the number of events
    AnalysisTools::Processing(jentry,m_nentries,(int)m_nentries/100);

    //--> Compute Run,lb,evt,npv
    ComputeBasicEventInfo();

    
    //======================= PU Weight for MC ONLY ===============================//
    if(!m_data)
      m_PUweight = m_pileupTool->GetCombinedWeight((int)m_Rd->GetVariable("RunNumber"),
						   (int)m_Rd->GetVariable("mc_channel_number"),
						   m_Rd->GetVariable("averageIntPerXing") );
    //___________________________________________________________________________//
 
     
 
    //================================ Diphotons selection ====================================//
    n[0]+=1;n_w[0]+=m_PUweight; 
    if(!EventTrigOK() ) continue;//--> Trigger
    n[1]+=1;n_w[1]+=m_PUweight;
    if ( !EventGRLOK() )      continue;//--> GRL 
    n[2]+=1;n_w[2]+=m_PUweight;   
    if( !EventPVOK() )    continue;//--> PV 
    n[3]+=1;n_w[3]+=m_PUweight;

    //----------------------- PRESELECTION CUTS  -----------------------------------
    int ilead=-99999; int isublead=-99999;
    double pt_lead=-99999; double pt_sublead=-99999;
    bool passPreSel = EventPreSelectionOK(&ilead,&isublead,&pt_lead,&pt_sublead);
    if( !passPreSel ) continue;//--> Preselection
    //--> Compute MC event weight
    ComputeIDScaleFactors(ilead,isublead,pt_lead,pt_sublead);
    weight = m_PUweight*m_SF_ID_lead*m_SF_ID_sublead;
    n[4]+=1;n_w[4]+=weight;
    hmgg[4]->Fill(mgg);hmgg_w[4]->Fill(mgg,weight);
    //----------------- PT CUT ------------------------
    if ( !PhotonPtOK(ilead,50) ) continue;//--> lead pt cut 
    if ( !PhotonPtOK(isublead,50) ) continue;//--> sublead pt cut 
    //--> Compute (most of the) output variables
    ComputeKinematicsOutput(ilead,isublead);
    n[5]+=1;n_w[5]+=weight;
    hmgg[5]->Fill(mgg);hmgg_w[5]->Fill(mgg,weight);
    //----------------- TIGHT CUT ------------------------
    if( !PhotonIsTightOK(ilead) )    continue;//--> Tight lead
    if( !PhotonIsTightOK(isublead) ) continue;//--> Tight sublead
    n[6]+=1;n_w[6]+=weight;
    hmgg[6]->Fill(mgg);hmgg_w[6]->Fill(mgg,weight);
    //----------------- ISOLATION CUT --------------------------------
    if( !PhotonIsolation_toolOK(ilead,8) )    continue;//--> iso lead
    if( !PhotonIsolation_toolOK(isublead,8) ) continue;//--> iso sublead
    n[7]+=1;n_w[7]+=weight;
    hmgg[7]->Fill(mgg);hmgg_w[7]->Fill(mgg,weight);
    //----------------- LARERROR CUT -------------------------------------
    if((int)m_Rd->GetVariable("larError") == 2 ) continue;//--> larError
    n[8]+=1;n_w[8]+=weight;
    hmgg[8]->Fill(mgg);hmgg_w[8]->Fill(mgg,weight);
    //----------------- TILEERROR CUT -------------------------------------
    if((int)m_Rd->GetVariable("tileError") == 2 ) continue;//--> tileError
    n[9]+=1;n_w[9]+=weight;
    hmgg[9]->Fill(mgg);hmgg_w[9]->Fill(mgg,weight);
    //------------------ OVERLAP REMOVAL -------------------------------------------------
    if( EventInZprimme(m_ee_events,m_RunNumber,m_EventNumber) ) continue;//--> Gee removal
    n[10]+=1;n_w[10]+=weight;
    hmgg[10]->Fill(mgg);hmgg_w[10]->Fill(mgg,weight);
    hmgg_final->Fill(mgg);hmgg_final_w->Fill(mgg,weight);
    //--------------------------------------------

    //--> Fill output tree
    tree_out->Fill();
    //--> truth mc_m branch for MC RS G* only 
    if( !m_data ){
      unsigned int mcsize=(unsigned int)m_Rd->GetVariable( "@mc_m.size()");
      for(unsigned int imc=0 ; imc<mcsize ;imc++ ){
	int pdgId=(int)m_Rd->GetVariable( Form("mc_pdgId[%d]",imc) );
	if(pdgId==5100039)
	  mymc_m=m_Rd->GetVariable(Form("mc_m[%d]",imc))/1000.;	
      }
    }

    //--------------------------------------------------------------------------
    if(m_data)
      PrintEventInfo(m_RunNumber,m_EventNumber,ilead,isublead,mgg,1000.);
    //-------------------------------------------------------------------------
  }
  //========================================================================//
  //========== End of the loop over all entries ============================//
  //========================================================================//


  //-----------------------------------------------------
  std::cout << "Initial    = " << (int)n[0] << std::endl;
  std::cout << "Trigger    = " << (int)n[1] << std::endl;
  std::cout << "GRL        = " << (int)n[2] << std::endl;
  std::cout << "PV         = " << (int)n[3] << std::endl;
  std::cout << "Preselec   = " << (int)n[4] << std::endl;
  std::cout << "Pt cut     = " << (int)n[5] << std::endl;
  std::cout << "Tight      = " << (int)n[6] << std::endl;
  std::cout << "Iso        = " << (int)n[7] << std::endl;
  std::cout << "LAr        = " << (int)n[8] << std::endl;
  std::cout << "Tile       = " << (int)n[9] << std::endl;
  std::cout << "ee removal = " << (int)n[10] << std::endl;
  if(!m_data){
    std::cout << "-----Weighted------"  << std::endl;
    std::cout << "Initial    = " << n_w[0] << std::endl;
    std::cout << "Trigger    = " << n_w[1] << std::endl;
    std::cout << "GRL        = " << n_w[1] << std::endl;
    std::cout << "PV         = " << n_w[3] << std::endl;
    std::cout << "Preselec   = " << n_w[4] << std::endl;
    std::cout << "Pt cut     = " << n_w[5] << std::endl;
    std::cout << "Tight      = " << n_w[6] << std::endl;
    std::cout << "Iso        = " << n_w[7] << std::endl;
    std::cout << "LAr        = " << n_w[8] << std::endl;
    std::cout << "Tile       = " << n_w[9] << std::endl;
    std::cout << "ee removal = " << n_w[10] << std::endl;
  }

  //------------------------------
  for(int ic=0;ic<nchecks;ic++){
    hCutFlow->SetBinContent(ic+1,n[ic]);
    hCutFlow_w->SetBinContent(ic+1,n_w[ic]);
    hist_mgg->Add(hmgg[ic]);
    hist_mgg_w->Add(hmgg_w[ic]);
  }
  //---------------------------------
  //  FillOutputFile();
}

///////////////////////////////////////
void GravitonSelector::FillOutputFile()
///////////////////////////////////////
{
  // TFile fout(output,"RECREATE");
  file_out->cd();
  hist_mgg->Write("mgg_cutflow",1);
  hist_mgg_w->Write("mgg_cutflow_w",1);
  hmgg_final->Write();
  hmgg_final_w->Write();
  hCutFlow->Write();
  hCutFlow_w->Write();
  //  tree_out->Write();
  file_out->Close();
}

/////////////////////////////////////////////////////////////////////////////
void GravitonSelector::ComputeIDScaleFactors( int ilead, int isublead,
					       double ptlead,double ptsublead)
/////////////////////////////////////////////////////////////////////////////
{
  if(!m_data){
    // double etas2_lead    = m_Rd->GetVariable( Form("ph_etas2[%d]",ilead) );
    // double etas2_sublead = m_Rd->GetVariable( Form("ph_etas2[%d]",isublead) );
    // int isconv_lead      = (int)m_Rd->GetVariable( Form("ph_isConv[%d]",ilead) );
    // int isconv_sublead   = (int)m_Rd->GetVariable( Form("ph_isConv[%d]",isublead) );
    // m_SF_ID_lead    = scaleForFFUnconvPhoton(ptlead,etas2_lead,isconv_lead);
    // m_SF_ID_sublead = scaleForFFUnconvPhoton(ptsublead,etas2_sublead,isconv_sublead);
    m_SF_ID_lead    = 1;
    m_SF_ID_sublead = 1;
  }else{
    m_SF_ID_lead    = 1;
    m_SF_ID_sublead = 1;
  }
}
///////////////////////////////////////////////
void GravitonSelector::ComputeBasicEventInfo()
//////////////////////////////////////////////
{
  m_RunNumber   = (int)m_Rd->GetVariable("RunNumber");
  m_EventNumber = (int)m_Rd->GetVariable("EventNumber");
  m_LumiBlock   = (int)m_Rd->GetVariable("lbn");
  m_NPV         = (int)m_Rd->GetVariable("@PV_nTracks.size()");
}
////////////////////////////////////////////////////////////////////////
void GravitonSelector::ComputeKinematicsOutput(int ilead, int isublead)
///////////////////////////////////////////////////////////////////////
{
  if(m_data){
    Lead_pT = PhotonRescaledPt(ilead);
    SubLead_pT = PhotonRescaledPt(isublead);
  }else{// Apply smearing for MC
    Lead_pT    = PhotonSmearedPt(ilead);
    SubLead_pT = PhotonSmearedPt(isublead);
    Lead_pT_smearedup    = PhotonSmearedPt(ilead,2);
    SubLead_pT_smearedup = PhotonSmearedPt(isublead,2);
    Lead_pT_smeareddown    = PhotonSmearedPt(ilead,1);
    SubLead_pT_smeareddown = PhotonSmearedPt(isublead,1);
  }

  Lead_eta       = m_Rd->GetVariable( Form("ph_etas2[%d]",ilead) );
  Lead_phi       = m_Rd->GetVariable( Form("ph_phi[%d]",ilead) );
  Lead_IsConv    = (int)m_Rd->GetVariable( Form("ph_isConv[%d]",ilead) );

  SubLead_eta    = m_Rd->GetVariable( Form("ph_etas2[%d]",isublead) );
  SubLead_phi    = m_Rd->GetVariable( Form("ph_phi[%d]",isublead) );
  SubLead_IsConv = (int)m_Rd->GetVariable( Form("ph_isConv[%d]",isublead) );
  
  double Lead_etaS1    = m_Rd->GetVariable( Form("ph_etas1[%d]",ilead) );
  double SubLead_etaS1 = m_Rd->GetVariable( Form("ph_etas1[%d]",isublead) );
  double PV_ID         = m_Rd->GetVariable("PV_z[0]");
  Lead_eta_PV          = InvMassCalc::EtaS1PVCorrected(Lead_etaS1,PV_ID);
  SubLead_eta_PV       = InvMassCalc::EtaS1PVCorrected(SubLead_etaS1,PV_ID);

  //-------------------------- Nominal ---------------------------------------//
  Lead_lv.SetPtEtaPhiM( Lead_pT,Lead_eta_PV,Lead_phi,0 );
  SubLead_lv.SetPtEtaPhiM( SubLead_pT,SubLead_eta_PV,SubLead_phi,0 );
  gamgam_lv    = Lead_lv + SubLead_lv;
  mgg          = gamgam_lv.M();
  ptgg         = gamgam_lv.Pt();
  ygg          = gamgam_lv.Rapidity();
  deltaphi     = Lead_lv.DeltaPhi(SubLead_lv);
  costhetastar = CosThetaStar_CS(Lead_lv,SubLead_lv);
  //-------------------------- Smeared Up ---------------------------------------//
  Lead_lv_smearedup.SetPtEtaPhiM( Lead_pT_smearedup,Lead_eta_PV,Lead_phi,0 );
  SubLead_lv_smearedup.SetPtEtaPhiM( SubLead_pT_smearedup,SubLead_eta_PV,SubLead_phi,0 );
  gamgam_lv_smearedup    = Lead_lv_smearedup + SubLead_lv_smearedup;
  mgg_smearedup          = gamgam_lv_smearedup.M();
  ptgg_smearedup         = gamgam_lv_smearedup.Pt();
  ygg_smearedup          = gamgam_lv_smearedup.Rapidity();
  costhetastar_smearedup = CosThetaStar_CS(Lead_lv_smearedup,SubLead_lv_smearedup);
  //-------------------------- Smeared Down ---------------------------------------//
  Lead_lv_smeareddown.SetPtEtaPhiM( Lead_pT_smeareddown,Lead_eta_PV,Lead_phi,0 );
  SubLead_lv_smeareddown.SetPtEtaPhiM( SubLead_pT_smeareddown,SubLead_eta_PV,SubLead_phi,0 );
  gamgam_lv_smeareddown    = Lead_lv_smeareddown + SubLead_lv_smeareddown;
  mgg_smeareddown          = gamgam_lv_smeareddown.M();
  ptgg_smeareddown         = gamgam_lv_smeareddown.Pt();
  ygg_smeareddown          = gamgam_lv_smeareddown.Rapidity();
  costhetastar_smeareddown = CosThetaStar_CS(Lead_lv_smeareddown,SubLead_lv_smeareddown);
  //-------------------------------------------------------------------------------------//  

}

/////////////////////////////////////////////////////////////////////
void GravitonSelector::PrintEventInfo( int RunNumber,int EventNumber,
					int LumiBlock,int ilead,
					int isublead,double mggval,
					double mggcut )
////////////////////////////////////////////////////////////////////
{

  if(mggval>mggcut){
    std::cout << "************************************************" << std::endl;
    std::cout <<  "RunNumber   = " << RunNumber << std::endl;
    std::cout <<  "EventNumber = " << EventNumber << std::endl;
    std::cout <<  "mgg         = " << mggval << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "Leader Photon" << std::endl;
    std::cout << "pT        = " << PhotonRescaledPt(ilead) << std::endl;
    std::cout << "eta       = " << m_Rd->GetVariable( Form("ph_etas2[%d]",ilead) ) << std::endl;
    std::cout << "phi       = " << m_Rd->GetVariable( Form("ph_phi[%d]",ilead) ) << std::endl;
    std::cout << "OQ        = " << (int)m_Rd->GetVariable( Form("ph_OQ[%d]",ilead) ) << std::endl;
    std::cout << "cone iso(tool) = " << PhotonIsolation_tool(ilead)/1000. << std::endl;
    std::cout << "cone iso(D3PD) = " << m_Rd->GetVariable( Form("ph_Etcone40_corrected[%d]",ilead) )/1000.
	      << std::endl;
    std::cout << "topo iso(tool) = " << PhotonTopoIsolation_tool(ilead)/1000. << std::endl;
    std::cout << "rhad1     = " << m_Rd->GetVariable(Form("ph_Ethad1[%d]",ilead))/PhotonRescaledPt(ilead) 
	      << std::endl;
    std::cout << "rhad      = " << m_Rd->GetVariable(Form("ph_Ethad[%d]",ilead))/PhotonRescaledPt(ilead)
	      << std::endl;
    std::cout << "e277      = " << m_Rd->GetVariable(Form("ph_E277[%d]",ilead)) <<std::endl; 
    std::cout << "reta      = " << m_Rd->GetVariable(Form("ph_reta[%d]",ilead)) <<std::endl; 
    std::cout << "rphi      = " << m_Rd->GetVariable(Form("ph_rphi[%d]",ilead)) <<std::endl; 
    std::cout << "weta2     = " << m_Rd->GetVariable(Form("ph_weta2[%d]",ilead))<< std::endl;
    std::cout << "f1        = " << m_Rd->GetVariable(Form("ph_f1[%d]",ilead)) << std::endl; 
    std::cout << "fside     = " << m_Rd->GetVariable(Form("ph_fside[%d]",ilead)) << std::endl;
    std::cout << "wtot      = " << m_Rd->GetVariable(Form("ph_wstot[%d]",ilead)) << std::endl;
    std::cout << "ws3       = " << m_Rd->GetVariable(Form("ph_ws3[%d]",ilead)) << std::endl; 
    std::cout << "DeltaE    = " 
	      << ( m_Rd->GetVariable(Form("ph_Emax2[%d]",ilead)) - m_Rd->GetVariable(Form("ph_Emins1[%d]",ilead)) )
	      << std::endl; 
    std::cout << "eratio    = " 
	      << ( m_Rd->GetVariable(Form("ph_emaxs1[%d]",ilead)) - m_Rd->GetVariable(Form("ph_Emax2[%d]",ilead)) )/( m_Rd->GetVariable(Form("ph_emaxs1[%d]",ilead)) + m_Rd->GetVariable(Form("ph_Emax2[%d]",ilead)) ) 
	      << std::endl; 
    std::cout << "isconv    = " << (int)m_Rd->GetVariable(Form("ph_isConv[%d]",ilead)) << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "SubLeader Photon" << std::endl;
    std::cout << "pT        = " << PhotonRescaledPt(isublead) << std::endl;
    std::cout << "eta       = " << m_Rd->GetVariable( Form("ph_etas2[%d]",isublead) ) << std::endl;
    std::cout << "phi       = " << m_Rd->GetVariable( Form("ph_phi[%d]",isublead) ) << std::endl;
    std::cout << "OQ        = " << (int)m_Rd->GetVariable( Form("ph_OQ[%d]",isublead) ) << std::endl;
    std::cout << "cone iso(tool) = " << PhotonIsolation_tool(isublead)/1000. << std::endl;
    std::cout << "cone iso(D3PD) = " <<m_Rd->GetVariable( Form("ph_Etcone40_corrected[%d]",isublead) )/1000.
	      << std::endl;
    std::cout << "topo iso(tool) = " << PhotonTopoIsolation_tool(isublead)/1000. << std::endl;
    std::cout << "rhad1     = " << m_Rd->GetVariable(Form("ph_Ethad1[%d]",isublead))/PhotonRescaledPt(isublead) 
	      << std::endl;
    std::cout << "rhad      = " << m_Rd->GetVariable(Form("ph_Ethad[%d]",isublead))/PhotonRescaledPt(isublead)
	      << std::endl;
    std::cout << "e277      = " << m_Rd->GetVariable(Form("ph_E277[%d]",isublead)) <<std::endl; 
    std::cout << "reta      = " << m_Rd->GetVariable(Form("ph_reta[%d]",isublead)) <<std::endl; 
    std::cout << "rphi      = " << m_Rd->GetVariable(Form("ph_rphi[%d]",isublead)) <<std::endl; 
    std::cout << "weta2     = " << m_Rd->GetVariable(Form("ph_weta2[%d]",isublead))<< std::endl;
    std::cout << "f1        = " << m_Rd->GetVariable(Form("ph_f1[%d]",isublead)) << std::endl; 
    std::cout << "fside     = " << m_Rd->GetVariable(Form("ph_fside[%d]",isublead)) << std::endl;
    std::cout << "wtot      = " << m_Rd->GetVariable(Form("ph_wstot[%d]",isublead)) << std::endl;
    std::cout << "ws3       = " << m_Rd->GetVariable(Form("ph_ws3[%d]",isublead)) << std::endl; 
    std::cout << "DeltaE    = " 
	      << ( m_Rd->GetVariable(Form("ph_Emax2[%d]",isublead)) - m_Rd->GetVariable(Form("ph_Emins1[%d]",isublead)) )
	      << std::endl; 
    std::cout << "eratio  = " 
	      << ( m_Rd->GetVariable(Form("ph_emaxs1[%d]",isublead)) - m_Rd->GetVariable(Form("ph_Emax2[%d]",isublead)) )/( m_Rd->GetVariable(Form("ph_emaxs1[%d]",isublead)) + m_Rd->GetVariable(Form("ph_Emax2[%d]",isublead)) ) 
	      << std::endl; 
    std::cout << "isconv  = " << (int)m_Rd->GetVariable(Form("ph_isConv[%d]",isublead)) << std::endl;
  }
}

/////////////////////////////////////////////////////////////////////////////
double GravitonSelector::CosThetaStar_CS(TLorentzVector v1,TLorentzVector v2)
////////////////////////////////////////////////////////////////////////////
{
  //---->  cosThetaStarCS ------------------------
  double l1_plus  = v1.E()+v1.Pz();
  double l1_minus = v1.E()-v1.Pz();
  double l2_plus  = v2.E()+v2.Pz();
  double l2_minus = v2.E()-v2.Pz();
  double num      = -(l2_plus*l1_minus-l1_plus*l2_minus);

  TLorentzVector v_sum = v1+v2;
  double denom    =  v_sum.M()*sqrt(pow(v_sum.M(),2) + 
					pow(v_sum.Pt(),2));
  return num/denom;
}
