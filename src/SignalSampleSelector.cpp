#include "SignalSampleSelector.h"
#include "SignalTemplateCreator.h"

#include "ToolsUtilities.h"
#include "ToolsCommons.h"
#include "InvMassCalculator.h"
#include <TError.h>
//////////////////////////////////////////////////////////////////////////////////
SignalSampleSelector::SignalSampleSelector(TTree* tree) : GravitonSelector(tree)
//////////////////////////////////////////////////////////////////////////////////
{
   if (tree == 0) {
     std::cout << "You need to specify an input tree !!!"<<std::endl;
   }
}
////////////////////////////////////////////////////////////////////
SignalSampleSelector::SignalSampleSelector() : GravitonSelector()
///////////////////////////////////////////////////////////////////
{
  //Default constructor
}
///////////////////////////////////////////////
SignalSampleSelector::~SignalSampleSelector()
///////////////////////////////////////////////
{
  //Default destructor
}
//////////////////////////////////////////////
void SignalSampleSelector::Init()
/////////////////////////////////////////////
{


}
////////////////////////////////////////
void SignalSampleSelector::InitOutTree()
////////////////////////////////////////
{
    //-----> Output variables  
  file_out   = new TFile(m_output,"RECREATE");
  tree_out   = new TTree("tree","tree");

  //------> Event informations
  tree_out->Branch("RunNumber",&m_RunNumber,"RunNumber/I");
  tree_out->Branch("EventNumber",&m_EventNumber,"EventNumber/I");
  tree_out->Branch("LumiBlock",&m_LumiBlock,"LumiBlock/I");
  tree_out->Branch("NPV",&m_NPV,"NPV/I");
  //------> Leader photon informations
  tree_out->Branch("pT_L",&Lead_pT,"pT_L/D");
  tree_out->Branch("phi_L",&Lead_phi,"phi_L/D");
  tree_out->Branch("eta_L",&Lead_eta,"eta_L/D");
  tree_out->Branch("Iso_L",&Lead_iso,"Iso_L/D");
  tree_out->Branch("IsTight_L",&Lead_IsTight,"IsTight_L/I");
  tree_out->Branch("IsConv_L",&Lead_IsConv,"IsConv_L/I");
  tree_out->Branch("weight_L",&Lead_weight,"weight_L/I");
  //------> Subleader photon informations
  tree_out->Branch("pT_SL",&SubLead_pT,"pT_SL/D");
  tree_out->Branch("phi_SL",&SubLead_phi,"phi_SL/D");
  tree_out->Branch("eta_SL",&SubLead_eta,"eta_SL/D");
  tree_out->Branch("Iso_SL",&SubLead_iso,"Iso_SL/D");
  tree_out->Branch("IsTight_SL",&SubLead_IsTight,"IsTight_SL/I");
  tree_out->Branch("IsConv_SL",&SubLead_IsConv,"IsConv_SL/I");
  tree_out->Branch("weight_SL",&SubLead_weight,"weight_SL/I");
  //------> Invariant mass calculation
  tree_out->Branch("mgg",&mgg,"mgg/D");
  //------> ptgg calculation
  tree_out->Branch("ptgg",&ptgg,"ptgg/D");
  //------> ygg calculation
  tree_out->Branch("ygg",&ygg,"ygg/D");
  //------> Cos(theta*) calculation
  tree_out->Branch("costhetastar",&costhetastar,"costhetastar/D");
  //------> Deltaphi calculation
  tree_out->Branch("deltaphi",&deltaphi,"deltaphi/D");

    //----> Event weight
  tree_out->Branch("weight",&weight,"weight/D");
  tree_out->Branch("gen_weight",&gen_weight,"gen_weight/D");
  tree_out->Branch("Lead_truth_eta",&Lead_truth_eta,"Lead_truth_eta/D");
  tree_out->Branch("Lead_truth_phi",&Lead_truth_phi,"Lead_truth_phi/D");
  tree_out->Branch("Lead_truth_pT",&Lead_truth_pT,"Lead_truth_pT/D");
  tree_out->Branch("Lead_truth_status",&Lead_truth_status,"Lead_truth_status/I");
  tree_out->Branch("SubLead_truth_eta",&SubLead_truth_eta,"SubLead_truth_eta/D");
  tree_out->Branch("SubLead_truth_phi",&SubLead_truth_phi,"SubLead_truth_phi/D");
  tree_out->Branch("SubLead_truth_pT",&SubLead_truth_pT,"SubLead_truth_pT/D");
  tree_out->Branch("SubLead_truth_status",&SubLead_truth_status,"SubLead_truth_status/I");
  //----> True mass (MC only)
  tree_out->Branch("mc_m",&mymc_m,"mc_m/D");
  tree_out->Branch("truth_mgg",&truth_mgg,"truth_mgg/D");
  tree_out->Branch("truth_ptgg",&truth_ptgg,"truth_ptgg/D");
  tree_out->Branch("truth_ygg",&truth_ygg,"truth_ygg/D");
  tree_out->Branch("truth_dphigg",&truth_dphigg,"truth_dphigg/D");
  tree_out->Branch("truth_costhetastar",&truth_costhetastar,"truth_costhetastar/D");


}
/////////////////////////////////////////////////////
void SignalSampleSelector::EventLoop(TString outfile)
//////////////////////////////////////////////////////
{

  m_output = outfile;
  InitOutTree();



  //==============================================================================//  
  //==================== Start the Loop over all entries =========================//
  //==============================================================================//  
  for(int jentry=0;jentry<m_nentries;jentry++){
    m_Rd->GetEntry(jentry);
    //--> Print the number of events
    AnalysisTools::Processing(jentry,m_nentries,(int)m_nentries/100);

    bool truth_info = EventTruthSelector(jentry,m_data);
    if( !m_data && !truth_info) continue; //--> Skip MC events without truth matching



    unsigned int mcsize=(unsigned int)m_Rd->GetVariable( "@mc_m.size()");
    for(unsigned int imc=0 ; imc<mcsize ;imc++ ){
      int pdgId=(int)m_Rd->GetVariable( Form("mc_pdgId[%d]",imc) );
      if(pdgId==5100039)
	mymc_m=m_Rd->GetVariable(Form("mc_m[%d]",imc))/1000.; 
    }





    //--> Compute run,lb,evt,npv
    ComputeBasicEventInfo();
    
    //======================= PU Weight for MC ONLY ===============================//
    if(!m_data)
      m_PUweight = m_pileupTool->GetCombinedWeight( (int)m_Rd->GetVariable("RunNumber"),
                                                    (int)m_Rd->GetVariable("mc_channel_number"),
                                                    m_Rd->GetVariable("averageIntPerXing") );
    //___________________________________________________________________________//

    //======================= Generator Weight for MC ONLY =======================//
    if(m_data)
      gen_weight = 1;
    else 
      gen_weight = SignalTemplateCreator::GetGravitonWeight(mymc_m,2000,0.10 );
    //___________________________________________________________________________//
 

    //================================ Diphotons selection ====================================//
    if( !EventTrigOK() ) continue;//--> Trigger
    if( !EventGRLOK() )  continue;//--> GRL 
    if( !EventPVOK() )   continue;//--> PV 
    //----------------------- PRESELECTION CUTS  --------------------------------------
    int ilead=-99999; int isublead=-99999;
    double ptlead=-99999; double ptsublead=-99999;
    bool passPreSel = EventPreSelectionOK(&ilead,&isublead,&ptlead,&ptsublead);
    if( !passPreSel ) continue;//--> Preselection

    //--> Compute MC event weight
    ComputeIDScaleFactors(ilead,isublead,ptlead,ptsublead);
    weight = m_PUweight*m_SF_ID_lead*m_SF_ID_sublead;
    Lead_weight    = m_PUweight*m_SF_ID_lead;
    SubLead_weight = m_PUweight*m_SF_ID_sublead;

    //----------------- PT CUT ----------------------------
    if ( !PhotonPtOK(ilead,40) ) continue;//--> lead pt cut 
    if ( !PhotonPtOK(isublead,30) ) continue;//--> sublead pt cut 
    //--> Compute (most of the) output variables
    ComputeKinematicsOutput(ilead,isublead);
    //----------------- LARERROR CUT --------------------------------
    if((int)m_Rd->GetVariable("larError") ==2 ) continue;//--> larError
    //----------------- TILEERROR CUT --------------------------------
    if((int)m_Rd->GetVariable("tileError") ==2 ) continue;//--> tile  Error
    //------------------ OVERLAP REMOVAL -------------------------------------------------
    if( EventInZprimme(m_ee_events,m_RunNumber,m_EventNumber) ) continue;//--> Gee removal
    //-----------------------------------------------------------------------------------


    //----> Tight computation
    Lead_IsTight    = (int)PhotonIsTightOK(ilead);
    SubLead_IsTight = (int)PhotonIsTightOK(isublead);
    if(m_isotype == "CONE"){
      Lead_iso    = PhotonIsolation_tool(ilead)/1000.;
      SubLead_iso = PhotonIsolation_tool(isublead)/1000.;
    }else if(m_isotype == "TOPO"){
      Lead_iso    = PhotonTopoIsolation_tool(ilead)/1000.;
      SubLead_iso = PhotonTopoIsolation_tool(isublead)/1000.;
    }else if(m_isotype == "TOPO_PTDEP"){
      Lead_iso    = PhotonTopoIsolationPtdep_tool(ilead)/1000.;
      SubLead_iso = PhotonTopoIsolationPtdep_tool(isublead)/1000.;
    }else Fatal("SignalSampleSelector::EventLoop","Wrong isolation type !!!");

    //--> Fill output tree
    tree_out->Fill();
    //------------------
  }
  //========================================================================//
  //========== End of the loop over all entries ============================//
  //========================================================================//
  FillOutputFile();

}

///////////////////////////////////////////////////////////////
bool SignalSampleSelector::EventTruthSelector(int entry,bool data)
//////////////////////////////////////////////////////////////
{

 if(data) return false;


  m_Ts->GetEntry(entry);
  std::vector<int> index_truth = m_Ts->GetRSGravChildIndex();
  double pt1 = m_Rd->GetVariable(Form("mc_pt[%d]",index_truth[0]));
  double pt2 = m_Rd->GetVariable(Form("mc_pt[%d]",index_truth[1]));
  int ilead_truth =-1; int isublead_truth=-1;
  if( index_truth.size() !=2 ) return false;
  if(pt1>pt2) {
    ilead_truth = index_truth[0];
    isublead_truth = index_truth[1];
  }else{
    ilead_truth = index_truth[1];
    isublead_truth = index_truth[0];
  }

  Lead_truth_eta     = m_Rd->GetVariable(Form("mc_eta[%d]",ilead_truth));
  Lead_truth_phi     = m_Rd->GetVariable(Form("mc_phi[%d]",ilead_truth));
  Lead_truth_pT      = m_Rd->GetVariable(Form("mc_pt[%d]",ilead_truth))/1000.;
  Lead_truth_status  = (int)m_Rd->GetVariable(Form("mc_status[%d]",ilead_truth));

  SubLead_truth_eta     = m_Rd->GetVariable(Form("mc_eta[%d]",isublead_truth));
  SubLead_truth_phi     = m_Rd->GetVariable(Form("mc_phi[%d]",isublead_truth));
  SubLead_truth_pT      = m_Rd->GetVariable(Form("mc_pt[%d]",isublead_truth))/1000.;
  SubLead_truth_status  = (int)m_Rd->GetVariable(Form("mc_status[%d]",ilead_truth));

  double truth_PV_ID = m_Rd->GetVariable("mc_vx_z[0]"); 
  TLorentzVector Lead_truth_lv;TLorentzVector SubLead_truth_lv;
  TLorentzVector gamgam_truth_lv;

  Lead_truth_lv.SetPtEtaPhiM( Lead_truth_pT,
                              InvMassCalc::EtaS1PVCorrected(Lead_truth_eta,truth_PV_ID),
                              Lead_truth_phi,0 );
  SubLead_truth_lv.SetPtEtaPhiM( SubLead_truth_pT,
                                 InvMassCalc::EtaS1PVCorrected(SubLead_truth_eta,truth_PV_ID),
                                 SubLead_truth_phi,0 );
  gamgam_truth_lv = Lead_truth_lv + SubLead_truth_lv;
  truth_mgg      = gamgam_truth_lv.M();
  truth_ptgg     = gamgam_truth_lv.Pt();
  truth_ygg      = gamgam_truth_lv.Rapidity();
  truth_dphigg = Lead_lv.DeltaPhi(SubLead_lv);

  //---->  cosThetaStarCS ------------------------
  double l1_plus  = Lead_truth_lv.E()+Lead_truth_lv.Pz();
  double l1_minus = Lead_truth_lv.E()-Lead_truth_lv.Pz();
  double l2_plus  = SubLead_truth_lv.E()+SubLead_truth_lv.Pz();
  double l2_minus = SubLead_truth_lv.E()-SubLead_truth_lv.Pz();
  double num      = -(l2_plus*l1_minus-l1_plus*l2_minus);
  double denom    =  gamgam_truth_lv.M()*sqrt(pow(gamgam_truth_lv.M(),2) + 
                                              pow(gamgam_truth_lv.Pt(),2));
  truth_costhetastar = num/denom;
  
  return true;
}
////////////////////////////////////////////
void SignalSampleSelector::FillOutputFile()
///////////////////////////////////////////
{
  // TFile file_out(output,"RECREATE");
  file_out->cd();
  tree_out->Write();
  file_out->Close();
}
