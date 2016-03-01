#include "BkgSampleSelector.h"
#include "ToolsUtilities.h"
#include "ToolsCommons.h"
#include <TError.h>
#include <TF1.h>

////////////////////////////////////////////////////////////////////////////
BkgSampleSelector::BkgSampleSelector(TTree* tree) : GravitonSelector(tree)
///////////////////////////////////////////////////////////////////////////
{
  if (tree == 0){
    std::cout << "You need to specify an input tree !!!"<<std::endl;
  }
  Init();
}
/////////////////////////////////////////////////////////////
BkgSampleSelector::BkgSampleSelector() : GravitonSelector()
////////////////////////////////////////////////////////////
{
  //Default constructor
}
//////////////////////////////////////
BkgSampleSelector::~BkgSampleSelector()
/////////////////////////////////////
{
  //Default destructor

}
//////////////////////////////
void BkgSampleSelector::Init()
//////////////////////////////
{
  looseprime2    = 0x26fc01;
  looseprime3    = 0x24fc01;
  looseprime4    = 0x04fc01;
  looseprime5    = 0x00fc01;
  // looseprime2    = 2620417;
  // looseprime3    = 2489345;
  // looseprime4    = 392193;
  // looseprime5    = 130049;
  m_truthoutput  = "toto_truth.root";
}
/////////////////////////////////////
void BkgSampleSelector::InitOutTree()
/////////////////////////////////////
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
  tree_out->Branch("IsLoosePrime2_L",&Lead_IsLoosePrime2,"IsLoosePrime2_L/I");
  tree_out->Branch("IsLoosePrime3_L",&Lead_IsLoosePrime3,"IsLoosePrime3_L/I");
  tree_out->Branch("IsLoosePrime4_L",&Lead_IsLoosePrime4,"IsLoosePrime4_L/I");
  tree_out->Branch("IsLoosePrime5_L",&Lead_IsLoosePrime5,"IsLoosePrime5_L/I");
  tree_out->Branch("IsConv_L",&Lead_IsConv,"IsConv_L/I");
  tree_out->Branch("weight_L",&Lead_weight,"weight_L/I");
  //------> Subleader photon informations
  tree_out->Branch("pT_SL",&SubLead_pT,"pT_SL/D");
  tree_out->Branch("phi_SL",&SubLead_phi,"phi_SL/D");
  tree_out->Branch("eta_SL",&SubLead_eta,"eta_SL/D");
  tree_out->Branch("Iso_SL",&SubLead_iso,"Iso_SL/D");
  tree_out->Branch("IsTight_SL",&SubLead_IsTight,"IsTight_SL/I");
  tree_out->Branch("IsLoosePrime2_SL",&SubLead_IsLoosePrime2,"IsLoosePrime2_SL/I");
  tree_out->Branch("IsLoosePrime3_SL",&SubLead_IsLoosePrime3,"IsLoosePrime3_SL/I");
  tree_out->Branch("IsLoosePrime4_SL",&SubLead_IsLoosePrime4,"IsLoosePrime4_SL/I");
  tree_out->Branch("IsLoosePrime5_SL",&SubLead_IsLoosePrime5,"IsLoosePrime5_SL/I");
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

  if( !m_data ){
    //----> Event weight
    tree_out->Branch("weight",&weight,"weight/D");
    tree_out->Branch("gen_weight",&gen_weight,"gen_weight/D");
    tree_out->Branch("nlo_weight",&nlo_weight,"nlo_weight/D");
    tree_out->Branch("nlo_weight_xabier",&nlo_weight_xabier,"nlo_weight_xabier/D");
    tree_out->Branch("Lead_truth_eta",&Lead_truth_eta,"Lead_truth_eta/D");
    tree_out->Branch("Lead_truth_phi",&Lead_truth_phi,"Lead_truth_phi/D");
    tree_out->Branch("Lead_truth_pT",&Lead_truth_pT,"Lead_truth_pT/D");
    tree_out->Branch("Lead_truth_status",&Lead_truth_status,"Lead_truth_status/I");
    tree_out->Branch("SubLead_truth_eta",&SubLead_truth_eta,"SubLead_truth_eta/D");
    tree_out->Branch("SubLead_truth_phi",&SubLead_truth_phi,"SubLead_truth_phi/D");
    tree_out->Branch("SubLead_truth_pT",&SubLead_truth_pT,"SubLead_truth_pT/D");
    tree_out->Branch("SubLead_truth_status",&SubLead_truth_status,"SubLead_truth_status/I");
    tree_out->Branch("truth_mgg",&truth_mgg,"truth_mgg/D");
    tree_out->Branch("truth_ptgg",&truth_ptgg,"truth_ptgg/D");
    tree_out->Branch("truth_ygg",&truth_ygg,"truth_ygg/D");
    tree_out->Branch("truth_dphigg",&truth_dphigg,"truth_dphigg/D");
    tree_out->Branch("truth_costhetastar",&truth_costhetastar,"truth_costhetastar/D");
  }

}
/////////////////////////////////////
void BkgSampleSelector::InitTruthTree()
/////////////////////////////////////
{
  file_truth = new TFile(m_truthoutput,"RECREATE");
  tree_truth = new TTree("tree_truth","tree_truth");
  tree_truth->Branch("gen_weight",&gen_weight,"gen_weight/D");
  tree_truth->Branch("nlo_weight",&nlo_weight,"nlo_weight/D");
  tree_truth->Branch("nlo_weight_xabier",&nlo_weight_xabier,"nlo_weight_xabier/D");
  tree_truth->Branch("Lead_truth_eta",&Lead_truth_eta,"Lead_truth_eta/D");
  tree_truth->Branch("Lead_truth_phi",&Lead_truth_phi,"Lead_truth_phi/D");
  tree_truth->Branch("Lead_truth_pT",&Lead_truth_pT,"Lead_truth_pT/D");
  tree_truth->Branch("Lead_truth_status",&Lead_truth_status,"Lead_truth_status/I");
  tree_truth->Branch("SubLead_truth_eta",&SubLead_truth_eta,"SubLead_truth_eta/D");
  tree_truth->Branch("SubLead_truth_phi",&SubLead_truth_phi,"SubLead_truth_phi/D");
  tree_truth->Branch("SubLead_truth_pT",&SubLead_truth_pT,"SubLead_truth_pT/D");
  tree_truth->Branch("SubLead_truth_status",&SubLead_truth_status,"SubLead_truth_status/I");
  tree_truth->Branch("truth_mgg",&truth_mgg,"truth_mgg/D");
  tree_truth->Branch("truth_ptgg",&truth_ptgg,"truth_ptgg/D");
  tree_truth->Branch("truth_ygg",&truth_ygg,"truth_ygg/D");
  tree_truth->Branch("truth_dphigg",&truth_dphigg,"truth_dphigg/D");
  tree_truth->Branch("truth_costhetastar",&truth_costhetastar,"truth_costhetastar/D");
}


/////////////////////////////////////////////////////////////////
void BkgSampleSelector::EventLoop(TString outfile,int Nrelaxedcut)
////////////////////////////////////////////////////////////////
{
  m_output = outfile;
  InitOutTree();


  unsigned int looseprimemask=0;
  if(Nrelaxedcut==2)
    looseprimemask=looseprime2;
  else if(Nrelaxedcut==3)
    looseprimemask=looseprime3;
  else if(Nrelaxedcut==4)
    looseprimemask=looseprime4;
  else if(Nrelaxedcut==5)
    looseprimemask=looseprime5;
  else
    looseprimemask=0;


  //==============================================================================//  
  //==================== Start the Loop over all entries =========================//
  //==============================================================================//  
  for(int jentry=0;jentry<m_nentries;jentry++){
    m_Rd->GetEntry(jentry);
    //--> Print the number of events
    AnalysisTools::Processing(jentry,m_nentries);

    bool truth_info = EventTruthSelector(jentry,m_data);
    if( !m_data && !truth_info) continue; //--> Skip MC events without truth matching

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
      gen_weight = Commons::GetMCWeight( (int)m_Rd->GetVariable("mc_channel_number") );
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


    // if( ((unsigned int)m_Rd->GetVariable(Form("ph_isEM[%d]",ilead))&looseprimemask)!=0 )    continue; 
    // if( ((unsigned int)m_Rd->GetVariable(Form("ph_isEM[%d]",isublead))&looseprimemask)!=0 ) continue; 

    // if( (PhotonIsEM(ilead)&looseprimemask)!=0 )    continue; 
    // if( (PhotonIsEM(isublead)&looseprimemask)!=0 ) continue; 


    //----> LoosePrime computation
    Lead_IsLoosePrime2    = (int)( (PhotonIsEM(ilead)&looseprime2) != 0 ? 1:0 );
    Lead_IsLoosePrime3    = (int)( (PhotonIsEM(ilead)&looseprime3) != 0 ? 1:0 );
    Lead_IsLoosePrime4    = (int)( (PhotonIsEM(ilead)&looseprime4) != 0 ? 1:0 );
    Lead_IsLoosePrime5    = (int)( (PhotonIsEM(ilead)&looseprime5) != 0 ? 1:0 );
    SubLead_IsLoosePrime2  = (int)( (PhotonIsEM(isublead)&looseprime2) != 0 ? 1:0 );
    SubLead_IsLoosePrime3  = (int)( (PhotonIsEM(isublead)&looseprime3) != 0 ? 1:0 );
    SubLead_IsLoosePrime4  = (int)( (PhotonIsEM(isublead)&looseprime4) != 0 ? 1:0 );
    SubLead_IsLoosePrime5  = (int)( (PhotonIsEM(isublead)&looseprime5) != 0 ? 1:0 );
    //----> Tight computation
    Lead_IsTight    = (int)PhotonIsTightOK(ilead);
    SubLead_IsTight = (int)PhotonIsTightOK(isublead);
    if(m_isotype == "CONE"){
      Lead_iso    = PhotonIsolation_tool(ilead)/1000.;
      SubLead_iso = PhotonIsolation_tool(isublead)/1000.;
    }else if(m_isotype == "TOPO"){
      Lead_iso    = PhotonTopoIsolation_tool(ilead)/1000.;
      SubLead_iso = PhotonTopoIsolation_tool(isublead)/1000.;
    }else Fatal("BkgSampleSelector::EventLoop","Wrong isolation type !!!");

 
    //--> Fill output tree
    tree_out->Fill();
    //------------------
  }
  //========================================================================//
  //========== End of the loop over all entries ============================//
  //========================================================================//

  FillOutputFile();
}
///////////////////////////////////////////////////////
void BkgSampleSelector::TruthEventLoop(TString outfile)
///////////////////////////////////////////////////////
{
  m_truthoutput = outfile;
  InitTruthTree();

  //==============================================================================//  
  //==================== Start the Loop over all entries =========================//
  //==============================================================================//  
  for(int jentry=0;jentry<m_nentries;jentry++){
    m_Rd->GetEntry(jentry);
    //--> Print the number of events
    AnalysisTools::Processing(jentry,m_nentries);
    gen_weight = Commons::GetMCWeight( (int)m_Rd->GetVariable("mc_channel_number") );
    bool truth_info = EventTruthSelector(jentry,m_data);
    if( truth_info )
      tree_truth->Fill();
  }
  //========================================================================//
  //========== End of the loop over all entries ============================//
  //========================================================================//
  FillOutputTruthFile();
}
////////////////////////////////////////
void BkgSampleSelector::FillOutputFile()
////////////////////////////////////////
{
  // TFile file_out(output,"RECREATE");
  file_out->cd();
  tree_out->Write();
  file_out->Close();
}
////////////////////////////////////////////
void BkgSampleSelector::FillOutputTruthFile()
////////////////////////////////////////////
{
  file_truth->cd();
  tree_truth->Write();
  file_truth->Close();
}
///////////////////////////////////////////////////////////////
bool BkgSampleSelector::EventTruthSelector(int entry,bool data)
//////////////////////////////////////////////////////////////
{

 if(data) return false;

  // int n_truth = (int)m_Rd->GetVariable("egammaTruth_n");
  // std::vector<int> index_truth_old;
  // for ( int itruth=0 ; itruth<n_truth;itruth++ ){
  //   int isHP = (int)m_Rd->GetVariable( Form("egammaTruth_isHardProcPhoton[%d]",itruth) );
  //   if( isHP != 1 ) continue;
  //   index_truth_old.push_back( (int)m_Rd->GetVariable( Form("egammaTruth_mc_index[%d]",itruth) ) );
  // }


  m_Ts->GetEntry(entry);
  std::vector<int> index_truth = m_Ts->GetDirectGamGamChildIndex();
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
  truth_costhetastar = CosThetaStar_CS(Lead_truth_lv,SubLead_truth_lv);

  nlo_weight = GetkFactor_40_30(truth_mgg);
  nlo_weight_xabier = GetkFactor_Xabier(truth_mgg);
  
  return true;
}
////////////////////////////////////////////////////////
double BkgSampleSelector::GetkFactor_Xabier(double mass)
////////////////////////////////////////////////////////
{
  TString filename = "rootfiles/kFactor_NLO_LO_rel16.root";
  TFile *f = new TFile(filename,"read");
  TVectorT<float>* vv_mass       = (TVectorT<float>*)f->Get("vv_mass");
  TVectorT<float>* vv_mass_er    = (TVectorT<float>*)f->Get("vv_mass_er");
  TVectorT<float>* vv_kFactorEff = (TVectorT<float>*)f->Get("vv_kFactorEff");

  int Nrow = vv_mass->GetNrows();
  double kFac = -999.;
  for( int irow = 0 ;irow < Nrow ;irow++){
    double temp_mass    = (*vv_mass)[irow];
    double temp_mass_er = (*vv_mass_er)[irow];
    double low = temp_mass - temp_mass_er;
    double up  = temp_mass + temp_mass_er;
    if( mass < up && mass >= low){
      kFac = (*vv_kFactorEff)[irow];
      break;
    }
  }
  f->Close();
  if( mass > 1300 )
    kFac=1;
  if ( kFac==-999.)
      Fatal("BkgSampleSelector::GetIrrkFactor", "No k factor computed !");

  return kFac;
}
///////////////////////////////////////////////////////
double BkgSampleSelector::GetkFactor_40_30(double mass)
///////////////////////////////////////////////////////
{
  TF1 fk("fk","exp([0]+[1]*x)+[2]",120,40000);
  fk.SetParameters(2.87873e-01,-1.16827e-03,7.74223e-01);//Fit 140-1500
  // fk.SetParameters(4.16025e-01,-9.12814e-04,5.61667e-01);//Fit 120-1500
  // fk.SetParameters(1.02548e-01,-1.49308e-03,1);//Fit 120-1500 with [2] fixed to 1

  if( mass < 120.) return 1.;
  else return fk.Eval(mass);

}
///////////////////////////////////////////////////////////////////////
double BkgSampleSelector::GetkFactor_50_50(double mass,TString config)
//////////////////////////////////////////////////////////////////////
{
  // fit with a 2nd order polynom up to 250
  TF1 fk_low("fk_low","[0]*x*x+[1]*x+[2]",120,250);
  TF1 fk_high("fk_high","exp([0]+[1]*x)+[2]",250,4000000);

  if(config == "nominal"){
    fk_low.SetParameters(-1.13633e-05,0.00494087,1.11658);
    fk_high.SetParameters(0.0213382,-0.00143431,0.942951);
  }else if( config == "var_pdfset" ){
    fk_low.SetParameters(-1.13753e-05,0.00467289,1.4204);
    fk_high.SetParameters(-0.063475,-0.00271677,1.41882);
  }else if( config == "var_pdfeig" ){
  fk_low.SetParameters(-1.13411e-05,0.0049442,1.15587);
  fk_high.SetParameters(-0.107886,-0.00188615,1.14368);
  }else if( config == "var_scale" ){
    fk_low.SetParameters(-1.31102e-05,0.00559751,1.17581);
    fk_high.SetParameters(0.0609006,-0.00166052,1.07523);
  }else if ( config == "var_iso" ){
    fk_low.SetParameters(-1.40246e-05,0.00601151,1.0313);
    fk_high.SetParameters(0.0237086,-0.0014787,0.975151);
  }else Fatal("BkgSampleSelector::GetkFactor_50_50()","Wrong config TString : "+config);

  if( mass < 120.) return 1.;
  else if(mass < 250.) return fk_low.Eval(mass);
  else return fk_high.Eval(mass);

}



