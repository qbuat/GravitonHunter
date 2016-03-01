#include "PersoCleaningStudies.h"
#include "ToolsUtilities.h"
#include "ToolsCommons.h"
#include <TError.h>
#include <TVectorT.h>
#include "GoodRunsLists/DQHelperFunctions.h"
//////////////////////////////////////////////////////////////
CleaningStudies::CleaningStudies(TTree* tree) : GravitonAnalysis(tree)
//////////////////////////////////////////////////////////////
{
   if (tree == 0) {
     std::cout << "You need to specify an input tree !!!"<<std::endl;
   }
   Init(tree);
}
///////////////////////////////////////////////
CleaningStudies::CleaningStudies() : GravitonAnalysis()
///////////////////////////////////////////////
{
  //Default constructor
}
///////////////////////////
CleaningStudies::~CleaningStudies()
///////////////////////////
{
  //Default destructor
  delete m_IMC;
}
///////////////////////////////////
void CleaningStudies::Init(TTree* tree)
//////////////////////////////////
{
  m_passTrigger  = false;
  m_inGRL        = false;
  m_passPV       = false;
  m_passPreSel   = false;

  m_output       = "toto.root";
  m_PU_mcfile    = Commons::PileupMCFile;
  m_PU_datafile  = Commons::PileupDataFile;
  m_IMC          = new InvMassCalc(tree);

  // SetPileUpTool(m_PU_mcfile,m_PU_datafile);

}

///////////////////////////////////////////////////////////
void CleaningStudies::InvMassCleaning(bool data,TString output)
///////////////////////////////////////////////////////////
{
  m_data=data;

  //=====Declare output variables====//
  TTree *tree=new TTree("tree","tree");
  int myRunNumber=0;
  int myEventNumber=0;
  int myLumiBlock=0;
  int myLArError=0;
  double mggGEV=0;
  double LeadpT      = 0;
  double Leadetas2   = 0;
  double Leadphi     = 0;
  double Leadrphi    = 0;
  double Leadreta    = 0;
  double Leadtime    = 0;
  int    LeadPassLC  = 0;
  int    LeadPassLC2 = 0;
  int    LeadPassJC  = 0;

  double SubLeadpT      = 0;
  double SubLeadetas2   = 0;
  double SubLeadphi     = 0;
  double SubLeadrphi    = 0;
  double SubLeadreta    = 0;
  double SubLeadtime    = 0;
  int    SubLeadPassLC  = 0;
  int    SubLeadPassLC2 = 0;
  int    SubLeadPassJC  = 0;
  //__________________//
  tree->Branch("mgg",&mggGEV,"mgg/D");
  tree->Branch("RunNumber",&myRunNumber,"RunNumber/I");
  tree->Branch("EventNumber",&myEventNumber,"EventNumber/I");
  tree->Branch("LumiBlock",&myLumiBlock,"LumiBlock/I");
  tree->Branch("larError",&myLArError,"larError/I");
  tree->Branch("Lead_pT",&LeadpT,"Lead_pT/D");
  tree->Branch("Lead_etas2",&Leadetas2,"Lead_etas2/D");
  tree->Branch("Lead_phi",&Leadphi,"Lead_phi/D");
  tree->Branch("Lead_PassLC",&LeadPassLC,"Lead_PassLC/I");
  tree->Branch("Lead_PassLC2",&LeadPassLC2,"Lead_PassLC2/I");
  tree->Branch("Lead_PassJC",&LeadPassJC,"Lead_PassJC/I");
  tree->Branch("SubLead_pT",&SubLeadpT,"SubLead_pT/D");
  tree->Branch("SubLead_etas2",&SubLeadetas2,"SubLead_etas2/D");
  tree->Branch("SubLead_phi",&SubLeadphi,"SubLead_phi/D");
  tree->Branch("SubLead_PassLC",&SubLeadPassLC,"SubLead_PassLC/I");
  tree->Branch("SubLead_PassLC2",&SubLeadPassLC2,"SubLead_PassLC2/I");
  tree->Branch("SubLead_PassJC",&SubLeadPassJC,"SubLead_PassJC/I");
  //________________________________//
  
  //==Declare the GRL object==//
  DQ::SetXMLFile(m_GRL);
  //__________________________//
  
  
    //======================================================================//  
  //====== Start the Loop over all entries ===============================//
  //======================================================================//  
  for(int jentry=0;jentry<m_nentries;jentry++){
    m_Rd->GetEntry(jentry);
  
   // ===clean all flags=======//     
    m_passTrigger       = false ;
    m_inGRL             = false ;
    m_passPV            = false ;
    m_passPreSel        = false ;
    //__________________________//
 
   //============ GRL Selection ==================//
    int Run=(int)m_Rd->GetVariable("RunNumber");
    int lbn=(int)m_Rd->GetVariable("lbn");
    if(m_data){
      m_inGRL = DQ::PassRunLB(Run ,lbn);
    }else{
      m_inGRL=true;
    }
    //___________________________________________//
    
 

    //=================== Trigger Requirement ====================//
    //For Photon D3PDs
    m_passTrigger    = (bool)m_Rd->GetVariable("EF_2g20_loose");
    //__________________________________________________________//
 
    
    //==============Primary Vertex Requirement==========//
    unsigned int n_PV=(unsigned int)m_Rd->GetVariable("@PV_nTracks.size()");
    for(unsigned int i=0;i<n_PV;i++){
      if( (int)m_Rd->GetVariable(Form("PV_nTracks[%d]",i))>2){
	  m_passPV=true;
	  break;
      }
    }
    //_________________________________________________//
    

    //=============================== Diphotons selection ==========================//
    if(!m_passTrigger) continue;
    if (!m_inGRL)    continue; 
    if( !m_passPV )    continue; 
    double pt_lead=-99999; double pt_sublead=-99999;
    int ilead=-99999; int isublead=-99999;
    m_passPreSel = EventPreSelectionOK(&ilead,&isublead,&pt_lead,&pt_sublead);
    if( !m_passPreSel ) continue;
    if( !PhotonIsTightOK(ilead) ) continue;
    if( !PhotonIsTightOK(isublead) ) continue;
    if( !PhotonIsolation_toolOK(ilead) )    continue;
    if( !PhotonIsolation_toolOK(isublead) ) continue;


    mggGEV=m_IMC->GetInvMass(jentry,ilead,isublead)/1000.;
    if(mggGEV>120.){
      LeadPassJC=1;
      if( !PassNoiseCleaning(ilead) ) LeadPassJC=0;
      SubLeadPassJC=1;
      if( !PassNoiseCleaning(isublead) ) SubLeadPassJC=0;
      LeadPassLC=1;
      if( !PassLArCleaning(ilead) ) LeadPassLC=0;
      SubLeadPassLC=1;
      if( !PassLArCleaning(isublead) ) SubLeadPassLC=0;
      LeadPassLC2=1;
      if( !PassPhotonCleaning(ilead) ) LeadPassLC2=0;
      SubLeadPassLC=1;
      if( !PassPhotonCleaning(isublead) ) SubLeadPassLC2=0;
      if(m_data)
	LeadpT = PhotonRescaledPt(ilead);
      else
	LeadpT = PhotonSmearedPt(ilead);
      Leadetas2     = m_Rd->GetVariable( Form("ph_etas2[%d]",ilead) );      
      Leadphi       = m_Rd->GetVariable( Form("ph_phi[%d]",ilead) );      
      Leadrphi      = m_Rd->GetVariable( Form("ph_rphi[%d]",ilead) );      
      Leadreta      = m_Rd->GetVariable( Form("ph_reta[%d]",ilead) );      
      Leadtime      = m_Rd->GetVariable( Form("ph_time[%d]",ilead) );      
      if(m_data)
	SubLeadpT = PhotonRescaledPt(isublead);
      else
	SubLeadpT = PhotonSmearedPt(isublead);
      SubLeadetas2  = m_Rd->GetVariable( Form("ph_etas2[%d]",isublead) );      
      SubLeadphi    = m_Rd->GetVariable( Form("ph_phi[%d]",isublead) );      
      SubLeadrphi   = m_Rd->GetVariable( Form("ph_rphi[%d]",isublead) );      
      SubLeadreta   = m_Rd->GetVariable( Form("ph_reta[%d]",isublead) );      
      SubLeadtime   = m_Rd->GetVariable( Form("ph_time[%d]",isublead) );      
      myRunNumber   = (int)m_Rd->GetVariable("RunNumber");
      myEventNumber = (int)m_Rd->GetVariable("EventNumber");
      myLumiBlock   = (int)m_Rd->GetVariable("lbn");
      myLArError    = (int)m_Rd->GetVariable("larError");
      tree->Fill();
    }
    
  //______________________________________________________________//
    
  }//End of the loop over all entries
  //========================================================================//
  //========== End of the loop over all entries ============================//
  //========================================================================//


  if(output !="DEFAULT") m_output=output;
  TFile *f1=new TFile(m_output,"RECREATE");  
  tree->Write();
  f1->Close();

  delete f1;
  delete tree;


}
//////////////////////////////////////////////////////////
void CleaningStudies::PhotonStudies(bool data, TString output)
//////////////////////////////////////////////////////////
{
  

  //=====Declare output variables====//
  TTree *tree=new TTree("photons","My photons tree");
  int myRunNumber=0;
  int myEventNumber=0;
  int myLumiBlock;
  int myLArError;
  double mymet=0;
  unsigned int phOQ=0;
  double phclpt=0;
  double phetas2=0;
  double phphi=0;
  double phreta=0;
  double phrphi=0;
  int phPassJC=0;
  int phPassLC=0;
  int phPassLC2=0;
  double phJC=0;
  int phistight=0;
  tree->Branch("RunNumber",&myRunNumber,"RunNumber/I");
  tree->Branch("EventNumber",&myEventNumber,"EventNumber/I");
  tree->Branch("LumiBlock",&myLumiBlock,"LumiBlock/I");
  tree->Branch("larError",&myLArError,"larError/I");
  tree->Branch("MET",&mymet,"MET/D");
  tree->Branch("phclpt",&phclpt,"phclpt/D");
  tree->Branch("phetas2",&phetas2,"phetas2/D");
  tree->Branch("phphi",&phphi,"phphi/D");
  tree->Branch("phreta",&phreta,"phreta/D");
  tree->Branch("phrphi",&phrphi,"phrphi/D");
  tree->Branch("phPassJC",&phPassJC,"phPassJC/I");
  tree->Branch("phPassLC",&phPassLC,"phPassLC/I");
  tree->Branch("phPassLC2",&phPassLC2,"phPassLC2/I");
  tree->Branch("phOQ",&phOQ,"phOQ/i");
  tree->Branch("phJC",&phJC,"phJC/D");
  tree->Branch("phistight",&phistight,"phistight/I");
  //________________________________//
  
  //==Declare the GRL object==//
  DQ::SetXMLFile(m_GRL);
  //__________________________//

  //======================================================================//  
  //====== Start the Loop over all entries ===============================//
  //======================================================================//  
  for(int jentry=0;jentry<m_nentries;jentry++){
    m_Rd->GetEntry(jentry);
    
    // ===clean all flags=======//     
    m_passTrigger       = false ;
    m_inGRL             = false ;
    m_passPV            = false ;
    m_passPreSel        = false; 
    //__________________________//
 
   //============ GRL Selection ==================//
    int Run = (int)m_Rd->GetVariable("RunNumber");
    int lbn = (int)m_Rd->GetVariable("lbn");
    if(data){
      m_inGRL = DQ::PassRunLB(Run ,lbn);
    }else{
      m_inGRL=true;
    }
    //___________________________________________//
    
 
    //=================== Trigger Requirement ====================//
    //For Photon D3PDs
    m_passTrigger    = (bool)m_Rd->GetVariable("EF_2g20_loose");
    //__________________________________________________________//

    //==============Primary Vertex Requirement==========//
    unsigned int n_PV=(unsigned int)m_Rd->GetVariable("@PV_nTracks.size()");
    for(unsigned int i=0;i<n_PV;i++){
      if( (int)m_Rd->GetVariable(Form("PV_nTracks[%d]",i))>2){
	  m_passPV=true;
	  break;
      }
    }
    //_________________________________________________//
    

    //=============================== Photons selection ==========================//
    if( !m_passTrigger) continue;
    if( !m_inGRL)       continue; 
    if( !m_passPV )     continue; 

 
    unsigned int phsize=(unsigned int)m_Rd->GetVariable("@ph_cl_pt.size()");
    for ( unsigned int iph=0 ; iph < phsize ; iph++ ){//Loop over reco photons
      if ( PhotonPtOK(iph) && PhotonIsLooseOK(iph) && //Select pt>25GeV,loose
	   PhotonOQOK(iph) && PhotonEtaOK(iph) ){//Select OQ and cracks
	if( (int)m_Rd->GetVariable(Form("ph_author[%d]",iph) )==128) continue;
	myRunNumber   = Run;
	myEventNumber = (int)m_Rd->GetVariable("EventNumber");
	myLumiBlock   = (int)m_Rd->GetVariable("lbn");
	myLArError    = (int)m_Rd->GetVariable("larError");
	mymet=m_Rd->GetVariable("MET_RefFinal_et");
	double phclE=m_Rd->GetVariable( Form("ph_cl_E[%d]",iph) );
	phetas2=m_Rd->GetVariable( Form("ph_etas2[%d]",iph) );
	phclpt=phclE/TMath::CosH(phetas2);
	phphi=m_Rd->GetVariable( Form("ph_phi[%d]",iph) );
	phreta=m_Rd->GetVariable( Form("ph_reta[%d]",iph) );
	phrphi=m_Rd->GetVariable( Form("ph_rphi[%d]",iph) );
	phOQ=(unsigned int) m_Rd->GetVariable( Form("ph_OQ[%d]",iph) );
	phPassJC=1;
	if( !PassNoiseCleaning(iph) ) phPassJC=0;
	phPassLC=1; 
	if( !PassLArCleaning(iph) ) phPassLC=0;
	phPassLC2=1;
	if( !PassPhotonCleaning(iph) ) phPassLC2=0;
	if( HasJetAssociated(iph) ){
	  int indexJ=(int)m_Rd->GetVariable( Form("ph_jet_AntiKt4TopoEMJets_index[%d]",iph) );
	  phJC=m_Rd->GetVariable( Form("jet_AntiKt4TopoEMJets_LArQuality[%d]",indexJ) );
	}
	phistight=0;
	if(PhotonIsTightOK(iph) ) phistight=1;
	tree->Fill();	
      }
    } 
    
  //_____________________________________________________________________________________//
    AnalysisTools myAT;
    myAT.ShowNumberOfProcessedEvents(jentry,m_nentries);
 }//End of the loop over all entries
  //========================================================================//
  //========== End of the loop over all entries ============================//
  //========================================================================//


  if(output !="DEFAULT") m_output=output;
  TFile *f1=new TFile(m_output,"RECREATE");  
  tree->Write();
  f1->Close();

  delete f1;
  delete tree;

}



///////////////////////////////////////////
bool CleaningStudies::HasJetAssociated(int iph)
///////////////////////////////////////////
{  
  float jet_emfrac=0;
  int is_jet_matched=(int)m_Rd->GetVariable( Form("ph_jet_AntiKt4TopoEMJets_matched[%d]",iph) );
  if( is_jet_matched>0 ){
    int index_jet=(int)m_Rd->GetVariable( Form("ph_jet_AntiKt4TopoEMJets_index[%d]",iph) );
    jet_emfrac= (float)m_Rd->GetVariable( Form("jet_AntiKt4TopoEMJets_emfrac[%d]",index_jet) );
  }
  bool okjetq = true;
  if( is_jet_matched <0 || (is_jet_matched>=0 &&jet_emfrac<0.95) ){
    okjetq = false;
  }
  return okjetq;
}
