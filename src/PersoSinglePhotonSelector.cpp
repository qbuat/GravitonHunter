#include "PersoSinglePhotonSelector.h"
#include "ToolsUtilities.h"
#include "ToolsCommons.h"
#include "GoodRunsLists/DQHelperFunctions.h"
//////////////////////////////////////////////////////////////
SinglePhotonSelector::SinglePhotonSelector(TTree* tree) : GravitonAnalysis(tree)
//////////////////////////////////////////////////////////////
{
   if (tree == 0) {
     std::cout << "You need to specify an input tree !!!"<<std::endl;
   }
   Init(tree);
}
///////////////////////////////////////////////
SinglePhotonSelector::SinglePhotonSelector() : GravitonAnalysis()
///////////////////////////////////////////////
{
  //Default constructor
}
///////////////////////////
SinglePhotonSelector::~SinglePhotonSelector()
///////////////////////////
{
  //Default destructor
}
///////////////////////////////////
void SinglePhotonSelector::Init(TTree* tree)
//////////////////////////////////
{
  m_passTrigger  = false;
  m_inGRL        = false;
  m_passPV       = false;
  m_passPreSel   = false;

  m_output       = "toto.root";
  m_PU_mcfile    = Commons::PileupMCFile;
  m_PU_datafile  = Commons::PileupDataFile;

  looseprime2    = 2620417;
  looseprime3    = 2489345;
  looseprime4    = 392193;
  looseprime5    = 130049;
  SetPileUpTool(m_PU_mcfile,m_PU_datafile);

}


////////////////////////////////////////////////////////////////////////////////////////////////////
void SinglePhotonSelector::EventLoop(bool data,TString output,int Nreversedcut,TString mctype)
///////////////////////////////////////////////////////////////////////////////////////////////////
{

  unsigned int looseprimemask=0;
  if(Nreversedcut==0)
    looseprimemask=0;
  else if(Nreversedcut==2)
    looseprimemask=looseprime2;
  else if(Nreversedcut==3)
    looseprimemask=looseprime3;
  else if(Nreversedcut==4)
    looseprimemask=looseprime4;
  else if(Nreversedcut==5)
    looseprimemask=looseprime5;

  std::cout << "NENTRIES " << m_nentries << std::endl;
  m_data  = data;//Choose between data and mc
  m_mctype= mctype;// Choose mc type "mc11a"/"mc11b"

  //--------------------------
  double Iso_g       = -9999.;
  double weight_g    = -9999.;
  int IsConv_g       = -9999;
  int NPV            = -9999;
  double pT_g        = -9999.;
  double eta_g        = -9999.;
  double phi_g        = -9999.;
  //------------------------------------------
  TString treename              = "tree";
  TString Branchname_g          = "Iso_g";
  TString Branchname_IsConv_g   = "IsConv_g";
  TString Branchname_pT_g       = "pT_g";
  TString Branchname_eta_g      = "eta_g";
  TString Branchname_phi_g      = "phi_g";
  TString Branchname_NPV        = "NPV";
  TString Branchname_W_g        = "weight_g";
  TString Branchname_Run        = "RunNumber";
  TString Branchname_Event      = "EventNumber";
  TString Branchname_LB         = "LumiBlock";
  //----------------------------------------------------------------------
  TTree *tree = new TTree(treename,"Isolation tree");
  tree->Branch(Branchname_g,&Iso_g,Branchname_g+"/D");
  tree->Branch(Branchname_IsConv_g,&IsConv_g,Branchname_IsConv_g+"/I");
  tree->Branch(Branchname_pT_g,&pT_g,Branchname_pT_g+"/D");
  tree->Branch(Branchname_NPV,&NPV,Branchname_NPV+"/I");
  tree->Branch(Branchname_W_g,&weight_g,Branchname_W_g+"/D");
  tree->Branch(Branchname_Run,&m_RunNumber,Branchname_Run+"/I");
  tree->Branch(Branchname_Event,&m_EventNumber,Branchname_Event+"/I");
  tree->Branch(Branchname_LB,&m_LumiBlock,Branchname_LB+"/I");
  //---------------------------------------------------------------------  

  //======weight ======//
  double weight   = 1;
  double PUweight = 1;
  //==================//

  //== Declare the GRL object ==//
  DQ::SetXMLFile(m_GRL);
  //____________________________//
  
  
  //======================================================================//  
  //================ Start the Loop over all entries =====================//
  //======================================================================//  
  for(int jentry=0;jentry<m_nentries;jentry++){
    m_Rd->GetEntry(jentry);
    
    //======== PU Weight && enrgy smearing for MC ONLY ==================//
    if(!m_data){
      PUweight = m_pileupTool->GetCombinedWeight((int)m_Rd->GetVariable("RunNumber"),
						 (int)m_Rd->GetVariable("mc_channel_number"),
						 m_Rd->GetVariable("averageIntPerXing") );
    }
    //__________________________________________________________________//
 
   // ==== clean all flags =====//     
    m_passTrigger       = false ;
    m_inGRL             = false ;
    m_passPV            = false ;
    m_passPreSel        = false ;
    //__________________________//
 
   //================= GRL Selection ====================//
    m_RunNumber   = (int)m_Rd->GetVariable("RunNumber");
    m_EventNumber = (int)m_Rd->GetVariable("EventNumber");
    m_LumiBlock   = (int)m_Rd->GetVariable("lbn");
    if(m_data)
      m_inGRL = DQ::PassRunLB(m_RunNumber ,m_LumiBlock);
    else
      m_inGRL = true;
    //__________________________________________________//
    
    //================== Trigger Requirement ===================//
    //For Photon D3PDs
    m_passTrigger    = ( (bool)m_Rd->GetVariable("EF_g20_loose") ||
			 (bool)m_Rd->GetVariable("EF_g40_loose") ||
			 (bool)m_Rd->GetVariable("EF_g60_loose") ||
			 (bool)m_Rd->GetVariable("EF_g80_loose")  
			 );
    //__________________________________________________________//
 
    
    //======================== Primary Vertex Requirement ======================//
    unsigned int n_PV=(unsigned int)m_Rd->GetVariable("@PV_nTracks.size()");
    NPV=(int)n_PV;
    if( !m_data )
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

    pT_g=-99999; 
    int i_g=-99999; 
    weight=PUweight;
    for( int iph=0;iph<(int)m_Rd->GetVariable("@ph_cl_pt.size()");iph++){
      if( !PhotonEtaOK(iph) ) continue;
      if( !PhotonPtOK(iph) ) continue;
      if( !PhotonOQOK(iph) ) continue;
      if( !PassPhotonCleaning(iph) ) continue;
      bool isem; 
      if ( Nreversedcut==0 )
	isem = PhotonIsTightOK(iph);
      else{
	isem = ( ((unsigned int)m_Rd->GetVariable(Form("ph_isEM[%d]",iph))&looseprimemask)!=0 &&
		 !PhotonIsTightOK(iph)
		 );
      }
      if( !m_data)
	isem = PhotonIsEM_MCOK(iph,true);
      if( !isem ) continue;
      PhotonRescaledPt(iph);
      if (PhotonRescaledPt(iph)> pT_g){
	pT_g = PhotonRescaledPt(iph);
	i_g = iph;
      }
    }
    if( i_g==-99999) continue;

    if(m_isotype == "CONE")
      Iso_g  = PhotonIsolation_tool(i_g)/1000.;
    else if(m_isotype == "TOPO")
      Iso_g  = PhotonTopoIsolation_tool(i_g)/1000.;
    else Fatal(" SinglePhotonSelector::EventLoop","Wrong isolation type !!!");
    if( (int)m_Rd->GetVariable("larError")>1 ) continue;


    IsConv_g  = (int)m_Rd->GetVariable(Form("ph_isConv[%d]",i_g));
    pT_g  = PhotonRescaledPt(i_g);
    eta_g = m_Rd->GetVariable(Form("ph_etas2[%d]",i_g));
    phi_g = m_Rd->GetVariable(Form("ph_phi[%d]",i_g));
    weight_g  = PUweight;
    tree->Fill();

    AnalysisTools myAT;
    myAT.ShowNumberOfProcessedEvents(jentry,m_nentries);
    
  }//End of the loop over all entries
  std::cout << "End of the loop over all entries " << std::endl;
  //========================================================================//
  //========== End of the loop over all entries ============================//
  //========================================================================//

  if(output !="DEFAULT") m_output=output;
  TFile f1(output,"RECREATE");  
  tree->Write();
  f1.Close();
  delete tree;
}

