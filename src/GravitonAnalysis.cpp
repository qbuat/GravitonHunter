#include "GravitonAnalysis.h"
#include "ToolsCommons.h"
#include "GoodRunsLists/DQHelperFunctions.h"
#include "egammaAnalysisUtils/CaloIsoCorrection.h"
#include "egammaAnalysisUtils/PhotonIDTool.h"
#include "egammaAnalysisUtils/FudgeMCTool.h"

///////////////////////////////////////////////
GravitonAnalysis::GravitonAnalysis(TTree* tree)
///////////////////////////////////////////////
{
   if (tree == 0)
     std::cout << "You need to specify an input tree !!!"<<std::endl;
   Init(tree);
}
///////////////////////////////////
GravitonAnalysis::GravitonAnalysis()
///////////////////////////////////
{
  //Default constructor
  m_nentries      = 0;
  m_mctype        = "";
  m_isotype       = "";
  m_GRL           = "";
  m_data          = true;
  m_Rd            = 0;
  m_pileupTool    = 0;
}
/////////////////////////////////////
GravitonAnalysis::~GravitonAnalysis()
////////////////////////////////////
{
  //Default destructor
  std::cout << "clean memory : GravitonAnalysis" << std::endl;
  delete m_Rd;
  std::cout<<"m_Rd deleted"<<std::endl;
  delete m_pileupTool;
  std::cout<<"m_pileupTool deleted"<<std::endl;
  delete m_Ts;
  std::cout<<"m_Ts deleted"<<std::endl;
}
////////////////////////////////////////////
void GravitonAnalysis::InitTree(TTree* tree)
////////////////////////////////////////////
{ m_Rd = new TreeReader(tree); }

/////////////////////////////////////////
void GravitonAnalysis::Init(TTree* tree)
/////////////////////////////////////////
{
  InitTree(tree);
  m_nentries      = (int)tree->GetEntriesFast();

  //=========== Energy Rescaler ===========//
  m_Rescale.Init("./rootfiles/EnergyRescalerData.root","2012","es2010");
  m_SmearFac  = 1;

  m_Rescale_geo21.setESModel(egEnergyCorr::es2012c);
  m_Rescale_geo21.initialize();


  CaloIsoCorrection::SetPtLeakageCorrectionsFile("./rootfiles/isolation_leakage_corrections.root");
  //=========== PileUp Reweighting ========//
  m_pileupTool = 0;
  //=========== Truth Selector ========//
  m_Ts = 0;

  //=========== Basic settings ==========//

  // SetMCType("");
  // SetIsoType("CONE");

  // SetGRL(Commons::GrlFile);
  // SetGeeFile(Commons::GeeFile);
  // SetPileUpTool(Commons::PileupMCFile,Commons::PileupDataFile);
  // SetTruthSelector(tree);
}

//////////////////////////////////////////
void GravitonAnalysis::SetGRL(TString GRL)
//////////////////////////////////////////
{
  m_GRL = GRL;
  DQ::SetXMLFile(m_GRL);
  std::cout << "A GRL has been set: " 
	    << GRL << std::endl;
}
//////////////////////////////////////////////
void GravitonAnalysis::SetEntries(int entries)
//////////////////////////////////////////////
{
  //--> Use this method to run on a subset of the sample
  //--> To be called after Init(tree);
  m_nentries = entries;
}
//////////////////////////////////////////////////
void GravitonAnalysis::SetMCType(TString mctype)
//////////////////////////////////////////////////
{ 
  m_mctype = mctype; 
  std::cout << "A MonteCarlo type has been set: " 
	    << mctype << std::endl;
}
//////////////////////////////////////////////////
void GravitonAnalysis::SetIsoType(TString st)
//////////////////////////////////////////////////
{ 
  m_isotype = st; 
  std::cout << "An isolation type has been set: "
	    << st << std::endl;
}
///////////////////////////////////////////////
void GravitonAnalysis::SetStreamType(bool data)
///////////////////////////////////////////////
{
  //--> data: true, mc:false
  m_data = data;
  if(m_data) std::cout << "The chosen stream is data" << std::endl;
  else std::cout << "The chosen stream is mc" << std::endl;
}
/////////////////////////////////////////////////////////////////////
void GravitonAnalysis::SetPileUpTool(TString mcfile,TString datafile)
//////////////////////////////////////////////////////////////////////
{
  if(m_data) m_pileupTool = 0;
  else {
    m_pileupTool = new Root::TPileupReweighting("PileUpTool");
    m_pileupTool->AddConfigFile(mcfile);
    m_pileupTool->AddLumiCalcFile(datafile);
    m_pileupTool->SetUnrepresentedDataAction(2);
    m_pileupTool->initialize();
  }
  std::cout << "The Pileup Reweighting tool has been initialized with: "
	    << mcfile << ", and " << datafile 
	    << std::endl;
}
//////////////////////////////////////////////////
void GravitonAnalysis::SetGeeFile(TString Geefile)
//////////////////////////////////////////////////
{
  m_GeeFile   = Geefile ;
  m_ee_events = build_eeset(m_GeeFile);
  std::cout << "A Z'->ee event list has been set: " << Geefile << std::endl;
}
////////////////////////////////////////////////////
void GravitonAnalysis::SetTruthSelector(TTree* tree)
///////////////////////////////////////////////////
{
  if(m_data) {
    m_Ts = 0;
    std::cout << "The TruthSelector IS NOT initialized" << std::endl;
  }
  else { 
    m_Ts = new TruthSelector(tree);
  std::cout << "The TruthSelector has been initialized" << std::endl;
  }
}
///////////////////////////////////
bool GravitonAnalysis::EventGRLOK()
///////////////////////////////////
{
  bool inGRL = false;
  int run = (int)m_Rd->GetVariable("RunNumber");
  int lbn = (int)m_Rd->GetVariable("lbn");
  if(m_data)
    inGRL = DQ::PassRunLB(run,lbn);
  else
   inGRL = true;
  return inGRL;
}
/////////////////////////////////////////
bool GravitonAnalysis::EventCompletedOK()
////////////////////////////////////////
{
  // Incomplete events: In 2012 data-taking the TTC restart was developed to recover 
  // certain detector busy conditions without a run-restart. In the lumi-block after a 
  // TTC restart there can be events with incomplete events (where some detector information 
  // is missing from the event). 
  // People should check for and remove such events from their analysis.
  bool isEventCompleted = true;
  if( ((unsigned int)m_Rd->GetVariable("coreFlags")&0x40000)!=0 ) 
    isEventCompleted = false;
  else
    isEventCompleted = true;
  return isEventCompleted;
}
///////////////////////////////////
bool GravitonAnalysis::EventPVOK()
///////////////////////////////////
{
  //--> Keep an event if at least one vertex
  //--> with 3 or more associated tracks
  bool passPV = false;
  int NPV = (int)m_Rd->GetVariable("@PV_nTracks.size()");
  for(int ipv=0;ipv<NPV;ipv++){
    if( (int)m_Rd->GetVariable(Form("PV_nTracks[%d]",ipv) )>2 ){
      passPV=true;
      break;
    }
  }
  return passPV;
}
////////////////////////////////////
bool GravitonAnalysis::EventTrigOK()
///////////////////////////////////
{
  //--> Trigger requirement
  bool passTrigger = (bool)m_Rd->GetVariable("EF_g35_loose_g25_loose");
  return passTrigger;
}
/////////////////////////////////////////////////////////////
bool GravitonAnalysis::PhotonPtOK(int iph,double CutValueGeV)
////////////////////////////////////////////////////////////
{
  //--> pt requirement (default value is 25 GeV)
  //--> Apply on the rescaled pT for data
  //--> Apply on the smeared pT for MC
  //--> Correction applied with egammaAnalysisUtils
  bool temp=false;
  float pt=0;
  if(m_data)
    pt = PhotonRescaledPt_geo21(iph);
  else
    pt = PhotonSmearedPt(iph);
 
  if( pt > CutValueGeV )
    temp = true;

  return temp;
}
///////////////////////////////////////////
bool GravitonAnalysis::PhotonEtaOK(int iph)
///////////////////////////////////////////
{
  //--> pseudorapidity requirement
  //--> Use the LAr cal sampling 2 eta value
  float eta=(float)m_Rd->GetVariable( Form("ph_etas2[%d]",iph) );
  bool temp=false;
  if( fabs(eta) < 1.37 || 
      (fabs(eta) > 1.56 && fabs(eta) < 2.37) ){
    temp = true;
  }
  return temp;
}
/////////////////////////////////////////
bool GravitonAnalysis::PhotonOQOK(int iph)
/////////////////////////////////////////
{
  //--> Object Quality requirement cut
  //--> See https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/LArCleaningAndObjectQuality
  bool temp=false;
  unsigned int OQbits=(unsigned int)m_Rd->GetVariable( Form("ph_OQ[%d]",iph) );
  unsigned int mask=34214;
  if( (OQbits&mask)==0)
    temp=true;
  
  return temp;
}
/////////////////////////////////////////////////////////////////
bool GravitonAnalysis::PhotonAmbiguityResolverOK(const int & iph)
/////////////////////////////////////////////////////////////////
{
  unsigned int isEM_word = (unsigned int)m_Rd->GetVariable(Form("ph_isEM[%d]", iph));
  if ((isEM_word & 0x800000) != 0)
    return false;
  else
    return true;
}

///////////////////////////////////////////////
bool GravitonAnalysis::PassLArCleaning(int iph)
//////////////////////////////////////////////
{
  //--> LArCleaning>0.8 requirement cut
  //--> See https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/LArCleaningAndObjectQuality
  bool temp=false;
  unsigned int OQbits=(unsigned int)m_Rd->GetVariable( Form("ph_OQ[%d]",iph) );
  unsigned int mask=134217728;
  if( (OQbits&mask)==0){
    temp=true;
  }
  return temp;
}
//////////////////////////////////////////////////
bool GravitonAnalysis::PassPhotonCleaning(int iph)
//////////////////////////////////////////////////
{
  //--> Complete PhotonCleaning requirement cut
  //--> See https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/LArCleaningAndObjectQuality
  bool temp=true;
  unsigned int OQbits=(unsigned int)m_Rd->GetVariable( Form("ph_OQ[%d]",iph) );
  unsigned int maskLC=134217728;
  unsigned int maskTIME=67108864;
  double reta=m_Rd->GetVariable( Form("ph_reta[%d]",iph) );
  double rphi=m_Rd->GetVariable( Form("ph_rphi[%d]",iph) );
  if( (OQbits&maskLC)!=0 &&
      (reta>0.98 || rphi>1.0 || (OQbits&maskTIME)!=0)//rphi,reta,time cut
      ){
    temp=false;
  }
  return temp;
}
/////////////////////////////////////////////////
bool GravitonAnalysis::PassNoiseCleaning(int iph)
////////////////////////////////////////////////
{  
  //gamma-'jet' cleaning
  // Old cleaning based on the jet cleaning
  // applied on the jet matched to the photon cluster
  // Not used since 2011
  bool okjetq = true;
  int is_jet_matched=(int)m_Rd->GetVariable( Form("ph_jet_AntiKt4TopoEMJets_matched[%d]",iph) );
  float jet_emFrac=0;
  float jet_qual=0;

  if(is_jet_matched){
    int index_jet=(int)m_Rd->GetVariable( Form("ph_jet_AntiKt4TopoEMJets_index[%d]",iph) );
    jet_emFrac= (float)m_Rd->GetVariable( Form("jet_AntiKt4TopoEMJets_emfrac[%d]",index_jet) );
    jet_qual= (float)m_Rd->GetVariable( Form("jet_AntiKt4TopoEMJets_LArQuality[%d]",index_jet) );
    if( jet_emFrac>0.95 && fabs(jet_qual)>0.8 ){
      okjetq=false;
    }
  }
  return okjetq;
}
///////////////////////////////////////////////
bool GravitonAnalysis::PhotonIsTightOK(int iph)
///////////////////////////////////////////////
{
  //--> Tight photon decision
  //--> The ID menu is recomputed from the D3PD
  //--> with the PhotonIDTool. 
  //--> Different for data and mc because of the use
  //--> of the pT of the photon (rescaled/smeared)
  bool temp=false;
  if(m_data)
    temp = PhotonIsEM_OK(iph, true);
    // temp = PhotonIsTight_D3PDOK(iph);
  else
    temp = PhotonIsEM_MCOK(iph,true);
  return temp;
}
///////////////////////////////////////////////
bool GravitonAnalysis::PhotonIsLooseOK(int iph)
///////////////////////////////////////////////
{
  //--> Loose photon decision
  //--> The ID menu is recomputed from the D3PD
  //--> with the PhotonIDTool. 
  //--> Different for data and mc because of the use
  //--> of the pT of the photon (rescaled/smeared)
  bool temp=false;
  if(m_data)
    temp = PhotonIsEM_OK(iph, false);
    // temp = PhotonIsLoose_D3PDOK(iph);
  else
    temp = PhotonIsEM_MCOK(iph,false);
  return temp;
}
////////////////////////////////////////////////////
bool GravitonAnalysis::PhotonIsTight_D3PDOK(int iph)
///////////////////////////////////////////////////
{
  //--> Tight ID decision stored in the D3PD (computed on the AOD)
  bool temp=false;
  int ph_isPhotonTight=(int)m_Rd->GetVariable( Form("ph_tight[%d]",iph) );
  if(ph_isPhotonTight==1)
    temp=true;
  return temp;
}
////////////////////////////////////////////////////
bool GravitonAnalysis::PhotonIsLoose_D3PDOK(int iph)
////////////////////////////////////////////////////
{
  //--> Loose ID decision stored in the D3PD (computed on the AOD)
  bool temp=false;
  int ph_isPhotonLoose=(int)m_Rd->GetVariable( Form("ph_loose[%d]",iph) );
  if(ph_isPhotonLoose==1)
    temp=true;
  return temp;
}
/////////////////////////////////////////////////////////////////
bool GravitonAnalysis::PhotonIsolationOK(int iph,double CutValueGeV)
////////////////////////////////////////////////////////////////
{
  //--> Isolation requirement from the
  //--> Cone40 isolation stored in the D3PD (computed on the AOD)
  //--> corrected for Energy Density and pT leakage
  bool temp=false;
  double Isol=m_Rd->GetVariable( Form("ph_Etcone40_corrected[%d]",iph) )/1000.;
  if( Isol < CutValueGeV )
    temp=true;
  return temp;
}
/////////////////////////////////////////////////////////////////////////
bool GravitonAnalysis::PhotonIsolation_toolOK(int iph,double CutValueGeV)
////////////////////////////////////////////////////////////////////////
{
  //--> Isolation requirement from the
  //--> isolation computed with CaloIsoCorrection tool
  //--> corrected for Energy Density and pT leakage
  //Default CutValue=5;
  bool temp=false;
  float Isol;
  if( m_isotype == "CONE" )
    Isol = PhotonIsolation_tool(iph)/1000.;
  else if( m_isotype == "TOPO")
    Isol = PhotonTopoIsolation_tool_geo21(iph)/1000.;
  else if( m_isotype == "TOPO_PTDEP")
    Isol = PhotonTopoIsolationPtdep_tool(iph)/1000.;

  else Fatal("GravitonAnalysis::PhotonIsolation_toolOK", "Wrong isolation type !!!");
  if(Isol < CutValueGeV )
    temp=true;
  return temp;
}
/////////////////////////////////////////////////////////////////
bool GravitonAnalysis::ElectronPtOK(int iele, double CutValueGeV)
////////////////////////////////////////////////////////////////
{
  //--> Electron pT requirement 
  //--> Applied on rescaled/smeared pT for data/mc
  //--> pT computed with the cluster energy and cluster eta
  bool temp=false;
  float E   =(float)m_Rd->GetVariable( Form("el_cl_E[%d]",iele) )/1000.;
  float eta =(float)m_Rd->GetVariable( Form( "el_cl_eta[%d]",iele) );
  float pt=0;
  if(m_data)
    pt = ElectronRescaledE(iele)/TMath::CosH(eta);
  else
    pt = E/TMath::CosH(eta);

  if( pt>CutValueGeV )
    temp = true;
  return temp;
}
//////////////////////////////////////////////
bool GravitonAnalysis::ElectronEtaOK(int iele)
//////////////////////////////////////////////
{
  //--> Electron pseudorapidity requirement
  //--> use the eta of the ECAL cluster
  float eta=(float)m_Rd->GetVariable( Form("el_cl_eta[%d]",iele) );
  bool temp=false;
  if( fabs(eta)<1.37 || (fabs(eta)>1.52 && fabs(eta)<2.37) ){
    temp=true;
  }
  return temp;
}
/////////////////////////////////////////////
bool GravitonAnalysis::ElectronOQOK(int iele)
////////////////////////////////////////////
{
  //--> Object Quality requirement cut
  //--> See https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/LArCleaningAndObjectQuality
  bool temp=false;
  unsigned int OQbits=(unsigned int)m_Rd->GetVariable( Form("el_OQ[%d]",iele) );
  unsigned int mask=1446;
  if( (OQbits&mask)==0){
    temp=true;
  }
  return temp;
}///////////////////////////////////////////////////
bool GravitonAnalysis::ElectronTrackIsoOK(int iele)
///////////////////////////////////////////////////
{
  //--> Electron track isolation requirement

  bool temp=false;
  float iso =  m_Rd->GetVariable( Form("el_ptcone20[%d]",iele) );
  if( iso<0.15 )
    temp=true;
  return temp;
}
///////////////////////////////////////////////////
bool GravitonAnalysis::ElectronIsMediumOK(int iele)
///////////////////////////////////////////////////
{
  //--> Electron Medium ID requirement 
  //--> Use the value stored in the D3PDs.
  bool temp=false;
  int isem=(int)m_Rd->GetVariable( Form("el_medium[%d]",iele) );
  if(isem==1)
    temp=true;
  return temp;
}
//////////////////////////////////////////////////////////////////////////////
bool GravitonAnalysis::EventPreSelectionOK(int *ilead,int *isublead,
					   double *pt_lead,double *pt_sublead)
//////////////////////////////////////////////////////////////////////////////
{
  //--> Diphoton Event preselection algorithm
  //--> Return the indices of the two selected photons and also the value of their
  //--> transverse impulsion in GeV.
  //--> The logic consists in selecting the two highest pT photons satisfying
  //--> Pt cut (30 GeV), Eta cut, Object Quality, Photon cleaning and loose ID.

  int nPhotPtLooseEtaOQCleaning=0;
  unsigned int phsize=(unsigned int)m_Rd->GetVariable("@ph_cl_pt.size()");
  for ( unsigned int iph=0 ; iph < phsize ; iph++ ){//Loop over reco photons

    if( PhotonPtOK(iph,30.) && PhotonEtaOK(iph) &&
	PhotonOQOK(iph) && PassPhotonCleaning(iph) &&
	PhotonIsLooseOK(iph) ) {
      ++nPhotPtLooseEtaOQCleaning;
      double pt_current_photon=0;
      if(m_data)
	pt_current_photon=PhotonRescaledPt_geo21(iph);
      else// Apply smearing for MC
	pt_current_photon=PhotonSmearedPt(iph);
    
      if ( pt_current_photon > (*pt_lead) ) { //found a new photon above the leading
	*pt_sublead=*pt_lead; //the leading so-far becomes the new subleading
	*pt_lead=pt_current_photon;
	*isublead=(*ilead);
	*ilead=iph;
      }else if ( pt_current_photon > (*pt_sublead) ) { //found a new photon above the subleading
	*pt_sublead=pt_current_photon;
	*isublead=iph;
      }
    }
  }// End of loop over reco photons
  if ((*ilead)<0 || (*isublead)<0)
    return 0;


//  float E1=(float)m_Rd->GetVariable( Form("ph_cl_E[%d]",*ilead) );
//  float etas21=(float)m_Rd->GetVariable( Form("ph_etas2[%d]",*ilead) );
//  float pt1=E1/TMath::CosH(etas21);
//  
//  float E2=(float)m_Rd->GetVariable( Form("ph_cl_E[%d]",*isublead) );
//  float etas22=(float)m_Rd->GetVariable( Form("ph_etas2[%d]",*isublead) );
//  float pt2=E2/TMath::CosH(etas22);
//  
//  if(pt1<pt2) {
//    std::cout<<"--------------------------------------------"<<std::endl;
//    std::cout <<"EventNumber "<< m_Rd->GetVariable("EventNumber") <<" RunNumber "<< m_Rd->GetVariable("RunNumber") << " i1 "<< *ilead <<" i2 "<<*isublead<<std::endl;
//    std::cout<<"reco photon : Lead_pT "<<pt1*PhotonSmearingFactor(*ilead,0)<<" SubLead_pT "<<pt2*PhotonSmearingFactor(*isublead,0)<< std::endl;
//    std::cout<<"ph_cl_pt : "<<pt1/1000. <<" "<< pt2/1000.<<std::endl;
//    std::cout<<"photon smearing factors : "<<PhotonSmearingFactor(*ilead,0)<<" "<<PhotonSmearingFactor(*isublead,0)<<std::endl;
//  }

  return 1;
}



//----------------------------------------
//////////////////////////////////////////////////////////////////////////////
std::vector<int> GravitonAnalysis::EventPreSelectionOK()
//////////////////////////////////////////////////////////////////////////////
{
  //--> Diphoton Event preselection algorithm
  //--> Return a vector of indices of the two selected photons and also the value of their
  //--> transverse impulsion in GeV.
  //--> The logic consists in selecting the two highest pT photons satisfying
  //--> Pt cut (30 GeV), Eta cut, Object Quality, Photon cleaning and loose ID.


  // vector of indices    
  std::multimap<float,int> pt_index;
  
  // sorted vector of indices
  std::vector<int> v_index_pt;
  
  int nPhotPtLooseEtaOQCleaning=0;
  unsigned int phsize=(unsigned int)m_Rd->GetVariable("@ph_cl_pt.size()");
  //  std::cout<<"phsize "<<phsize<<std::endl;
  if(phsize < 2) return v_index_pt;
  for ( unsigned int iph=0 ; iph < phsize ; iph++ ){//Loop over reco photons

    if( PhotonPtOK(iph,30.) && PhotonEtaOK(iph) &&
	PhotonOQOK(iph) && PassPhotonCleaning(iph) && PhotonAmbiguityResolverOK(iph) && 
	PhotonIsLooseOK(iph) ) {
      ++nPhotPtLooseEtaOQCleaning;
      double pt_current_photon=0;
      if(m_data)
	pt_current_photon=PhotonRescaledPt_geo21(iph);
      else// Apply smearing for MC
	pt_current_photon=PhotonSmearedPt(iph);
      
      pt_index.insert(std::pair<float,int>(pt_current_photon,iph));

//      if ( pt_current_photon > (*pt_lead) ) { //found a new photon above the leading
//	*pt_sublead=*pt_lead; //the leading so-far becomes the new subleading
//	*pt_lead=pt_current_photon;
//	*isublead=(*ilead);
//	*ilead=iph;
//      }else if ( pt_current_photon > (*pt_sublead) ) { //found a new photon above the subleading
//	*pt_sublead=pt_current_photon;
//	*isublead=iph;
//      }



    }
  }// End of loop over reco photons

  if(pt_index.size() < 2) return v_index_pt;
  for( std::multimap<float,int>::reverse_iterator ii=pt_index.rbegin(); ii!=pt_index.rend() && (*ii).first != 0; ++ii) {
    //    std::cout << (*ii).first << " ";
    //    std::cout << (*ii).second << " ";
    v_index_pt.push_back((*ii).second);
  }
  
  
  





//  float E1=(float)m_Rd->GetVariable( Form("ph_cl_E[%d]",v_index_pt[0]) );
//  float etas21=(float)m_Rd->GetVariable( Form("ph_etas2[%d]",v_index_pt[0]) );
//  float pt1=E1/TMath::CosH(etas21);
//  
//  float E2=(float)m_Rd->GetVariable( Form("ph_cl_E[%d]",v_index_pt[1]) );
//  float etas22=(float)m_Rd->GetVariable( Form("ph_etas2[%d]",v_index_pt[1]) );
//  float pt2=E2/TMath::CosH(etas22);
//  
//  if(pt1<pt2) {
//    std::cout<<"---- In presel step (GravitonAnalysis) ------ PT1 < PT2 -------------"<<std::endl;
//    std::cout <<"EventNumber "<< m_Rd->GetVariable("EventNumber") <<" RunNumber "<< m_Rd->GetVariable("RunNumber") << " i1 "<< v_index_pt[0] <<" i2 "<<v_index_pt[1]<<std::endl;
////    std::cout<<"reco photon : Lead_pT "<<pt1*PhotonSmearingFactor(v_index_pt[0],0)<<" SubLead_pT "<<pt2*PhotonSmearingFactor(v_index_pt[1],0)<< std::endl;
//    std::cout<<"ph_cl_pt : "<<pt1/1000. <<" "<< pt2/1000.<<std::endl;
////    std::cout<<"photon smearing factors : "<<PhotonSmearingFactor(v_index_pt[0],0)<<" "<<PhotonSmearingFactor(v_index_pt[1],0)<<std::endl;
//  }

  return v_index_pt;
}

//----------------------------------------




/////////////////////////////////////////////////////////
float GravitonAnalysis::PhotonTopoIsolation_tool(int iph)
/////////////////////////////////////////////////////////
{
  //--> Return the topological cone isolation value in GeV
  //--> The cone size is DeltaR = 0.4 and the isolation is 
  //--> corrected for Ambiant energy and cluster leakage

  float ED_median = (float)m_Rd->GetVariable(Form("ph_ED_median[%d]",iph));
  float energy;
  if(m_data)
    energy = 1000*PhotonRescaledE(iph);
  else
    energy = 1000*PhotonSmearedE(iph);

  float etaS2        = (float)m_Rd->GetVariable( Form("ph_etas2[%d]",iph) );
  float etaPointing  = (float)m_Rd->GetVariable( Form("ph_etap[%d]",iph) );
  float etaCluster   = (float)m_Rd->GetVariable( Form("ph_cl_eta[%d]",iph) );
  float radius       = 40;
  bool  is_mc        = (!m_data);
  float topoEtcone40 = (float)m_Rd->GetVariable( Form("ph_topoEtcone40[%d]",iph) );
  int   isConv       = (int)m_Rd->GetVariable( Form("ph_isConv[%d]",iph) );
  float iso = CaloIsoCorrection::GetPtEDCorrectedTopoIsolation( ED_median,
								energy,
								etaS2,
								etaPointing,
								etaCluster,
								radius,
								is_mc,
								topoEtcone40,
								isConv,
								CaloIsoCorrection::PHOTON );
  return iso;
}

////////////////////////////////////////////////////////////
float GravitonAnalysis::PhotonTopoIsolation_tool_geo21(int iph)
///////////////////////////////////////////////////////////
{
  //--> Return the cone isolation value in GeV
  //--> The cone size is DeltaR = 0.4 and the isolation is 
  //--> corrected for Ambiant energy and cluster leakage


  float radius = 40.;
  bool is_mc = (!m_data);

  float Isol = CaloIsoCorrection::GetPtEDCorrectedTopoIsolation(
								m_Rd->GetVariable(Form("ph_ED_median[%d]", iph)),
								m_Rd->GetVariable(Form("ph_cl_E[%d]", iph)),
								m_Rd->GetVariable(Form("ph_etas2[%d]", iph)),
								m_Rd->GetVariable(Form("ph_etap[%d]", iph)),
								m_Rd->GetVariable(Form("ph_cl_eta[%d]", iph)),
								radius,
								is_mc,
								m_Rd->GetVariable(Form("ph_topoEtcone40[%d]", iph)),
								(int)m_Rd->GetVariable(Form("ph_isConv[%d]", iph)),
								CaloIsoCorrection::PHOTON,
								CaloIsoCorrection::REL17);
  return Isol;
}

/////////////////////////////////////////////////////
float GravitonAnalysis::PhotonIsolation_tool(int iph)
/////////////////////////////////////////////////////
{
  //--> Return the cone isolation value in GeV
  //--> The cone size is DeltaR = 0.4 and the isolation is 
  //--> corrected for Ambiant energy and cluster leakage
  float Etcone40              = (float)m_Rd->GetVariable( Form("ph_Etcone40[%d]",iph) );
  float Etcone40_ED_corrected = (float)m_Rd->GetVariable( Form("ph_Etcone40_ED_corrected[%d]",iph) );
  float energy;
  if(m_data)
    energy = 1000*PhotonRescaledE(iph);
  else
    energy = 1000*PhotonSmearedE(iph);

  float etaS2                 = (float)m_Rd->GetVariable( Form("ph_etas2[%d]",iph) );
  float etaPointing           = (float)m_Rd->GetVariable( Form("ph_etap[%d]",iph) );
  float etaCluster            = (float)m_Rd->GetVariable( Form("ph_cl_eta[%d]",iph) );
  float radius                = 40;
  bool is_mc                  = (!m_data);
  float Etcone_value          = Etcone40;
  bool isConversion           = false;
  if( (int)m_Rd->GetVariable( Form("ph_isConv[%d]",iph) )==1){
    isConversion=true;
  }
  float Isol = CaloIsoCorrection::GetPtEDCorrectedIsolation( Etcone40,
							     Etcone40_ED_corrected,
							     energy,
							     etaS2,
							     etaPointing,
							     etaCluster,
							     radius,
							     is_mc,
							     Etcone_value,
							     isConversion,
							     CaloIsoCorrection::PHOTON );
  return Isol;
}
/////////////////////////////////////////////////////
float GravitonAnalysis::PhotonTopoIsolationPtdep_tool(int iph)
/////////////////////////////////////////////////////
{
  //--> Return the topological cone isolation value in GeV
  //--> The cone size is DeltaR = 0.4 and the isolation is 
  //--> corrected for Ambiant energy and cluster leakage
  //--> An additional correction available in the Commons class
  //--> is applied to mitigate the residual pT dependency at very high pT.
  float pt;
  if(m_data)
    pt = PhotonRescaledPt_geo21(iph);
  else
    pt = PhotonSmearedPt(iph);

// pt is in GeV, Commons::GetTopoIsoPtcorr is in GeV, PhotonTopoIsolation_tool is in MeV
  double Isol = PhotonTopoIsolation_tool_geo21(iph)/1000. - Commons::GetTopoIsoPtcorr(pt);
  return 1000*Isol;

}
///////////////////////////////////////////////////////
float GravitonAnalysis::PhotonIsolationError_tool(int iph)
///////////////////////////////////////////////////////
{
  //This method calculate the error as the difference
  //between the central value corrected by the tool
  // and the shifted value obtained by the tool
  // float energy=(float)m_Rd->GetVariable( Form("ph_cl_E[%d]",iph) );
  float energy      = 1000*PhotonRescaledE(iph);
  float etaS2       = (float)m_Rd->GetVariable( Form("ph_etas2[%d]",iph) );
  float etaPointing = (float)m_Rd->GetVariable( Form("ph_etap[%d]",iph) );
  float etaCluster  = (float)m_Rd->GetVariable( Form("ph_cl_eta[%d]",iph) );
  float radius      = 40;
  bool is_mc        = false;
  // the correction for the central value is only for data
  // for Etcone40_corrected
  float Isol_error = CaloIsoCorrection::GetPtCorrectedIsolationError( energy,
								      etaS2,
								      etaPointing,
								      etaCluster,
								      radius,
								      is_mc, 
								      CaloIsoCorrection::PHOTON );
  return Isol_error;
}
//////////////////////////////////////////////////
float GravitonAnalysis::PhotonRescaledPt(int iph)
/////////////////////////////////////////////////
{
  //--> Return the rescaled pT from egammaAnalysisUtils
  double cl_eta  = m_Rd->GetVariable( Form("ph_cl_eta[%d]",iph) );
  // double cl_phi  = m_Rd->GetVariable( Form("ph_cl_phi[%d]",iph) );
  // double cl_et   = m_Rd->GetVariable( Form("ph_cl_pt[%d]",iph))/1000.;
  double unCorrE = m_Rd->GetVariable( Form("ph_cl_E[%d]",iph) );
  double CorrE = -999;;
  // std::string part_type;
  if( (int)m_Rd->GetVariable( Form("ph_isConv[%d]",iph) )>0)
    CorrE = m_Rescale.applyEnergyCorrection( cl_eta,unCorrE,
					     egRescaler::EnergyRescalerUpgrade::Converted,
					     egRescaler::EnergyRescalerUpgrade::Nominal);
    // part_type="CONVERTED_PHOTON";
  else
    CorrE = m_Rescale.applyEnergyCorrection( cl_eta,unCorrE,
					     egRescaler::EnergyRescalerUpgrade::Unconverted,
					     egRescaler::EnergyRescalerUpgrade::Nominal);
    // part_type="UNCONVERTED_PHOTON";
  
  // double cl_e_corr = m_Rescale.applyEnergyCorrectionGeV(cl_eta,cl_phi,unCorrE,cl_et,0,part_type);

  double etaS2   = m_Rd->GetVariable( Form("ph_etas2[%d]",iph) );
  float pt=(float)CorrE/TMath::CosH(etaS2)/1000.;
  return pt;
}
//////////////////////////////////////////////////
float GravitonAnalysis::PhotonRescaledE(int iph)
/////////////////////////////////////////////////
{
  //--> Returns the rescaled Energy from egammaAnalysisUtils
  double cl_eta  = m_Rd->GetVariable( Form("ph_cl_eta[%d]",iph) );
  // double cl_phi  = m_Rd->GetVariable( Form("ph_cl_phi[%d]",iph) );
  double unCorrE = m_Rd->GetVariable( Form("ph_cl_E[%d]",iph) );
  // double cl_et   = m_Rd->GetVariable( Form("ph_cl_pt[%d]",iph))/1000.;
  double CorrE   = -999;
  // std::string part_type;
  if( (int)m_Rd->GetVariable( Form("ph_isConv[%d]",iph) )>0){
    CorrE = m_Rescale.applyEnergyCorrection( cl_eta,unCorrE,
					     egRescaler::EnergyRescalerUpgrade::Converted,
					     egRescaler::EnergyRescalerUpgrade::Nominal);
    // part_type="CONVERTED_PHOTON";
  }else{
    CorrE = m_Rescale.applyEnergyCorrection( cl_eta,unCorrE,
					     egRescaler::EnergyRescalerUpgrade::Unconverted,
					     egRescaler::EnergyRescalerUpgrade::Nominal);
    // part_type="UNCONVERTED_PHOTON";
  }

  // float cl_e_corr = m_Rescale.applyEnergyCorrectionGeV(cl_eta,cl_phi,unCorrE,cl_et,0,part_type);
  return CorrE/1000.;
}

//////////////////////////////////////////////////
float GravitonAnalysis::PhotonRescaledE_geo21(int iph)
/////////////////////////////////////////////////
{
  // float local_eta = (float)m_Rd->GetVariable(Form("ph_cl_eta[%d]", iph));
  
  egRescaler::EnergyRescalerUpgrade::ParticleType ptype;
  if ((int)m_Rd->GetVariable(Form("ph_convFlag[%d]", iph)) % 10 == 0)
    ptype  = egRescaler::EnergyRescalerUpgrade::Unconverted;
  else
    ptype = egRescaler::EnergyRescalerUpgrade::Converted;


  PATCore::ParticleDataType::DataType dataType = PATCore::ParticleDataType::Data;
  // Scale::Variation scaleVar = Scale::Nominal;
  // Resolution::Variation resVar = Resolution::None;

  double varSF = 1.;

  if (m_Rd->GetVariable(Form("ph_cl_E[%d]", iph)) <= 0.)
    return 0.;

  double corrected_energy = m_Rescale_geo21.getCorrectedEnergyPhoton(
								    (unsigned int)m_Rd->GetVariable("RunNumber"),
								    PATCore::ParticleDataType::Data,
								    m_Rd->GetVariable(Form("ph_rawcl_Es0[%d]", iph)),
								    m_Rd->GetVariable(Form("ph_rawcl_Es1[%d]", iph)),
								    m_Rd->GetVariable(Form("ph_rawcl_Es2[%d]", iph)),
								    m_Rd->GetVariable(Form("ph_rawcl_Es3[%d]", iph)),
								    m_Rd->GetVariable(Form("ph_cl_eta[%d]", iph)),
								    m_Rd->GetVariable(Form("ph_cl_phi[%d]", iph)),
								    m_Rd->GetVariable(Form("ph_cl_E[%d]", iph)),
								    m_Rd->GetVariable(Form("ph_cl_etaCalo[%d]", iph)),
								    m_Rd->GetVariable(Form("ph_cl_phiCalo[%d]", iph)),
								    m_Rd->GetVariable(Form("ph_ptconv[%d]", iph)),
								    m_Rd->GetVariable(Form("ph_pt1conv[%d]", iph)),
 								    m_Rd->GetVariable(Form("ph_pt2conv[%d]", iph)),
								    (int)m_Rd->GetVariable(Form("ph_convtrk1nPixHits[%d]", iph)),
								    (int)m_Rd->GetVariable(Form("ph_convtrk1nSCTHits[%d]", iph)),
								    (int)m_Rd->GetVariable(Form("ph_convtrk2nPixHits[%d]", iph)),
								    (int)m_Rd->GetVariable(Form("ph_convtrk2nSCTHits[%d]", iph)),
								    m_Rd->GetVariable(Form("ph_Rconv[%d]", iph)),
								    egEnergyCorr::Scale::Nominal,
								    egEnergyCorr::Resolution::None,
								    egEnergyCorr::Resolution::SigmaEff90,
								    varSF);
  								    
  return (float)corrected_energy / 1000.;
}


//////////////////////////////////////////////////
float GravitonAnalysis::PhotonRescaledPt_geo21(int iph)
/////////////////////////////////////////////////
{
  float corr_e = PhotonRescaledE_geo21(iph);
  float corr_pt = corr_e / TMath::CosH(m_Rd->GetVariable(Form("ph_etas2[%d]", iph)));
  return corr_pt;
}

//////////////////////////////////////////////////
float GravitonAnalysis::ElectronRescaledE(int iele)
/////////////////////////////////////////////////
{
  //--> Returns the rescale E from egammaAnalysisUtils
  double cl_eta = m_Rd->GetVariable( Form("el_cl_eta[%d]",iele) );
  // double cl_phi = m_Rd->GetVariable( Form("el_cl_phi[%d]",iele) );
  double unCorrE= m_Rd->GetVariable( Form("el_cl_E[%d]",iele) );

  // double cl_et  = 0;
  // int nPIX=(int)m_Rd->GetVariable( Form("el_nPixHits[%d]",iele) );
  // int nSCT=(int)m_Rd->GetVariable( Form("el_nSCTHits[%d]",iele) );

  // if(nPIX+nSCT>=4){
  //   cl_et=unCorrE/TMath::Cos(m_Rd->GetVariable( Form("el_tracketa[%d]",iele) ));
  // }else{
  //   cl_et=m_Rd->GetVariable( Form("el_cl_pt[%d]",iele) )/1000.;
  // }
  // std::string part_type="ELECTRON";
  // double cl_e_corr = m_Rescale.applyEnergyCorrectionGeV(cl_eta,cl_phi,unCorrE,cl_et,0,part_type);


  float CorrE = m_Rescale.applyEnergyCorrection( cl_eta,unCorrE,
						 egRescaler::EnergyRescalerUpgrade::Electron,
						 egRescaler::EnergyRescalerUpgrade::Nominal);

  return CorrE/1000.;
}
/////////////////////////////////////////////////////////
float GravitonAnalysis::PhotonSmearedPt(int iph,int type)
////////////////////////////////////////////////////////
{
  //--> Returns the smeared pT from egammaAnalysisUtils
  float E=(float)m_Rd->GetVariable( Form("ph_cl_E[%d]",iph) )/1000.;
  float smearcorr = PhotonSmearingFactor(iph,type);
  float corrE=smearcorr*E;
  float etas2=(float)m_Rd->GetVariable( Form("ph_etas2[%d]",iph) );
  float pt=corrE/TMath::CosH(etas2);
  return pt;
}
////////////////////////////////////////////////////////
float GravitonAnalysis::PhotonSmearedE(int iph,int type)
///////////////////////////////////////////////////////
{
  //--> Returns the smeared E from egammaAnalysisUtils
  float E=(float)m_Rd->GetVariable( Form("ph_cl_E[%d]",iph) )/1000.;
  float smearcorr = PhotonSmearingFactor(iph,type);
  float corrE=smearcorr*E;
  return corrE;
}
//////////////////////////////////////////////////////////////
float GravitonAnalysis::PhotonSmearingFactor(int iph,int type)
/////////////////////////////////////////////////////////////
{
  //--> Returns the smearign factor 
  // int evt=(int)m_Rd->GetVariable("EventNumber");
  // m_Rescale.SetRandomSeed(1771561+evt+(iph*10));
  // bool mcWithConstantTerm=false;// if you use a MC without constant term,set true 
  float Ene  = (float)m_Rd->GetVariable( Form("ph_cl_E[%d]",iph) );
  float eta  =(float)m_Rd->GetVariable(Form("ph_cl_eta[%d]",iph));
  // float smearcorr = m_Rescale.getSmearingCorrectionGeV(eta, E, type, mcWithConstantTerm,"2012");
  float smearcorr = 1;
  if( type == 0)
    smearcorr = m_Rescale.getSmearingCorrection( eta, Ene, 
						 egRescaler::EnergyRescalerUpgrade::NOMINAL);
  else if( type==1)
    smearcorr = m_Rescale.getSmearingCorrection(eta, Ene, 
						egRescaler::EnergyRescalerUpgrade::ERR_DOWN);
  else if( type==2) 
    smearcorr = m_Rescale.getSmearingCorrection(eta, Ene, 
						egRescaler::EnergyRescalerUpgrade::ERR_UP);
  else Fatal("GravitonAnalysis::PhotonSmearingFactor","Wrong type !");

  return smearcorr;
}
//////////////////////////////////////////////////
unsigned int GravitonAnalysis::PhotonIsEM(int iph)
//////////////////////////////////////////////////
{ 
  //-->Returns the isem word recomputed from in D3PD for data events
  double et    = 1000*PhotonRescaledPt(iph);
  // double et    = m_Rd->GetVariable(Form("ph_cl_pt[%d]",iph));
  double etaS2 = m_Rd->GetVariable( Form("ph_etas2[%d]",iph) );
  //--------------------------------------------------------------------------------------
  double dv_rhad1  = m_Rd->GetVariable(Form("ph_Ethad1[%d]",iph))/et;
  double dv_rhad   = m_Rd->GetVariable(Form("ph_Ethad[%d]",iph))/et;
  double dv_e277   = m_Rd->GetVariable(Form("ph_E277[%d]",iph)); 
  double dv_reta   = m_Rd->GetVariable(Form("ph_E277[%d]",iph)) != 0. ? m_Rd->GetVariable(Form("ph_E237[%d]",iph))/m_Rd->GetVariable(Form("ph_E277[%d]",iph)) : 0. ;
  double dv_rphi   = m_Rd->GetVariable(Form("ph_E237[%d]",iph)) != 0. ? m_Rd->GetVariable(Form("ph_E233[%d]",iph))/m_Rd->GetVariable(Form("ph_E237[%d]",iph)) : 0. ;
  double dv_weta2  = m_Rd->GetVariable(Form("ph_weta2[%d]",iph));
  double dv_f1     = m_Rd->GetVariable(Form("ph_f1[%d]",iph)); 
  double dv_fside  = m_Rd->GetVariable(Form("ph_fside[%d]",iph));
  double dv_wtot   = m_Rd->GetVariable(Form("ph_wstot[%d]",iph));
  double dv_w1     = m_Rd->GetVariable(Form("ph_ws3[%d]",iph)); 
  double dv_deltae = ( m_Rd->GetVariable(Form("ph_Emax2[%d]",iph)) - m_Rd->GetVariable(Form("ph_Emins1[%d]",iph)) ); 
  double dv_eratio = ( m_Rd->GetVariable(Form("ph_emaxs1[%d]",iph)) + m_Rd->GetVariable(Form("ph_Emax2[%d]",iph)) ) != 0. ? ( m_Rd->GetVariable(Form("ph_emaxs1[%d]",iph)) - m_Rd->GetVariable(Form("ph_Emax2[%d]",iph)) )/( m_Rd->GetVariable(Form("ph_emaxs1[%d]",iph)) + m_Rd->GetVariable(Form("ph_Emax2[%d]",iph)) ) : 0.; 
  int dv_conv      = (int)m_Rd->GetVariable(Form("ph_isConv[%d]",iph));
  //--------------------------------------------------------------------------------------
  PhotonIDTool robust_isEM_nominal=PhotonIDTool( et,etaS2,
						 dv_rhad1,dv_rhad,dv_e277,
						 dv_reta,dv_rphi,dv_weta2,
						 dv_f1,dv_fside,dv_wtot,
						 dv_w1,dv_deltae,dv_eratio,
						 dv_conv ); 
  //----------------------------------------------------------------------  
  return robust_isEM_nominal.isEM(4,2012);
}
/////////////////////////////////////////////////////////////////
bool GravitonAnalysis::PhotonIsEM_OK(int iph,bool loose0_tight1)
/////////////////////////////////////////////////////////////////
{ 
  //--> Returns the Photon ID decision for loose or tight
  //--> Recomputed on the D3PDs for data events
  double et    = 1000. * PhotonRescaledPt_geo21(iph);
  // double et    = m_Rd->GetVariable(Form("ph_cl_pt[%d]",iph));
  double etaS2 = m_Rd->GetVariable( Form("ph_etas2[%d]",iph) );
  //--------------------------------------------------------------------------------------
  double dv_rhad1  = m_Rd->GetVariable(Form("ph_Ethad1[%d]",iph))/et;
  double dv_rhad   = m_Rd->GetVariable(Form("ph_Ethad[%d]",iph))/et;
  double dv_e277   = m_Rd->GetVariable(Form("ph_E277[%d]",iph)); 
  double dv_reta   = m_Rd->GetVariable(Form("ph_E277[%d]",iph)) != 0. ? m_Rd->GetVariable(Form("ph_E237[%d]",iph))/m_Rd->GetVariable(Form("ph_E277[%d]",iph)) : 0. ;
  double dv_rphi   = m_Rd->GetVariable(Form("ph_E237[%d]",iph)) != 0. ? m_Rd->GetVariable(Form("ph_E233[%d]",iph))/m_Rd->GetVariable(Form("ph_E237[%d]",iph)) : 0. ;
  double dv_weta2  = m_Rd->GetVariable(Form("ph_weta2[%d]",iph));
  double dv_f1     = m_Rd->GetVariable(Form("ph_f1[%d]",iph)); 
  double dv_fside  = m_Rd->GetVariable(Form("ph_fside[%d]",iph));
  double dv_wtot   = m_Rd->GetVariable(Form("ph_wstot[%d]",iph));
  double dv_w1     = m_Rd->GetVariable(Form("ph_ws3[%d]",iph)); 
  double dv_deltae = ( m_Rd->GetVariable(Form("ph_Emax2[%d]",iph)) - m_Rd->GetVariable(Form("ph_Emins1[%d]",iph)) ); 
  double dv_eratio = ( m_Rd->GetVariable(Form("ph_emaxs1[%d]",iph)) + m_Rd->GetVariable(Form("ph_Emax2[%d]",iph)) ) != 0. ? ( m_Rd->GetVariable(Form("ph_emaxs1[%d]",iph)) - m_Rd->GetVariable(Form("ph_Emax2[%d]",iph)) )/( m_Rd->GetVariable(Form("ph_emaxs1[%d]",iph)) + m_Rd->GetVariable(Form("ph_Emax2[%d]",iph)) ) : 0.; 
  int dv_conv      = (int)m_Rd->GetVariable(Form("ph_isConv[%d]",iph));
  //--------------------------------------------------------------------------------------
  PhotonIDTool robust_isEM_nominal=PhotonIDTool( et,etaS2,
						 dv_rhad1,dv_rhad,dv_e277,
						 dv_reta,dv_rphi,dv_weta2,
						 dv_f1,dv_fside,dv_wtot,
						 dv_w1,dv_deltae,dv_eratio,
						 dv_conv ); 
  //----------------------------------------------------------------------  
  if(!loose0_tight1)
    return robust_isEM_nominal.PhotonCutsLoose(4);
  else
    return robust_isEM_nominal.PhotonCutsTight(2012);
}
//////////////////////////////////////////////////////////////////////////////
bool GravitonAnalysis::ElectronIsEM_OK(int iele,int isConv,bool loose0_tight1)
//////////////////////////////////////////////////////////////////////////////
{ 
  //--> Returns the Photon ID decision for loose or tight
  //--> applied on electrons
  //--> Recomputed on the D3PDs for data events
  
  double cl_eta = m_Rd->GetVariable( Form("el_cl_eta[%d]",iele) );
  double et     = 1000*ElectronRescaledE(iele)/TMath::Cos(cl_eta);
  //----------------------------------------------------------------------
  double dv_rhad1  = m_Rd->GetVariable(Form("el_Ethad1[%d]",iele))/et;
  double dv_rhad   = m_Rd->GetVariable(Form("el_Ethad[%d]",iele))/et;
  double dv_e277   = m_Rd->GetVariable(Form("el_E277[%d]",iele)); 
  double dv_reta   = m_Rd->GetVariable(Form("el_E277[%d]",iele)) != 0. ? m_Rd->GetVariable(Form("el_E237[%d]",iele))/m_Rd->GetVariable(Form("el_E277[%d]",iele)) : 0. ;
  double dv_rphi   = m_Rd->GetVariable(Form("el_E237[%d]",iele)) != 0. ? m_Rd->GetVariable(Form("el_E233[%d]",iele))/m_Rd->GetVariable(Form("el_E237[%d]",iele)) : 0. ;
  double dv_weta2  = m_Rd->GetVariable(Form("el_weta2[%d]",iele));
  double dv_f1     = m_Rd->GetVariable(Form("el_f1[%d]",iele)); 
  double dv_fside  = m_Rd->GetVariable(Form("el_fside[%d]",iele));
  double dv_wtot   = m_Rd->GetVariable(Form("el_wstot[%d]",iele));
  double dv_w1     = m_Rd->GetVariable(Form("el_ws3[%d]",iele)); 
  double dv_deltae = ( m_Rd->GetVariable(Form("el_Emax2[%d]",iele)) - m_Rd->GetVariable(Form("el_Emins1[%d]",iele)) ); 
  double dv_eratio = ( m_Rd->GetVariable(Form("el_emaxs1[%d]",iele)) + m_Rd->GetVariable(Form("el_Emax2[%d]",iele)) ) != 0. ? ( m_Rd->GetVariable(Form("el_emaxs1[%d]",iele)) - m_Rd->GetVariable(Form("el_Emax2[%d]",iele)) )/( m_Rd->GetVariable(Form("el_emaxs1[%d]",iele)) + m_Rd->GetVariable(Form("el_Emax2[%d]",iele)) ) : 0.; 
  int dv_conv      = isConv;
  //---------------------------------------------------------------------------------
  PhotonIDTool robust_isEM=PhotonIDTool( et,cl_eta,
						 dv_rhad1,dv_rhad,dv_e277,
						 dv_reta,dv_rphi,dv_weta2,
						 dv_f1,dv_fside,dv_wtot,
						 dv_w1,dv_deltae,dv_eratio,
						 dv_conv ); 
  //----------------------------------------------------------------------------------  
  if(!loose0_tight1)
    return robust_isEM.PhotonCutsLoose(4);
  else
    return robust_isEM.PhotonCutsTight(2012);
}
///////////////////////////////////////////////////////////////////////////////
bool GravitonAnalysis::PhotonIsEM_MCOK(int iph,bool loose0_tight1,bool VERBOSE)
//////////////////////////////////////////////////////////////////////////////
{ 
 
  //--> Returns the Photon ID decision for loose or tight
  //--> Recomputed on the D3PDs for mc events
  //--> The shower shape variables are corrected using the Fudge Factors

  float E          = m_Rd->GetVariable( Form("ph_cl_E[%d]",iph) )/1000.;
  float smearcorr  = PhotonSmearingFactor(iph,0);
  float corrE      = smearcorr*E;
  double et        = 1000*corrE/cosh(m_Rd->GetVariable( Form("ph_etas2[%d]",iph) ));
  double etaS2     = m_Rd->GetVariable( Form("ph_etas2[%d]",iph) );
  //------------------------------------------------------------------
  double dv_rhad1  = m_Rd->GetVariable(Form("ph_Ethad1[%d]",iph))/et;
  double dv_rhad   = m_Rd->GetVariable(Form("ph_Ethad[%d]",iph))/et;
  double dv_e277   = m_Rd->GetVariable(Form("ph_E277[%d]",iph)); 
  double dv_reta   = m_Rd->GetVariable(Form("ph_E277[%d]",iph)) != 0. ? m_Rd->GetVariable(Form("ph_E237[%d]",iph))/m_Rd->GetVariable(Form("ph_E277[%d]",iph)) : 0. ;
  double dv_rphi   = m_Rd->GetVariable(Form("ph_E237[%d]",iph)) != 0. ? m_Rd->GetVariable(Form("ph_E233[%d]",iph))/m_Rd->GetVariable(Form("ph_E237[%d]",iph)) : 0. ;
  double dv_weta2  = m_Rd->GetVariable(Form("ph_weta2[%d]",iph));
  double dv_f1     = m_Rd->GetVariable(Form("ph_f1[%d]",iph)); 
  double dv_fside  = m_Rd->GetVariable(Form("ph_fside[%d]",iph));
  double dv_wtot   = m_Rd->GetVariable(Form("ph_wstot[%d]",iph));
  double dv_w1     = m_Rd->GetVariable(Form("ph_ws3[%d]",iph)); 
  double dv_deltae = ( m_Rd->GetVariable(Form("ph_Emax2[%d]",iph)) - m_Rd->GetVariable(Form("ph_Emins1[%d]",iph)) ); 
  double dv_eratio = ( m_Rd->GetVariable(Form("ph_emaxs1[%d]",iph)) + m_Rd->GetVariable(Form("ph_Emax2[%d]",iph)) ) != 0. ? ( m_Rd->GetVariable(Form("ph_emaxs1[%d]",iph)) - m_Rd->GetVariable(Form("ph_Emax2[%d]",iph)) )/( m_Rd->GetVariable(Form("ph_emaxs1[%d]",iph)) + m_Rd->GetVariable(Form("ph_Emax2[%d]",iph)) ) : 0.; 
  int dv_conv      = (int)m_Rd->GetVariable(Form("ph_isConv[%d]",iph));
  //-----------------------------------------------------------------------  

  FudgeMCTool myFFTool( et, etaS2, dv_conv, 0 ); 
  int mySelection = 14; 
  myFFTool.SetPreselection( mySelection );	
  myFFTool.FudgeShowers( et,etaS2,
			 dv_rhad1,dv_rhad,dv_e277,
			 dv_reta,dv_rphi,dv_weta2,
			 dv_f1,dv_fside,dv_wtot,
			 dv_w1,dv_deltae,dv_eratio,
			 dv_conv ); 
  //----------------------------------------------------------------------
  PhotonIDTool robust_isEM = PhotonIDTool( et,etaS2,
					   dv_rhad1,dv_rhad,dv_e277,
					   dv_reta,dv_rphi,dv_weta2,
					   dv_f1,dv_fside,dv_wtot,
					   dv_w1,dv_deltae,dv_eratio,
					   dv_conv ); 
  //----------------------------------------------------------------------  
  if(VERBOSE){
    std::cout << "et    : " << et       << std::endl;
    std::cout << "etaS2 : " << etaS2    << std::endl;
    std::cout << "rhad1 : " << dv_rhad1 << std::endl;
    std::cout << "rhad  : " << dv_rhad  << std::endl;
    std::cout << "e277  : " << dv_e277  << std::endl;
    std::cout << "reta  : " << dv_reta  << std::endl;
    std::cout << "rphi  : " << dv_rphi  << std::endl;
    std::cout << "weta2 : " << dv_weta2 << std::endl;
    std::cout << "f1    : " << dv_f1    << std::endl;
    std::cout << "fside : " << dv_fside << std::endl;
    std::cout << "wtot  : " << dv_wtot  << std::endl;
    std::cout << "w1    : " << dv_w1    << std::endl;
    std::cout << "deltae: " << dv_deltae<< std::endl;
    std::cout << "eratio: " << dv_eratio<< std::endl;
    std::cout << "conv  : " << dv_conv  << std::endl;  
  }
  if(!loose0_tight1)
    return robust_isEM.PhotonCutsLoose(4);
  else
    return robust_isEM.PhotonCutsTight(2012);


}
////////////////////////////////////////////////////////////////////////////////
double GravitonAnalysis::scaleForFFUnconvPhoton(double pT,double eta,int isconv)
////////////////////////////////////////////////////////////////////////////////
{
  // temporary obtained by comparing FFMC to electron extrapolation
  // apply this scale factor for FFed unconverted photon
  // in 1.81 < |eta| < 2.37 (use ph_eta2 in DPD)
  // pT in GeV
  if( isconv ==1 ){
    return 1;
  }else{
    if( fabs(eta)<1.81 ){
      return 1;
    }else{
      if (pT<25.) return 0.86;
      else if (pT<30.) return 0.89;
      else if (pT<35.) return 0.94;
      else if (pT<40.) return 0.92;
      else if (pT<45.) return 0.98;
      else return 0.97;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////////////
unsigned long long GravitonAnalysis::make_runevent(unsigned int run, unsigned int event)
////////////////////////////////////////////////////////////////////////////////////////
{ return static_cast<unsigned long long>(run) << 32 | event; }
///////////////////////////////////////////////////////////////////////////////
std::set<unsigned long long> GravitonAnalysis::build_eeset(const char* filename)
///////////////////////////////////////////////////////////////////////////////
{
  std::set<unsigned long long> contents;
  std::ifstream eventfile(filename);
  if (!eventfile.is_open())
    throw;
  int run = 0, event = 0, lb=0;
  while (!eventfile.eof()){
    eventfile >> run >> event>> lb;
    if (!eventfile.good())
      break;
    contents.insert(make_runevent(run, event));
  }
  return contents;
}
////////////////////////////////////////////////////////////////////////////////////
bool GravitonAnalysis::EventInZprimme(const std::set<unsigned long long>& container, 
				      unsigned int run, unsigned int event)
////////////////////////////////////////////////////////////////////////////////////
{ return container.find(make_runevent(run, event)) != container.end(); }
