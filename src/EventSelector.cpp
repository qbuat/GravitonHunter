#include "EventSelector.h"
#include "ToolsUtilities.h"
#include "ToolsCommons.h"
#include "InvMassCalculator.h"
#include "egammaAnalysisUtils/CaloIsoCorrection.h"
#include <TError.h>

#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <cstring>
#include <assert.h>

////////////////////////////////////////////////////////////////////////////
EventSelector::EventSelector(TTree* tree) : GravitonAnalysis(tree)
///////////////////////////////////////////////////////////////////////////
{
  if (tree == 0){
    std::cout << "You need to specify an input tree !!!"<<std::endl;
  }
}
///////////////////////////////////////////////////
EventSelector::EventSelector() : GravitonAnalysis()
///////////////////////////////////////////////////
{
  //Default constructor
}
//////////////////////////////////////
EventSelector::~EventSelector()
/////////////////////////////////////
{
  //Default destructor
}
/////////////////////////////////
bool EventSelector::InitOutTree()
/////////////////////////////////
{
  //-----> Output variables  
  file_out   = new TFile(m_output,"RECREATE");
  tree_out   = new TTree("tree","tree");
  hist_mgg   = new TObjArray();
  hist_mgg_w = new TObjArray();


  //------> Event informations
  tree_out->Branch("RunNumber",&m_RunNumber,"RunNumber/I");
  tree_out->Branch("EventNumber",&m_EventNumber,"EventNumber/I");
  tree_out->Branch("LumiBlock",&m_LumiBlock,"LumiBlock/I");
  tree_out->Branch("NPV",&m_NPV,"NPV/I");
  tree_out->Branch("mu",&m_mu,"mu/D");
  tree_out->Branch("PV_z",&PV_ID,"PV_z/D");


  //------> MET variables
  tree_out->Branch("MET_Topo_et",&m_MET_Topo_et,"MET_Topo_et/D");
  tree_out->Branch("MET_Topo_sumet",&m_MET_Topo_sumet,"MET_Topo_sumet/D");
  tree_out->Branch("MET_Topo_sumet_EMB",&m_MET_Topo_sumet_EMB,"MET_Topo_sumet_EMB/D");
  tree_out->Branch("MET_Topo_sumet_EME",&m_MET_Topo_sumet_EME,"MET_Topo_sumet_EME/D");
  tree_out->Branch("MET_Topo_sumet_CentralReg",&m_MET_Topo_sumet_CentralReg,"MET_Topo_sumet_CentralReg/D");
  tree_out->Branch("MET_Topo_sumet_EndcapRegion",&m_MET_Topo_sumet_EndcapRegion,"MET_Topo_sumet_EndcapRegion/D");
  tree_out->Branch("MET_RefFinal_et",&m_MET_RefFinal_et,"MET_RefFinal_et/D");
  tree_out->Branch("MET_RefFinal_sumet",&m_MET_RefFinal_sumet,"MET_RefFinal_sumet/D");
  tree_out->Branch("MET_RefFinal_sumet_CentralReg",&m_MET_RefFinal_sumet_CentralReg,"MET_RefFinal_sumet_CentralReg/D");
  tree_out->Branch("MET_RefFinal_sumet_EndcapRegion",&m_MET_RefFinal_sumet_EndcapRegion,"MET_RefFinal_sumet_EndcapRegion/D");

  //------> Variables for the diphoton system
  tree_out->Branch("mgg",&mgg,"mgg/D");
  tree_out->Branch("ptgg",&ptgg,"ptgg/D");
  tree_out->Branch("ygg",&ygg,"ygg/D");
  tree_out->Branch("costhetastar",&costhetastar,"costhetastar/D");
  tree_out->Branch("deltaphi",&deltaphi,"deltaphi/D");
  //---
  tree_out->Branch("mgg_smeareddown",&mgg_smeareddown,"mgg_smeareddown/D");
  tree_out->Branch("ptgg_smeareddown",&ptgg_smeareddown,"ptgg_smeareddown/D");
  tree_out->Branch("ygg_smeareddown",&ygg_smeareddown,"ygg_smeareddown/D");
  tree_out->Branch("costhetastar_smeareddown",&costhetastar_smeareddown,"costhetastar_smeareddown/D");
  tree_out->Branch("deltaphi_smeareddown",&deltaphi_smeareddown,"deltaphi_smeareddown/D");
  //---
  tree_out->Branch("mgg_smearedup",&mgg_smearedup,"mgg_smearedup/D");
  tree_out->Branch("ptgg_smearedup",&ptgg_smearedup,"ptgg_smearedup/D");
  tree_out->Branch("ygg_smearedup",&ygg_smearedup,"ygg_smearedup/D");
  tree_out->Branch("costhetastar_smearedup",&costhetastar_smearedup,"costhetastar_smearedup/D");
  tree_out->Branch("deltaphi_smearedup",&deltaphi_smearedup,"deltaphi_smearedup/D");

  //------> Variables for the leading photon
  tree_out->Branch("pT_smearedup_L",&Lead_pT_smearedup,"pT_smearedup_L/D");
  tree_out->Branch("pT_smeareddown_L",&Lead_pT_smeareddown,"pT_smeareddown_L/D");
  tree_out->Branch("pT_L",&Lead_pT,"pT_L/D");
  tree_out->Branch("phi_L",&Lead_phi,"phi_L/D");
  tree_out->Branch("eta_L",&Lead_eta,"eta_L/D");
  tree_out->Branch("eta_PV_L",&Lead_eta_PV,"eta_PV_L/D");
  tree_out->Branch("ED_med_L", &Lead_ed_med, "ED_med_L/D");
  tree_out->Branch("Iso_uncor_L", &Lead_iso_uncor, "Iso_uncor_L/D");
  tree_out->Branch("Iso_ptcorr_L", &Lead_iso_ptcorr, "Iso_ptcorr_L/D");
  tree_out->Branch("Iso_L",&Lead_iso,"Iso_L/D");
  if(m_isotype == "TOPO_PTDEP")
    tree_out->Branch("Iso_L_mod",&Lead_iso_mod,"Iso_L_mod/D");
  tree_out->Branch("IsTight_L",&Lead_IsTight,"IsTight_L/I");
  tree_out->Branch("IsLoosePrime2_L",&Lead_IsLoosePrime2,"IsLoosePrime2_L/I");
  tree_out->Branch("IsLoosePrime3_L",&Lead_IsLoosePrime3,"IsLoosePrime3_L/I");
  tree_out->Branch("IsLoosePrime4_L",&Lead_IsLoosePrime4,"IsLoosePrime4_L/I");
  tree_out->Branch("IsLoosePrime5_L",&Lead_IsLoosePrime5,"IsLoosePrime5_L/I");
  tree_out->Branch("IsConv_L",&Lead_IsConv,"IsConv_L/I");
  //------> Variables for the subleading photon
  tree_out->Branch("pT_SL",&SubLead_pT,"pT_SL/D");
  tree_out->Branch("pT_smearedup_SL",&SubLead_pT_smearedup,"pT_smearedup_SL/D");
  tree_out->Branch("pT_smeareddown_SL",&SubLead_pT_smeareddown,"pT_smeareddown_SL/D");
  tree_out->Branch("phi_SL",&SubLead_phi,"phi_SL/D");
  tree_out->Branch("eta_SL",&SubLead_eta,"eta_SL/D");
  tree_out->Branch("eta_PV_SL",&SubLead_eta_PV,"eta_PV_SL/D");
  tree_out->Branch("ED_med_SL", &SubLead_ed_med, "ED_med_SL/D");
  tree_out->Branch("Iso_uncor_SL", &SubLead_iso_uncor, "Iso_uncor_SL/D");
  tree_out->Branch("Iso_ptcorr_SL", &SubLead_iso_ptcorr, "Iso_ptcorr_SL/D");
  tree_out->Branch("Iso_SL",&SubLead_iso,"Iso_SL/D");
  if(m_isotype == "TOPO_PTDEP")
    tree_out->Branch("Iso_SL_mod",&SubLead_iso_mod,"Iso_SL_mod/D");
  tree_out->Branch("IsTight_SL",&SubLead_IsTight,"IsTight_SL/I");
  tree_out->Branch("IsLoosePrime2_SL",&SubLead_IsLoosePrime2,"IsLoosePrime2_SL/I");
  tree_out->Branch("IsLoosePrime3_SL",&SubLead_IsLoosePrime3,"IsLoosePrime3_SL/I");
  tree_out->Branch("IsLoosePrime4_SL",&SubLead_IsLoosePrime4,"IsLoosePrime4_SL/I");
  tree_out->Branch("IsLoosePrime5_SL",&SubLead_IsLoosePrime5,"IsLoosePrime5_SL/I");
  tree_out->Branch("IsConv_SL",&SubLead_IsConv,"IsConv_SL/I");

  //--> Photons cluster variables
  tree_out->Branch("ph_cl_E",       &m_ph_cl_E     );  
  tree_out->Branch("ph_cl_pt",      &m_ph_cl_pt    );  
  tree_out->Branch("ph_cl_eta",     &m_ph_cl_eta   );  
  tree_out->Branch("ph_cl_phi",     &m_ph_cl_phi   );  
  tree_out->Branch("ph_rawcl_E",    &m_ph_rawcl_E  );
  tree_out->Branch("ph_rawcl_pt",   &m_ph_rawcl_pt );
  tree_out->Branch("ph_rawcl_eta",  &m_ph_rawcl_eta);
  tree_out->Branch("ph_rawcl_phi",  &m_ph_rawcl_phi);
  tree_out->Branch("ph_ED_median",  &m_ph_ED_median);



  if( !m_data ){
    //------> Weighting variables
    tree_out->Branch("weight",&weight,"weight/D");
    tree_out->Branch("gen_weight",&gen_weight,"gen_weight/D");
    tree_out->Branch("weight_L",&Lead_weight,"weight_L/I");
    tree_out->Branch("weight_SL",&SubLead_weight,"weight_SL/I");
    //------> Variables at truth level
    tree_out->Branch("truth_eta_L",&Lead_truth_eta,"truth_eta_L/D");
    tree_out->Branch("truth_phi_L",&Lead_truth_phi,"truth_phi_L/D");
    tree_out->Branch("truth_pT_L",&Lead_truth_pT,"truth_pT_L/D");
    tree_out->Branch("truth_status_L",&Lead_truth_status,"truth_status_L/I");
    //---
    tree_out->Branch("truth_eta_SL",&SubLead_truth_eta,"truth_eta_SL/D");
    tree_out->Branch("truth_phi_SL",&SubLead_truth_phi,"truth_phi_SL/D");
    tree_out->Branch("truth_pT_SL",&SubLead_truth_pT,"truth_pT_SL/D");
    tree_out->Branch("truth_status_SL",&SubLead_truth_status,"truth_status_SL/I");
    //---
    tree_out->Branch("truth_mgg",&truth_mgg,"truth_mgg/D");
    tree_out->Branch("truth_ptgg",&truth_ptgg,"truth_ptgg/D");
    tree_out->Branch("truth_ygg",&truth_ygg,"truth_ygg/D");
    tree_out->Branch("truth_dphigg",&truth_dphigg,"truth_dphigg/D");
    tree_out->Branch("truth_costhetastar",&truth_costhetastar,"truth_costhetastar/D");
  }
  if( !m_data && m_mctype == "pythia_rs"){
    //----> True mass (MC only)
    tree_out->Branch("mc_m",&mymc_m,"mc_m/D");
  }
  return true;
}














/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







///////////////////////////////////
void EventSelector::PrintParticleInfo(int index) {

  std::cout<<"--------> "<<index<<std::endl;
  std::cout<<"pt "<<m_Rd->GetVariable(Form("mc_pt[%d]",index))<<std::endl;
  std::cout<<"eta "<<m_Rd->GetVariable(Form("mc_eta[%d]",index))<<std::endl;
  std::cout<<"phi "<<m_Rd->GetVariable(Form("mc_phi[%d]",index))<<std::endl;
  std::cout<<"status "<<m_Rd->GetVariable(Form("mc_status[%d]",index))<<std::endl;
  std::cout<<"pdgid "<<m_Rd->GetVariable(Form("mc_pdgId[%d]",index))<<std::endl;
  std::cout<<"child index "<<m_Rd->GetVariable(Form("mc_child_index[%d]",index))<<std::endl;
  std::cout<<"parent index "<<m_Rd->GetVariable(Form("mc_parent_index[%d]",index))<<std::endl;
  std::cout<<"barcode "<<m_Rd->GetVariable(Form("mc_barcode[%d]",index))<<std::endl;

}


//////////////////////////////////////////
void EventSelector::PrintEventInfo(int limit) {
  
  std::cout<<"======================== eventnumber "<<m_Rd->GetVariable("EventNumber")<<"============================"<<std::endl;
  //  std::cout<<"in file "<<file->GetName()<<std::endl;
  std::cout<<"index \t pt \t\t eta \t\t phi \t\t status \t\t pdgid \t\t child index \t\t parent index"<<std::endl;
  for(int p = 0; p < m_Rd->GetVariable("@mc_pt.size()") && p <= limit; p++) {
    std::cout<<p<<" \t ";
    std::cout<<m_Rd->GetVariable(Form("mc_pt[%d]",p))<<" \t ";
    std::cout<<m_Rd->GetVariable(Form("mc_eta[%d]",p))<<" \t\t ";
    std::cout<<m_Rd->GetVariable(Form("mc_phi[%d]",p))<<" \t\t ";
    std::cout<<m_Rd->GetVariable(Form("mc_status[%d]",p))<<" \t\t ";
    std::cout<<m_Rd->GetVariable(Form("mc_pdgId[%d]",p))<<" \t\t ";
    std::cout<<m_Rd->GetVariable(Form("mc_child_index[%d]",p))<<" \t\t ";
    std::cout<<m_Rd->GetVariable(Form("mc_parent_index[%d]",p))<<" \t\t "<<std::endl;
  }

}




///////////////////////////////////////
int EventSelector::GetOriginalPhotonIndex(int index) {
  
  int index_test = index;
  
  while( true  ) {
    int parent_index = m_Rd->GetVariable(Form("mc_parent_index[%d]",index_test));
    if(parent_index == 0) break;
    if(m_Rd->GetVariable(Form("mc_pdgId[%d]",parent_index)) != 22) break;
    index_test = parent_index;
  }
  
  return index_test;

}

/////////////////////////////////////////////////////////////////////
std::vector<int> EventSelector::HardScatterPhotonIndices(std::vector<int> v_ph_id) {
    
  std::vector<int> v_id;
  for(unsigned int iph = 0; iph < v_ph_id.size(); iph++) {
    int index_test = v_ph_id[iph];

    while( true ) {
      int parent_index = m_Rd->GetVariable(Form("mc_parent_index[%d]",index_test));
      if(parent_index == 0) break;
      if(m_Rd->GetVariable(Form("mc_pdgId[%d]",parent_index)) != 22) break;
      index_test = parent_index;
    }
    v_id.push_back(index_test);
    
  }

  return v_id;

}





//////////////////////////////////////
bool EventSelector::IsPhotonFromHardProc(int index) {
  
  int iph = GetOriginalPhotonIndex(index);

  //  std::cout<<"in IsPhotonFromHardProc"<<iph<<" "<<index<<std::endl;
  if( m_Rd->GetVariable(Form("mc_status[%d]",iph)) == 23 &&
      m_Rd->GetVariable(Form("mc_pdgId[%d]",iph)) == 22 ) 
    return true;
  else return false;
  
}

////////////////////////////////
bool EventSelector::IsPhotonFromUE(int index) {

  int iph = GetOriginalPhotonIndex(index);
  int parent_index = m_Rd->GetVariable(Form("mc_parent_index[%d]",iph));
  if( m_Rd->GetVariable(Form("mc_status[%d]",iph)) == 1 &&
      m_Rd->GetVariable(Form("mc_pdgId[%d]",iph)) == 22 &&
      m_Rd->GetVariable(Form("mc_child_index[%d]",iph)) == 0 &&
      parent_index == 0 ) 
    return true;
  else return false;
    
}




///////////////////////////////////////
bool EventSelector::IsPhotonFromRadiation(int index) {

  if( IsPhotonFromHardProc(index) || IsPhotonFromUE(index) ) return false;

  int iph = GetOriginalPhotonIndex(index);

  //  std::cout<<"in IsPhotonFromRadiation"<<index<<" is associated to mother w/ index "<<iph<<std::endl;
  int parent_index = m_Rd->GetVariable(Form("mc_parent_index[%d]",iph));

  if( parent_index < 2 || m_Rd->GetVariable(Form("mc_status[%d]",parent_index)) == 21 ) return true;
  
  //  PrintParticleInfo(parent_index);
  TLorentzVector lv_p;
  lv_p.SetPtEtaPhiM( m_Rd->GetVariable(Form("mc_pt[%d]",parent_index)),
		     m_Rd->GetVariable(Form("mc_eta[%d]",parent_index)),
		     m_Rd->GetVariable(Form("mc_phi[%d]",parent_index)),
		     m_Rd->GetVariable(Form("mc_m[%d]",parent_index)) );

  if(lv_p.Pt() < 1 || fabs(lv_p.Eta()) > 10 ) return true;

  return false;
  
}


//////////////////////////////////
bool EventSelector::IsPhotonFromFrag(int index) {

  int iph = GetOriginalPhotonIndex(index);
  int parent_index = m_Rd->GetVariable(Form("mc_parent_index[%d]",iph));
  if(parent_index == 0) return false;

  if( m_Rd->GetVariable(Form("mc_pdgId[%d]",parent_index)) != 22 && 
      m_Rd->GetVariable(Form("mc_pt[%d]",parent_index)) > 1 && 
      m_Rd->GetVariable(Form("mc_status[%d]",iph)) != 23 )
    return true;
  else return false;

}



//////////////////////////////////////////////////////////////////////////////
float EventSelector::GetPartonETIso(int index, float cone_size, float parton_minpt_threshold) {

  //  std::cout<<"Isolation for photon index "<<index <<std::endl;
  
  float Eiso = 0;

  TLorentzVector g;
  g.SetPtEtaPhiM( m_Rd->GetVariable(Form("mc_pt[%d]",index)),
		  m_Rd->GetVariable(Form("mc_eta[%d]",index)),
		  m_Rd->GetVariable(Form("mc_phi[%d]",index)),
		  m_Rd->GetVariable(Form("mc_m[%d]",index)) );

  for(int i = 0; i < (int)m_Rd->GetVariable("@mc_pdgId.size()"); i++) {
    if( m_Rd->GetVariable(Form("mc_status[%d]",i)) != 1 ) continue;
    int mypdgid = m_Rd->GetVariable(Form("mc_pdgId[%d]",i));
    if( abs(mypdgid) == 12 || abs(mypdgid) == 14 || abs(mypdgid) == 16 ) continue;  // exclude neutrino
    if( i == index ) continue;
    
    TLorentzVector p;
    p.SetPtEtaPhiM(m_Rd->GetVariable(Form("mc_pt[%d]",i)),
		   m_Rd->GetVariable(Form("mc_eta[%d]",i)),
		   m_Rd->GetVariable(Form("mc_phi[%d]",i)),
		   m_Rd->GetVariable(Form("mc_m[%d]",i)) );
    
    if(p.Pt() < parton_minpt_threshold) continue;
    
    if(g.DeltaR(p) < cone_size) {
      float k = 1;
      if( abs(mypdgid) == 15 ) k=0.1783 * (1+0.1);
      else if( abs(mypdgid) == 13 ) k=0.1;
      Eiso += k*p.Pt();
    }
    
  }


//  if(Eiso > 100000) {
//    std::cout<<"--------------- photon with index "<< index <<" has etiso="<<Eiso<<" with the following particle in the cone"<<std::endl;
//    std::cout<<m_Rd->GetVariable("@mc_pdgId.size()")<<" MC particle in this event"<<std::endl;
//    PrintEventInfo();
//    for(int i = 0; i < (int)m_Rd->GetVariable("@mc_pdgId.size()"); i++) {
//      if( m_Rd->GetVariable(Form("mc_status[%d]",i)) != 1 ) continue;
//      int mypdgid = m_Rd->GetVariable(Form("mc_pdgId[%d]",i));
//      if( abs(mypdgid) == 12 || abs(mypdgid) == 14 || abs(mypdgid) == 16 ) continue;  // exclude neutrino
//      if( i == index ) continue;
//      
//      TLorentzVector p;
//      p.SetPtEtaPhiM(m_Rd->GetVariable(Form("mc_pt[%d]",i)),
//		     m_Rd->GetVariable(Form("mc_eta[%d]",i)),
//		     m_Rd->GetVariable(Form("mc_phi[%d]",i)),
//		     m_Rd->GetVariable(Form("mc_m[%d]",i)) );
//      
//      if(p.Pt() < parton_minpt_threshold) continue;
//      
//      if(g.DeltaR(p) < cone_size) {
//	float k = 1;
//	if( abs(mypdgid) == 15 ) k=0.1783 * (1+0.1);
//	else if( abs(mypdgid) == 13 ) k=0.1;
//	//      Eiso += k*p.Pt();
//	PrintParticleInfo(i);
//	std::cout<<"this particle is";
//	if(IsPhotonFromHardProc(i))  std::cout<<"from hardproc"<<std::endl;
//	if(IsPhotonFromFrag(i))      std::cout<<"from frag"<<std::endl;
//	if(IsPhotonFromUE(i))        std::cout<<"from ue"<<std::endl;
//	if(IsPhotonFromRadiation(i)) std::cout<<"from radiation"<<std::endl;
//
//      }
//      
//    }
//
//  }





  return Eiso;
  
}




/////////////////////////
float EventSelector::ComputeX(int iph) {
  
  if(IsPhotonFromFrag(iph)) {
    TLorentzVector lv_g, lv_q; 
    lv_g.SetPtEtaPhiM(m_Rd->GetVariable(Form("mc_pt[%d]",iph)),
		      m_Rd->GetVariable(Form("mc_eta[%d]",iph)),
		      m_Rd->GetVariable(Form("mc_phi[%d]",iph)),
		      m_Rd->GetVariable(Form("mc_m[%d]",iph)));
    
    int orig_ph_index = GetOriginalPhotonIndex(iph);
    
    //    int iq = m_Rd->GetVariable(Form("mc_parent_index[%d]",iph));
    int iq = m_Rd->GetVariable(Form("mc_parent_index[%d]",orig_ph_index));
    lv_q.SetPtEtaPhiM(m_Rd->GetVariable(Form("mc_pt[%d]",iq)),
		      m_Rd->GetVariable(Form("mc_eta[%d]",iq)),
		      m_Rd->GetVariable(Form("mc_phi[%d]",iq)),
		      m_Rd->GetVariable(Form("mc_m[%d]",iq)));

    //    std::cout<<"frag photon : X = "<<lv_g.E()<<" / "<< lv_q.E()<<" = "<<lv_g.E()/lv_q.E()<<std::endl;
    //    if(lv_g.E()/lv_q.E() > 1 || lv_g.E()/lv_q.E() < 0.1) {
    //      std::cout<<"gamma: pt = "<< lv_g.Pt() <<std::endl;
    //      std::cout<<"       e = "<< lv_g.E() <<std::endl;
    //      std::cout<<"       eta = "<< lv_g.Eta() <<std::endl;
    //      std::cout<<"       phi = "<< lv_g.Phi() <<std::endl;
    //      std::cout<<"       m = "<< lv_g.M() <<std::endl;
    //      PrintParticleInfo(iph);
    //    
    //      std::cout<<"quark: pt = "<< lv_q.Pt() <<std::endl;
    //      std::cout<<"       e = "<< lv_q.E() <<std::endl;
    //      std::cout<<"       eta = "<< lv_q.Eta() <<std::endl;
    //      std::cout<<"       phi = "<< lv_q.Phi() <<std::endl;
    //      std::cout<<"       m = "<< lv_q.M() <<std::endl;
    //      PrintParticleInfo(iq);
    //    
    //    }

    return lv_g.E()/lv_q.E();

  }
  else return -999;
  
}




/////////////////////////////////////////////
////////////////// MAIN /////////////////////
/////////////////////////////////////////////

void EventSelector::TestEventRecord(TString filename) {

  
  TFile * fout = new TFile("testeventrecorttree.root","recreate");
  TTree * out_tree = new TTree("truth_tree","truth_tree");

  int isfromhardproc[3] = {-999,-999,-999};
  int isfromfrag[3] = {-999,-999,-999};
  int isfromradiation[3] = {-999,-999,-999};
  int isfromue[3] = {-999,-999,-999};
  double x[3] = {-999,-999,-999};
  double pt[3] = {-999,-999,-999};
  double eta[3] = {-999,-999,-999};
  double phi[3] = {-999,-999,-999};
  double iso_part[3] = {-999,-999,-999};
  double mgg = 0;

  out_tree->Branch("isfromhardproc",    &isfromhardproc,"isfromhardproc[3]/I");	
  out_tree->Branch("isfromfrag",        &isfromfrag,"isfromfrag[3]/I");	
  out_tree->Branch("isfromradiation",   &isfromradiation,"isfromradiation[3]/I");	
  out_tree->Branch("isfromue",		&isfromue,"isfromue[3]/I");
  out_tree->Branch("x",		        &x,"x[3]/D");
  out_tree->Branch("pt",		&pt,"pt[3]/D");		
  out_tree->Branch("eta",		&eta,"eta[3]/D");
  out_tree->Branch("phi",		&phi,"phi[3]/D");
  out_tree->Branch("iso_part",		&iso_part,"iso_part[3]/D");
  out_tree->Branch("mgg",               &mgg,"mgg/D");
													      
													      




       
  // gammagamma
  //  file = new TFile("/sps/atlas/j/jbrown/DATASETS/mc12a/mc12_8TeV.158339.Pythia8_AU2CTEQ6L1_2DP20_Mass_120_200.merge.NTUP_PHOTON.e1367_s1499_s1504_r3658_r3549_p1032_tid00920673_00/NTUP_PHOTON.00920673._000001.root.1","READ");
  
  // gamma+jet
  //  file = new TFile("/sps/atlas/j/jbrown/DATASETS/mc12a/mc12_8TeV.180189.Pythia8_AU2CTEQ6L1_Photos_gamma_plus_brehm_2DP20_120_200.merge.NTUP_PHOTON.e1713_s1581_s1586_r3658_r3549_p1344_tid01189039_00/NTUP_PHOTON.01189039._000001.root.1","READ");



/////  TTree * t = (TTree*)file->Get("photon");
/////  m_Rd = new TreeReader(t);

  int nentries = m_Rd->GetEntries();  
  //  nentries = 50;

  std::vector<int> all_photons_indices;
  std::vector<int> HSphotons_indices;
  std::map<int,int> stable_to_mother_association;    

  //  return;

  for(int i = 0; i < nentries; i++) {

    m_Rd->GetEntry(i);

    if(i%100==0)    std::cout<<"========== entry #"<<i<<" "<<double(i)/double(nentries)*100. <<"% ==========="<<std::endl;
    std::cout<<"========== entry #"<<i<<" ==========="<<std::endl;

    for(int p = 0; p < m_Rd->GetVariable("@mc_pt.size()"); p++) {
      if( m_Rd->GetVariable(Form("mc_pt[%d]",p)) < 20000 ) continue;
      if( m_Rd->GetVariable(Form("mc_status[%d]",p)) != 1 ) continue;
      if( m_Rd->GetVariable(Form("mc_pdgId[%d]",p)) != 22 ) continue;
      if( fabs(m_Rd->GetVariable(Form("mc_eta[%d]",p))) > 2.7 ) continue;
      all_photons_indices.push_back(p);
      //      std::cout<<"Found photon w/ index "<<p<<std::endl;
    }


    if(all_photons_indices.size() < 1) continue;
    //    PrintEventInfo(all_photons_indices[all_photons_indices.size()-1]);    

    //    HSphotons_indices = HardScatterPhotonIndices(all_photons_indices);
    for(int ii = 0; ii<(int)all_photons_indices.size(); ii++) {
      stable_to_mother_association.insert(std::pair<int,int>(all_photons_indices[ii],GetOriginalPhotonIndex(all_photons_indices[ii])));
      //      std::cout<<"stable_to_mother_association["<<all_photons_indices[ii] <<"] = " << stable_to_mother_association[all_photons_indices[ii]] << std::endl;
    }



    std::multimap<float,int> pt_index;    
//    for(int j = 0; j<(int)HSphotons_indices.size(); j++) 
//      pt_index.insert(std::pair<float,int>(m_Rd->GetVariable(Form("mc_pt[%d]",HSphotons_indices[j])),HSphotons_indices[j]));
    for(int j = 0; j<(int)all_photons_indices.size(); j++) 
      pt_index.insert(std::pair<float,int>(m_Rd->GetVariable(Form("mc_pt[%d]",all_photons_indices[j])),all_photons_indices[j]));


    std::vector<int> v_index;
    for( std::multimap<float,int>::reverse_iterator ii=pt_index.rbegin(); ii!=pt_index.rend() && (*ii).first != 0; ++ii)
      v_index.push_back((*ii).second);
    
    

    for(unsigned int v = 0; v < v_index.size(); v++) {

      if(v>2) break;
//      std::cout<<v_index[v]<<" is originating from "<< stable_to_mother_association[v_index[v]]<<std::endl;

      bool HardProcPhoton  = IsPhotonFromHardProc (v_index[v]);
      bool RadiationPhoton = IsPhotonFromRadiation(v_index[v]);
      bool FragPhoton      = IsPhotonFromFrag     (v_index[v]);
      bool UEPhoton        = IsPhotonFromUE       (v_index[v]);

      if( HardProcPhoton + RadiationPhoton + FragPhoton + UEPhoton != 1 ) 
	std::cout<<"XXXXXXXXXXXXXXXXX photon is badly flagged "<<std::endl;

      if(FragPhoton)  x[v] = ComputeX(v_index[v]);

      isfromhardproc[v]  = HardProcPhoton;
      isfromfrag[v]      = FragPhoton;
      isfromradiation[v] = RadiationPhoton;
      isfromue[v]        = UEPhoton;

      pt[v] = m_Rd->GetVariable(Form("mc_pt[%d]",v_index[v]));
      eta[v] = m_Rd->GetVariable(Form("mc_eta[%d]",v_index[v]));
      phi[v] = m_Rd->GetVariable(Form("mc_phi[%d]",v_index[v]));

      iso_part[v] = GetPartonETIso(v_index[v], 0.4, 0);

    }


    double truth_PV_ID = m_Rd->GetVariable("mc_vx_z[0]"); 
    TLorentzVector Lead_truth_lv, SubLead_truth_lv, gamgam_truth_lv;
    
    Lead_truth_lv.SetPtEtaPhiM   ( pt[0], InvMassCalc::EtaS1PVCorrected(eta[0], truth_PV_ID), phi[0],0 );
    SubLead_truth_lv.SetPtEtaPhiM( pt[1], InvMassCalc::EtaS1PVCorrected(eta[1], truth_PV_ID), phi[1],0 );

    gamgam_truth_lv = Lead_truth_lv + SubLead_truth_lv;
    mgg = gamgam_truth_lv.M();


    out_tree->Fill();

    all_photons_indices.clear();
    HSphotons_indices.clear();
    v_index.clear();
    pt_index.clear();
    stable_to_mother_association.clear();
    
    for(int z = 0; z < 3; z++) {
      isfromhardproc[z] = -999;
      isfromfrag[z] = -999;
      isfromradiation[z] = -999;
      isfromue[z] = -999;
      x[z] = -999;
      pt[z] = -999;
      eta[z] = -999;
      phi[z] = -999;
      iso_part[z] = -999;
    }
    mgg = 0;


  }



//  fout->cd();
//  out_tree->Write();
//  fout->Close();
  
  
}














/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////













///////////////////////////////
void EventSelector::EventLoop()
///////////////////////////////
{


  //============== Declare Cutflow variables ===============//
  const int nchecks = 10;
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
    AnalysisTools::Processing(jentry,m_nentries);

    bool truth_info = EventTruthSelector(jentry,m_data);
    if( !m_data && !truth_info) continue; //--> Skip MC events without truth matching

    //--> Compute run,lb,evt,npv
    ComputeBasicEventInfo();
    
    //======================= PU Weight for MC ONLY ===============================//
    if(m_data) m_PUweight = 1;
    else
      m_PUweight = m_pileupTool->GetCombinedWeight( (int)m_Rd->GetVariable("RunNumber"),
						    (int)m_Rd->GetVariable("mc_channel_number"),
						    m_Rd->GetVariable("averageIntPerXing") );
    //___________________________________________________________________________//

    //======================= Generator Weight for MC ONLY =======================//
    if(m_data) gen_weight = 1;
    else gen_weight = Commons::GetMCWeight( (int)m_Rd->GetVariable("mc_channel_number") );
    //___________________________________________________________________________//
 

    //================================ Diphotons selection ====================================//
    n[0]+=1;n_w[0]+=m_PUweight; 

    // ------------------- TRIGGER -----------------
    if(!EventTrigOK() ) continue;//--> Trigger
    n[1]+=1;n_w[1]+=m_PUweight;
    // ------------------  GRL -------------------
    if ( !EventGRLOK() )      continue;//--> GRL 
    n[2]+=1;n_w[2]+=m_PUweight;   
    //----------------- LARERROR CUT -------------------------------------
    if((int)m_Rd->GetVariable("larError") == 2 ) continue;//--> larError
    //----------------- TILEERROR CUT -------------------------------------
    if((int)m_Rd->GetVariable("tileError") == 2 ) continue;//--> tileError
    //------------------ EVENT COMPLETED -------------------------------------------------
    if( !EventCompletedOK() ) continue;//--> Event completed
    n[3]+=1;n_w[3]+=weight;

    //------------------- PV -----------------
    if( !EventPVOK() )    continue;//--> PV 
    n[4]+=1;n_w[4]+=m_PUweight;
    //----------------------- PRESELECTION CUTS  --------------------------------------
    int ilead=-99999; int isublead=-99999;
    double ptlead=-99999; double ptsublead=-99999;
    //    bool passPreSel = EventPreSelectionOK(&ilead,&isublead,&ptlead,&ptsublead);
    //    if( !passPreSel ) continue;//--> Preselection

    std::vector<int> indices = EventPreSelectionOK();
    if(indices.size()<2) continue;
    ilead=indices[0];
    isublead=indices[1];

    //--> Compute MC event weight
    weight = m_PUweight;
    Lead_weight    = m_PUweight;
    SubLead_weight = m_PUweight;
    //--> Compute Kinematics variables
    ComputeKinematics(ilead,isublead);
    n[5]+=1;n_w[5]+=weight;
    hmgg[5]->Fill(mgg);hmgg_w[5]->Fill(mgg,weight);


//    if(Lead_pT<SubLead_pT) {
//      float E1=(float)m_Rd->GetVariable( Form("ph_cl_E[%d]",ilead) );
//      float etas21=(float)m_Rd->GetVariable( Form("ph_etas2[%d]",ilead) );
//      float pt1=E1/TMath::CosH(etas21);
//
//      float E2=(float)m_Rd->GetVariable( Form("ph_cl_E[%d]",isublead) );
//      float etas22=(float)m_Rd->GetVariable( Form("ph_etas2[%d]",isublead) );
//      float pt2=E2/TMath::CosH(etas22);
//
//      std::cout<<"--------------------------------------------"<<std::endl;
//      std::cout<<"AFTER PTCUT : Lead_pT "<<Lead_pT<<" SubLead_pT "<<SubLead_pT << std::endl;
//      std::cout<<"ph_cl_pt : "<<pt1/1000. <<" "<< pt2/1000.<<std::endl;
//      std::cout<<"photon smearing factors : "<<GravitonAnalysis::PhotonSmearingFactor(ilead,0)<<" "<<GravitonAnalysis::PhotonSmearingFactor(isublead,0)<<std::endl;
//    }


    //----------------- PT CUT ----------------------------
    if ( !PhotonPtOK(ilead, 50) ) continue;//--> lead pt cut 
    if ( !PhotonPtOK(isublead, 50) ) continue;//--> sublead pt cut 
    n[6]+=1;n_w[6]+=weight;
    hmgg[6]->Fill(mgg);hmgg_w[6]->Fill(mgg,weight);

    //--> Compute Identification decision variables
    ComputeIdDecisions(ilead,isublead);
    //--> Compute Isolation variables
    ComputeIsolations(ilead,isublead);


    //--> Fill output tree if the event pass the cleaning and overlap removal criteria
    if( (int)m_Rd->GetVariable("larError") != 2 && 
	(int)m_Rd->GetVariable("tileError") != 2 && 
	EventCompletedOK() &&
	!EventInZprimme(m_ee_events,m_RunNumber,m_EventNumber) )
      tree_out->Fill();
    //------------------

    //----------------- TIGHT CUT ------------------------
    if( !PhotonIsTightOK(ilead) )    continue;//--> Tight lead
    if( !PhotonIsTightOK(isublead) ) continue;//--> Tight sublead
    n[7]+=1;n_w[7]+=weight;
    hmgg[7]->Fill(mgg);hmgg_w[7]->Fill(mgg,weight);
  
  //----------------- ISOLATION CUT --------------------------------
    if( !PhotonIsolation_toolOK(ilead, 8) )    continue;//--> iso lead
    if( !PhotonIsolation_toolOK(isublead, 8) ) continue;//--> iso sublead
    n[8]+=1;n_w[8]+=weight;
    hmgg[8]->Fill(mgg);hmgg_w[8]->Fill(mgg,weight);

    //------------------ OVERLAP REMOVAL -------------------------------------------------
    // if( EventInZprimme(m_ee_events,m_RunNumber,m_EventNumber) ) continue;//--> Gee removal
    // n[9]+=1;n_w[9]+=weight;
    // hmgg[9]->Fill(mgg);hmgg_w[9]->Fill(mgg,weight);
    // hmgg_final->Fill(mgg);hmgg_final_w->Fill(mgg,weight);
    //--------------------------------------------

    if( mgg < 179.1034055) continue;
    n[9]+=1;n_w[9]+=weight;
    hmgg[9]->Fill(mgg);hmgg_w[9]->Fill(mgg,weight);
    hmgg_final->Fill(mgg);hmgg_final_w->Fill(mgg,weight);

  }
  //========================================================================//
  //========== End of the loop over all entries ============================//
  //========================================================================//

  //-----------------------------------------------------
  std::cout << "Initial         = " << (int)n[0] << std::endl;
  std::cout << "Trigger         = " << (int)n[1] << std::endl;
  std::cout << "GRL             = " << (int)n[2] << std::endl;
  std::cout << "LAr             = " << (int)n[3] << std::endl;
  std::cout << "PV              = " << (int)n[4] << std::endl;
  std::cout << "Preselec        = " << (int)n[5] << std::endl;
  std::cout << "Pt cut          = " << (int)n[6] << std::endl;
  std::cout << "Tight           = " << (int)n[7] << std::endl;
  std::cout << "Iso             = " << (int)n[8] << std::endl;
  std::cout << "mgg             = " << (int)n[9] << std::endl;
  if(!m_data){
    std::cout << "-----Weighted------"  << std::endl;
    std::cout << "Initial         = " << (int)n_w[0] << std::endl;
    std::cout << "Trigger         = " << (int)n_w[1] << std::endl;
    std::cout << "GRL             = " << (int)n_w[2] << std::endl;
    std::cout << "LAr             = " << (int)n_w[3] << std::endl;
    std::cout << "PV              = " << (int)n_w[4] << std::endl;
    std::cout << "Preselec        = " << (int)n_w[5] << std::endl;
    std::cout << "Pt cut          = " << (int)n_w[6] << std::endl;
    std::cout << "Tight           = " << (int)n_w[7] << std::endl;
    std::cout << "Iso             = " << (int)n_w[8] << std::endl;
    std::cout << "mgg             = " << (int)n_w[9] << std::endl;
  }


  //------------------------------
  for(int ic=0;ic<nchecks;ic++){
    hCutFlow->SetBinContent(ic+1,n[ic]);
    hCutFlow_w->SetBinContent(ic+1,n_w[ic]);
    hist_mgg->Add(hmgg[ic]);
    hist_mgg_w->Add(hmgg_w[ic]);
  }
  //---------------------------------


  FillOutputFile();
}

////////////////////////////////////////
void EventSelector::FillOutputFile()
////////////////////////////////////////
{
  // TFile file_out(output,"RECREATE");
  file_out->cd();
  hist_mgg->Write("mgg_cutflow",1);
  hist_mgg_w->Write("mgg_cutflow_w",1);
  hmgg_final->Write();
  hmgg_final_w->Write();
  hCutFlow->Write();
  hCutFlow_w->Write();
  tree_out->Write();
  file_out->Close();
}




///////////////////////////////////////////////
void EventSelector::ComputeBasicEventInfo()
//////////////////////////////////////////////
{
  m_RunNumber   = (int)m_Rd->GetVariable("RunNumber");
  m_EventNumber = (int)m_Rd->GetVariable("EventNumber");
  m_LumiBlock   = (int)m_Rd->GetVariable("lbn");
  m_NPV         = (int)m_Rd->GetVariable("@PV_nTracks.size()");
  m_mu          = m_Rd->GetVariable("averageIntPerXing");
  PV_ID         = m_Rd->GetVariable("PV_z[0]");

}
//////////////////////////////////////////////////////////////
void EventSelector::ComputeKinematics(int ilead, int isublead)
///////////////////////////////////////////////////////////////
{



  m_ph_cl_E.clear();
  m_ph_cl_pt.clear();
  m_ph_cl_eta.clear();
  m_ph_cl_phi.clear();
  m_ph_rawcl_E.clear();
  m_ph_rawcl_pt.clear();
  m_ph_rawcl_eta.clear();
  m_ph_rawcl_phi.clear();
  m_ph_ED_median.clear();



  if(m_data){
    Lead_pT = PhotonRescaledPt_geo21(ilead);
    SubLead_pT = PhotonRescaledPt_geo21(isublead);
  }else{// Apply smearing for MC
    Lead_pT    = PhotonSmearedPt(ilead);
    SubLead_pT = PhotonSmearedPt(isublead);
    Lead_pT_smearedup    = PhotonSmearedPt(ilead,2);
    SubLead_pT_smearedup = PhotonSmearedPt(isublead,2);
    Lead_pT_smeareddown    = PhotonSmearedPt(ilead,1);
    SubLead_pT_smeareddown = PhotonSmearedPt(isublead,1);
  }


  m_MET_Topo_et			  = m_Rd->GetVariable("MET_Topo_et")/1000.;                  
  m_MET_Topo_sumet		  = m_Rd->GetVariable("MET_Topo_sumet")/1000.;               
  m_MET_Topo_sumet_EMB		  = m_Rd->GetVariable("MET_Topo_sumet_EMB")/1000.;           
  m_MET_Topo_sumet_EME		  = m_Rd->GetVariable("MET_Topo_sumet_EME")/1000.;           
  m_MET_Topo_sumet_CentralReg	  = m_Rd->GetVariable("MET_Topo_sumet_CentralReg")/1000.;    
  m_MET_Topo_sumet_EndcapRegion	  = m_Rd->GetVariable("MET_Topo_sumet_EndcapRegion")/1000.;     
  m_MET_RefFinal_et		  = m_Rd->GetVariable("MET_RefFinal_et")/1000.;              
  m_MET_RefFinal_sumet		  = m_Rd->GetVariable("MET_RefFinal_sumet")/1000.;           
  m_MET_RefFinal_sumet_CentralReg = m_Rd->GetVariable("MET_RefFinal_sumet_CentralReg")/1000.;
  m_MET_RefFinal_sumet_EndcapRegion  = m_Rd->GetVariable("MET_RefFinal_sumet_EndcapRegion")/1000.; 




  Lead_eta       = m_Rd->GetVariable( Form("ph_etas2[%d]",ilead) );
  Lead_phi       = m_Rd->GetVariable( Form("ph_phi[%d]",ilead) );
  Lead_IsConv    = (int)m_Rd->GetVariable( Form("ph_isConv[%d]",ilead) );

  SubLead_eta    = m_Rd->GetVariable( Form("ph_etas2[%d]",isublead) );
  SubLead_phi    = m_Rd->GetVariable( Form("ph_phi[%d]",isublead) );
  SubLead_IsConv = (int)m_Rd->GetVariable( Form("ph_isConv[%d]",isublead) );

  for(int i = 0; i < m_Rd->GetVariable("@ph_rawcl_E.size()"); i++) {

    m_ph_cl_E.push_back(m_Rd->GetVariable(Form("ph_cl_E[%d]",i))/1000.);
    m_ph_cl_pt.push_back(m_Rd->GetVariable(Form("ph_cl_pt[%d]",i))/1000.);
    m_ph_cl_eta.push_back(m_Rd->GetVariable(Form("ph_cl_eta[%d]",i)));
    m_ph_cl_phi.push_back(m_Rd->GetVariable(Form("ph_cl_phi[%d]",i)));
    m_ph_rawcl_E.push_back(m_Rd->GetVariable(Form("ph_rawcl_E[%d]",i))/1000.);
    m_ph_rawcl_pt.push_back(m_Rd->GetVariable(Form("ph_rawcl_pt[%d]",i))/1000.);
    m_ph_rawcl_eta.push_back(m_Rd->GetVariable(Form("ph_rawcl_eta[%d]",i)));
    m_ph_rawcl_phi.push_back(m_Rd->GetVariable(Form("ph_rawcl_phi[%d]",i)));
    m_ph_ED_median.push_back(m_Rd->GetVariable(Form("ph_ED_median[%d]",i))/1000.);

  }






  
  double Lead_etaS1    = m_Rd->GetVariable( Form("ph_etas1[%d]",ilead) );
  double SubLead_etaS1 = m_Rd->GetVariable( Form("ph_etas1[%d]",isublead) );
  //  double PV_ID         = m_Rd->GetVariable("PV_z[0]");
  PV_ID         = m_Rd->GetVariable("PV_z[0]");
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
  deltaphi_smearedup     = Lead_lv_smearedup.DeltaPhi(SubLead_lv_smearedup);
  costhetastar_smearedup = CosThetaStar_CS(Lead_lv_smearedup,SubLead_lv_smearedup);
  //-------------------------- Smeared Down ---------------------------------------//
  Lead_lv_smeareddown.SetPtEtaPhiM( Lead_pT_smeareddown,Lead_eta_PV,Lead_phi,0 );
  SubLead_lv_smeareddown.SetPtEtaPhiM( SubLead_pT_smeareddown,SubLead_eta_PV,SubLead_phi,0 );
  gamgam_lv_smeareddown    = Lead_lv_smeareddown + SubLead_lv_smeareddown;
  mgg_smeareddown          = gamgam_lv_smeareddown.M();
  ptgg_smeareddown         = gamgam_lv_smeareddown.Pt();
  ygg_smeareddown          = gamgam_lv_smeareddown.Rapidity();
  deltaphi_smeareddown     = Lead_lv_smeareddown.DeltaPhi(SubLead_lv_smeareddown);
  costhetastar_smeareddown = CosThetaStar_CS(Lead_lv_smeareddown,SubLead_lv_smeareddown);
  //-------------------------------------------------------------------------------------//  

}
/////////////////////////////////////////////////////////////
void EventSelector::ComputeIdDecisions(int ilead,int isublead) 
/////////////////////////////////////////////////////////////
{
  //----> LoosePrime computation
  Lead_IsLoosePrime2    = (int)( (PhotonIsEM(ilead)&0x26fc01) == 0 ? 1:0 );
  Lead_IsLoosePrime3    = (int)( (PhotonIsEM(ilead)&0x24fc01) == 0 ? 1:0 );
  Lead_IsLoosePrime4    = (int)( (PhotonIsEM(ilead)&0x04fc01) == 0 ? 1:0 );
  Lead_IsLoosePrime5    = (int)( (PhotonIsEM(ilead)&0x00fc01) == 0 ? 1:0 );
  SubLead_IsLoosePrime2 = (int)( (PhotonIsEM(isublead)&0x26fc01) == 0 ? 1:0 );
  SubLead_IsLoosePrime3 = (int)( (PhotonIsEM(isublead)&0x24fc01) == 0 ? 1:0 );
  SubLead_IsLoosePrime4 = (int)( (PhotonIsEM(isublead)&0x04fc01) == 0 ? 1:0 );
  SubLead_IsLoosePrime5 = (int)( (PhotonIsEM(isublead)&0x00fc01) == 0 ? 1:0 );
  //----> Tight computation
  Lead_IsTight    = (int)PhotonIsTightOK(ilead);
  SubLead_IsTight = (int)PhotonIsTightOK(isublead);
}
/////////////////////////////////////////////////////////////
void EventSelector::ComputeIsolations(int ilead,int isublead) 
/////////////////////////////////////////////////////////////
{
  Lead_ed_med = m_Rd->GetVariable(Form("ph_ED_median[%d]", ilead))/ 1000. ;
  SubLead_ed_med = m_Rd->GetVariable(Form("ph_ED_median[%d]", isublead))/ 1000. ;
  Lead_iso_uncor = m_Rd->GetVariable(Form("ph_topoEtcone40[%d]", ilead)) / 1000.;
  SubLead_iso_uncor = m_Rd->GetVariable(Form("ph_topoEtcone40[%d]", isublead)) / 1000.;

  float radius = 40.;
  bool isMC = (!m_data);
  Lead_iso_ptcorr = m_Rd->GetVariable(Form("ph_topoEtcone40[%d]", ilead));
  Lead_iso_ptcorr -= CaloIsoCorrection::GetPtCorrectedTopoIsolation(m_Rd->GetVariable(Form("ph_cl_E[%d]", ilead)),
								    m_Rd->GetVariable(Form("ph_etas2[%d]", ilead)),
								   m_Rd->GetVariable(Form("ph_etap[%d]", ilead)),
								   m_Rd->GetVariable(Form("ph_cl_eta[%d]", ilead)),
								   radius,
								   isMC,
								   (bool)m_Rd->GetVariable(Form("ph_isConv[%d]", ilead)),
								   CaloIsoCorrection::PHOTON);
  Lead_iso_ptcorr /= 1000.;

  SubLead_iso_ptcorr = m_Rd->GetVariable(Form("ph_topoEtcone40[%d]", isublead));
  SubLead_iso_ptcorr -= CaloIsoCorrection::GetPtCorrectedTopoIsolation(m_Rd->GetVariable(Form("ph_cl_E[%d]", isublead)),
								      m_Rd->GetVariable(Form("ph_etas2[%d]", isublead)),
								      m_Rd->GetVariable(Form("ph_etap[%d]", isublead)),
								      m_Rd->GetVariable(Form("ph_cl_eta[%d]", isublead)),
								      radius,
								      isMC,
								      (bool)m_Rd->GetVariable(Form("ph_isConv[%d]", isublead)),
								      CaloIsoCorrection::PHOTON);
  SubLead_iso_ptcorr /= 1000.;


  if(m_isotype == "CONE"){
    Lead_iso    = PhotonIsolation_tool(ilead)/1000.;
    SubLead_iso = PhotonIsolation_tool(isublead)/1000.;
  }else if(m_isotype == "TOPO"){
    Lead_iso    = PhotonTopoIsolation_tool_geo21(ilead)/1000.;
    SubLead_iso = PhotonTopoIsolation_tool_geo21(isublead)/1000.;
  }else if(m_isotype == "TOPO_PTDEP"){
    Lead_iso        = PhotonTopoIsolation_tool_geo21(ilead)/1000.;
    SubLead_iso     = PhotonTopoIsolation_tool_geo21(isublead)/1000.;
    Lead_iso_mod    = PhotonTopoIsolationPtdep_tool(ilead)/1000.;
    SubLead_iso_mod = PhotonTopoIsolationPtdep_tool(isublead)/1000.;
  }else Fatal("EventSelector::ComputeIsolations()","Wrong isolation type !!!");  
}

/////////////////////////////////////////////////////////////////////////////
double EventSelector::CosThetaStar_CS(TLorentzVector v1,TLorentzVector v2)
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
///////////////////////////////////////////////////////////////
bool EventSelector::EventTruthSelector(int entry,bool data)
//////////////////////////////////////////////////////////////
{

  if(data) return false;
  
  m_Ts->GetEntry(entry);

  std::vector<int> index_truth;
  if( m_mctype == "pythia_gg" )
    index_truth = m_Ts->GetDirectGamGamChildIndex(); 
  else if( m_mctype == "pythia_gj" )
    index_truth = m_Ts->GetGamJetChildIndex();
  else if( m_mctype == "pythia_rs" )
    index_truth = m_Ts->GetRSGravChildIndex();
  else Fatal("EventSelector::EventTruthSelector()", "Wrong MC sample type !");

  //jb, 08/10/2013, removed this part to not reject events with number of mc photons from HS != 2
//  if( index_truth.size() !=2 ) return false;
//
//  double pt1 = m_Rd->GetVariable(Form("mc_pt[%d]",index_truth[0]));
//  double pt2 = m_Rd->GetVariable(Form("mc_pt[%d]",index_truth[1]));
//  int ilead_truth =-1; int isublead_truth=-1;
//
//  if(pt1>pt2) {
//    ilead_truth = index_truth[0];
//    isublead_truth = index_truth[1];
//  }else{
//    ilead_truth = index_truth[1];
//    isublead_truth = index_truth[0];
//  }



  int ilead_truth =-1; int isublead_truth=-1;
  
  
  if(index_truth.size() >= 2) { 
    // more than 2 photons 
    std::multimap<float,int> pt_index;    
    for(int i = 0; i<(int)index_truth.size(); i++) 
      pt_index.insert(std::pair<float,int>(m_Rd->GetVariable(Form("mc_pt[%d]",index_truth[i]))/1000.,index_truth[i]));

    std::vector<int> v_index;
    for( std::multimap<float,int>::reverse_iterator ii=pt_index.rbegin(); ii!=pt_index.rend() && (*ii).first != 0; ++ii) {
      //      std::cout << "pt= "<<(*ii).first << " index= "<<(*ii).second<<std::endl;
      v_index.push_back((*ii).second);
    }
    //    std::cout<<std::endl;

    //fill variables with indices v_index[0] and [1]
    ilead_truth = v_index[0];
    isublead_truth = v_index[1];

    
    Lead_truth_eta     = m_Rd->GetVariable(Form("mc_eta[%d]",ilead_truth));
    Lead_truth_phi     = m_Rd->GetVariable(Form("mc_phi[%d]",ilead_truth));
    Lead_truth_pT      = m_Rd->GetVariable(Form("mc_pt[%d]",ilead_truth))/1000.;
    Lead_truth_status  = (int)m_Rd->GetVariable(Form("mc_status[%d]",ilead_truth));
    
    SubLead_truth_eta     = m_Rd->GetVariable(Form("mc_eta[%d]",isublead_truth));
    SubLead_truth_phi     = m_Rd->GetVariable(Form("mc_phi[%d]",isublead_truth));
    SubLead_truth_pT      = m_Rd->GetVariable(Form("mc_pt[%d]",isublead_truth))/1000.;
    SubLead_truth_status  = (int)m_Rd->GetVariable(Form("mc_status[%d]",isublead_truth));
    
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
    truth_costhetastar = CosThetaStar_CS(Lead_truth_lv,SubLead_truth_lv);
    
  }
  
  
  else if(index_truth.size() == 0) {
    Lead_truth_eta    = Lead_truth_phi    = Lead_truth_pT    = Lead_truth_status    =
    SubLead_truth_eta = SubLead_truth_phi = SubLead_truth_pT = SubLead_truth_status =
    truth_mgg         = truth_ptgg        = truth_ygg        = truth_dphigg         = truth_costhetastar = -999;
  }
  
  else if(index_truth.size() == 1) { 

    ilead_truth = index_truth[0];    
    
    Lead_truth_eta     = m_Rd->GetVariable(Form("mc_eta[%d]",ilead_truth));
    Lead_truth_phi     = m_Rd->GetVariable(Form("mc_phi[%d]",ilead_truth));
    Lead_truth_pT      = m_Rd->GetVariable(Form("mc_pt[%d]",ilead_truth))/1000.;
    Lead_truth_status  = (int)m_Rd->GetVariable(Form("mc_status[%d]",ilead_truth));
    
    SubLead_truth_eta = SubLead_truth_phi = SubLead_truth_pT = SubLead_truth_status =
    truth_mgg         = truth_ptgg        = truth_ygg        = truth_dphigg         = truth_costhetastar = -999;

  }

  
  if( m_mctype == "pythia_rs" ){
    unsigned int mcsize=(unsigned int)m_Rd->GetVariable( "@mc_m.size()");
    for(unsigned int imc=0 ; imc<mcsize ;imc++ ){
      int pdgId=(int)m_Rd->GetVariable( Form("mc_pdgId[%d]",imc) );
      if(pdgId==5100039){
	mymc_m=m_Rd->GetVariable(Form("mc_m[%d]",imc))/1000.; 
	break;
      }
    }
  }

  
  return true;
}
