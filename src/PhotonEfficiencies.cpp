#define PhotonEfficiencies_cxx
#include "PhotonEfficiencies.h"
#include "ToolsCommons.h"
#include <TAxis.h>
#include <TMath.h>

////////////////////////////////////////////////////////////////////////////
PhotonEfficiencies::PhotonEfficiencies(TTree* tree) : GravitonAnalysis(tree)
////////////////////////////////////////////////////////////////////////////
{
   if (tree == 0) {
     std::cout << "You need to specify an input tree !!!"<<std::endl;
   }
   Init();

}
//////////////////////////////////////////////////////////////
PhotonEfficiencies::PhotonEfficiencies() : GravitonAnalysis()
//////////////////////////////////////////////////////////////
{
  //Default constructor
}
/////////////////////////////////////////
PhotonEfficiencies::~PhotonEfficiencies()
/////////////////////////////////////////
{
  //destructor
  delete m_map_lead_truth;
  delete m_map_lead_reco;
  delete m_map_lead_recotight;
  delete m_map_lead_recotightiso;
  delete m_map_lead_reco_isoapplied;
  delete m_map_lead_recotight_isoapplied;
  delete m_map_sublead_truth;
  delete m_map_sublead_reco;
  delete m_map_sublead_recotight;
  delete m_map_sublead_recotightiso;
  delete m_map_sublead_reco_isoapplied;
  delete m_map_sublead_recotight_isoapplied;

}

//////////////////////////////
void PhotonEfficiencies::Init()
//////////////////////////////
{
  SetEtaBins();
  SetPtBins();
  m_map_lead_truth = new TH2F( "map_lead_truth", "map_lead_truth",
			       (int)m_bins_pt.size()-1, &m_bins_pt[0],
			       (int)m_bins_eta.size()-1, &m_bins_eta[0] );
  m_map_lead_reco = new TH2F( "map_lead_reco", "map_lead_reco",
			      (int)m_bins_pt.size()-1, &m_bins_pt[0] ,
			      (int)m_bins_eta.size()-1, &m_bins_eta[0] );
  m_map_lead_recotight = new TH2F( "map_lead_recotight", "map_lead_recotight",
				   (int)m_bins_pt.size()-1,  &m_bins_pt[0],
				   (int)m_bins_eta.size()-1, &m_bins_eta[0] );
  m_map_lead_recotightiso = new TH2F( "map_lead_recotightiso", "map_lead_recotightiso",
				     (int)m_bins_pt.size()-1,  &m_bins_pt[0],
				      (int)m_bins_eta.size()-1, &m_bins_eta[0] );
  m_map_lead_reco_isoapplied = new TH2F( "map_lead_reco_isoapplied", "map_lead_reco_isoapplied",
					 (int)m_bins_pt.size()-1,&m_bins_pt[0],
					 (int)m_bins_eta.size()-1,&m_bins_eta[0]);
  m_map_lead_recotight_isoapplied = new TH2F( "map_lead_recotight_isoapplied",
					      "map_lead_recotight_isoapplied",
					      (int)m_bins_pt.size()-1,&m_bins_pt[0] ,
					      (int)m_bins_eta.size()-1,&m_bins_eta[0]);
  m_map_sublead_truth = new TH2F( "map_sublead_truth", "map_sublead_truth",
			       (int)m_bins_pt.size()-1, &m_bins_pt[0],
			       (int)m_bins_eta.size()-1, &m_bins_eta[0] );
  m_map_sublead_reco = new TH2F( "map_sublead_reco", "map_sublead_reco",
			      (int)m_bins_pt.size()-1, &m_bins_pt[0] ,
			      (int)m_bins_eta.size()-1, &m_bins_eta[0] );
  m_map_sublead_recotight = new TH2F( "map_sublead_recotight", "map_sublead_recotight",
				   (int)m_bins_pt.size()-1,  &m_bins_pt[0],
				   (int)m_bins_eta.size()-1, &m_bins_eta[0] );
  m_map_sublead_recotightiso = new TH2F( "map_sublead_recotightiso", "map_sublead_recotightiso",
				     (int)m_bins_pt.size()-1,  &m_bins_pt[0],
				      (int)m_bins_eta.size()-1, &m_bins_eta[0] );
  m_map_sublead_reco_isoapplied = new TH2F( "map_sublead_reco_isoapplied", "map_sublead_reco_isoapplied",
					 (int)m_bins_pt.size()-1,&m_bins_pt[0],
					 (int)m_bins_eta.size()-1,&m_bins_eta[0]);
  m_map_sublead_recotight_isoapplied = new TH2F( "map_sublead_recotight_isoapplied",
					      "map_sublead_recotight_isoapplied",
					      (int)m_bins_pt.size()-1,&m_bins_pt[0] ,
					      (int)m_bins_eta.size()-1,&m_bins_eta[0]);


}

/////////////////////////////////////
void PhotonEfficiencies::SetEtaBins()
////////////////////////////////////
{
  // double bins[]= { -2.37,-1.81,-1.52,-1.37,-0.6,
  // 		   0.,0.6,1.37,1.52,1.81,2.37 };
  // m_bins_eta.assign(bins,bins+11);
  double bins[]= { 0.,0.6,1.37,1.52,1.81,2.37 };
  m_bins_eta.assign(bins,bins+6);
}
/////////////////////////////////////
void PhotonEfficiencies::SetPtBins()
////////////////////////////////////
{
  int n = 500;
  for(int i=0;i<n;i++)
    m_bins_pt.push_back(5*i);
}
///////////////////////////////////////////////////////
void PhotonEfficiencies::SetMonteCarloSample(TString s)
///////////////////////////////////////////////////////
{
  if(s=="SMGG") m_mcsample = s;
  else if(s=="SMGJ") m_mcsample = s;
  else if(s=="RSGG") m_mcsample = s;
  else if(s=="ADDGG") m_mcsample = s;
  else Fatal("PhotonEfficiencies::SetMonteCarloSample", "Wrong TString !");
}


///////////////////////////////////////////////////
void PhotonEfficiencies::EventLoop(TString MCTYPE)
//////////////////////////////////////////////////
{
  if(m_data) Fatal("PhotonEfficiencies::EventLoop","This class is not meant to run over data !");

  //==============================================================================//  
  //==================== Start the Loop over all entries =========================//
  //==============================================================================//  
  for(int entry=0;entry<m_nentries;entry++){
    //--> Print the number of events
    AnalysisTools::Processing(entry,m_nentries);
    m_Ts->GetEntry(entry);
    std::vector<int> truth_index;
    if(m_mcsample=="SMGG") truth_index= m_Ts->GetDirectGamGamChildIndex();
    else if(m_mcsample=="SMGJ") truth_index= m_Ts->GetDirectGamGamChildIndex();
    else if(m_mcsample=="RSGG") truth_index= m_Ts->GetRSGravChildIndex();
    else if(m_mcsample=="ADDGG") truth_index= m_Ts->GetDirectGamGamChildIndex();
    else Fatal("PhotonEfficiencies::EventLoop", "Wrong TString !");

    if(truth_index.size() !=2 ) continue;
    // need to rearrange the vector as a function of mc_pt here or in the TruthSelector class;
    m_Rd->GetEntry(entry);
    double gen_weight = 1;//Commons::GetMCWeight( (int)m_Rd->GetVariable("mc_channel_number") );
    double PU_weight = m_pileupTool->GetCombinedWeight( (int)m_Rd->GetVariable("RunNumber"),
							(int)m_Rd->GetVariable("mc_channel_number"),
							m_Rd->GetVariable("averageIntPerXing") );



    //    std::cout<<gen_weight<<" "<<PU_weight<<std::endl;

    m_map_lead_truth->Fill( m_Rd->GetVariable(Form("mc_pt[%d]",truth_index[0]))/1000.,
			    m_Rd->GetVariable(Form("mc_eta[%d]",truth_index[0])),
			    gen_weight );
    m_map_sublead_truth->Fill( m_Rd->GetVariable(Form("mc_pt[%d]",truth_index[1]))/1000.,
			       m_Rd->GetVariable(Form("mc_eta[%d]",truth_index[1])),
			       gen_weight );

    int lead_reco_index = -9999; int sublead_reco_index = -9999;
    for(int iph=0; iph<(int)m_Rd->GetVariable("@ph_truth_index.size()");iph++){
      if( (int)m_Rd->GetVariable(Form("ph_truth_index[%d]",iph)) < 0) continue;
      if( (int)m_Rd->GetVariable(Form("ph_truth_index[%d]",iph)) == truth_index[0])
	lead_reco_index = iph;
      if( (int)m_Rd->GetVariable(Form("ph_truth_index[%d]",iph)) == truth_index[1])
	sublead_reco_index = iph;
    }      

    //--------- LEADING PHOTON ------------------------------------------------------
    if( lead_reco_index>=0 && PhotonOQOK(lead_reco_index)){//reco
      m_map_lead_reco->Fill( m_Rd->GetVariable(Form("mc_pt[%d]",truth_index[0]))/1000.,
			     m_Rd->GetVariable(Form("mc_eta[%d]",truth_index[0])),
			     gen_weight*PU_weight);
      if( PhotonIsTightOK(lead_reco_index) ){//recotight
	m_map_lead_recotight->Fill( m_Rd->GetVariable(Form("mc_pt[%d]",truth_index[0]))/1000.,
				    m_Rd->GetVariable(Form("mc_eta[%d]",truth_index[0])),
				    gen_weight*PU_weight );
	if( PhotonIsolation_toolOK(lead_reco_index) ){//recotightiso
	  m_map_lead_recotightiso->Fill( m_Rd->GetVariable(Form("mc_pt[%d]",truth_index[0]))/1000.,
					 m_Rd->GetVariable(Form("mc_eta[%d]",truth_index[0])),
					 gen_weight*PU_weight );

	}//end of recotightiso
      }//end of recotight
      if( PhotonIsolation_toolOK(lead_reco_index) ){//reco_isoapplied
	m_map_lead_reco_isoapplied->Fill( m_Rd->GetVariable(Form("mc_pt[%d]",truth_index[0]))/1000.,
					  m_Rd->GetVariable(Form("mc_eta[%d]",truth_index[0])),
					  gen_weight*PU_weight );
	if( PhotonIsTightOK(lead_reco_index) ){
	  m_map_lead_recotight_isoapplied->Fill( m_Rd->GetVariable(Form("mc_pt[%d]",truth_index[0]))/1000.,
						 m_Rd->GetVariable(Form("mc_eta[%d]",truth_index[0])),
						 gen_weight*PU_weight );
	}//end of tightreco_isapplied
      }// end of reco_isoapplied
    }//end of reco
    //--------- SUBLEADING PHOTON ------------------------------------------------------
    if( sublead_reco_index>=0 && PhotonOQOK(sublead_reco_index)){//reco
      m_map_sublead_reco->Fill( m_Rd->GetVariable(Form("mc_pt[%d]",truth_index[1]))/1000.,
				m_Rd->GetVariable(Form("mc_eta[%d]",truth_index[1])),
				gen_weight*PU_weight );
      if( PhotonIsTightOK(sublead_reco_index) ){//recotight
	m_map_sublead_recotight->Fill( m_Rd->GetVariable(Form("mc_pt[%d]",truth_index[1]))/1000.,
				       m_Rd->GetVariable(Form("mc_eta[%d]",truth_index[1])),
				       gen_weight*PU_weight );
	
	if( PhotonIsolation_toolOK(sublead_reco_index) ){//recotightiso
	  m_map_sublead_recotightiso->Fill( m_Rd->GetVariable(Form("mc_pt[%d]",truth_index[1]))/1000.,
					    m_Rd->GetVariable(Form("mc_eta[%d]",truth_index[1])),
					    gen_weight*PU_weight );
	
	}//end of recotightiso
      }//end of recotight
      if( PhotonIsolation_toolOK(sublead_reco_index) ){//reco_isoapplied
	m_map_sublead_reco_isoapplied->Fill( m_Rd->GetVariable(Form("mc_pt[%d]",truth_index[1]))/1000.,
					     m_Rd->GetVariable(Form("mc_eta[%d]",truth_index[1])),
					     gen_weight*PU_weight );
	if( PhotonIsTightOK(sublead_reco_index) ){
	  m_map_sublead_recotight_isoapplied->Fill( m_Rd->GetVariable(Form("mc_pt[%d]",truth_index[1]))/1000.,
						    m_Rd->GetVariable(Form("mc_eta[%d]",truth_index[1])),
						    gen_weight*PU_weight );
	}//end of tightreco_isapplied
      }// end of reco_isoapplied
    }//end of reco

  }//End of entry loop	
  //========================================================================//
  //========== End of the loop over all entries ============================//
  //========================================================================//
}
////////////////////////////////////////////////////
void PhotonEfficiencies::CreateOutputFile(TString st)
///////////////////////////////////////////////////
{
  TFile f(st,"RECREATE");
  f.Add(m_map_lead_truth);
  f.Add(m_map_lead_reco);
  f.Add(m_map_lead_recotight);
  f.Add(m_map_lead_recotightiso);
  f.Add(m_map_lead_reco_isoapplied);
  f.Add(m_map_lead_recotight_isoapplied);
  f.Add(m_map_sublead_truth);
  f.Add(m_map_sublead_reco);
  f.Add(m_map_sublead_recotight);
  f.Add(m_map_sublead_recotightiso);
  f.Add(m_map_sublead_reco_isoapplied);
  f.Add(m_map_sublead_recotight_isoapplied);
  f.Write();
  f.Close();
}
