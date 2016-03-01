#include "SignalTemplateCreator.h"
#include "ToolsUtilities.h"
#include "ToolsCommons.h"
#include <TError.h>
//////////////////////////////////////////////////////////////////////////////////
SignalTemplateCreator::SignalTemplateCreator(TTree* tree) : EventSelector(tree)
//////////////////////////////////////////////////////////////////////////////////
{
   if (tree == 0) {
     std::cout << "You need to specify an input tree !!!"<<std::endl;
   }
   Init(tree);
}
////////////////////////////////////////////////////////////////////
SignalTemplateCreator::SignalTemplateCreator() : EventSelector()
///////////////////////////////////////////////////////////////////
{
  //Default constructor
}
///////////////////////////////////////////////
SignalTemplateCreator::~SignalTemplateCreator()
///////////////////////////////////////////////
{
  //Default destructor
  delete m_template;
  delete m_template_truth;
  delete m_hmgg_array;
  delete m_hmgg_w_array;
  delete m_cutflow_array;
  delete m_cutflow_w_array;
  delete m_cutflowsum2_w_array;
}
//////////////////////////////////////////////
void SignalTemplateCreator::Init(TTree* tree)
/////////////////////////////////////////////
{

  m_data         = false;
  //--> In GeV : Set mass = LowMass + i*MassSpacing (0<i<Nmasses)
  SetNmasses(25);
  SetLowMass(500);
  SetMassSpacing(100);
  SetCoupling(0.1);  
  InitListOfMasses();

  SetLeadPtCut(40.);
  SetSubleadPtCut(30.);
  SetLeadIsoCut(5.);
  SetSubleadIsoCut(5.);
}
//////////////////////////////////////////////
void SignalTemplateCreator::InitListOfMasses()
//////////////////////////////////////////////
{
  //---- PREVIOUS VERSION TO KEEP THE DEPENDENCY (DO NOT UNCOMMENT)
  // ///CAUTION IN TEV
  // int Nmass = 72;

  // for(int iMass=0;iMass<Nmass;iMass++){   
  //   // double mass =  0.45+ (iMass*0.025);
  //   double mass =  0.04 * iMass + 0.13;
  //   m_vec_mass.push_back(mass);
  // }
  ///CAUTION IN TEV
  m_vec_mass.clear();
  for(int iMass=0;iMass<m_Nmasses;iMass++)   
    m_vec_mass.push_back(m_MassSpacing * iMass + m_LowMass);

}
/////////////////////////////////////////////////////////////////////////////////////////////////
double SignalTemplateCreator::GetGravitonWeight(double truemass,double polemass,double coupling)
////////////////////////////////////////////////////////////////////////////////////////////////
{
  //------------------------------------------------------------------------------------------------
  //Implementation from Daniel Hayden, See:
  //https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HowToWeightMC2012#RS_Graviton_Flat_Sample_Reweight
  //------------------------------------------------------------------------------------------------


  // Desired Graviton Mass (GeV)
  Double_t mRes = polemass;
  // Desired Graviton Width (GeV) given k/Mpl Coupling i.e. 0.01, 0.03, 0.05, 0.10.
  Double_t wRes = 1.43813*mRes*coupling*coupling;

  Double_t m2Res = mRes*mRes;
  Double_t GamMRat = wRes/mRes;
  Double_t sHat = (truemass)*(truemass);
  Double_t weightBW = 1.0/(pow(sHat - m2Res,2) + pow(sHat * GamMRat,2)); 

  return weightBW;
}

//////////////////////////////////////////////////////////////
void SignalTemplateCreator::CreateOutputFile(TString output)
//////////////////////////////////////////////////////////////
{
  TFile *fout = new TFile(output,"RECREATE");
  m_template         = new TObjArray();
  m_template_truth   = new TObjArray();
  m_hmgg_array       = new TObjArray();
  m_hmgg_w_array       = new TObjArray();
  m_cutflow_array    = new TObjArray();
  m_cutflow_w_array    = new TObjArray();
  m_cutflowsum2_w_array    = new TObjArray();
  for (int imass=0; imass < m_Nmasses; imass++){
    double I = vec_sighist[imass]->Integral();
    vec_sighist[imass]->Scale(1./I);
    m_template->Add(vec_sighist[imass]);

    double I_truth = vec_sighist_truth[imass]->Integral();
    vec_sighist_truth[imass]->Scale(1./I_truth);
    m_template_truth->Add(vec_sighist_truth[imass]);
    m_hmgg_array->Add(vec_hmgg_final[imass]);
    m_hmgg_w_array->Add(vec_hmgg_final[imass]);
    m_cutflow_array->Add(vec_hCutFlow[imass]);
    m_cutflow_w_array->Add(vec_hCutFlow_w[imass]);
    m_cutflowsum2_w_array->Add(vec_hCutFlowSum2_w[imass]);

  }
  m_template->Write("template",1);
  m_template_truth->Write("template_truth",1);
  m_hmgg_array->Write("hmgg_full",1);
  m_hmgg_w_array->Write("hmgg_full_weighted",1);
  m_cutflow_array->Write("cutflow_unweighted",1);
  m_cutflow_w_array->Write("cutflow_weighted",1);
  m_cutflowsum2_w_array->Write("cutflow_sum2_weighted",1);
  fout->Close();
}

///////////////////////////////////////////////////////
void SignalTemplateCreator::EventLoop(double coupling)
///////////////////////////////////////////////////////
{
  SetCoupling(coupling);
  double generator_weight[m_Nmasses];
  vec_hCutFlow        = new TH1D* [m_Nmasses];
  vec_hCutFlow_w      = new TH1D* [m_Nmasses];
  vec_hCutFlowSum2_w  = new TH1D* [m_Nmasses];
  vec_hmgg_final      = new TH1D* [m_Nmasses];
  vec_hmgg_final_w    = new TH1D* [m_Nmasses];
  vec_sighist         = new TH1D* [m_Nmasses];
  vec_sighist_truth   = new TH1D* [m_Nmasses];
    
  //============== Declare output variables ===============//
  int nchecks = 12;
  double n[m_Nmasses][nchecks];
  double n_w[m_Nmasses][nchecks];
  double n2_w[m_Nmasses][nchecks];

  for(int im=0;im<m_Nmasses;im++){

    for(int ic=0 ;ic<nchecks;ic++){
      n[im][ic]=0;n_w[im][ic]=0;n2_w[im][ic]=0;
    }
    TString hname = Form("CutFlow_%d_%1.2f",(int)m_vec_mass[im],m_coupling);
    vec_hCutFlow[im]   = new TH1D(hname,hname,nchecks,0,nchecks);
    hname += "_w";
    vec_hCutFlow_w[im] = new TH1D(hname,hname,nchecks,0,nchecks);
    hname += "_Sum2";
    vec_hCutFlowSum2_w[im] = new TH1D(hname,hname,nchecks,0,nchecks);
    
    //-----------------------------------------------------------------
    TString hmgg_name = Form("hmgg_%d_%1.2f",(int)m_vec_mass[im],m_coupling);
    vec_hmgg_final[im]   = new TH1D(hmgg_name,hmgg_name,Commons::nBins,Commons::binning);
    hmgg_name += "_w";
    vec_hmgg_final_w[im] = new TH1D(hmgg_name,hmgg_name,Commons::nBins,Commons::binning);
    //--------------------------------------------------------------------------------------//
    TString hsearch_name = Form("%d",im);
    TString htruth_name  = Form("truth_%d",im);
    double center = m_vec_mass[im];
    vec_sighist[im] = new TH1D(hsearch_name,hsearch_name,Commons::nBins_search,Commons::binning_search);
    vec_sighist_truth[im] = new TH1D(htruth_name,htruth_name,200,center-100,center+100);
  }
  
  
  
  //==============================================================================//  
  //==================== Start the Loop over all entries =========================//
  //==============================================================================//  
  for(int jentry=0;jentry<m_nentries;jentry++){
    m_Rd->GetEntry(jentry);
    //--> Print the number of events
    AnalysisTools::Processing(jentry,m_nentries,(int)m_nentries/100);
    
    
    //======================= PU Weight for MC ONLY ===============================//
    if(!m_data)
      m_PUweight = m_pileupTool->GetCombinedWeight((int)m_Rd->GetVariable("RunNumber"),
						   (int)m_Rd->GetVariable("mc_channel_number"),
						   m_Rd->GetVariable("averageIntPerXing") );
    //___________________________________________________________________________//
    //--> truth mc_m branch for MC RS G* only 
    unsigned int mcsize=(unsigned int)m_Rd->GetVariable( "@mc_m.size()");
    for(unsigned int imc=0 ; imc<mcsize ;imc++ ){
      int pdgId=(int)m_Rd->GetVariable( Form("mc_pdgId[%d]",imc) );
      if(pdgId==5100039){
	mymc_m=m_Rd->GetVariable(Form("mc_m[%d]",imc))/1000.;	
	break;
      }
    }
    
    for(int imass=0;imass<m_Nmasses;imass++){
      generator_weight[imass] = GetGravitonWeight(mymc_m,m_vec_mass[imass],m_coupling);
      n[imass][0]+=1;n_w[imass][0]+=m_PUweight*generator_weight[imass]; 
      n2_w[imass][0]+= m_PUweight*generator_weight[imass]*m_PUweight*generator_weight[imass]; 
    }
    
    
    //================================ Diphotons selection ====================================//


    //----------------------- TRIGGER CUT  ---------------------------------------
    if(!EventTrigOK() ) continue;//--> Trigger
    for(int imass=0;imass<m_Nmasses;imass++){
      n[imass][1]+=1;
      n_w[imass][1]+=m_PUweight*generator_weight[imass]; 
      n2_w[imass][1]+=m_PUweight*generator_weight[imass]*m_PUweight*generator_weight[imass]; 
    }
    //------------------------- GRL CUT  ---------------------------------------
    if ( !EventGRLOK() )      continue;//--> GRL 
    for(int imass=0;imass<m_Nmasses;imass++){
      n[imass][2]+=1;
      n_w[imass][2] +=m_PUweight*generator_weight[imass]; 
      n2_w[imass][2]+=m_PUweight*generator_weight[imass]*m_PUweight*generator_weight[imass]; 
      
    }
    //------------------------- PRIMARY VERTEX CUT  ---------------------------------------
    if( !EventPVOK() )    continue;//--> PV 
    for(int imass=0;imass<m_Nmasses;imass++){
      n[imass][3]+=1;n_w[imass][3]+=m_PUweight*generator_weight[imass]; 
      n2_w[imass][3]+=m_PUweight*generator_weight[imass]*m_PUweight*generator_weight[imass]; 
    }
    //----------------------- PRESELECTION CUTS  ---------------------------------------
    int ilead=-99999; int isublead=-99999;
    double pt_lead=-99999; double pt_sublead=-99999;
    bool passPreSel = EventPreSelectionOK(&ilead,&isublead,&pt_lead,&pt_sublead);
    if( !passPreSel ) continue;//--> Preselection
    for(int imass=0;imass<m_Nmasses;imass++){
      n[imass][4]+=1;n_w[imass][4]+=m_PUweight*generator_weight[imass]; 
      n2_w[imass][4]+=m_PUweight*generator_weight[imass]*m_PUweight*generator_weight[imass]; 
    }
    ComputeKinematics(ilead,isublead);


    //----------------------------- PT CUT ---------------------------------------------
    if ( !PhotonPtOK(ilead,40) ) continue;//--> lead pt cut 
    if ( !PhotonPtOK(isublead,30) ) continue;//--> sublead pt cut 
    for(int imass=0;imass<m_Nmasses;imass++){
      n[imass][5]+=1;n_w[imass][5]+=m_PUweight*generator_weight[imass]; 
      n2_w[imass][5]+=m_PUweight*generator_weight[imass]*m_PUweight*generator_weight[imass]; 
    }
    //----------------- TIGHT CUT ------------------------
    if( !PhotonIsTightOK(ilead) )    continue;//--> Tight lead
    if( !PhotonIsTightOK(isublead) ) continue;//--> Tight sublead
    for(int imass=0;imass<m_Nmasses;imass++){
      n[imass][6]+=1;n_w[imass][6]+=m_PUweight*generator_weight[imass]; 
      n2_w[imass][6]+=m_PUweight*generator_weight[imass]*m_PUweight*generator_weight[imass]; 
    }
    //----------------- ISOLATION CUT --------------------------------
    bool IsoLeadOK;
    bool IsoSubLeadOK;
    IsoLeadOK    = PhotonIsolation_toolOK(ilead,m_isol_cut);
    IsoSubLeadOK = PhotonIsolation_toolOK(isublead,m_isosl_cut);

    if( !IsoLeadOK )    continue;//--> iso lead
    if( !IsoSubLeadOK ) continue;//--> iso sublead
    for(int imass=0;imass<m_Nmasses;imass++){
      n[imass][7]+=1;n_w[imass][7]+=m_PUweight*generator_weight[imass]; 
      n2_w[imass][7]+=m_PUweight*generator_weight[imass]*m_PUweight*generator_weight[imass]; 
    }
    //----------------- ADDITIONAL PT CUT ------------------------
    if ( !PhotonPtOK(ilead,m_ptl_cut) ) continue;//--> lead pt cut 
    if ( !PhotonPtOK(isublead,m_ptsl_cut) ) continue;//--> sublead pt cut 
    for(int imass=0;imass<m_Nmasses;imass++){
      n[imass][8]+=1;n_w[imass][8]+=m_PUweight*generator_weight[imass]; 
      n2_w[imass][8]+=m_PUweight*generator_weight[imass]*m_PUweight*generator_weight[imass]; 
    }
    //--------------- MGG CUT -------------------------------
    if( mgg < Commons::binning_search[0] ) continue;
    for(int imass=0;imass<m_Nmasses;imass++){
      n[imass][9]+=1;n_w[imass][9]+=m_PUweight*generator_weight[imass]; 
      n2_w[imass][9]+=m_PUweight*generator_weight[imass]*m_PUweight*generator_weight[imass]; 
    }
    //----------------- LARERROR CUT -------------------------------------
    if((int)m_Rd->GetVariable("larError") == 2 ) continue;//--> larError
    //----------------- TILEERROR CUT -------------------------------------
    if((int)m_Rd->GetVariable("tileError") == 2 ) continue;//--> tileError
    //------------------ EVENT COMPLETED -------------------------------------------------
    if( !EventCompletedOK() ) continue;//--> Event completed
    for(int imass=0;imass<m_Nmasses;imass++){
      n[imass][10]+=1;n_w[imass][10]+=m_PUweight*generator_weight[imass]; 
      n2_w[imass][10]+=m_PUweight*generator_weight[imass]*m_PUweight*generator_weight[imass]; 
      vec_hmgg_final[imass]->Fill(mgg,generator_weight[imass]);
      vec_hmgg_final_w[imass]->Fill(mgg,m_PUweight*generator_weight[imass]);
      vec_sighist[imass]->Fill(mgg,m_PUweight*generator_weight[imass]);
      vec_sighist_truth[imass]->Fill(mymc_m,generator_weight[imass]);
    }
    //--------------------------------------------

  }
  //========================================================================//
  //========== End of the loop over all entries ============================//
  //========================================================================//


  //------------------------------
  for(int imass=0;imass<m_Nmasses;imass++){
    for(int ic=0;ic<nchecks;ic++){
      vec_hCutFlow[imass]->SetBinContent(ic+1,n[imass][ic]);
      vec_hCutFlow_w[imass]->SetBinContent(ic+1,n_w[imass][ic]);
      vec_hCutFlowSum2_w[imass]->SetBinContent(ic+1,n2_w[imass][ic]);
    }
  //---------------------------------
  }
}


