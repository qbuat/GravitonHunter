#include "InvMassCalculator.h"


////////////////////////////////////////////////////////////
InvMassCalc::InvMassCalc() : GravitonAnalysis()
///////////////////////////////////////////////////////////
{
  //Default constructor
}
////////////////////////////////////////////////////////////
InvMassCalc::InvMassCalc(TTree* tree) : GravitonAnalysis()
///////////////////////////////////////////////////////////
{
  //Default constructor
  InitTree(tree);
}
///////////////////////////////////////
InvMassCalc::~InvMassCalc()
///////////////////////////////////////
{
  //Default destructor
}

/////////////////////////////////////////////////////////////////
double InvMassCalc::GetInvMass(int entry,int ilead,int isublead)
////////////////////////////////////////////////////////////////
{
  m_Rd->GetEntry(entry);
  double temp=-9999.;
  if( ilead==-1 || isublead==-1){
     std::cout << "Missing leader or subleader !!" << std::endl;
  }
  int ndiph=(int)m_Rd->GetVariable("diphoton_n");
  int ncand=0;
  for(int idiph=0; idiph<ndiph; idiph++){
    int diph_ilead=(int)m_Rd->GetVariable( Form("diphoton_ph_index[%d][0]",idiph) );
    int diph_isublead=(int)m_Rd->GetVariable( Form("diphoton_ph_index[%d][1]",idiph) );
    if( (diph_ilead==ilead && diph_isublead==isublead) 
	||(diph_ilead==isublead && diph_isublead==ilead) ){
      temp=m_Rd->GetVariable( Form("diphoton_HPV_mgg[%d]",idiph) );
      ++ncand;
    }
  }
  if(ncand==0){
    std::cout << "No matching ph_ and diphoton_ " << std::endl;
  }
  
  return temp;
}
///////////////////////////////////////////////////////////////////////////////
double InvMassCalc::GetInvMass2(int entry,int ilead,int isublead,int smearing)
//////////////////////////////////////////////////////////////////////////////
{
  m_Rd->GetEntry(entry);
  double DiPhoton_zcommon=-999;
  float energy_corrected_EMscale_leading=-99999;
  float energy_corrected_EMscale_subleading=-99999;
  float eta_corrected_leading=-99999;
  float eta_corrected_subleading=-99999;
  int ndiph=(int)m_Rd->GetVariable("diphoton_n");
  int ncand=0;
  for (int idiph=0;idiph<ndiph;idiph++) { //for each diphoton object
    int diph_ilead=(int)m_Rd->GetVariable( Form("diphoton_ph_index[%d][0]",idiph) );
    int diph_isublead=(int)m_Rd->GetVariable( Form("diphoton_ph_index[%d][1]",idiph) );
    if( (diph_ilead==ilead && diph_isublead==isublead) 
	||(diph_ilead==isublead && diph_isublead==ilead) ){
      DiPhoton_zcommon=m_Rd->GetVariable( Form("diphoton_HPV_zcommon[%d]",idiph) );
      ncand++;
    }
  }
  if(ncand==0){
    std::cout << "No matching ph_ and diphoton_ " << std::endl;
  }

  float ph_etas1_lead=(float)m_Rd->GetVariable( Form("ph_etas1[%d]",ilead) );
  float ph_etas1_sublead=(float)m_Rd->GetVariable( Form("ph_etas1[%d]",isublead) );
 
  CorrectEta(DiPhoton_zcommon,
	     ph_etas1_lead,ph_etas1_sublead,
	     &eta_corrected_leading,&eta_corrected_subleading);

  // End correction to take mass "A La CaloPointing"
  //--------------------------------------------------------------------------
  
  float energy_correction_leading=1;
  float energy_correction_subleading=1;
  if(smearing>=0){
    energy_correction_leading=PhotonSmearingFactor(ilead,smearing);
    energy_correction_subleading=PhotonSmearingFactor(isublead,smearing);
  }
  float ph_E_lead=(float)m_Rd->GetVariable( Form("ph_E[%d]",ilead) );
  float ph_E_sublead=(float)m_Rd->GetVariable( Form("ph_E[%d]",isublead) );
  energy_corrected_EMscale_leading=energy_correction_leading*ph_E_lead;
  energy_corrected_EMscale_subleading=energy_correction_subleading*ph_E_sublead;

  //================================
  //final correction
  float ph_pt_corrected_leading=energy_corrected_EMscale_leading/cosh(eta_corrected_leading); 
  //final
  float ph_pt_corrected_subleading=energy_corrected_EMscale_subleading/cosh(eta_corrected_subleading); 
  float phi_lead=(float)m_Rd->GetVariable( Form("ph_phi[%d]",ilead) );
  float phi_sublead=(float)m_Rd->GetVariable( Form("ph_phi[%d]",isublead) );

  TLorentzVector Vp1_corrected;
  TLorentzVector Vp2_corrected;
  Vp1_corrected.SetPtEtaPhiE(ph_pt_corrected_leading,
			     eta_corrected_leading,
			     phi_lead,
			     energy_corrected_EMscale_leading); //final corrections
  Vp2_corrected.SetPtEtaPhiE(ph_pt_corrected_subleading,
			     eta_corrected_subleading,
			     phi_sublead,
			     energy_corrected_EMscale_subleading); //final
  TLorentzVector Vdiphoton_corrected=Vp1_corrected+Vp2_corrected;
  double mass_gg_corrected=Vdiphoton_corrected.M();
  return mass_gg_corrected;
}
//////////////////////////////////////////////////////////////////////////////
double InvMassCalc::GetInvMass3(int entry,int ilead,int isublead,int smearing)
//////////////////////////////////////////////////////////////////////////////
{
  m_Rd->GetEntry(entry);
  double DiPhoton_zcommon=-999;
  float energy_corrected_EMscale_leading=-99999;
  float energy_corrected_EMscale_subleading=-99999;
  float eta_corrected_leading=-99999;
  float eta_corrected_subleading=-99999;
  int ndiph=(int)m_Rd->GetVariable("diphoton_n");
  int ncand=0;
  for (int idiph=0;idiph<ndiph;idiph++) { //for each diphoton object
    int diph_ilead=(int)m_Rd->GetVariable( Form("diphoton_ph_index[%d][0]",idiph) );
    int diph_isublead=(int)m_Rd->GetVariable( Form("diphoton_ph_index[%d][1]",idiph) );
    if( (diph_ilead==ilead && diph_isublead==isublead) 
	||(diph_ilead==isublead && diph_isublead==ilead) ){
      DiPhoton_zcommon=m_Rd->GetVariable( Form("diphoton_HPV_zcommon_primVtxLH[%d]",idiph) );
      ncand++;
    }
  }
  if(ncand==0){
    std::cout << "No matching ph_ and diphoton_ " << std::endl;
  }

  float ph_etas1_lead=(float)m_Rd->GetVariable( Form("ph_etas1[%d]",ilead) );
  float ph_etas1_sublead=(float)m_Rd->GetVariable( Form("ph_etas1[%d]",isublead) );
 
  CorrectEta(DiPhoton_zcommon,
	     ph_etas1_lead,
	     ph_etas1_sublead,
	     &eta_corrected_leading,
	     &eta_corrected_subleading);
  // End correction to take mass "A La CaloPointing"
  //--------------------------------------------------------------------------
  
  float energy_correction_leading=1;
  float energy_correction_subleading=1;
  if(smearing>=0){
    energy_correction_leading=PhotonSmearingFactor(ilead,smearing);
    energy_correction_subleading=PhotonSmearingFactor(isublead,smearing);
  }
  float ph_E_lead=(float)m_Rd->GetVariable( Form("ph_E[%d]",ilead) );
  float ph_E_sublead=(float)m_Rd->GetVariable( Form("ph_E[%d]",isublead) );
  energy_corrected_EMscale_leading=energy_correction_leading*ph_E_lead;
  energy_corrected_EMscale_subleading=energy_correction_subleading*ph_E_sublead;

  //================================
  //final correction
  float ph_pt_corrected_leading=energy_corrected_EMscale_leading/cosh(eta_corrected_leading); 
  //final
  float ph_pt_corrected_subleading=energy_corrected_EMscale_subleading/cosh(eta_corrected_subleading); 
  float phi_lead=(float)m_Rd->GetVariable( Form("ph_phi[%d]",ilead) );
  float phi_sublead=(float)m_Rd->GetVariable( Form("ph_phi[%d]",isublead) );

  TLorentzVector Vp1_corrected;
  TLorentzVector Vp2_corrected;
  Vp1_corrected.SetPtEtaPhiE(ph_pt_corrected_leading,
			     eta_corrected_leading,
			     phi_lead,
			     energy_corrected_EMscale_leading);
  //final corrections
  Vp2_corrected.SetPtEtaPhiE(ph_pt_corrected_subleading,
			     eta_corrected_subleading,
			     phi_sublead,
			     energy_corrected_EMscale_subleading); //final
  TLorentzVector Vdiphoton_corrected=Vp1_corrected+Vp2_corrected;
  double mass_gg_corrected=Vdiphoton_corrected.M();

  return mass_gg_corrected;
}
//////////////////////////////////////////////////////////////////
void InvMassCalc::CorrectEta(float PV_z,
			     float initial_etas1_leading,
			     float initial_etas1_subleading,
			     float *corrected_etas1_leading,
			     float *corrected_etas1_subleading)
////////////////////////////////////////////////////////////////////
{
  //####  PATCH v0.05 : more pretty code
  // TO GET THE CORRECTED INVARIANT MASS WITH PAU NTUPLES
  // without to have to reproduce PAU ntuples (for which the selection is different)
  
  double R_photon_front;
  double Z_photon_front;
    
  if (fabs(initial_etas1_leading) < 1.5) {      // barrel
    R_photon_front=ReturnRZ_1stSampling_cscopt2(initial_etas1_leading);
    Z_photon_front=R_photon_front*sinh(initial_etas1_leading);
  }
  else {                    // endcap
    Z_photon_front=ReturnRZ_1stSampling_cscopt2(initial_etas1_leading);
    R_photon_front=Z_photon_front/sinh(initial_etas1_leading);
  }
    
  *corrected_etas1_leading=asinh((Z_photon_front-PV_z)/R_photon_front);
    
  if (fabs(initial_etas1_subleading) < 1.5) {      // barrel
    R_photon_front=ReturnRZ_1stSampling_cscopt2(initial_etas1_subleading);
    Z_photon_front=R_photon_front*sinh(initial_etas1_subleading);
  }
  else {                    // endcap
    Z_photon_front=ReturnRZ_1stSampling_cscopt2(initial_etas1_subleading);
    R_photon_front=Z_photon_front/sinh(initial_etas1_subleading);
  }
  
  *corrected_etas1_subleading=asinh((Z_photon_front-PV_z)/R_photon_front);
}


////////////////////////////////////////////////////////////////////////////////////////////////
double InvMassCalc::GetCorrectedInvMass(double Elead,double etaS1lead,double philead,
					      double Esublead,double etaS1sublead,double phisublead,
					      double PV_ID)
///////////////////////////////////////////////////////////////////////////////////////////////
{

  double R_photon_front;
  double Z_photon_front;
  //leading photon
  if (fabs(etaS1lead) < 1.5) {                                    // barrel
    R_photon_front=ReturnRZ_1stSampling_cscopt2(etaS1lead);
    Z_photon_front=R_photon_front*sinh(etaS1lead);
  }
  else {                                                          // endcap
    Z_photon_front=ReturnRZ_1stSampling_cscopt2(etaS1lead);
    R_photon_front=Z_photon_front/sinh(etaS1lead);
  }
  
  double eta_corrected_leading=asinh((Z_photon_front- PV_ID)/R_photon_front);
  //subleading photon
  if (fabs(etaS1sublead) < 1.5) {      // barrel
    R_photon_front=ReturnRZ_1stSampling_cscopt2(etaS1sublead);
    Z_photon_front=R_photon_front*sinh(etaS1sublead);
 }
  else {                    // endcap
    Z_photon_front=ReturnRZ_1stSampling_cscopt2(etaS1sublead);
    R_photon_front=Z_photon_front/sinh(etaS1sublead);
  }
  
  double eta_corrected_subleading=asinh((Z_photon_front- PV_ID)/R_photon_front);
  
  //leading photon : 
  float energy_corrected_EMscale_leading = Elead;  //apply energy corrections due to EMscale
  float ph_pt_corrected_leading= energy_corrected_EMscale_leading / cosh(eta_corrected_leading);
  TLorentzVector Vp1_corrected;
  Vp1_corrected.SetPtEtaPhiE(ph_pt_corrected_leading,
			     eta_corrected_leading, 
			     philead ,
			     energy_corrected_EMscale_leading);
  
  //subleading photon :
  float energy_corrected_EMscale_subleading= Esublead;
  float ph_pt_corrected_subleading=energy_corrected_EMscale_subleading/cosh(eta_corrected_subleading);
  TLorentzVector Vp2_corrected;
  Vp2_corrected.SetPtEtaPhiE(ph_pt_corrected_subleading,
			     eta_corrected_subleading,
			     phisublead,
			     energy_corrected_EMscale_subleading);
  
  TLorentzVector Vdiphoton_corrected=Vp1_corrected+Vp2_corrected;
  double mass_gg_corrected = Vdiphoton_corrected.M();
 
  return mass_gg_corrected  ;

}

////////////////////////////////////////////////////////////////////////
double InvMassCalc::ReturnRZ_1stSampling_cscopt2(double eta_1st_sampling)
////////////////////////////////////////////////////////////////////////
{
  //adapted from CaloDepthTool.cxx, double CaloDepthTool::cscopt2_parametrized(const CaloCell_ID::CaloSample sample,
  //  const double eta, const double /*phi*/ )
  // No warranty !!!
  
  double millimeter=1;
  
  double radius = -99999;
  
  float aeta_1st_sampling = fabs(eta_1st_sampling);
  
  if (aeta_1st_sampling<1.5) { //barrel
    if (aeta_1st_sampling < 0.8)
      radius = (1558.859292 - 4.990838*aeta_1st_sampling - 21.144279*aeta_1st_sampling*aeta_1st_sampling)*millimeter;
    else
      radius = (1522.775373 + 27.970192*aeta_1st_sampling - 21.104108*aeta_1st_sampling*aeta_1st_sampling)*millimeter;
  }
  else { //endcap
    if (aeta_1st_sampling < 1.5)
      radius = (12453.297448 - 5735.787116*aeta_1st_sampling)*millimeter;
    else
      radius = 3790.671754*millimeter;
    if (eta_1st_sampling < 0.) radius = -radius;
  }
  
  return radius;
}
///////////////////////////////////////////////////////////
double InvMassCalc::EtaS1PVCorrected(double etas1, double pv)
//////////////////////////////////////////////////////////
{
  
  double R_photon_front;
  double Z_photon_front;
  
  if (fabs(etas1) < 1.5) {      // barrel
    R_photon_front=ReturnRZ_1stSampling_cscopt2(etas1);
    Z_photon_front=R_photon_front*sinh(etas1);
  }
  else {                    // endcap
    Z_photon_front=ReturnRZ_1stSampling_cscopt2(etas1);
    R_photon_front=Z_photon_front/sinh(etas1);
  }
  
  return asinh((Z_photon_front-pv)/R_photon_front);
}  
