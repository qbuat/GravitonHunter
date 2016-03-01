#include "BkgFitFramework.h"
#include "ToolsCommons.h"
#include <TError.h>
#include <RooMsgService.h>
#include <TIterator.h>


////////////////////////////////////////////////////////////////////////////////////////////
BkgFitFramework::BkgFitFramework(TTree* tree, std::pair<double,double > mggbin,TString FitParFile)
////////////////////////////////////////////////////////////////////////////////////////////
{
  // Constructor: list all the method to call to setup 
  // the criteria on the different variables.
  // It also runs the Init method with mggbin and FitParFile as input)

  if (tree == 0)
    Fatal("BkgFitFramework::BkgFitFramework", "No input tree specified !!");
  if (FitParFile == "")
    std::cout << "No fit parameters file specified !" << std::endl;

  SetTree(tree);
  SetMggBounds(0.,5000);
  SetIsoBounds(-100.,3000);
  SetPhotonsIsoCut(5.,5.);
  SetPhotonsPtCut(40.,30.);
  SetEtaCategory("NONE");
  SetLoosePrimeType(0);
  SetPtDependentIsocut(false);
  Init(mggbin,FitParFile);
}
//////////////////////////////////
BkgFitFramework::BkgFitFramework()
////////////////////////////////
{
  // Default constructor
  // Dummy setup
  m_tree           = 0;
  m_mggbin         = std::make_pair(0.,10000.);
  m_fitparfile     = "";
  SetMggBounds(0.,5000);
  SetIsoBounds(-100.,3000);
  SetPhotonsPtCut(40,30);
  SetPhotonsIsoCut(5,5);
  SetLoosePrimeType(0);
  SetPtDependentIsocut(false);
}
////////////////////////////////////////////////////////////
void BkgFitFramework::SetTree(TTree* tree,int entriesforfit)
///////////////////////////////////////////////////////////
{
  // This method sets the tree used to run the fits (internal m_tree variable)
  // The second argument can be used to run only on a subset a the sample

  if( entriesforfit != -1)
    m_tree = tree->CloneTree(entriesforfit);
  else
    m_tree = tree;

  if (tree == 0)
    Fatal("BkgFitFramework::SetTree", "No dataset tree specified !!");

}
/////////////////////////////////////////////////////////////////////
void BkgFitFramework::SetPhotonsIsoCut(double l_cut,double sl_cut)
////////////////////////////////////////////////////////////////////
{  
  m_IsoCut_L  = l_cut;
  m_IsoCut_SL = sl_cut;
}
/////////////////////////////////////////////////////////////////
void BkgFitFramework::SetPhotonsPtCut(double l_cut,double sl_cut)
////////////////////////////////////////////////////////////////
{
  m_PtCut_L = l_cut;
  m_PtCut_SL = sl_cut;
}
/////////////////////////////////////////////////////////
void BkgFitFramework::SetIsoBounds(double min,double max)
/////////////////////////////////////////////////////////
{  m_Isomin = min; m_Isomax = max; }
/////////////////////////////////////////////////////////
void BkgFitFramework::SetMggBounds(double min,double max)
/////////////////////////////////////////////////////////
{  m_Mggmin = min; m_Mggmax = max; }

//////////////////////////////////////////////////////////////////////////
void BkgFitFramework::Init(std::pair<double,double> mggbin, TString FitParFile)
//////////////////////////////////////////////////////////////////////////
{
  // Initialize the fit framework
  // Construct the different variables (RooRealVar) declared in the include file
  // Create the dataset (with application of the criteria defined in the by the 
  // SetXXXX() methods.
  // It also take as input the mgg range and the name of the files with parameters 
  // initial values and ranges.

  m_mggbin     = mggbin;
  SetParametersFile(FitParFile);


  // If using the modified isolation variable, the isolation variable as an other name
  // than the default one
  if( m_ptdepisocut ){
    m_Iso_L  = new RooRealVar("Iso_L_mod","E^{iso,mod}_{T,1} [GeV]",0,m_Isomin,m_Isomax);
    m_Iso_SL = new RooRealVar("Iso_SL_mod","E^{iso,mod}_{T,2} [GeV]",0,m_Isomin,m_Isomax);
  }else{
    m_Iso_L  = new RooRealVar("Iso_L","E^{iso}_{T,1} [GeV]",0,m_Isomin,m_Isomax);
    m_Iso_SL = new RooRealVar("Iso_SL","E^{iso}_{T,2} [GeV]",0,m_Isomin,m_Isomax);
  }


  m_eta_L  = new RooRealVar("eta_L","#eta_{1}",0,-4,4);
  m_eta_SL = new RooRealVar("eta_SL","#eta_{2}",0,-4,4);
  m_pt_L   = new RooRealVar("pT_L","#p_{T,1}",60,0,7000);
  m_pt_SL  = new RooRealVar("pT_SL","#p_{T,2}",60,0,7000);
  m_mgg    = new RooRealVar("mgg","m_{#gamma#gamma} [GeV]",100,m_Mggmin,m_Mggmax);
  //  m_mgg    = new RooRealVar("mgg","m_{#gamma#gamma} [GeV]",100, 180, 5000);
  //---------------------------------------------------

  // Define the categories based on the identification (tight and looseprime) menu
  // decision for each photon.
  m_IsTight_L = new RooCategory("IsTight_L","IsTight_L");
  m_IsTight_L->defineType("Tight",1);
  m_IsTight_L->defineType("NotTight",0);
  m_IsTight_SL = new RooCategory("IsTight_SL","IsTight_SL");
  m_IsTight_SL->defineType("Tight",1);
  m_IsTight_SL->defineType("NotTight",0);

  m_IsLoosePrime_L = new RooCategory(Form("IsLoosePrime%d_L",m_looseprime),
				     Form("IsLoosePrime%d_L",m_looseprime) );
  m_IsLoosePrime_L->defineType("LoosePrime",1);
  m_IsLoosePrime_L->defineType("NotLoosePrime",0);
  m_IsLoosePrime_SL = new RooCategory(Form("IsLoosePrime%d_SL",m_looseprime),
				      Form("IsLoosePrime%d_SL",m_looseprime) );
  m_IsLoosePrime_SL->defineType("LoosePrime",1);
  m_IsLoosePrime_SL->defineType("NotLoosePrime",0);

  //--------------------------------------------------------
  // Arrange the variables in RooArgSet to construct the RooDataSet
  RooArgSet Variables_L(*m_Iso_L,*m_IsTight_L,*m_eta_L,*m_pt_L,"leading");
  RooArgSet Variables_SL(*m_Iso_SL,*m_IsTight_SL,*m_eta_SL,*m_pt_SL,"subleading");
  if( m_looseprime != 0 ){
    Variables_L.add(*m_IsLoosePrime_L);
    Variables_SL.add(*m_IsLoosePrime_SL);
  }
  RooArgSet  Variables(Variables_L,Variables_SL,"variables");
  Variables.add(*m_mgg);


  //Construct the TString used to apply the selection to the RooDataSet.

  m_etacut_st = Commons::EtaCategory(m_etacat);
  m_ptcut_st  = Form("pT_L>%f && pT_SL>%f",m_PtCut_L,m_PtCut_SL);
  m_mggcut_st = Form("mgg>%f && mgg<%f",m_mggbin.first,m_mggbin.second);
  m_isocut_st = Form("Iso_L<%f && Iso_SL<%f",m_IsoCut_L,m_IsoCut_SL);



  if( m_ptdepisocut )
    m_isocut_st = Form("Iso_L_mod<%f && Iso_SL_mod<%f",m_IsoCut_L,m_IsoCut_SL);

  TString looseprimecut =  Form( "IsLoosePrime%d_L == 1 && IsLoosePrime%d_SL == 1",
				 m_looseprime,m_looseprime );
    
  //--------------------------------------------------------
  std::cout << looseprimecut << std::endl;
  std::cout << m_etacut_st << std::endl;
  std::cout << m_ptcut_st  << std::endl;
  std::cout << m_mggcut_st << std::endl;

  TString cut = Form("(%s) && (%s) && (%s) && (%s)", m_etacut_st.Data(), m_ptcut_st.Data(), m_mggcut_st.Data(), looseprimecut.Data());

  if (m_looseprime == 0)
    cut = Form("(%s) && (%s) && (%s)", 
	       m_etacut_st.Data(), m_ptcut_st.Data(), m_mggcut_st.Data());
  std::cout << "Selection = \n " << cut << std::endl;
  //-------------------------------------------------------------
  // Apply the different criteria .
  // After each step, the preceeding dataset is deleted to prevent memory overload
  // Construct the initial (no selection) RooDataSet 

  m_tree = m_tree->CopyTree(cut.Data());

  std::cout << "Building the dataset" << std::endl;
  m_DataSet_mggcut = new RooDataSet( "DataSet","DataSet",m_tree,Variables );
  std::cout << "The dataset has been constructed" << std::endl;

  // if( m_looseprime != 0 ) m_DataSet_looseprimecut = (RooDataSet*)m_DataSet->reduce(looseprimecut);
  // else m_DataSet_looseprimecut = (RooDataSet*)m_DataSet->Clone();
  delete m_DataSet;
  // m_DataSet_etacut      = (RooDataSet*)m_DataSet_looseprimecut->reduce(m_etacut_st);
  delete m_DataSet_looseprimecut;
  // m_DataSet_ptcut       = (RooDataSet*)m_DataSet_etacut->reduce(m_ptcut_st);
  delete m_DataSet_etacut;
  // m_DataSet_mggcut      = (RooDataSet*)m_DataSet_ptcut->reduce(m_mggcut_st);
  delete m_DataSet_ptcut;

  //------------------------------------------------------------------
}

/////////////////////////////////////////////////////////////////
void BkgFitFramework::SetConstantParameters(const RooArgSet& set)
/////////////////////////////////////////////////////////////////
{
  //--> This method set constant
  //    all RooRealVar parameters in a given RooArgSet
  TIterator *iterat= set.createIterator();
  RooRealVar *next = 0;
  while( (0 != (next= (RooRealVar*)iterat->Next())) )
    next->setConstant();
  delete iterat;  
}
///////////////////////////////////
BkgFitFramework::~BkgFitFramework()
///////////////////////////////////
{
  // //destructor
  delete m_Iso_L;
  delete m_Iso_SL;
  delete m_pt_L;
  delete m_pt_SL;
  delete m_eta_L;
  delete m_eta_SL;
  delete m_mgg; 
  delete m_IsTight_L;
  delete m_IsTight_SL;
  delete m_IsLoosePrime_L;
  delete m_IsLoosePrime_SL;
  delete m_DataSet_mggcut;
}
