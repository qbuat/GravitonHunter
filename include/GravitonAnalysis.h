// In this class ,the G->gammmagamma event preselection and all 
// the selection methods are defined. 
// The other classes which need those 
// methods inherits from GravitonAnalysis .
// See for example EventSelector and SignalTemplateCreator
// To use this class one needs to initialize 
// the constructor with an input tree . 
// Then an initialisation is needed :
// - GRL file (with SetGRL() )
// - G*->ee overlap removal file with SetGeeFile()
// - PileUp tool with SetPileUpTool()
// - Isolation type with SetIsoType()
// - Choose between data and mc with SetStreamType()
// To run on a subset of the dataset, one can use
// the SetEntries() method
#ifndef GravitonAnalysis_h
#define GravitonAnalysis_h

#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include "TString.h"
#include "TMath.h"
#include "ToolsTreeReader.h"
#include "ToolsTruthSelector.h"
#include "ElectronPhotonFourMomentumCorrection/egammaEnergyCorrectionTool.h"
#include "egammaAnalysisUtils/EnergyRescalerUpgrade.h"
#include "PileupReweighting/TPileupReweighting.h"
#include <stdlib.h>


class GravitonAnalysis{

 protected:

  int            m_nentries; //nb of tree entries
  TString        m_mctype;   // mc type 
  TString        m_isotype;   // mc type 
  TString        m_GRL;      // GRL string name
  TString        m_GeeFile;  // Gee string name
  bool           m_data;     // data=true,MC=false

  TreeReader    *m_Rd;
  TruthSelector *m_Ts;
  //== Energy Rescaling/Smearing ==//
  egRescaler::EnergyRescalerUpgrade m_Rescale;

  //== Energy Rescaling/Smearing GEO21 ==//
  AtlasRoot::egammaEnergyCorrectionTool m_Rescale_geo21;

  float          m_SmearFac;
  //== PileUp reweighting tool  ==//
  Root::TPileupReweighting* m_pileupTool;
  //===ee overlap removal ==//
  std::set<unsigned long long> m_ee_events;


 public:
  GravitonAnalysis();
  GravitonAnalysis(TTree*);
  virtual ~GravitonAnalysis();
  void     InitTree(TTree*);
  void     Init(TTree*);
  void     SetGRL(TString GRL);
  void     SetEntries(int entries);
  void     SetMCType(TString mctype = "" );
  void     SetIsoType(TString st= "CONE");
  void     SetStreamType(bool data = true);
  void     SetPileUpTool(TString mcfile,TString datafile);
  void     SetGeeFile(TString Geefile );
  void     SetTruthSelector(TTree* tree);

 protected:
  //--> Selection method (cannot be used outside the
  //--> GravitonAnalysis and daughter classes

  //========== Event Selection methods ==============//
  bool     EventGRLOK();
  bool     EventPVOK();
  bool     EventTrigOK();
  bool     EventCompletedOK();
  //=========== Photon Selection methods ============//
  bool     PhotonPtOK(int iph,double CutValueGeV=25.);
  bool     PhotonEtaOK(int iph);
  bool     PhotonOQOK(int iph);
  bool     PhotonAmbiguityResolverOK(const int & iph);
  bool     PhotonIsTightOK(int iph);
  bool     PhotonIsLooseOK(int iph);
  bool     PhotonIsolationOK(int iph,double CutValueGeV=5);
  bool     PhotonIsolation_toolOK(int iph,double CutValueGeV=5.);
  bool     PassPhotonCleaning(int iph);
  bool     PassNoiseCleaning(int iph);
  bool     PassLArCleaning(int iph);
  bool     PhotonIsTight_D3PDOK(int iph);
  bool     PhotonIsLoose_D3PDOK(int iph);
  //=========== Electron Selection methods ==============//
  bool     ElectronPtOK(int iele,double CutValueGeV=25.);
  bool     ElectronEtaOK(int iele);
  bool     ElectronOQOK(int iele);
  bool     ElectronTrackIsoOK(int iele);
  bool     ElectronIsMediumOK(int iele);
  //=========== Event level Selection method ============//
  bool     EventPreSelectionOK(int*,int*,
			       double*,double*);
  std::vector<int>     EventPreSelectionOK();

  //=== Modified variable calculation ===//
  float    PhotonTopoIsolation_tool(int iph);
  float    PhotonTopoIsolation_tool_geo21(int iph);
  float    PhotonIsolation_tool(int iph);
  float    PhotonTopoIsolationPtdep_tool(int iph);
  float    PhotonIsolationError_tool(int iph);
  float    PhotonRescaledPt(int iph);
  float    PhotonRescaledE(int iph);

  float    PhotonRescaledPt_geo21(int iph);
  float    PhotonRescaledE_geo21(int iph);

  float    ElectronRescaledE(int iph);
  float    PhotonSmearedPt(int iph,//type: 0=nom,1=down,2=up
			   int type=0);
  float    PhotonSmearedE(int iph,//type: 0=nom,1=down,2=up
			  int type=0);
  float    PhotonSmearingFactor(int iph,//type: 0=nom,1=down,2=up
				int type=0);

  //======= IsEM selection with PhotonIDTool ================//
  bool     PhotonIsEM_OK(int,bool);//data with rescaled pT
  bool     ElectronIsEM_OK(int,int,bool);//data with rescaled pT
  bool     PhotonIsEM_MCOK(int,bool,//MC with smeared pT+FF
			   bool VERBOSE =false);
  unsigned int PhotonIsEM(int);


 public:
  //--> Methods where the tree is not needed
  //--> Can be used for other purposes than
  //--> GravitonAnalysis and its daughters
  double   scaleForFFUnconvPhoton(double pT,//Apply a SF to FF
				  double eta,//for 1.81<|eta|<2.37
				  int isconv);


  //================ Remove ee overlaping events ==================//
  std::set<unsigned long long> build_eeset(const char* filename);    
  unsigned long long make_runevent( unsigned int run, 
				    unsigned int event );
  bool EventInZprimme( const std::set<unsigned long long>& container, 
		       unsigned int run, 
		       unsigned int event );

};

#endif // #ifdef GravitonAnalysis_cxx
