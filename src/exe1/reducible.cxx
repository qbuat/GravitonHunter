#include <iostream>
#include <map>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TChain.h>
#include <TBranch.h>
#include <TH1.h>
#include <TSystem.h>
#include <TF1.h> 
#include <TError.h>
#include "ToolsCommons.h"
#include "ToolsTreeReader.h"
#include "ToolsUtilities.h"

#include "BkgExtrapolation.h"
#include "Bkg2DFit.h"
#include "BkgEstimatorUpgrade.h"
#include "GoodRunsLists/DQHelperFunctions.h"

using std::cerr;

int main(int argc, char ** argv){


  //------------- PROTECTION AGAINST WRONG NUMBER OF PARAMETERS ------------
  if( argc<3 ){
    std::cerr << "Wrong usage ! "
	      << " data_x etacat filename looseprime \n" ;
    exit(1);
  }
  //------------------------------------------------------------------------



  //------------- INPUT PARAMETERS VALUES ------------------------------------
  TString etacat    = argv[1];// eta category : --> "NONE", "CC", "CE", "EC", "EE", "EE_S", "EE_O"
  TString filename  = argv[2];
  int looseprime = atoi(argv[3]);
  TString iso_type  = "TOPO_PTDEP"; //argv[1]; // iso_type --> "TOPO" or "CONE"
  double isocut_l   = 8.; //atof(argv[2]); // lead isocut value
  double isocut_sl  = 8.; //atof(argv[2]); // sublead isocut value
  double leadcut    = 50; //atoi(argv[4]); //leading photon cut
  double subleadcut = 50.; //atoi(argv[5]); //subleading photon cut
  bool doDATA       = false;
  bool doIRR        = false;
  bool doRED        = true;
  bool doPUR        = false;
  bool doSYS        = false;
  bool doFINAL      = false;
  //--------------------------------------------------------------------------
  std::cerr << "Processing : " << argv[0] << " " << iso_type << " " 
	    << Form("%1.1f",isocut_l) << " " << etacat << " " 
	    << Form("%d %d",(int)leadcut,(int)subleadcut) << " "
	    << Form("%d %d %d %d %d %d %s \n",(int)doDATA,(int)doIRR,(int)doRED,(int)doPUR,(int)doSYS,(int)doFINAL, filename.Data());
  
  std::cout << "looseprime = " << looseprime << std::endl;
  Commons::Setup();

  // Commons::SetPurityUncert(0.02);
  //----> Normalize between 180 and 409 (bin 26 to 47) for 50-50 config
  Commons::SetNormBin(26,47);

  std::cout << "--------> RUN PROCESSING <-------" << std::endl;
  if(doDATA)  std::cout << "---> Process DATA HISTOGRAMS <---" << std::endl;
  if(doIRR)   std::cout << "---> Process IRREDUCIBLE BKG <---" << std::endl;
  if(doRED)   std::cout << "---> Process REDUCIBLE BKG <-----" << std::endl;
  if(doPUR)   std::cout << "---> Process PURITY ESTIMATE <---" << std::endl;
  if(doSYS)   std::cout << "---> Process UNCERTAINTIES <-----" << std::endl;
  if(doFINAL) std::cout << "---> Process FINAL ESTIMATE <----" << std::endl;
  std::cout << "----------------><---------------" << std::endl;

  std::cout << "-----------> RUN CONFIGURATION <----------" << std::endl;
  std::cout << Form( "--> Normalization region : [%1.2f,%1.2f]",
		     Commons::binning[Commons::norm_bin.first-1],
		     Commons::binning[Commons::norm_bin.second]) << std::endl;
  std::cout << Form( "--> Pt cuts : (%d,%d)",(int)leadcut,(int)subleadcut) << std::endl; 
  std::cout << "--> Type of isolation : " << iso_type << std::endl;
  std::cout << Form( "--> Isolation cut value = %1.1f GeV",isocut_l) << std::endl; 
  if( iso_type.Contains("PTDEP" ))
    std::cout << "--> Use pT dependent isocut" << std::endl; 
  std::cout << "--> Eta category : " << etacat << std::endl;
  std::cout << "-----------------><-----------------" << std::endl;

  
  TString iso_type_name;
  if( iso_type.Contains("TOPO") )
    iso_type_name = "TOPO";
  else if(iso_type.Contains("CONE") )
    iso_type_name = "CONE";
  else Fatal("main()","Wrong iso type !");

  //------------------ Data Handling -------------------------------------------------------
  TChain *trD = new TChain("tree");
  trD->Add("/sps/atlas/q/qbuat/DIPHOTON_2012_REF/fuerSarahEtAl/Look/datafile_TOPO_PTDEP_run*.root");
  // int nentriesD = (int)trD->GetEntriesFast();


  //------------------------- Fit Parameters Files -------------------------------
  TString fitpar_ext = etacat + "_";
  fitpar_ext += iso_type_name + "_";
  fitpar_ext += Form("%d_%d", (int)leadcut, (int)subleadcut);
  TString FitParFile = "FitParameters/FitStudy_"+fitpar_ext+ Form("_L%d.config", looseprime);
  std::cout << "Using the param file " << FitParFile << std::endl;
  //------------------------------------------------------------------------------

  //------------------ Output Naming ------------------------------------------------
  TString extension = etacat + "_" + iso_type + "_" + Form("%1.1f",isocut_l);
  extension += Form("_lead%d_sublead%d", (int)leadcut, (int)subleadcut);
  extension += Form("_norm%d_%d_", Commons::norm_bin.first, Commons::norm_bin.second);
  extension += filename;


  //------------------------------------------------------------------------
  TString redfile      = "./rootfiles/reducible_hist_" + extension + Form("_L%d.root", looseprime);
  //-------------------------------------------------------------------------



  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //
  //                    REDUCIBLE REDUCIBLE REDUCIBLE

  if( doRED ){
    std::cout << "********************  REDUCIBLE *******************" << std::endl;

    //--------------- Create Reducible background mgg file --------------------------------------
    std::pair<double,double> fitrange( Commons::binning[Commons::norm_bin.first-1],1000 );
    std::pair<double,double> mggbin( Commons::binning[Commons::norm_bin.first-1],1000);

    //    pair<double,double> fitrange( Commons::binning[Commons::norm_bin.first-1],1000 );
    //    pair<double,double> mggbin( Commons::binning[Commons::norm_bin.first-1],3000);
    
    BkgExtrapolation * BE = new BkgExtrapolation();
    BE->SetTree(trD);
    BE->SetFitRange();
    BE->SetNPE(1000);
    BE->SetMggBounds(Commons::binning[0],Commons::binning[Commons::nBins]);
    BE->SetPhotonsIsoCut(isocut_l,isocut_sl);
    BE->SetPhotonsPtCut(leadcut,subleadcut);
    BE->SetEtaCategory(etacat);
    if( iso_type == "TOPO_PTDEP" )
      BE->SetPtDependentIsocut(true);

    if (looseprime != 0)
      BE->SetLoosePrimeType(looseprime);
    BE->Init(mggbin,FitParFile);
    BE->SetFitRange(fitrange.first,fitrange.second);
    BE->Fitter(true);
    BE->Extrapolate(true);
    BE->StoreToRootFile(redfile);
    //    BE->StoreToRootFile("./rootfiles/reducible_hist_"+extension+"data_extended_mgg_range_redfit.root");

    delete BE;
    //----------------------------------------------------------------------------------------------
  }    
  std::cout << "***** MEMORY CLEANING ***** "<< std::endl;
  delete trD;
  std::cout << "***** END OF THE RUN ***** "<< std::endl;

  return 0;

}


