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
#include "BkgCategoryMerging.h"

using std::cerr;

int main(int argc, char ** argv){



  //------------- PROTECTION AGAINST WRONG NUMBER OF PARAMETERS ------------
  if( argc<2 ){
    std::cerr << "Wrong usage ! \n";
    exit(1);
  }
  //------------------------------------------------------------------------

  //------------- INPUT PARAMETERS VALUES ------------------------------------
  TString iso_type  = "TOPO_PTDEP";//argv[1]; // iso_type --> "TOPO" or "CONE"
  double isocut_l   = 8.;//atof(argv[2]); // lead isocut value
  double isocut_sl  = 8.;//atof(argv[2]); // sublead isocut value
  TString etacat    = argv[1];//argv[3];// eta category : --> "NONE", "CC", "CE", "EC", "EE", "EE_S", "EE_O"
  double leadcut    = 50.;//atoi(argv[4]); //leading photon cut
  double subleadcut = 50.;//atoi(argv[5]); //subleading photon cut
  bool doDATA       = false;
  bool doIRR        = false;
  bool doRED        = false;
  bool doPUR        = false;
  bool doSYS        = true;
  bool doFINAL      = true;
  TString filename  = argv[2];
  //--------------------------------------------------------------------------
  std::cerr << "Processing : " << argv[0] << " " << iso_type << " " 
	    << Form("%1.1f",isocut_l) << " " << etacat << " " 
	    << Form("%d %d",(int)leadcut,(int)subleadcut) << " "
	    << Form("%d %d %d %d %d %d \n",(int)doDATA,(int)doIRR,(int)doRED,(int)doPUR,(int)doSYS,(int)doFINAL)<< " " 
	    << filename;


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


  std::vector<TString> categories;
  categories.push_back("CC");
  categories.push_back("CE");
  categories.push_back("EC");
  categories.push_back("EE_S");
  categories.push_back("EE_O");

  for (unsigned int icat = 0; icat < categories.size(); icat++) {
    
    TString etacat = categories[icat];

    //------------------ Output Naming ------------------------------------------------
    TString extension = etacat+"_"+iso_type+"_"+Form("%1.1f",isocut_l);
    extension += Form("_lead%d_sublead%d",(int)leadcut,(int)subleadcut);
    extension += Form("_norm%d_%d_",Commons::norm_bin.first,Commons::norm_bin.second);
    extension += filename;



    TString datafile     = "./rootfiles/data_hist_"+extension+".root";
    TString irrfile      = "./rootfiles/irreducible_hist_"+extension+".root";
    //------------------------------------------------------------------------
    TString redfile      = "./rootfiles/reducible_hist_"+extension+"_L0.root";
    TString redfile_L4   = "./rootfiles/reducible_hist_"+extension+"_L4.root";
    TString redfile_L5   = "./rootfiles/reducible_hist_"+extension+"_L5.root";
    //------------------------------------------------------------------------
    TString yieldfile_L2 = "./rootfiles/2DFitResult_"+extension+"_L2.root";
    TString yieldfile_L3 = "./rootfiles/2DFitResult_"+extension+"_L3.root";
    TString yieldfile_L4 = "./rootfiles/2DFitResult_"+extension+"_L4.root";
    TString yieldfile_L5 = "./rootfiles/2DFitResult_"+extension+"_L5.root";
    //-------------------------------------------------------------------------
    
  

    if( doFINAL ){
      std::cout << "********************  FINAL ESTIMATE *******************" << std::endl;
      //----------- Final Background/data plot and output ---------------------------------------------
      BkgEstimatorUpgrade * BES = new BkgEstimatorUpgrade();
      BES->SetVerbose(true);
      BES->SetUncertFlag(doSYS);
      BES->SetDataFile(datafile);
      BES->SetIrrFile(irrfile) ;
      BES->SetRedFile(redfile) ;
      BES->SetYieldFile(yieldfile_L4); 
      if(doSYS){
	BES->SetRedSystFiles(redfile_L4,redfile_L5);
	BES->SetYieldSystFiles(yieldfile_L2,yieldfile_L3,yieldfile_L5);
      }
      BES->Init();
      BES->CreateBATInput("./rootfiles/Bkg_Template_"+extension+".root");
      BES->CreateBumpHunterInput("./rootfiles/Bkg_Total_"+extension+".root");  
      BES->GetFinalYieldsPerBin_NoUncert();
      delete BES;
      std::cout << "***************** END OF FINAL ESTIMATE ****************" << std::endl;
      //------------------------------------------------------------------------------------------------
    }
  }



    //------------------ Output Naming ------------------------------------------------
  TString final_extension = iso_type +"_" + Form("%1.1f",isocut_l);
  final_extension += Form("_lead%d_sublead%d",(int)leadcut,(int)subleadcut);
  final_extension += Form("_norm%d_%d_", Commons::norm_bin.first,Commons::norm_bin.second);
  final_extension += filename;

       // // //////--> Bkg total
     TString ftot_CC   =   "./rootfiles/Bkg_Total_CC_"+final_extension+".root";
     TString ftot_CE   =   "./rootfiles/Bkg_Total_CE_"+final_extension+".root";
     TString ftot_EC   =   "./rootfiles/Bkg_Total_EC_"+final_extension+".root";
     TString ftot_EE_O =   "./rootfiles/Bkg_Total_EE_O_"+final_extension+".root";
     TString ftot_EE_S =   "./rootfiles/Bkg_Total_EE_S_"+final_extension+".root";
     TString fout_tot  =   "./rootfiles/Bkg_Total_allcat_"+final_extension+".root";
     BkgCategoryMerging BCM_0("full",ftot_CC,ftot_CE,ftot_EC,ftot_EE_O,ftot_EE_S);
     BCM_0.StoreToRootFile(fout_tot);
   
   
     // // --> Bkg search region
     TString fsearch_CC   =   "./rootfiles/Bkg_Template_CC_"+final_extension+".root";
     TString fsearch_CE   =   "./rootfiles/Bkg_Template_CE_"+final_extension+".root";
     TString fsearch_EC   =   "./rootfiles/Bkg_Template_EC_"+final_extension+".root";
     TString fsearch_EE_O =   "./rootfiles/Bkg_Template_EE_O_"+final_extension+".root";
     TString fsearch_EE_S =   "./rootfiles/Bkg_Template_EE_S_"+final_extension+".root";
     TString fout_search  =   "./rootfiles/Bkg_Template_allcat_"+final_extension+".root";
     BkgCategoryMerging BCM_1("search",fsearch_CC,fsearch_CE,fsearch_EC,fsearch_EE_O,fsearch_EE_S);
     BCM_1.StoreToRootFile_Search(fout_search);



  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  return 0;

}


