#include <iostream>
#include <stdlib.h>
#include <map>

#include <TSystem.h>
#include <TError.h>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TBranch.h>
#include <TGraphErrors.h>

#include "ToolsCommons.h"
#include "ToolsTreeReader.h"
#include "SingleParticleIsolationFitter.h"
#include "IsolationGrapher.h"

using std::cerr;

int main(int argc, char ** argv){
  

  //------------- PROTECTION AGAINST WRONG NUMBER OF PARAMETERS ------------

  std::cout<<argc<<std::endl;
  if( argc<8 ){
    std::cerr << "Wrong usage ! "
	      << argv[0] << " "
	      << "sample absetamin absetamax ptmin ptmax isconv outfile \n" ;
    exit(1);
  }
  //------------------------------------------------------------------------

  //------------- INPUT PARAMETERS VALUES ------------------------------------
  TString sample    = argv[1];// input file

  double  absetamin = atof(argv[2]); // lower |eta| bound
  double  absetamax = atof(argv[3]);// upper |eta| bound
  double  ptmin     = atof(argv[4]);// liwer pt bound
  double  ptmax     = atof(argv[5]);// upper pt bound
  int     isconv    = atof(argv[6]);// upper pt bound

  TString outfile   = argv[7]; // output file name
  //--------------------------------------------------------------------------
  std::cerr << "Processing : " << argv[0]  << " " 
	    << sample <<" " 
	    << absetamin << " " 
	    << absetamax << " " 
	    << ptmin << " "    
	    << ptmax << " "    
	    << isconv << " "  
	    << outfile  
	    << "\n";

  Commons::Setup();
  

  TChain *tree = new TChain("tree");
  tree->Add(sample);
  std::cout<<tree->GetEntries()<<std::endl;
  std::cout << "TOTO" << std::endl;

  SingleParticleIsolationFitter* IF = new SingleParticleIsolationFitter();
  IF->SetEntries(tree);
  IF->SetTree(tree);
  
  //  IF->SetParFile("FitParameters/Iso_PtNPVeta.config");
  IF->SetParFile("FitParameters/Iso_singlephoton.config");

  IF->SetPtEtaIsoTightConvNames("_pt","_eta","_iso40","_istight","_isconv");
  IF->SetPtBounds(ptmin,ptmax);
  IF->SetEtaBounds(absetamin,absetamax);
  IF->SetIsConv(isconv);

  if(sample.Contains("noPU")) IF->SetIsoBounds(-2000,10000);
  else                        IF->SetIsoBounds(-5000,25000);
  
  
  IF->Init_Vars();
  IF->Init_DataSets();
  IF->Fitter(true);
  //  IF->Fitter(false);
  IF->StoreToRootFile(outfile);

}


