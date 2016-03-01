#include <iostream>
#include <map>
#include <utility>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TChain.h>
#include <TBranch.h>
#include <TSystem.h>
#include <TF1.h> 
#include <TError.h>
#include "ToolsCommons.h"
#include "ToolsTreeReader.h"
#include "ToolsUtilities.h"

#include "BkgExtrapolation.h"


using std::cerr;

int main(int argc, char ** argv){
  

  //------------- PROTECTION AGAINST WRONG NUMBER OF PARAMETERS ------------
  if( argc<2 ){
    std::cerr << "Wrong usage ! "
	      << " analysis_new_x iso_type runnumber \n" ;
    exit(1);
  }
  //------------------------------------------------------------------------

  //------------- INPUT PARAMETERS VALUES ------------------------------------
  TString iso_type  = argv[1]; // iso_type --> "TOPO" or "CONE"
  int runnumber     = atoi(argv[2]); // runnumber
  //--------------------------------------------------------------------------
  std::cerr << "Processing : " << argv[0] << " " << iso_type << " " 
	    << Form("%d",runnumber) << "\n";


  Commons::Setup();
  
  TString iso_type_name;
  if( iso_type.Contains("TOPO") )
    iso_type_name = "TOPO";
  else if(iso_type.Contains("CONE") ){
    Fatal("main()", "No pt dependent cut on iso has been derived yet");
    iso_type_name = "CONE";
  }else Fatal("main()","Wrong iso type !");

  // TFile* fin = new TFile( Commons::outdir+"runs_19_06_13/datafile_"+iso_type_name+Form("_run%d.root",runnumber),"read");
  TFile* fin = new TFile( Commons::outdir+"pythiagamgam_TOPO.root","read");
  TTree *trD = (TTree*)fin->Get("tree");
  TreeReader Rd_D(trD);
  int nentriesD = (int)trD->GetEntries();

  // TFile* fout = new TFile( Commons::outdir+"runs_19_06_13/datafile_"+iso_type_name+Form("_PTDEP_run%d.root",runnumber),"RECREATE");
  TFile* fout = new TFile( Commons::outdir+"pythiagamgam_TOPO_PTDEP.root","RECREATE");
  TTree* trD_new = (TTree*)trD->CloneTree();

  double Iso_L_mod;
  double Iso_SL_mod;
  TBranch* biso_l; 
  TBranch* biso_sl;

  if (nentriesD != 0){
    biso_l  = trD_new->Branch("Iso_L_mod",&Iso_L_mod,"Iso_L_mod/D");   
    biso_sl = trD_new->Branch("Iso_SL_mod",&Iso_SL_mod,"Iso_SL_mod/D");
    for(int entry=0;entry<nentriesD;entry++){
      Rd_D.GetEntry(entry);
      Iso_L_mod  = Rd_D.GetVariable("Iso_L") - Commons::GetTopoIsoPtcorr(Rd_D.GetVariable("pT_L"));
      Iso_SL_mod = Rd_D.GetVariable("Iso_SL") - Commons::GetTopoIsoPtcorr(Rd_D.GetVariable("pT_SL"));
      biso_l->Fill();
      biso_sl->Fill();
    }
  }
  std::cout << "The entries have been processed" << std::endl;
  trD_new->Write();
  fout->Close();

  // std::cout << "***** MEMORY CLEANING ***** "<< std::endl;
  // delete trD;
  // std::cout << "***** END OF THE RUN ***** "<< std::endl;

  return 0;

}


