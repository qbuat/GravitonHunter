#include <iostream>
#include <map>
#include <TString.h>
#include <TSystem.h>
#include <TError.h>

#include "ToolsCommons.h"
#include "ToolsChainMaker.h"
#include "SignalTemplateCreator.h"

using std::cerr;

int main(int argc, char ** argv){
  

  //------------- PROTECTION AGAINST WRONG NUMBER OF PARAMETERS ------------
  if( argc<8 ){
    std::cerr << "Wrong usage ! "
	      << argv[0]
	      << " input_file iso_type isocut etacat leadcut subleadcut coupling"
	      << " outfile_label \n" ;
    exit(1);
  }
  //------------------------------------------------------------------------

  //------------- INPUT PARAMETERS VALUES ------------------------------------
  TString input_file    = argv[1]; //Name of the input file
  TString iso_type      = argv[2]; // iso_type --> "TOPO" "TOPO_PTDEP" or "CONE"
  double isocut_l       = atof(argv[3]); // lead isocut value
  double isocut_sl      = atof(argv[3]); // sublead isocut value
  TString etacat        = argv[4];// eta category : --> "NONE", "CC", "CE", "EC", "EE", "EE_S", "EE_O"
  double leadcut        = atoi(argv[5]); //leading photon cut
  double subleadcut     = atoi(argv[6]); //subleading photon cut
  double coupling       = atof(argv[7]); //RS coupling value
  TString outfile_label = argv[8]; //label used to name the file
  //--------------------------------------------------------------------------
  std::cerr << "Processing : " << argv[0] << " " << input_file << " " 
	    << iso_type << " " 
	    << Form("%1.1f",isocut_l) << " " << etacat << " " 
	    << Form("%d %d",(int)leadcut,(int)subleadcut) << " "
	    << outfile_label << "\n";

  Commons::Setup();

  TString output = "./rootfiles/";
  output += "Signal_Total_"+etacat+"_"+iso_type;
  output += Form("_%1.1f_lead%d_sublead%d_coupling%1.2f",isocut_l,(int)leadcut,(int)subleadcut,coupling);
  output += outfile_label+".root";

  std::cout << "-----------> RUN CONFIGURATION <----------" << std::endl;
  std::cout << Form( "--> Pt cuts : (%d,%d)",(int)leadcut,(int)subleadcut) << std::endl; 
  std::cout << "--> Type of isolation : " << iso_type << std::endl;
  std::cout << Form( "--> Isolation cut value = %1.1f GeV",isocut_l) << std::endl; 
  if( iso_type.Contains("PTDEP" ))
    std::cout << "--> Use pT dependent isocut" << std::endl; 
  std::cout << "--> Eta category : " << etacat << std::endl;
  std::cout << "--> Coupling : " << coupling  << std::endl;
  std::cout << "-----------------><-----------------" << std::endl;

  
  TString iso_type_name;
  if( iso_type.Contains("TOPO") )
    iso_type_name = "TOPO";
  else if(iso_type.Contains("CONE") )
    iso_type_name = "CONE";
  else Fatal("main()","Wrong iso type !");

  //------------------ Data Handling -------------------------------------------------------
  std::ifstream ifs(input_file);
  //---------------------------------------
  std::string argStr;
  char buf[256+1];
  unsigned int delpos;
  while (true){
    ifs.read(buf,256);
    if (ifs.eof()){
      if (ifs.gcount() == 0) break;
      delpos = ifs.gcount()-1;
    }else{
      delpos = ifs.gcount();
    }
    buf[delpos] = 0x00;
    argStr += buf;
  }
  // cout<<"argStr  ="<<argStr<<"."<<endl;
  //---------------------------------------

  //======================================//
  TChain *chain=0;
  ChainMaker MyChainMaker("photon");
  chain = new TChain("photon");
  std::cout << "Build the files chain ..." ;
  MyChainMaker.SetTreeFiles(argStr);
  std::cout << "done !" << std::endl;
  chain = MyChainMaker.GetChain();
  //______________________________________//


  std::cout << "STEP 0" << std::endl;
  SignalTemplateCreator *STC = new SignalTemplateCreator(chain);
  std::cout << "STEP 1" << std::endl;
  STC->SetStreamType(false);//data=true
  std::cout << "STEP 2" << std::endl;
  STC->SetPileUpTool(Commons::PileupMCFile,Commons::PileupDataFile);
  std::cout << "STEP 3" << std::endl;
  STC->SetTruthSelector(chain);
  std::cout << "STEP 4" << std::endl;
  STC->SetIsoType(iso_type);
  std::cout << "STEP 5" << std::endl;
  STC->SetMCType("pythia_rs");
  std::cout << "STEP 6" << std::endl;
  STC->SetOutputName(output);
  std::cout << "STEP 7" << std::endl;

  //--> In GeV : Set mass = LowMass + i*MassSpacing (0<i<Nmasses)
  STC->SetNmasses(25);
  STC->SetLowMass(500);
  STC->SetMassSpacing(100);
  STC->SetCoupling(coupling);  

  std::cout << "STEP 8" << std::endl;
  STC->SetLeadPtCut(leadcut);
  STC->SetSubleadPtCut(subleadcut);
  STC->SetLeadIsoCut(isocut_l);
  STC->SetSubleadIsoCut(isocut_sl);

  std::cout << "STEP 9" << std::endl;
  STC->EventLoop();
  std::cout << "STEP 10" << std::endl;


  STC->CreateOutputFile(output);
  std::cout << "DONE !" << std::endl;
  delete STC;
  std::cout << "Memory has been freed" << std::endl;





}
