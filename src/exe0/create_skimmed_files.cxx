#include <iostream>
#include <map>
#include <TString.h>
#include <TSystem.h>
#include <TError.h>

#include "ToolsCommons.h"
#include "ToolsChainMaker.h"
#include "EventSelector.h"

using std::cerr;

int main(int argc, char ** argv){
  

  //------------- PROTECTION AGAINST WRONG NUMBER OF PARAMETERS ------------
  if( argc<5 ){
    std::cerr << "Wrong usage ! "
	      << argv[0]
	      << " iso_type data_or_mc stream_type input_file output_file"
	      << "\n" ;
    exit(1);
  }
  //------------------------------------------------------------------------

  //------------- INPUT PARAMETERS VALUES ------------------------------------
  TString iso_type      = argv[1]; // iso_type --> "TOPO" "TOPO_PTDEP" or "CONE"
  bool data             = atoi(argv[2]); //--> data=true, mc=false.
  TString mc_type       = argv[3]; // mc type--> data, "pythia_rs", "pythia_gg", "pythia_gj"
  TString input_file    = argv[4]; //Name of the input file
  TString output_file   = argv[5]; //label used to name the file
  //--------------------------------------------------------------------------
  std::cerr << "Processing : " 
	    << argv[0] << " " 
	    << iso_type << " " 
	    << data << " " 
	    << mc_type << " " 
	    << input_file << " " 
	    << output_file 
	    << "\n";
  //--------------------------------------------------------------------------

  if( iso_type != "CONE" && iso_type !="TOPO" && iso_type != "TOPO_PTDEP" )
    Fatal("create_skimmed_files", "WRONG ISOLATION TYPE");
  if( !data) {
    if(mc_type != "pythia_rs" && mc_type !="pythia_gg" && mc_type != "pythia_gj") 
      Fatal("create_skimmed_files", "WRONG MC TYPE");
  }

  std::cout << "-----------> RUN CONFIGURATION <----------" << std::endl;
  std::cout << "--> Type of isolation : " << iso_type << std::endl;
  if(data) 
    std::cout << "--> Type of stream : data" << std::endl;
  else {
    std::cout << "--> Type of stream : montecarlo" << std::endl;
    std::cout << "--> Type of montecarlo : " << mc_type << std::endl;
  }
  std::cout << "--> Input file  : " << input_file << std::endl;
  std::cout << "--> Output file : " << output_file << std::endl;
  std::cout << "-----------------><-----------------" << std::endl;



  Commons::Setup();

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

  std::cout << "START OF THE RUN" << std::endl;
  EventSelector *ES = new EventSelector(chain);

  std::cout << "SET THE STREAM TYPE (DATA/MC))" << std::endl;
  ES->SetStreamType(data);//data=true                                                                      

  if(data){
    std::cout << "DATA ONLY: SET THE GRL" << std::endl;
    ES->SetGRL(Commons::GrlFile);
    std::cout << "DATA ONLY: SET THE LIST OF G->ee EVENTS IN SR" << std::endl;
    ES->SetGeeFile(Commons::GeeFile);
  }else{
    std::cout << "MC ONLY: SET THE PILEUP TOOL" << std::endl;
    ES->SetPileUpTool(Commons::PileupMCFile,Commons::PileupDataFile);
    std::cout << "MC ONLY: SET THE TRUTH SELECTOR" << std::endl;
    ES->SetTruthSelector(chain);
  }
  std::cout << "SET THE ISOLATION TYPE" << std::endl;
  ES->SetIsoType(iso_type);
  if(!data){
    std::cout << "MC ONLY: SET THE MC_TYPE" << std::endl;
    ES->SetMCType(mc_type);
  }
  std::cout << "SET THE OUTPUT NAME" << std::endl;
  ES->SetOutputName(output_file);
  std::cout << "INITIALIZE THE OUTPUT TREE" << std::endl;
  ES->InitOutTree();
  std::cout << "START OF THE EVENT LOOP" << std::endl;
  ES->EventLoop();
  std::cout << "END OF THE EVENT LOOP" << std::endl;
  delete ES;
  std::cout << "Memory has been freed" << std::endl;


}
