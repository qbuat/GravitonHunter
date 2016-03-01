#include "EventProcessing.h"
#include "EventSelector.h"
#include "SignalTemplateCreator.h"
#include "ToolsChainMaker.h"
#include "ToolsCommons.h"

/////////////////////////////////////
EventProcessing::EventProcessing() {}
/////////////////////////////////////

//////////////////////////////////////
EventProcessing::~EventProcessing() {}
//////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
void EventProcessing::MethodSkeleton(TString treename,bool data,TString mctype)
/////////////////////////////////////////////////////////////////////////////
{
  //======================================//
  TChain *chain=0;
  ChainMaker MyChainMaker(treename);
  chain = new TChain(treename);
  MyChainMaker.SetTreeFiles(_DataSetList);
  chain = MyChainMaker.GetChain();
  //______________________________________//
}


/////////////////////////////////////////////////////////////////////
void EventProcessing::DataProcessing(TString output,TString iso_type)
/////////////////////////////////////////////////////////////////////
{
  //======================================//
  TChain *chain=0;
  ChainMaker MyChainMaker("photon");
  chain = new TChain("photon");
  std::cout << "Build the files chain ..." ;
  MyChainMaker.SetTreeFiles(_DataSetList);
  std::cout << "done !" << std::endl;
  chain = MyChainMaker.GetChain();
  //______________________________________//

  std::cout << "STEP 0" << std::endl;
  EventSelector *ES = new EventSelector(chain);
  std::cout << "STEP 1" << std::endl;
  ES->SetStreamType(true);//data=true
  std::cout << "STEP 2" << std::endl;
  ES->SetGRL(Commons::GrlFile);
  std::cout << "STEP 3" << std::endl;
  ES->SetGeeFile(Commons::GeeFile);
  std::cout << "STEP 3" << std::endl;
  ES->SetIsoType(iso_type);
  std::cout << "STEP 4" << std::endl;
  ES->SetOutputName(output);
  std::cout << "STEP 5" << std::endl;
  ES->InitOutTree();
  std::cout << "STEP 6" << std::endl;
  ES->EventLoop();
  std::cout << "DONE !" << std::endl;
  delete ES;
  std::cout << "Memory has been freed" << std::endl;

}

///////////////////////////////////////////////////////////////////////////////////////////
void EventProcessing::MonteCarloProcessing(TString output,TString iso_type,TString mc_type)
//////////////////////////////////////////////////////////////////////////////////////////
{
  //======================================//
  TChain *chain=0;
  ChainMaker MyChainMaker("photon");
  chain = new TChain("photon");
  std::cout << "Build the files chain ..." ;
  MyChainMaker.SetTreeFiles(_DataSetList);
  std::cout << "done !" << std::endl;
  chain = MyChainMaker.GetChain();
  //______________________________________//

  std::cout << "STEP 0" << std::endl;
  EventSelector *ES = new EventSelector(chain);
  std::cout << "STEP 1" << std::endl;
  ES->SetStreamType(false);//data=true
  std::cout << "STEP 2" << std::endl;
  ES->SetPileUpTool(Commons::PileupMCFile,Commons::PileupDataFile);
  std::cout << "STEP 3" << std::endl;
  ES->SetTruthSelector(chain);
  std::cout << "STEP 4" << std::endl;
  ES->SetIsoType(iso_type);
  std::cout << "STEP 5" << std::endl;
  ES->SetMCType(mc_type);
  std::cout << "STEP 6" << std::endl;
  ES->SetOutputName(output);
  std::cout << "STEP 7" << std::endl;
  ES->InitOutTree();
  std::cout << "STEP 8" << std::endl;
  ES->EventLoop();
  std::cout << "DONE !" << std::endl;
  delete ES;
  std::cout << "Memory has been freed" << std::endl;

}
////////////////////////////////////////////////////////////////////////////
void EventProcessing::RSTemplateProcessing( TString output,TString iso_type,
					    double ptl_c,double ptsl_c,
					    double isol_c,double isosl_c )
////////////////////////////////////////////////////////////////////////////
{
  //======================================//
  TChain *chain=0;
  ChainMaker MyChainMaker("photon");
  chain = new TChain("photon");
  std::cout << "Build the files chain ..." ;
  MyChainMaker.SetTreeFiles(_DataSetList);
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
  STC->SetCoupling(0.1);  

  std::cout << "STEP 8" << std::endl;
  STC->SetLeadPtCut(ptl_c);
  STC->SetSubleadPtCut(ptsl_c);
  STC->SetLeadIsoCut(isol_c);
  STC->SetSubleadIsoCut(isosl_c);

  std::cout << "STEP 9" << std::endl;
  STC->EventLoop();
  std::cout << "STEP 10" << std::endl;
  STC->CreateOutputFile(output);
  std::cout << "DONE !" << std::endl;
  delete STC;
  std::cout << "Memory has been freed" << std::endl;

}








