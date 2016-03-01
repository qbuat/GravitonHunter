#define AnalysisWraper_cxx
#include "AnalysisWraper.h"
#include "GravitonSelector.h"
#include "BkgSampleSelector.h"
#include "SignalSampleSelector.h"
#include "PhotonEfficiencies.h"
#include "SignalTemplateCreator.h"
#include "ToolsUtilities.h"
#include "ToolsCommons.h"

// #include "GravitonSelector_old.h"
// #include "PhotonEfficiencies_old.h"
#include <TFile.h>
#include <TGraphErrors.h>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////
void AnalysisWraper::MethodSkeleton(TString treename,bool data,TString mctype)
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











///////////////////////////////////////////////////////////////////////////////////////////////////////
void AnalysisWraper::Analysis(TString treename,bool data,TString output,TString iso_type,TString mctype) 
///////////////////////////////////////////////////////////////////////////////////////////////////////
{
  //======================================//
  TChain *chain=0;
  ChainMaker MyChainMaker(treename);
  chain = new TChain(treename);
  std::cout << "Build the files chain ..." ;
  MyChainMaker.SetTreeFiles(_DataSetList);
  std::cout << "done !" << std::endl;
  chain = MyChainMaker.GetChain();
  //______________________________________//

  //==== Invariant mass spectrum =============//
  // GravitonSelector_old *t =new GravitonSelector_old(chain);
  // if(treename=="physics"){
  //   t->InitWZD3PD();
  // }else{
  //   t->InitPhD3PD();
  // }
  // t->EventLoop(data,output,mctype);

  std::cout << "STEP 0" << std::endl;
  GravitonSelector *GS =new GravitonSelector(chain);
  std::cout << "STEP 1" << std::endl;
  GS->SetStreamType(data);
  std::cout << "STEP 2" << std::endl;
  GS->SetIsoType(iso_type);
  std::cout << "STEP 3" << std::endl;
  GS->EventLoop(output);
  std::cout << "DONE !" << std::endl;
  delete GS;
}

///////////////////////////////////////////////////////////////////////////////////////////
void AnalysisWraper::SignalFile(TString treename,bool data,TString output,TString iso_type)
///////////////////////////////////////////////////////////////////////////////////////////
{
  //======================================//
  TChain *chain=0;
  ChainMaker MyChainMaker(treename);
  chain = new TChain(treename);
  std::cout << "Build the files chain ..." ;
  MyChainMaker.SetTreeFiles(_DataSetList);
  std::cout << "done !" << std::endl;
  chain = MyChainMaker.GetChain();
  //______________________________________//

  std::cout << "STEP 0" << std::endl;
  SignalSampleSelector *SSS = new SignalSampleSelector(chain);
  std::cout << "STEP 1" << std::endl;
  SSS->SetStreamType(data);
  std::cout << "STEP 2" << std::endl;
  SSS->SetTruthSelector(chain);
  std::cout << "STEP 3" << std::endl;
  SSS->SetIsoType(iso_type);
  std::cout << "STEP 4" << std::endl;
  SSS->EventLoop(output);
  std::cout << "DONE !" << std::endl;
  delete SSS;
  std::cout << "TOTO" << std::endl;
}
////////////////////////////////////////////////////////////////////////////////
void AnalysisWraper::Signal_RS(TString treename,TString output,TString iso_type)
////////////////////////////////////////////////////////////////////////////////
{
  //======================================//
  TChain *chain=0;
  ChainMaker MyChainMaker(treename);
  chain = new TChain(treename);
  std::cout << "Build the files chain ..." ;
  MyChainMaker.SetTreeFiles(_DataSetList);
  std::cout << "done !" << std::endl;
  chain = MyChainMaker.GetChain();
  //______________________________________//

  std::cout << "STEP 0" << std::endl;
  SignalTemplateCreator *STC = new SignalTemplateCreator(chain);
  std::cout << "STEP 1" << std::endl;
  STC->SetStreamType(false);
  std::cout << "STEP 2" << std::endl;
  STC->SetTruthSelector(chain);
  std::cout << "STEP 3" << std::endl;
  STC->SetIsoType(iso_type);
  std::cout << "STEP 4" << std::endl;
  STC->EventLoop(0.1);
  std::cout << "STEP 5" << std::endl;
  STC->CreateOutputFile(output);
  std::cout << "DONE !" << std::endl;


  delete STC;
  std::cout << "TOTO" << std::endl;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void AnalysisWraper::BkgFile(TString treename,bool data,TString output,TString iso_type,int Nrelaxedcut)
////////////////////////////////////////////////////////////////////////////////////////////////////////
{
  //======================================//
  TChain *chain=0;
  ChainMaker MyChainMaker(treename);
  chain = new TChain(treename);
  std::cout << "Build the files chain ..." ;
  MyChainMaker.SetTreeFiles(_DataSetList);
  std::cout << "done !" << std::endl;
  chain = MyChainMaker.GetChain();
  //______________________________________//

  std::cout << "STEP 0" << std::endl;
  BkgSampleSelector *BS = new BkgSampleSelector(chain);
  std::cout << "STEP 1" << std::endl;
  BS->SetStreamType(data);
  std::cout << "STEP 2" << std::endl;
  BS->SetTruthSelector(chain);
  std::cout << "STEP 3" << std::endl;
  BS->SetIsoType(iso_type);
  std::cout << "STEP 4" << std::endl;
  BS->EventLoop(output,Nrelaxedcut);
  std::cout << "DONE !" << std::endl;
  delete BS;
}
//////////////////////////////////////////////////////////////////
void AnalysisWraper::BkgTruthFile(TString treename,TString output)
//////////////////////////////////////////////////////////////////
{
  //======================================//
  TChain *chain=0;
  ChainMaker MyChainMaker(treename);
  chain = new TChain(treename);
  std::cout << "Build the files chain ..." ;
  MyChainMaker.SetTreeFiles(_DataSetList);
  std::cout << "done !" << std::endl;
  chain = MyChainMaker.GetChain();
  //______________________________________//

  std::cout << "STEP 0" << std::endl;
  BkgSampleSelector *BS = new BkgSampleSelector(chain);
  std::cout << "STEP 1" << std::endl;
  BS->SetStreamType(false);
  std::cout << "STEP 2" << std::endl;
  BS->SetTruthSelector(chain);
  std::cout << "STEP 3" << std::endl;
  BS->SetIsoType("CONE");
  std::cout << "STEP 4" << std::endl;
  BS->TruthEventLoop(output);
  std::cout << "DONE !" << std::endl;
  delete BS;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
void AnalysisWraper::StoreEffMaps(TString treename,TString iso_type,TString mc_type,TString output)
///////////////////////////////////////////////////////////////////////////////////////////////////
{
  //======================================//
  TChain *chain=0;
  ChainMaker MyChainMaker(treename);
  chain = new TChain(treename);
  std::cout << "Build the files chain ..." ;
  MyChainMaker.SetTreeFiles(_DataSetList);
  std::cout << "done !" << std::endl;
  chain = MyChainMaker.GetChain();
  //______________________________________//



  std::cout << "STEP 0" << std::endl;
  PhotonEfficiencies *PE =new PhotonEfficiencies(chain);
  std::cout << "STEP 1" << std::endl;
  PE->SetStreamType(false);// Run over mc
  std::cout << "STEP 2" << std::endl;
  PE->SetIsoType(iso_type);
  std::cout << "STEP 3" << std::endl;
  PE->SetMonteCarloSample(mc_type); 
  PE->SetPileUpTool(Commons::PileupMCFile,Commons::PileupDataFile);
  PE->SetTruthSelector(chain);
  PE->EventLoop(mc_type); 
 std::cout << "STEP 4" << std::endl;
  PE->CreateOutputFile(output);
  std::cout << "DONE !" << std::endl;

}
/////////////////////////////////////////////////////////////////////////////////////
void AnalysisWraper::StorePhotonEff(TString treename,bool data,TString output,TString mctype)
////////////////////////////////////////////////////////////////////////////////////
{ 
  //======================================//
  TChain *chain=0;
  ChainMaker MyChainMaker(treename);
  chain = new TChain(treename);
  std::cout << "Build the files chain ..." ;
  MyChainMaker.SetTreeFiles(_DataSetList);
  std::cout << "done !" << std::endl;
  chain = MyChainMaker.GetChain();
  //______________________________________//


}
// ////////////////////////////////////////////////////////////////////////////
// void AnalysisWraper::CreateSignalTemplate(TString treename,TString output,
// 					  TString iso_type,double coupling)
// ///////////////////////////////////////////////////////////////////////////
// {
//   //======================================//
//   TChain *chain=0;
//   ChainMaker MyChainMaker(treename);
//   chain = new TChain(treename);
//   std::cout << "Build the files chain ..." ;
//   MyChainMaker.SetTreeFiles(_DataSetList);
//   std::cout << "done !" << std::endl;
//   chain = MyChainMaker.GetChain();
//   //______________________________________//

//   SignalTemplateCreator *STC;
//   STC = new SignalTemplateCreator(chain);
//   STC->SetIsoType(iso_type);
//   STC->EventLoop(coupling);
//   STC->CreateTemplateFile(output);
//   delete STC;

// }
