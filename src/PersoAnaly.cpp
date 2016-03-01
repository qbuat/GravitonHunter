#define AnalyPerso_cxx
#include "PersoAnaly.h"
#include "PersoIsolationStudies.h"
#include "PersoCleaningStudies.h"
#include "PersoSinglePhotonSelector.h"
#include "ToolsPileUpToolHandler.h"
#include "ToolsCommons.h"

#include <TFile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <cmath>

//////////////////////////////////////////////////////////////////////
void AnalyPerso::MethodSkeleton(TString treename,bool data,TString mctype)
/////////////////////////////////////////////////////////////////////
{
  //======================================//
  TChain *chain=0;
  ChainMaker MyChainMaker(treename);
  chain = new TChain(treename);
  MyChainMaker.SetTreeFiles(_DataSetList);
  chain = MyChainMaker.GetChain();
  //______________________________________//
}

///////////////////////////////////////////////////////////////////////////////
void AnalyPerso::Isolation(TString treename,bool data,TString output,TString mctype)
///////////////////////////////////////////////////////////////////////////////
{
  //======================================//
  TChain *chain=0;
  ChainMaker MyChainMaker(treename);
  chain = new TChain(treename);
  MyChainMaker.SetTreeFiles(_DataSetList);
  chain = MyChainMaker.GetChain();
  //______________________________________//

  //== Isolation Studies ==//
   // IsolationStudies *IS=new IsolationStudies(chain);

  //    IS->CutOptimisation(data);
  //    TGraphErrors *gr5=IS->IsovsTightEff(data,5.);
  //    TString graphname5  = gr5->GetName();
  //    TCanvas *c1=new TCanvas("c1","Iso(5GeV)/Tight",800,800);
  //    c1->cd();
  //    gr5->GetYaxis()->SetRangeUser(60.,90.);
  //    gr5->Draw("AP");
  //    c1->SaveAs("../plots/"+graphname5+".png");
  //    c1->SaveAs("../plots/"+graphname5+".eps");
  
  //    TGraphErrors *gr10=IS->IsovsTightEff(data,10.);
  //    TString graphname10 = gr10->GetName();
  //    TCanvas *c2=new TCanvas("c2","Iso(10GeV)/Tight",800,800);
  //    c2->cd();
  //    gr10->GetYaxis()->SetRangeUser(60.,90.);
  //    gr10->Draw("AP");
  //    c2->SaveAs("../plots/"+graphname10+".png");
  //    c2->SaveAs("../plots/"+graphname10+".eps");
  
  
}
///////////////////////////////////////////////////////////////////
void AnalyPerso::FudgeFactorsEffCalc(TString treename,TString graphname)
///////////////////////////////////////////////////////////////////
{
  //======================================//
  TChain *chain=0;
  ChainMaker MyChainMaker(treename);
  chain = new TChain(treename);
  MyChainMaker.SetTreeFiles(_DataSetList);
  chain = MyChainMaker.GetChain();
  //______________________________________//

  // FudgeFactorsEff *FFE=new FudgeFactorsEff(chain);
  // FFE->InitPhD3PD();
  // FFE->PhotonsCounting();
  // FFE->ShowEfficiencies();
  // TGraphErrors* gr=FFE->TightRecoEff_PU();
  // gr->SetName(graphname);
  // // gr->Draw("AP");
  // TObjArray Hlist(0);
  // Hlist.Add(gr);
  // TFile *f=new TFile("PileUpSys.root","update");
  // Hlist.Write();
  // f->Close();
  // delete f;
  // delete gr;
  // pair<double,double> eff_FF = FFE->EventLevelEff();
  // TFile *fIso=new TFile("IsolCutValue.root","update");
  // TObjArray IsoList(0);
  // TGraphErrors * grIso=FFE->EventLevelEff_Iso();
  // grIso->SetName(graphname);
  // IsoList.Add(grIso);
  // IsoList.Write();
  // fIso->Close();
  // delete fIso;
  // delete grIso;
  // TFile *fPU=new TFile("PileUpSys_ADD.root","update");
  // TObjArray PUList(0);
  // TGraphErrors * grPU=FFE->EventLevelEff_PU();
  // grPU->SetName(graphname);
  // PUList.Add(grPU);
  // PUList.Write();
  // fPU->Close();
  // delete fPU;
  // delete grPU;
  // //Efficiency without FudgeFactors
  // std::cout <<"Efficiency without FudgeFactors" << std::endl;
  // bool UseFF=false;
  // pair<double,double> eff_noFF = FFE->EventLevelEff(UseFF);
  // std::cout << " SYS FF (%) = " 
  // 	    << 100*(eff_FF.first-eff_noFF.first) 
  // 	    << " +/- " 
  // 	    << 100*sqrt(eff_FF.second*eff_FF.second+eff_noFF.second*eff_noFF.second) 
  // 	    << std::endl;

}

///////////////////////////////////////////////////////////////
void AnalyPerso::Cleaning(TString treename,bool data,TString output)
///////////////////////////////////////////////////////////////
{
  //==== Method to study photons 
  //==== from gravitons (high pt photons) ==========//

  //======================================//
  TChain *chain=0;
  ChainMaker MyChainMaker(treename);
  chain = new TChain(treename);
  MyChainMaker.SetTreeFiles(_DataSetList);
  chain = MyChainMaker.GetChain();
  //______________________________________//
  CleaningStudies *t =new CleaningStudies(chain);
  // t->InitPhD3PD();
   // t->InvMassCleaning(data,"invmasscleaning.root");
   t->PhotonStudies(data,output);
}

////////////////////////////////////////////////////////////////////////////////////////////
void AnalyPerso::SinglePhotonIsoFiles(TString treename,bool data,TString iso_type,TString mctype)
///////////////////////////////////////////////////////////////////////////////////////////
{
  //======================================//
  TChain *chain=0;
  ChainMaker MyChainMaker(treename);
  chain = new TChain(treename);
  MyChainMaker.SetTreeFiles(_DataSetList);
  chain = MyChainMaker.GetChain();
  //______________________________________//

   SinglePhotonSelector *IS=new SinglePhotonSelector(chain);
   // IS->InitPhD3PD();
   IS->SetIsoType(iso_type);
   if(data){
     IS->EventLoop(data,"/sps/atlas/q/qbuat/OUTDS/GravitonAnalysis/SinglePhotonTight.root",
		   0,mctype);
     IS->EventLoop(data,"/sps/atlas/q/qbuat/OUTDS/GravitonAnalysis/SinglePhotonLooseprime2.root",
		   2,mctype);
     IS->EventLoop(data,"/sps/atlas/q/qbuat/OUTDS/GravitonAnalysis/SinglePhotonLooseprime3.root",
		   3,mctype);
     IS->EventLoop(data,"/sps/atlas/q/qbuat/OUTDS/GravitonAnalysis/SinglePhotonLooseprime4.root",
		   4,mctype);
     IS->EventLoop(data,"/sps/atlas/q/qbuat/OUTDS/GravitonAnalysis/SinglePhotonLooseprime5.root",
		   5,mctype);
   }else{
     IS->EventLoop(data,"/sps/atlas/q/qbuat/OUTDS/GravitonAnalysis/SinglePhotonDiphoton.root",
		   0,mctype);
   }
   delete IS;
}

//////////////////////////////////////////////////////////////////////////
void AnalyPerso::GeneratePileUpMC(TString treename,TString MC_PU_filename)
/////////////////////////////////////////////////////////////////////////
{
  //======================================//
  TChain *chain=0;
  ChainMaker MyChainMaker(treename);
  chain = new TChain(treename);
  MyChainMaker.SetTreeFiles(_DataSetList);
  chain = MyChainMaker.GetChain();
  //______________________________________//

// // == PileUp Reweighting == //
  PileUpToolHandler* MPU=new PileUpToolHandler(chain);
  MPU->GenerateMCFile(MC_PU_filename);

}
//////////////////////////////////////////////////
TCanvas* AnalyPerso::CheckPileUpMC(TString treename)
//////////////////////////////////////////////////
{
  //======================================//
  TChain *chain=0;
  ChainMaker MyChainMaker(treename);
  chain = new TChain(treename);
  MyChainMaker.SetTreeFiles(_DataSetList);
  chain = MyChainMaker.GetChain();
  //______________________________________//

// // == PileUp Reweighting == //
  PileUpToolHandler* MPU=new PileUpToolHandler(chain);

  MPU->SetMCFile(Commons::PileupMCFile);
  MPU->SetDataFile(Commons::PileupDataFile);
  TCanvas* c = MPU->CheckPileUpWeights();
  return c;
}
