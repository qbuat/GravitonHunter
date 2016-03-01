// This class is used as a wrapper 
// to run the analysis.
// One need to declare it in 
// the main algorithm 
// The input variable are the :
// (1) treename, used to chain the dataset files.
// (2) a boolean true if running on data, false for MC
// (3) a name of output rootfile
// One must start by setting the input 
// string with the SetDataSetList method
// (See for example RunData.C)
#ifndef AnalysisWraper_h
#define AnalysisWraper_h

#include <TChain.h>
#include <TString.h>
#include <vector>
#include <utility>
#include <iostream>
#include <stdlib.h>

#include "ToolsChainMaker.h"
class AnalysisWraper {

 public :
  AnalysisWraper();
  virtual ~AnalysisWraper();
  virtual void MethodSkeleton( TString treename,
			       bool data,
			       TString mctype="" );
  virtual void Analysis( TString treename,
			 bool data,
			 TString output,
			 TString iso_type= "CONE",
			 TString mctype="" );

  virtual void SignalFile( TString treename,
			   bool data,
			   TString output,
			   TString iso_type );
  virtual void Signal_RS( TString treename,
			  TString output,
			  TString iso_type );

  virtual void BkgFile( TString treename,
			bool data,
			TString output,
			TString iso_type,
			int Nrelaxedcut);
  virtual void BkgTruthFile(TString treename,
			    TString output);
  virtual void StoreEffMaps( TString treename,
			     TString iso_type,
			     TString mc_type,
			     TString output );

  virtual void StorePhotonEff( TString treename,
			       bool data,
			       TString output,
			       TString mctype="");
  /* virtual void CreateSignalTemplate( TString treename, */
  /* 				     TString output="temp.root", */
  /* 				     TString iso_type ="CONE", */
  /* 				     double coupling=0.1); */
  
  void SetDataSetList(std::string list){_DataSetList=list;}

 private:
  std::string     _DataSetList;
  
};

#endif

#ifdef AnalysisWraper_cxx
AnalysisWraper::AnalysisWraper()
{
}
AnalysisWraper::~AnalysisWraper()
{
}
#endif
