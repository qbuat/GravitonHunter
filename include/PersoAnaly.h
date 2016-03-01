#ifndef AnalyPerso_h
#define AnalyPerso_h

#include <TChain.h>
#include <TString.h>
#include <TCanvas.h>
#include <vector>
#include <utility>
#include <iostream>
#include <stdlib.h>

#include "ToolsChainMaker.h"
class AnalyPerso {

 public :
  AnalyPerso();
  virtual ~AnalyPerso();

  virtual void MethodSkeleton( TString treename,
			       bool data,
			       TString mctype="" );

  virtual void GeneratePileUpMC( TString treename,
				 TString MC_PU_filename );
  TCanvas* CheckPileUpMC( TString treename );

  virtual void Isolation( TString treename,
			  bool data,
			  TString output,
			  TString mctype );

  virtual void FudgeFactorsEffCalc( TString treename,
				    TString graphname );

  virtual void Cleaning( TString treename,
			 bool data,
			 TString output );
  
  virtual void SinglePhotonIsoFiles( TString treename,
				     bool data,
				     TString iso_type,
				     TString mctype);

  void SetDataSetList(std::string list){_DataSetList=list;}

 private:
  std::string     _DataSetList;
  
};

#endif

#ifdef AnalyPerso_cxx
AnalyPerso::AnalyPerso()
{
}
AnalyPerso::~AnalyPerso()
{
}
#endif
