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
#ifndef EventProcessing_h
#define EventProcessing_h

#include <TChain.h>
#include <TString.h>
#include <vector>
#include <utility>
#include <iostream>
#include <stdlib.h>


class EventProcessing {

 public :
  EventProcessing();
  virtual ~EventProcessing();
  virtual void MethodSkeleton( TString treename,
			       bool data,
			       TString mctype="" );
  virtual void DataProcessing(TString output,TString iso_type);
  virtual void MonteCarloProcessing(TString output,TString iso_type,TString mc_type);
  virtual void RSTemplateProcessing(TString output,TString iso_type,
				    double ptl_c,double ptsl_c,
				    double isol_c,double isosl_c);

  void SetDataSetList(std::string list){_DataSetList=list;}

 private:
  std::string     _DataSetList;
  
};

#endif

