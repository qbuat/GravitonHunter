#ifndef IsolationGrapher_h
#define IsolationGrapher_h

#include <iostream>
#include <vector>
#include <TGraphErrors.h>
#include <TString.h>
#include <TCanvas.h>
#include <TObjArray.h>

#include "IsolationFitter.h"

class IsolationGrapher
{
 private:
  
  int m_Nbins;
  TString m_type;
  std::vector<double> bin_center;
  std::vector<double> bin_width;
  std::vector<double> iso_fitmean;
  std::vector<double> iso_fitmean_error;
  std::vector<double> iso_fitwidth;
  std::vector<double> iso_fitwidth_error;

  std::vector<double> iso_mean;
  std::vector<double> iso_mean_error;
  std::vector<double> iso_rms;
  std::vector<double> iso_rms_error;


  TObjArray*       m_canvas;
  IsolationFitter* m_IF;

 public :

  IsolationGrapher( IsolationFitter* IF, 
		    int Nbins,const double* bins,
		    TString type);
  virtual ~IsolationGrapher();

  void SetNbins(int N) {m_Nbins=N;}
  void SetType(TString t="pT") {m_type=t;}
  void FillGraphs();

  TGraphErrors* GetFitMeanGraph();
  TGraphErrors* GetFitWidthGraph();
  TGraphErrors* GetHistMeanGraph();
  TGraphErrors* GetHistRMSGraph();
  TObjArray*    GetPlots(){return m_canvas;} 

 private :
  void SetBounds(double min, double max);

};

#endif

