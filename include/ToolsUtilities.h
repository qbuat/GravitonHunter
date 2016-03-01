//header of OQFlagTranslator
#ifndef ANALYSISTOOLS_H
#define ANALYSISTOOLS_H

#include <TCanvas.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphAsymmErrors.h>

class AnalysisTools
{

 private:


 public:
  AnalysisTools();
  virtual ~AnalysisTools();
  double efficiency(double k, double n);
  double sigma_efficiency(double k, double n);
  double RoundUpDouble(double x,int pre = 1000);
  void   ShowNumberOfProcessedEvents(int jentry,
				     int nentries,
				     int period=10000);
  static void Processing( int jentry,
			  int nentries,
			  int period = 10000);

  static TGraphAsymmErrors* DataGraph(const TH1D& hd);
  static TGraphAsymmErrors* BkgGraph(const TH1D& hb);
  static TH1D* HistAddedUncert(const TH1D& h1,const TH1D& h2);
  static TH1D* GetTruncatedHist(const TH1D& hfull, int nbins,const double* bins);
  static TCanvas* Get2PadPlot(TString name,TString title);


  static TH1D* GetPtggHist(TString name);
  static TH1D* GetAbsDeltaphiHist(TString name);
  static TH1D* GetAbsCosthetastarHist(TString name);
  static TH1D* GetDeltaetaggHist(TString name);
  static TH1D* GetYggHist(TString name);
  static TH1D* GetDeltarHist(TString name);
  static TH1D* GetPtHist(TString name);
  static TH1D* GetEtaHist(TString name);
  static TH1D* GetPhiHist(TString name);
  static TH2D* GetEtaPhiMap(TString name);

};
#endif




