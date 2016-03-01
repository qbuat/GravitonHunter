#ifndef PlottingBackground_h
#define PlottingEfficiencies_h

#include <iostream>
#include <vector>
#include <TCanvas.h>
#include <TString.h>
#include <TFile.h>
#include <TH2.h>
#include <TLatex.h>
#include <TGraphErrors.h>

class PlottingEfficiencies{

 private:
  TFile* m_filemaps;
  std::vector<double> m_pt_binedges;
  TLatex* m_lat;

 public:
  
  PlottingEfficiencies( TString filemaps="NONE");
  virtual ~PlottingEfficiencies();
  
  void InitMaps(TString filemaps);
  //--------------------------------
  void SetPtBinning(int Nbins,double *edges);
  TCanvas * RecoTightEfficiencyPlot(int etabin);
 private:
  std::vector<TGraphErrors*> BuildPtEfficiencyGraphs( const TH2 & m_den,
						      const TH2 & m_num );

};
#endif // #ifdef PlottingEfficiencies_h
