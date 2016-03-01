#ifndef SignalTemplateFitter_h
#define SignalTemplateFitter_h

#include <iostream>
#include <TTree.h>
#include <TCanvas.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooBreitWigner.h>
#include <RooFitResult.h>
#include <RooPlot.h>

class SignalTemplateFitter 
{

 private:
  TTree          *m_tree;
  RooDataSet     *m_DataSet;
  RooBreitWigner *m_BW_pdf;
  RooRealVar     *m_mass;
  RooRealVar     *m_var_width;
  RooFitResult   *m_FitRes_BW;
  RooPlot        *m_frame_BW;
  double          m_Xmin;
  double          m_Xmax;
  double          m_polemass;
  double          m_coupling;

 public :

  SignalTemplateFitter();
  virtual ~SignalTemplateFitter();
  void SetTree(TTree* tree);
  void SetPoleMass(double m = 1000.){m_polemass = m;}
  void SetCoupling(double c = 0.1){m_coupling = c;}
  void SetFitRange(double Xmin=0,
		   double Xmax = 1000);
  void BW_Fitter(bool doFit=true);
  double GetFittedWidth(){return m_var_width->getVal();}

  TCanvas* PlotBWFit();

  void PlotBWFit_old(bool doSave=true,
		     TString name ="");
  
};

#endif

