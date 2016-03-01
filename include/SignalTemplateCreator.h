// This class is meant to create the various
// RS signal templates
// It inherits from the GravitonAnalysis class 
#ifndef SignalTemplateCreator_h
#define SignalTemplateCreator_h

#include <iostream>
#include <utility>
#include <math.h>
#include <vector>

#include <TTree.h>
#include <TString.h>
#include <TH2.h>
#include <TFile.h>
#include <TObjArray.h>

#include "EventSelector.h"

class SignalTemplateCreator : public EventSelector
{
 private:

  double              m_coupling;
  int                 m_Nmasses;
  double              m_LowMass;
  double              m_MassSpacing;

  //----------------
  double m_ptl_cut;
  double m_ptsl_cut;
  double m_isol_cut;
  double m_isosl_cut;


  TObjArray*          m_template;
  TObjArray*          m_template_truth;
  TObjArray*          m_hmgg_array;
  TObjArray*          m_hmgg_w_array;
  TObjArray*          m_cutflow_array;
  TObjArray*          m_cutflow_w_array;
  TObjArray*          m_cutflowsum2_w_array;


  std::vector<double> m_vec_mass;

  TH1D**     vec_hCutFlow;
  TH1D**     vec_hCutFlow_w;
  TH1D**     vec_hCutFlowSum2_w;
  TH1D**     vec_hmgg_final;
  TH1D**     vec_hmgg_final_w;
  TH1D**     vec_sighist;
  TH1D**     vec_sighist_truth;



 public :
  SignalTemplateCreator();
  SignalTemplateCreator(TTree*);
  virtual ~SignalTemplateCreator();
  void     Init(TTree*);
  static double   GetGravitonWeight(double truemass,
				    double polemass, 
				    double coupling=0.10 );

  void     EventLoop(double coupling=0.10);
  void     CreateOutputFile(TString output="DEFAULT");


  void SetCoupling(double c) {m_coupling=c;}
  void SetNmasses(int N) {m_Nmasses=N;}
  void SetLowMass(double m) {m_LowMass=m;}
  void SetMassSpacing(double dm) {m_MassSpacing=dm;};

  void SetLeadPtCut(double c) {m_ptl_cut=c;}
  void SetSubleadPtCut(double c) {m_ptsl_cut=c;}
  void SetLeadIsoCut(double c) {m_isol_cut=c;}
  void SetSubleadIsoCut(double c) {m_isosl_cut=c;}


 private :
  void     InitListOfMasses();

};

#endif

