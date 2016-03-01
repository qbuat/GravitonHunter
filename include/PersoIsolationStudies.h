// This class inherits from Graviton Analysis.
// It is meant to study the isolation cut optimisation
// on an asymetric way between leader and subleader photons.
// The systematic study of such a cut has not been studied yet.
#ifndef IsolationStudies_h
#define IsolationStudies_h

#include <iostream>
#include <vector>
#include <utility>
#include "TTree.h"
#include "TString.h"
#include "TGraphErrors.h"
#include "ToolsTreeReader.h"
#include "GravitonAnalysis.h"
class IsolationStudies : public GravitonAnalysis
{
 private:

  int m_Npt;
  std::vector<double> m_bin_pt;
  std::vector<double> m_err_bin_pt;
  std::vector<double> Nph_isotight_pt;
  std::vector<double> Nph_tight_pt;
  bool m_passTrigger;
  bool m_inGRL;
  bool m_passPV;
  bool m_passPreSel;

 public :
  IsolationStudies();
  IsolationStudies(TTree*);
  virtual ~IsolationStudies();
  void          Init(TTree*);
  void          CutOptimisation(bool data);
  void          IsoPerNPVbins(bool,TString output,
			      TString mctype="");

};

#endif

