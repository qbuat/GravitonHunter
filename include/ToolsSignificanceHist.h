//--> Compute the significance
// between two histograms
// using a (D-B)/errB method
// or following arXiv:1111.2062v2

#ifndef SignificanceHist_h
#define SignificanceHist_h

#include <TH2.h>
#include <RooCurve.h>
#include <RooHist.h>
#include <TGraphAsymmErrors.h>
#include <iostream>
class SignificanceHist
{
 private:
  TH1F* m_hh_obs;
  TH1F* m_hh_exp; 

 public :

  SignificanceHist( const TH1F& hh_obs, 
		    const TH1F& hh_exp );
  SignificanceHist( const RooHist& hh_obs,
		    const RooCurve& hh_exp );
  SignificanceHist( const RooHist& hh_obs,
		    const RooCurve& hh_exp,
		    const RooCurve& hh_exp_syst );

  virtual ~SignificanceHist();

  TH1F* GetSignificanceHist(double sigma=2.5);
  TH1F* GetSignificanceHist_withsyst(double sigma=2.5);
  TH1F* GetChiHist(double sigma=2.5);
  TH1F* GetRatioHist(double deviation=0.5);
  TGraphAsymmErrors* GetDiffHist();
};

#endif

