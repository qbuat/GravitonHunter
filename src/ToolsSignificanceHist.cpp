#include "ToolsSignificanceHist.h"
#include <TMath.h>
#include <Math/SpecFuncMathCore.h>
#include <TError.h>

//////////////////////////////////////////////////////////////////////////
SignificanceHist::SignificanceHist(const TH1F& hh_obs, const TH1F& hh_exp) 
//////////////////////////////////////////////////////////////////////////
{
  m_hh_obs = (TH1F*)hh_obs.Clone(hh_obs.GetName());
  m_hh_exp = (TH1F*)hh_exp.Clone(hh_exp.GetName());
}
/////////////////////////////////////////////////////////////////////////////////
SignificanceHist::SignificanceHist(const RooHist& hh_obs, const RooCurve& hh_exp) 
/////////////////////////////////////////////////////////////////////////////////
{
  int Nbins = hh_obs.GetN();
  double * x_obs_low = hh_obs.GetEXlow();
  double * x_obs = hh_obs.GetX();
  double * y_obs = hh_obs.GetY(); 
  double x_bins[Nbins+1];
  for(int ibin=0;ibin<Nbins;ibin++) {
    x_bins[ibin]= x_obs[ibin]-x_obs_low[ibin];
  }
  x_bins[Nbins] = x_obs[Nbins-1]+hh_obs.GetErrorXhigh(Nbins-1);

  m_hh_obs = new TH1F(hh_obs.GetName(),hh_obs.GetTitle(),Nbins,x_bins);
  m_hh_exp = new TH1F(hh_exp.GetName(),hh_exp.GetTitle(),Nbins,x_bins);
  for( int ibin=1;ibin<=Nbins;ibin++){
    m_hh_obs->SetBinContent(ibin,y_obs[ibin-1]);
    double low = m_hh_exp->GetBinLowEdge(ibin);
    double up  = m_hh_exp->GetBinLowEdge(ibin+1);
    double exp_int = hh_exp.average(low,up);
    m_hh_exp->SetBinContent(ibin,exp_int);
  }
}
/////////////////////////////////////////////////////////////////////////////////
SignificanceHist::SignificanceHist( const RooHist& hh_obs, const RooCurve& hh_exp, 
				    const RooCurve& hh_exp_syst ) 
/////////////////////////////////////////////////////////////////////////////////
{
  int Nbins = hh_obs.GetN();
  double * x_obs_low = hh_obs.GetEXlow();
  double * x_obs = hh_obs.GetX();
  double * y_obs = hh_obs.GetY(); 
  double x_bins[Nbins+1];
  for(int ibin=0;ibin<Nbins;ibin++)
    x_bins[ibin]= x_obs[ibin]-x_obs_low[ibin];

  x_bins[Nbins] = x_obs[Nbins-1]+hh_obs.GetErrorXhigh(Nbins-1);
  m_hh_obs = new TH1F(hh_obs.GetName(),hh_obs.GetTitle(),Nbins,x_bins);
  m_hh_exp = new TH1F(hh_exp.GetName(),hh_exp.GetTitle(),Nbins,x_bins);
  for( int ibin=1;ibin<=Nbins;ibin++){
    m_hh_obs->SetBinContent(ibin,y_obs[ibin-1]);
    double low = m_hh_exp->GetBinLowEdge(ibin);
    double up  = m_hh_exp->GetBinLowEdge(ibin+1);
    double exp_int = hh_exp.average(low,up);
    double exp_syst_int = hh_exp_syst.average(low,up);
    m_hh_exp->SetBinContent(ibin,exp_int);
    m_hh_exp->SetBinError(ibin,fabs(exp_int-exp_syst_int));
  }
}
/////////////////////////////////////
SignificanceHist::~SignificanceHist()
/////////////////////////////////////
{
  //Default destructor
  delete m_hh_obs;
  delete m_hh_exp;
}
/////////////////////////////////////////////////////////
TH1F* SignificanceHist::GetSignificanceHist(double sigma)
/////////////////////////////////////////////////////////
{

  //---> from arXiv:1111.2062v2 
  TH1F* hh_signi = (TH1F*)m_hh_obs->Clone("hh_signi");
  for(int ibin=1 ;ibin<= hh_signi->GetNbinsX();ibin++){
    double obs  = m_hh_obs->GetBinContent(ibin);
    double exp  = m_hh_exp->GetBinContent(ibin);
    double sig  = -9999999;
    if( obs>exp ){
      double Q = ROOT::Math::inc_gamma_c(obs,exp);
      double pval = 1-Q;
      if(pval>0.5)
  	sig = 0;
      else
  	sig = sqrt(2.)*TMath::ErfInverse(1.-2*pval);
    }else{
      double Q = ROOT::Math::inc_gamma_c(obs+1,exp);
      double pval = Q;
      if(pval>0.5)
  	sig = 0;
      else
  	sig = - sqrt(2.)*TMath::ErfInverse(1.-2*pval);
    }
    hh_signi->SetBinContent(ibin,sig);
    hh_signi->SetBinError(ibin,0);
  }
  hh_signi->GetYaxis()->SetTitle("Significance");
  hh_signi->GetYaxis()->SetRangeUser(-1*sigma,1*sigma);
  hh_signi->SetFillColor(2);//convention for the stat only error

  //--> Some drawing commands to fit a low panel plot
  //--> works fine with AnalysisTools::Get2PadPlot()
  hh_signi->SetTitleSize(0.11,"Y");
  hh_signi->SetTitleOffset(0.60,"Y");
  hh_signi->SetTitleSize(0.11,"X");
  hh_signi->SetLabelSize(0.11,"X");
  hh_signi->SetLabelSize(0.11,"Y");
  hh_signi->GetYaxis()->SetNdivisions(6);

  return hh_signi;
}
///////////////////////////////////////////////////////////////////
TH1F* SignificanceHist::GetSignificanceHist_withsyst(double sigma)
//////////////////////////////////////////////////////////////////
{
  //---> Need the systematic uncertainty quoted as the bin error.
  TH1F* hh_signi = (TH1F*)m_hh_obs->Clone("hh_signi");

  //---> from arXiv:1111.2062v2 
  for(int ibin=1;ibin<=hh_signi->GetNbinsX();ibin++){
    int    D = m_hh_obs->GetBinContent(ibin);
    double B = m_hh_exp->GetBinContent(ibin);
    double S = m_hh_exp->GetBinError(ibin);
    double a = B*B/(S*S);
    double b = B/(S*S);
    double sig = -9999999;
    if( D>(int)B ){
      double pval_temp = 0;
      for(int n=0;n<D;n++){
        double temp = TMath::Power(b,a)/ROOT::Math::tgamma(a);
	temp *= ROOT::Math::tgamma(n+a)/(TMath::Factorial(n)*TMath::Power(1+b,n+a));
        pval_temp += temp;
      }
      double pval = 1-pval_temp;
      if(pval>0.5)
        sig = 0;
      else
        sig = sqrt(2.)*TMath::ErfInverse(1.-2*pval);
    }else{
      double pval = 0;
      for(int n=0;n<D+1;n++){
        double temp = TMath::Power(b,a)/ROOT::Math::tgamma(a);
	temp *= ROOT::Math::tgamma(n+a)/(TMath::Factorial(n)*TMath::Power(1+b,n+a));
        pval += temp;
      }
      if(pval>0.5)
        sig = 0;
      else
        sig = -sqrt(2.)*TMath::ErfInverse(1.-2*pval);
    }
    hh_signi->SetBinContent(ibin,sig);
    hh_signi->SetBinError(ibin,0);
  }
    
  hh_signi->GetYaxis()->SetTitle("Significance");
  hh_signi->GetYaxis()->SetRangeUser(-1*sigma,1*sigma);
  hh_signi->SetFillColor(4);//convention for the syst uncert
  
  //--> Some drawing commands to fit a low panel plot
  //--> works fine with AnalysisTools::Get2PadPlot()
  hh_signi->SetTitleSize(0.11,"Y");
  hh_signi->SetTitleOffset(0.60,"Y");
  hh_signi->SetTitleSize(0.11,"X");
  hh_signi->SetLabelSize(0.11,"X");
  hh_signi->SetLabelSize(0.11,"Y");
  hh_signi->GetYaxis()->SetNdivisions(6);
  
  return hh_signi;
}

////////////////////////////////////////////////
TH1F* SignificanceHist::GetChiHist(double sigma)
////////////////////////////////////////////////
{
  TH1F* hh_chi =  (TH1F*)m_hh_obs->Clone("hh_chi");
  for(int ibin=1 ;ibin<= hh_chi->GetNbinsX();ibin++){
    double obs       = m_hh_obs->GetBinContent(ibin);
    double exp       = m_hh_exp->GetBinContent(ibin);
    double sigma_exp = sqrt(exp);
    double chi       = (obs-exp)/sigma_exp;
    if( sigma_exp !=0 )
      hh_chi->SetBinContent(ibin,chi);
    else 
      hh_chi->SetBinContent(ibin,0);
    hh_chi->SetBinError(ibin,0);
  }
  hh_chi->GetYaxis()->SetTitle("#chi");
  hh_chi->GetYaxis()->SetRangeUser(-1*sigma,1*sigma);
  //--> Some drawing commands to fit a low panel plot
  //--> works fine with AnalysisTools::Get2PadPlot()
  hh_chi->SetFillColor(2);
  hh_chi->SetTitleSize(0.11,"Y");
  hh_chi->SetTitleOffset(0.60,"Y");
  hh_chi->SetTitleSize(0.11,"X");
  hh_chi->SetLabelSize(0.11,"X");
  hh_chi->SetLabelSize(0.11,"Y");
  hh_chi->GetYaxis()->SetNdivisions(6);

  return hh_chi;
}

//////////////////////////////////////////////////////
TH1F* SignificanceHist::GetRatioHist(double deviation)
//////////////////////////////////////////////////////
{

  TH1F* hh_ratio =  (TH1F*)m_hh_obs->Clone("hh_ratio");
  hh_ratio->Sumw2();
  hh_ratio->Divide(m_hh_exp);
  hh_ratio->GetYaxis()->SetTitle("Data/Background");
  hh_ratio->GetYaxis()->SetRangeUser(1-deviation,1+deviation);
  //--> Some drawing commands to fit a low panel plot
  //--> works fine with AnalysisTools::Get2PadPlot()
  hh_ratio->SetTitleSize(0.11,"Y");
  hh_ratio->SetTitleOffset(0.60,"Y");
  hh_ratio->SetTitleSize(0.11,"X");
  hh_ratio->SetLabelSize(0.11,"X");
  hh_ratio->SetLabelSize(0.11,"Y");
  hh_ratio->GetYaxis()->SetNdivisions(6);
  return hh_ratio;

}

/////////////////////////////////////////////
TGraphAsymmErrors* SignificanceHist::GetDiffHist()
/////////////////////////////////////////////
{
  TGraphAsymmErrors* gr = new TGraphAsymmErrors();
  gr->SetName("hh_difference");
  for(int ibin=0 ;ibin<= m_hh_obs->GetNbinsX();ibin++){
    double obs       = m_hh_obs->GetBinContent(ibin);
    double exp       = m_hh_exp->GetBinContent(ibin);

    double bin_width = m_hh_exp->GetBinLowEdge(ibin+1) - m_hh_exp->GetBinLowEdge(ibin);
    double bin_center = m_hh_exp->GetBinLowEdge(ibin) + bin_width/2;
    gr->SetPoint(ibin-1,bin_center,obs-exp);
    gr->SetPointError(ibin,bin_center-m_hh_exp->GetBinLowEdge(ibin),
		      m_hh_exp->GetBinLowEdge(ibin+1)-bin_center,
		      sqrt(obs),sqrt(obs));
  }
  gr->GetXaxis()->SetRangeUser(m_hh_obs->GetBinLowEdge(1),m_hh_obs->GetBinLowEdge(m_hh_obs->GetNbinsX()+1));
  gr->GetYaxis()->SetTitle("Data-Expectation");
  gr->GetYaxis()->SetNdivisions(6);
  return gr;
}
