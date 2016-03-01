#include "BkgUncertaintyBuilder.h"
#include "ToolsCommons.h"

///////////////////////////////////////////////////////////////////////////////
BkgUncertaintyBuilder::BkgUncertaintyBuilder(const TH1D& hnom,const TH1D& hsys)
//////////////////////////////////////////////////////////////////////////////
{
  SetNormRange(Commons::norm_bin.first,Commons::norm_bin.second);
  SetHists(hnom,hsys);
  Builder();
}

////////////////////////////////////////////////
BkgUncertaintyBuilder::~BkgUncertaintyBuilder()
////////////////////////////////////////////////
{
  delete m_hsys;
  delete m_hnom;
  delete m_hratio;
  delete m_hbkgfrac;
}
///////////////////////////////////////////////////////////////////////
void BkgUncertaintyBuilder::SetHists(const TH1D& hnom,const TH1D& hsys)
///////////////////////////////////////////////////////////////////////
{
  m_hsys = (TH1D*)hsys.Clone(hsys.GetName());
  m_hnom = (TH1D*)hnom.Clone(hnom.GetName());
}
//////////////////////////////////////////////////////////
void BkgUncertaintyBuilder::SetNormRange(int min,int max)
/////////////////////////////////////////////////////////
{ m_binmin = min; m_binmax = max;}

/////////////////////////////////////
void BkgUncertaintyBuilder::Builder()
/////////////////////////////////////
{

  //--> Ratio histogram
  m_hratio = (TH1D*)m_hsys->Clone("hh_ratio");
  m_hratio->Divide(m_hnom);
  m_hratio->GetYaxis()->SetTitle("Ratio to nominal");

  int Nbins = m_hratio->GetNbinsX();

  //--> Relative uncertainty histogram
  m_hreluncert = (TH1D*)m_hsys->Clone("hh_reluncert");
  m_hreluncert->GetYaxis()->SetTitle("Relative Uncertainty");

  //--> |sys-nom|/nom histogram
  m_hbkgfrac = (TH1D*)m_hsys->Clone("hh_bkgfrac");
  m_hbkgfrac->GetYaxis()->SetTitle("Uncertainty as a Fraction of Nominal");
  for( int ibin=0; ibin<=Nbins+1;ibin++){
    double sys_nom = m_hratio->GetBinContent(ibin);
    double BC = fabs(sys_nom-1);// |sys-nom|/nom
    m_hbkgfrac->SetBinContent(ibin,BC);
    m_hreluncert->SetBinContent(ibin,sys_nom-1);
    //if( BC < 1e-13 ) {
    //  std::cout<<"M="<<m_hratio->GetBinLowEdge(ibin)<<" BC="<< BC <<" sys/nom="<<sys_nom<<" "<<
    //	m_hsys->GetName()<<" "<<m_hsys->GetBinContent(ibin)<<" "<<
    //	m_hnom->GetName()<<" "<<m_hnom->GetBinContent(ibin)<<std::endl;
    //}
  }

}




