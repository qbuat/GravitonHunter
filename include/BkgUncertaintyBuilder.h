// This class compute the uncertainty on the mgg
// bkg shape. But it can also be used to compute the uncertainty
// on any other shape. This uncertainty is computed as a simple ratio
// histogram of as |sys-nom|/sys

#ifndef BkgUncertaintyBuilder_h
#define BkgUncertaintyBuilder_h

#include <TH2.h>
#include <TMath.h>

class BkgUncertaintyBuilder
{
 private:

  TH1D* m_hsys;
  TH1D* m_hnom;
  TH1D* m_hratio;
  TH1D* m_hreluncert;
  TH1D* m_hbkgfrac;

  int m_binmin;
  int m_binmax;

 public :

  BkgUncertaintyBuilder(){;}
  BkgUncertaintyBuilder(const TH1D& hsys,const TH1D& hnom);
  virtual ~BkgUncertaintyBuilder();

  void SetNormRange(int min,int max);
  void SetHists(const TH1D& hsys,const TH1D& hnom);

  void Builder();

  TH1D* GetRatioHist(){return m_hratio;}
  TH1D* GetRelUncertHist(){return m_hreluncert;}
  TH1D* GetBkgFracHist(){return m_hbkgfrac;}



};

#endif
