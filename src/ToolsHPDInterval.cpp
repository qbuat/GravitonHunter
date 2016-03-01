#include "ToolsHPDInterval.h"

////////////////////////////////////////
HPDInterval::HPDInterval(const TH1D& hh) 
////////////////////////////////////////
{
  m_hh = *(TH1D*)hh.Clone("m_hh");

}
/////////////////////////////////////
HPDInterval::~HPDInterval()
/////////////////////////////////////
{
  //Default destructor
}
//////////////////////////////////////////////////////////////////
std::pair<double,double> HPDInterval::GetHPDInterval(bool VERBOSE)
/////////////////////////////////////////////////////////////////
{

  //The result will be more precise 
  //if the hist granularity is high
  double y = m_hh.GetBinContent(m_hh.GetMaximumBin());
  double iter = y/10000.;
  int i =0;
  double ydown = 0;  double xdown = 0;
  double xup   = 0;  double yup   = 0;
  int binup    = 0;  int bindown  = 0;
  double p = 0;
  while( p< 68){
    for (int ibin = m_hh.GetMaximumBin(); ibin< m_hh.GetNbinsX(); ibin++){
      double yiter = m_hh.GetBinContent(ibin);
      if(yiter < y-i*iter)
	break;
      else{
	yup   = yiter ;
	xup   = m_hh.GetBinCenter(ibin);
	binup = ibin;
      }
    }
    for (int ibin =  m_hh.GetMaximumBin(); ibin>0; ibin--){
      double yiter = m_hh.GetBinContent(ibin);
      if(yiter < y-i*iter)
	break;
      else{
	ydown   = yiter ;
	xdown   = m_hh.GetBinCenter(ibin);
	bindown = ibin;
      }
    }
    p = 100 * m_hh.Integral(bindown,binup)/m_hh.Integral();
    i++;
  }
 
  std::pair<double,double> temp(xdown,xup);

  if (VERBOSE) { 
    std::cout  << " Xmin : " << m_hh.GetBinCenter(bindown)
	       << " Xmax : " << m_hh.GetBinCenter(binup)
	       << " Inte : " << 100*m_hh.Integral(bindown,binup)/m_hh.Integral()
	       << " i : "    << i 
	       << " ydown : " << ydown << " yup : " << yup
	       << std::endl;    
  }
  return temp;
}
