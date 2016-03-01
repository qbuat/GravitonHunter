//--> Compute the highest 
// probability density interval
// @ 68% confidence level around a central value
#ifndef HPDInterval_h
#define HPDInterval_h

#include <iostream>
#include <utility>
#include <TH2.h>

class HPDInterval
{
 private:
  TH1D m_hh;



 public :

  HPDInterval(const TH1D& hh);
  virtual ~HPDInterval();

  std::pair<double,double> GetHPDInterval(bool VERBOSE=false);

};

#endif

