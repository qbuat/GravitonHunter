/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "BkgRooDijetFunction.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 

ClassImp(RooDijetFunction) 

 RooDijetFunction::RooDijetFunction( const char *name, const char *title, 
				     RooAbsReal& _x,
				     RooAbsReal& _k1,
				     RooAbsReal& _k2 ) :
RooAbsPdf(name,title), 
  x("x","x",this,_x),
  k1("k1","k1",this,_k1),
  k2("k2","k2",this,_k2)
 { 
 } 


 RooDijetFunction::RooDijetFunction(const RooDijetFunction& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   k1("k1",this,other.k1),
   k2("k2",this,other.k2)
 { 
 } 



 Double_t RooDijetFunction::evaluate() const 
 { 

   return pow(x,k1*(1-log(x)))*pow(x,k2*log(x));

   // return pow(x, k1) * pow(x, k2 * log(x));



 } 


