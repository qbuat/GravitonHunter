/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROODIJETFUNCTION
#define ROODIJETFUNCTION

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RooDijetFunction : public RooAbsPdf {
public:
  RooDijetFunction() {} ; 
  RooDijetFunction(const char *name, const char *title,
		   RooAbsReal& _x,
		   RooAbsReal& _k1,
		   RooAbsReal& _k2);
  RooDijetFunction(const RooDijetFunction& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooDijetFunction(*this,newname); }
  inline virtual ~RooDijetFunction() { }

protected:

  RooRealProxy x ;
  RooRealProxy k1 ;
  RooRealProxy k2 ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooDijetFunction,1) // Your description goes here...
};
 
#endif
