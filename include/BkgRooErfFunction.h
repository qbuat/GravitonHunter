/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOERFFUNCTION
#define ROOERFFUNCTION

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RooErfFunction : public RooAbsPdf {
public:
  RooErfFunction() {} ; 
  RooErfFunction(const char *name, const char *title,
		   RooAbsReal& _x,
		   RooAbsReal& _k1,
		   RooAbsReal& _k2);
  RooErfFunction(const RooErfFunction& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooErfFunction(*this,newname); }
  inline virtual ~RooErfFunction() { }

protected:

  RooRealProxy x ;
  RooRealProxy k1 ;
  RooRealProxy k2 ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooErfFunction,1) // Your description goes here...
};
 
#endif
