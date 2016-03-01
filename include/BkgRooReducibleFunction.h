
#ifndef ROOREDUCIBLEFUNCTION
#define ROOREDUCIBLEFUNCTION

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooCategoryProxy.h"

#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RooArgList;

class RooReducibleFunction : public RooAbsPdf {
public:
  RooReducibleFunction(); 
  RooReducibleFunction(const char *name, const char *title,RooAbsReal& x);
  RooReducibleFunction(const char *name, const char *title, 
		       RooAbsReal& _x,const RooArgList& _coefList,Int_t lowestOrder=1);
  RooReducibleFunction(const RooReducibleFunction& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooReducibleFunction(*this,newname); }
  virtual ~RooReducibleFunction();

protected:

  RooRealProxy _x ;
  RooListProxy _coefList;
  Int_t _lowestOrder;
  TIterator* _coefIter;  
  Double_t evaluate() const ;

private:

  ClassDef(RooReducibleFunction,1) // Your description goes here...
};
 
#endif
