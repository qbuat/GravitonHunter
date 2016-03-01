
#include "Riostream.h" 

#include "BkgRooReducibleFunction.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 

ClassImp(RooReducibleFunction) 


RooReducibleFunction::RooReducibleFunction()
{
  _coefIter = _coefList.createIterator();
}

 RooReducibleFunction::RooReducibleFunction( const char *name, const char *title, 
					     RooAbsReal& x,const RooArgList& coefList,
					     Int_t lowestOrder):
   RooAbsPdf(name,title), _x("x","x",this,x),
   _coefList("coefList","List of coefficients",this),
   _lowestOrder(lowestOrder)
 { 

   // Constructor
   _coefIter = _coefList.createIterator() ;
   
   // Check lowest order
   if (_lowestOrder<0) {
     cout << "RooReducibleFunction::ctor(" << GetName() 
	  << ") WARNING: lowestOrder must be >=0, setting value to 0" << endl ;
     _lowestOrder=0 ;
  }

  TIterator* coefIter = coefList.createIterator() ;
  RooAbsArg* coef ;
  while((coef = (RooAbsArg*)coefIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      cout << "RooReducibleFunction::ctor(" << GetName() << ") ERROR: coefficient " << coef->GetName() 
	   << " is not of type RooAbsReal" << endl ;
      assert(0) ;
    }
    _coefList.add(*coef) ;
  }
  delete coefIter;
 } 


RooReducibleFunction::RooReducibleFunction(const char* name, const char* title,
					   RooAbsReal& x) :
  RooAbsPdf(name, title),
  _x("x", "Dependent", this, x),
  _coefList("coefList","List of coefficients",this),
  _lowestOrder(1)
{
  _coefIter = _coefList.createIterator() ;
}

RooReducibleFunction::RooReducibleFunction(const RooReducibleFunction& other, const char* name) :
  RooAbsPdf(other, name), 
  _x("x", this, other._x), 
  _coefList("coefList",this,other._coefList),
  _lowestOrder(other._lowestOrder) 
{
  // Copy constructor
  _coefIter = _coefList.createIterator() ;
}

RooReducibleFunction::~RooReducibleFunction()
{
  // Destructor
  delete _coefIter ;

}


 Double_t RooReducibleFunction::evaluate() const 
 { 

   Int_t order(_lowestOrder) ;
   Double_t sum(order<1 ? 0 : 1) ;

   _coefIter->Reset() ;
   
   RooAbsReal* coef ;
   const RooArgSet* nset = _coefList.nset() ;
   while((coef=(RooAbsReal*)_coefIter->Next())) {
     sum += coef->getVal(nset)*TMath::Power(log(_x),order++) ;
   }

   return exp(sum);


 } 



