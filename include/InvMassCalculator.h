// This class inherits from GravitonAnalysis.
#ifndef InvMassCalc_h
#define InvMassCalc_h

#include <iostream>
#include "TTree.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "GravitonAnalysis.h"


class InvMassCalc : public GravitonAnalysis
{
 private:
  void   CorrectEta( float,float,
		     float,float*,float* );
  static double ReturnRZ_1stSampling_cscopt2(double);

 public :
  InvMassCalc();
  InvMassCalc(TTree* tree);
  virtual ~InvMassCalc();
  //=========Invariant Mass Calculation===========//
  static double   GetCorrectedInvMass( double,double,
				       double,double,
				       double,double,
				       double );
  static double   EtaS1PVCorrected(double,double);
  double   GetInvMass(int,int,int);
  double   GetInvMass2(int, int,int,int smearing=-1);
  //-1: no smearing,0:nominal,1down,2:up
  double   GetInvMass3(int,int,int,int smearing=-1);

};

#endif

