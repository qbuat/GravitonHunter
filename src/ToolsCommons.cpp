#include "ToolsCommons.h"
#include <TMath.h>
#include <TF1.h>
TString                  Commons::outdir;
TString                  Commons::GeeFile;
TString                  Commons::GrlFile;
TString                  Commons::PileupMCFile;
TString                  Commons::PileupDataFile;
std::pair<int,int>       Commons::norm_bin;
std::pair<double,double> Commons::norm_val;
double                   Commons::int_lumi_fb;
double                   Commons::int_lumi_pb;
std::vector<double>      Commons::yields;
std::vector<double>      Commons::yields_er;
double                   Commons::purity_relative_uncert;
int                      Commons::nBins;
double*                  Commons::binning;
int                      Commons::nBins_search;
double*                  Commons::binning_search;

/////////////////////////////////
void Commons::Setup(bool VERBOSE)
/////////////////////////////////
{
  //  Commons::outdir = "/sps/atlas/q/qbuat/OUTDS/GravitonAnalysis8TeV/";
  Commons::outdir = "/sps/atlas/j/jbrown/OUTDS/GravitonAnalysis8TeV/";
  //--> Useful files 
  Commons::GeeFile = "/sps/atlas/q/qbuat/OUTDS/GravitonAnalysis8TeV/Gee_events_8TeV.txt";
  Commons::GrlFile = "./grl/data12_8TeV.periodAllYear_DetStatus-v61-pro14-02_DQDefects-00-01-00_PHYS_StandardGRL_All_Good.xml";
  Commons::PileupMCFile   = "./rootfiles/pu_mc12a.prw.root";
  // Commons::PileupDataFile = "./rootfiles/ilumicalc_histograms_EF_g35_loose_g25_loose_200841-208354.root";
  Commons::PileupDataFile = "./rootfiles/ilumicalc_histograms_EF_g35_loose_g25_loose_200841-208354.root";

  //--> Analysis with full 2012 dataset 
  Commons::int_lumi_pb = 20339;// (2012) with e/g GRL + trigger EF_g35loose_g25loose
  Commons::int_lumi_fb = 20.3;

  // SetYields(24378,5380,1796,1887);// Yields for 9.4/fb (cone iso<5)
  // SetYieldsStatError(195,59,45,22);// Yields stat error for 9.4/fb (cone iso<5)

  // SetYields(23933,4601,1784,1261);// Yields for 9.4/fb (topo iso<4)
  // SetYieldsStatError(195,59,45,22);// Yields stat error for 9.4/fb (topo iso<4)

  SetYields(24109,5624,2205,1848);// Yields for 9.4/fb (topo iso<5)
  SetYieldsStatError(194,70,55,28);// Yields stat error for 9.4/fb (topo iso<5)

  SetNormVal(140,400);
  SetBinning();
  SetBinning_SearchRegion(400);

  SetPurityUncert(0.1);

  if(VERBOSE) {
    double totyield = 0;
    for(int i=0;i<(int)yields.size();i++) totyield += yields[i];
    std::cout << "\E[31;1m ----------------------------------------------- \E[m" << std::endl;
    std::cout << "\E[32;1m ------------ RUNNING CONFIGURATION ------------ \E[m" << std::endl;
    std::cout << "\E[31;1m ----------------------------------------------- \E[m" << std::endl;
    std::cout << "\E[32;1m Integrated Luminosity : "
	      << int_lumi_fb << "/fb \E[m" << std::endl;
    std::cout << "\E[34;1m GRL : "
	      << GrlFile << "\E[m" << std::endl;
    std::cout << "\E[35;1m Yields : gg = " 
	      << yields[0] << "(~"<<(int)(100*yields[0]/totyield)<<")%"
	      << ", gj = " << yields[1] << "(~"<<(int)(100*yields[1]/totyield)<<")%"
	      << ", jg = " << yields[2] << "(~"<<(int)(100*yields[2]/totyield)<<")%"
	      << ", jj = " << yields[3] << "(~"<<(int)(100*yields[3]/totyield)<<")%"
	      << "\E[m" << std::endl;
    std::cout << "\E[36;1m The purity uncertainty is set to " 
	      <<  purity_relative_uncert
	      << "\E[m" << std::endl;

    std::cout << "\E[31;1m ----------------------------------------------- \E[m" << std::endl;
    std::cout << "\E[31;1m ----------------------------------------------- \E[m" << std::endl;
  }






}

///////////////////////////////////////////////////
void Commons::SetPurityUncert(double uncert)
{ 
  Commons::purity_relative_uncert= uncert;
}
///////////////////////////////////////////////////

///////////////////////////////////////////////////
void Commons::SetNormVal(double Xmin, double Xmax)
//////////////////////////////////////////////////
{
  Commons::norm_val.first  = Xmin;
  Commons::norm_val.second = Xmax;
}
////////////////////////////////////////////////
void Commons::SetNormBin(int BinMin, int BinMax)
////////////////////////////////////////////////
{
  Commons::norm_bin.first  = BinMin;
  Commons::norm_bin.second = BinMax;
}
//////////////////////////////////////////////////////////////////////////////////
void Commons::SetYields(double Ngg,double Ngj,double Njg,double Njj,bool VERBOSE)
/////////////////////////////////////////////////////////////////////////////////
{
  Commons::yields.clear();
  Commons::yields.push_back(Ngg);
  Commons::yields.push_back(Ngj);
  Commons::yields.push_back(Njg);
  Commons::yields.push_back(Njj);

  double totyield = 0;
  for(int i=0;i<(int)yields.size();i++) totyield += yields[i];
  if( VERBOSE )
    std::cout << "\E[35;1m Yields : gg = " 
	      << yields[0] << "(~"<<(int)(100*yields[0]/totyield)<<")%"
	      << ", gj = " << yields[1] << "(~"<<(int)(100*yields[1]/totyield)<<")%"
	      << ", jg = " << yields[2] << "(~"<<(int)(100*yields[2]/totyield)<<")%"
	      << ", jj = " << yields[3] << "(~"<<(int)(100*yields[3]/totyield)<<")%"
	      << "\E[m" << std::endl;
}
/////////////////////////////////////////////////////////////////////////////////////////
void Commons::SetYieldsStatError(double Ngg_er,double Ngj_er,double Njg_er,double Njj_er)
/////////////////////////////////////////////////////////////////////////////////////////
{
  Commons::yields_er.clear();
  Commons::yields_er.push_back(Ngg_er);
  Commons::yields_er.push_back(Ngj_er);
  Commons::yields_er.push_back(Njg_er);
  Commons::yields_er.push_back(Njj_er);
}
//////////////////////////
void Commons::SetBinning()
//////////////////////////
{
  //----------> Logarithmic binning
  Commons::nBins      = 100;
  double xmin_GeV     = 70;
  double xmax_GeV     = 3000;
  double logxmin_GeV  = log10(xmin_GeV);
  double logxmax_GeV  = log10(xmax_GeV);
  double binwidth_GeV = (logxmax_GeV-logxmin_GeV)/double(nBins);
  Commons::binning = new double[nBins+1];
  Commons::binning[0] = xmin_GeV;
  for (int i=1 ;i<=Commons::nBins ;i++){
    Commons::binning[i] = pow(10,logxmin_GeV+i*binwidth_GeV);
  }

  SetNormBin(19,47);//log binning
  // std::cout << "[" << Commons::binning[19] << "," << Commons::binning[47] << "]" << std::endl;
  //---------------------------------------------------------

  //---------> Linear binning
  // Commons::nBins      = 300;
  // double xmin_GeV     = 0;
  // double xmax_GeV     = 3000;
  // Commons::binning = new double[nBins+1];
  // Commons::binning[0] = xmin_GeV;
  // for (int i=1 ;i<=nBins ;i++){
  //   Commons::binning[i] = Commons::binning[i-1]+10;
  // }
  // SetNormBin(15,40);//linear binning
  //---------------------------------------------------
}
/////////////////////////////////////////////////////
void Commons::SetBinning_SearchRegion(double mgg_cut)
/////////////////////////////////////////////////////
{
  int istart = -999;
  for(int i=0;i<=Commons::nBins;i++){
    if( Commons::binning[i] < mgg_cut) continue;
    else{
      istart = i;
      break;
    }
  }


  //----------> Logarithmic binning
  Commons::nBins_search      = Commons::nBins-istart;
  Commons::binning_search = new double[nBins_search+1];
  for (int i=0 ;i<=nBins_search ;i++)
    Commons::binning_search[i] = Commons::binning[istart+i];
}

///////////////////////////////////////////
double Commons::GetMCWeight(int mcchannel)
//////////////////////////////////////////
{
  double temp;
  double sigma;
  double filter_eff;
  double Ngen;
  if(mcchannel==158339){
    sigma = 6.2868e-03;//nb 
    filter_eff = 0.47908;
    Ngen = 300000;
  }else if(mcchannel==158340){
    sigma = 1.6436e-03;//nb 
    filter_eff = 0.46320;
    Ngen = 200000;
  }else if(mcchannel==158341){
    sigma = 1.8900e-04;//nb 
    filter_eff = 0.42919 ;
    Ngen = 200000;
  }else if(mcchannel==158342){
    sigma = 1.4265e-05 ;//nb 
    filter_eff = 0.38662 ;
    Ngen = 200000;
  }else if(mcchannel==158343){
    sigma = 1.0637e-06 ;//nb 
    filter_eff = 0.33949 ;
    Ngen = 100000;
  }else if(mcchannel==180189){
    sigma = 6.9192 ;
    filter_eff = 7.1429e-04  ;
    Ngen = 97600;
  }else if(mcchannel==180190){
    sigma = 1.3328 ;
    filter_eff = 8.7400e-04 ;
    Ngen = 79400 ;
  }else if(mcchannel==180191){
    sigma = 1.1139e-01 ;
    filter_eff = 7.9171e-04 ;
    Ngen = 49200;
  }else if(mcchannel==180192){
    sigma = 6.2127e-03 ;
    filter_eff = 5.1235e-04 ;
    Ngen = 48600;
  }else if(mcchannel==180193){
    sigma = 3.4891E-04 ;
    filter_eff = 2.8202E-04 ;
    Ngen = 19200;
  }
  //else Fatal("ToolsCommons::GetMCWeight()","Wrong mc channel" );
  temp = int_lumi_pb*sigma*1e3*filter_eff/Ngen;
  return temp;
}

//////////////////////////////////////////////////////////////////
bool Commons::EtaCategory(double eta_L,double eta_SL, TString cat)
//////////////////////////////////////////////////////////////////
{
  if( cat== "CC")//Central-Central
    if( fabs(eta_L)<1.37 && fabs(eta_SL)<1.37) return true;
    else return false;
  else if( cat=="CE")// Central-Endcap
    if( fabs(eta_L)<1.37 && fabs(eta_SL)>1.37 ) return true;
    else return false;
  else if( cat=="EC")// Endcap-Central
    if( fabs(eta_L)>1.37 && fabs(eta_SL)<1.37) return true;
    else return false;
  else if( cat=="EE" )//Endcap-Endcap
    if( fabs(eta_L)>1.37 && fabs(eta_SL)>1.37) return true;
    else return false;
  else if( cat=="EE_S" )//Endcap-Endcap same sign
    if( fabs(eta_L)>1.37 && fabs(eta_SL)>1.37 && eta_L*eta_SL>0 ) return true;
    else return false;
  else if( cat=="EE_O" )//Endcap-Endcap opposite sign
    if( fabs(eta_L)>1.37 && fabs(eta_SL)>1.37 && eta_L*eta_SL<0 ) return true;
    else return false;
  else if(cat=="NONE")
    return true;
  else{
    Fatal("Commons::EtaCategory", "WRONG CATEGORY !!");
    return false;
  }
}
/////////////////////////////////////////
TString Commons::EtaCategory(TString cat)
/////////////////////////////////////////
{
  TString st;
  if(cat=="CC") st="abs(eta_L)<1.37 && abs(eta_SL)<1.37";
  else if(cat=="CE") st="abs(eta_L)<1.37 && abs(eta_SL)>1.37";
  else if(cat=="EC") st="abs(eta_L)>1.37 && abs(eta_SL)<1.37";
  else if(cat=="EE") st="abs(eta_L)>1.37 && abs(eta_SL)>1.37";
  else if(cat=="EE_S") st="abs(eta_L)>1.37 && abs(eta_SL)>1.37 && eta_L*eta_SL>0";
  else if(cat=="EE_O") st="abs(eta_L)>1.37 && abs(eta_SL)>1.37 && eta_L*eta_SL<0";
  else if(cat=="NONE") st="abs(eta_L)<4 && abs(eta_SL)<4";
  else Fatal("Commons::EtaCategory","WRONG CATEGORY!!");
  return st;
}
////////////////////////////////////////////////////////////
bool Commons::CosThetaStarCategory(double costh,TString cat)
////////////////////////////////////////////////////////////
{
  if( cat== "CosThetaStar0")//0<|cos(theta*)|<0.2 
    if( fabs(costh)>=0 && fabs(costh)<0.2) return true;
    else return false;
  else if( cat== "CosThetaStar1")//0.2<|cos(theta*)|<0.4 
    if( fabs(costh)>=0.2 && fabs(costh)<0.4) return true;
    else return false;
  else if( cat== "CosThetaStar2")//0.4<|cos(theta*)|<0.6 
    if( fabs(costh)>=0.4 && fabs(costh)<0.6) return true;
    else return false;
  else if( cat== "CosThetaStar3")//0.6<|cos(theta*)|<0.8 
    if( fabs(costh)>=0.6 && fabs(costh)<0.8) return true;
    else return false;
  else if( cat== "CosThetaStar4")//0.8<|cos(theta*)|<1 
    if( fabs(costh)>=0.8 && fabs(costh)<=1) return true;
    else return false;
  else if(cat=="NONE") return true;
  else{
    Fatal("Commons::CosThetaStarCategory", "WRONG CATEGORY !!");
    return false;
  }
}
//////////////////////////////////////////////////
TString Commons::CosThetaStarCategory(TString cat)
//////////////////////////////////////////////////
{
  TString st;
  if(cat=="CosThetaStar0") st="abs(costhetastar)>0 && abs(costhetastar)<0.2";
  else if(cat=="CosThetaStar1") st="abs(costhetastar)>0.2 && abs(costhetastar)<0.4";
  else if(cat=="CosThetaStar2") st="abs(costhetastar)>0.4 && abs(costhetastar)<0.6";
  else if(cat=="CosThetaStar3") st="abs(costhetastar)>0.6 && abs(costhetastar)<0.8";
  else if(cat=="CosThetaStar4") st="abs(costhetastar)>0.8 && abs(costhetastar)<1";
  else if(cat=="NONE") st = "abs(costhetastar) >-10 && abas(costhetastar)<10";
  else Fatal("Commons::CosThetaStarCategory","WRONG CATEGORY !!");
  return st;
}

///////////////////////////////////////////////////////
double Commons::GetkFactor_40_30(double mass)
///////////////////////////////////////////////////////
{
  // TF1 fk("fk","exp([0]+[1]*x)+[2]",120,40000);
  // fk.SetParameters(2.87873e-01,-1.16827e-03,7.74223e-01);//Fit 140-1500
  // if( mass < 120.) return 1.;
  // else return fk.Eval(mass);

  TF1 fk_low("fk_low","[0]*x*x+[1]*x+[2]",120,250);
  TF1 fk_high("fk_high","exp([0]+[1]*x)+[2]",250,4000000);

  fk_low.SetParameters( -7.0206e-06, 0.00241007 , 1.57908    );
  fk_high.SetParameters( 0.236178   , -0.00167918, 0.928866   );

  if( mass < 120.) return 1.;
  else if(mass < 250.) return fk_low.Eval(mass);
  else return fk_high.Eval(mass);


}
//////////////////////////////////////////////////////
double Commons::GetkFactor_50_50(double mass)
//////////////////////////////////////////////////////
{
  // fit with a 2nd order polynom up to 250
  TF1 fk_low("fk_low","[0]*x*x+[1]*x+[2]",120,250);
  TF1 fk_high("fk_high","exp([0]+[1]*x)+[2]",250,4000000);

  fk_low.SetParameters(-1.13633e-05,0.00494087,1.11658);
  fk_high.SetParameters(0.0213382,-0.00143431,0.942951);
  
  if( mass < 120.) return 1.;
  else if(mass < 250.) return fk_low.Eval(mass);
  else return fk_high.Eval(mass);

}

////////////////////////////////////////
double Commons::GetTopoIsoCut(double pT)
///////////////////////////////////////
{
  //TopoIsoCut derived between 50 and 1800 gev with a RS flat sample.

  TF1 f("f","[0]+[1]*x+[2]*x*x",0,100000);
  f.SetParameters(4.9337,0.000477272,2.63033e-06);
  return (double)f.Eval(pT);

// p0                        =       4.9337   +/-   0.301794
// p1                        =  0.000477272   +/-   0.00079499
// p2                        =  2.63033e-06   +/-   4.49163e-07
}

///////////////////////////////////////////
double Commons::GetTopoIsoPtcorr(double pT)
///////////////////////////////////////////
{
  //TopoIsoCut derived between 50 and 1800 gev with a RS flat sample.

  TF1 f("f","[0]+[1]*x+[2]*x*x - 5",0,100000);
  f.SetParameters(4.9337,0.000477272,2.63033e-06);
  return (double)f.Eval(pT);

// p0                        =       4.9337   +/-   0.301794
// p1                        =  0.000477272   +/-   0.00079499
// p2                        =  2.63033e-06   +/-   4.49163e-07
}






///////////////////////////////////////////
double Commons::GetIsoDataMCDiff(double pT, TString L_or_SL, bool capped)
///////////////////////////////////////////
{
  //Difference of meanCB between data and MC
  // capped = 1 : correction is capped at max value above pt=500GeV
  // capped = 0 : correction is extrapolated above 500 GeV

  if(L_or_SL == "L") {
    if(pT>500 && capped) pT=500;
    return 0.001907*pT;
  }

  else {
    if(pT>500 && capped) pT=500;
    return 0.002042*pT;
  }

}


/////////////////////////////////////////////
//double Commons::GetIsoDataMCDiff(double pT, TString L_or_SL, bool capped)
/////////////////////////////////////////////
//{
//  //Difference of meanCB between data and MC
//  // capped = 1 : correction is capped at max value above pt=500GeV
//  // capped = 0 : correction is extrapolated above 500 GeV
//
//  if(L_or_SL == "L") {
//    if(pT>500 && capped) pT=500;
//    return -0.265799 + 0.00787452*pT -3.77404e-05*pT*pT + 6.78405e-08*pT*pT*pT;    
//  }
//
//  else {
//    if(pT>500 && capped) pT=500;
//    return -0.265799 + 0.00787452*pT -3.77404e-05*pT*pT + 6.78405e-08*pT*pT*pT;    
//  }
//
//}










// BEFORE MODIF 18 sep 
/////////////////////////////////////////////
//double Commons::GetIsoDataMCDiff(double pT, TString L_or_SL)
/////////////////////////////////////////////
//{
//
//  if(L_or_SL == "L") {
//    if(pT>=500) {
//      pT=500;
//      return -0.265799 + 0.00787452*pT -3.77404e-05*pT*pT + 6.78405e-08*pT*pT*pT;
//    }
//    else
//      return -0.265799 + 0.00787452*pT -3.77404e-05*pT*pT + 6.78405e-08*pT*pT*pT;
//  }
//
//  else {
//    if(pT>=500) {
//      pT=500;
//      return -0.14177 + 0.00395668*pT -2.24508e-06*pT*pT;
//    }
//    else
//      return -0.14177 + 0.00395668*pT -2.24508e-06*pT*pT;
//  }
//  
//
//}




