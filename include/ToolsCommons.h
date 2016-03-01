// In this class, the different initialisations
// of files, binning, lumi, yields, normalization
// is done.
// Feel free to modify the values !
#ifndef Commons_h
#define Common_h
#include <utility>
#include <iostream>
#include <vector>
#include <TString.h>

class Commons{

 public:
  static TString                  outdir;
  static TString                  GeeFile;
  static TString                  GrlFile;
  static TString                  PileupMCFile;
  static TString                  PileupDataFile;
  static std::pair<int,int>       norm_bin;
  static std::pair<double,double> norm_val;
  static double                   int_lumi_fb;
  static double                   int_lumi_pb;
  static std::vector<double>      yields;
  static std::vector<double>      yields_er;
  static double                   purity_relative_uncert;
  static int                      nBins;
  static double*                  binning;

  static int                      nBins_search;
  static double*                  binning_search;


  static void    Setup(bool VERBOSE=false);
  static void    SetNormVal( double Xmin=140.,
			     double Xmax=400. );
  static void    SetNormBin( int BinMin=20,
			     int BinMax=47 );
  static void    SetYields( double Ngg,
			    double Ngj,
			    double Njg,
			    double Njj,
			    bool VERBOSE=false );
  static void    SetYieldsStatError( double Ngg_er,
				     double Ngj_er,
				     double Njg_er,
				     double Njj_er );
  static void    SetBinning();
  static void    SetBinning_SearchRegion(double mgg_cut);
  
  static double  GetMCWeight(int mchannel);
  static bool    EtaCategory(double eta_L,double eta_SL,TString cat="NONE");
  static TString EtaCategory(TString cat="NONE");
  static bool    CosThetaStarCategory(double costh,TString cat="NONE");
  static TString CosThetaStarCategory(TString cat="NONE");

  static double GetkFactor_40_30(double mass);
  static double GetkFactor_50_50(double mass);
  static double GetTopoIsoCut(double pT);
  static double GetTopoIsoPtcorr(double pT);
  static double GetIsoDataMCDiff(double pT, TString L_or_SL, bool capped);

  static void SetPurityUncert(double uncert);


 private:
  //This class is not meant to be instanciate
  // Constructor unavailable (private)
  Commons(){;}
  ~Commons(){;}
};
#endif
