#include <iostream>
#include <fstream>
#include <map>
#include <utility>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TChain.h>
#include <TBranch.h>
#include <TSystem.h>
#include <TF1.h> 
#include <TError.h>
#include <TH1.h>
#include "ToolsCommons.h"
#include "ToolsTreeReader.h"
#include "ToolsUtilities.h"
#include "ToolsChainMaker.h"
#include "BkgExtrapolation.h"


using std::cerr;

int main(int argc, char ** argv){
  

  //------------- PROTECTION AGAINST WRONG NUMBER OF PARAMETERS ------------
  if( argc<2 ){
    std::cerr << "Wrong usage ! "
	      << argv[0] 
	      << "inputfile \n" ;
    exit(1);
  }
  //------------------------------------------------------------------------

  //------------- INPUT PARAMETERS VALUES ------------------------------------
  std::string file      = argv[1]; // 
  //--------------------------------------------------------------------------
  std::cerr << "Processing : " << argv[0] << " " << file << " " 
	    << "\n";


  Commons::Setup();
  

  std::ifstream ifs(file.c_str());
  //---------------------------------------
  std::string argStr;
  char buf[256+1];
  unsigned int delpos;
  while (true){
    ifs.read(buf,256);
    if (ifs.eof()){
      if (ifs.gcount() == 0) break;
      delpos = ifs.gcount()-1;
    }else{
      delpos = ifs.gcount();
    }
    buf[delpos] = 0x00;
    argStr += buf;
  }
  // cout<<"argStr  ="<<argStr<<"."<<endl;
  //---------------------------------------


  std::vector<double> mu_bins_v;
  mu_bins_v.push_back(5);
  mu_bins_v.push_back(10);
  mu_bins_v.push_back(15);
  mu_bins_v.push_back(20);
  mu_bins_v.push_back(25);
  mu_bins_v.push_back(30);
  mu_bins_v.push_back(35);
  mu_bins_v.push_back(40);

  //======================================//
  TChain *chain=0;
  ChainMaker MyChainMaker("photon");
  chain = new TChain("photon");
  MyChainMaker.SetTreeFiles(argStr);
  chain = MyChainMaker.GetChain();
  //______________________________________//

  std::cout << "Declare the tree reader" << std::endl;
  TreeReader Rd(chain);

  std::cout << "Declare the histogram" << std::endl;
  const int Nbins = 44;
  const double binning[]= {-300,-250,-225,-200,-175,-160,-150,-140,-130,-120,-110,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,-0.2,0,
			   0.2  ,10 ,20 ,30,40  ,50, 60, 70, 80,90,100,110,120,130,140,150,160,175,200,225,250,300};

  // TH1F* hz = new TH1F("hz","hz",100,-10,10);
  TH1F* hz = new TH1F("hz","hz",Nbins,binning);
  hz->Sumw2();
  hz->GetXaxis()->SetTitle("#Delta z=z^{reco}-z^{truth} [mm]");
  const int size = mu_bins_v.size()-1;
  TH1F* hz_mu[size];
  for(int i=0 ; i<size;i++){
    std::cout << "Declare the histogram " << i << std::endl;
    // hz_mu[i] = new TH1F( Form("hz_mu%d",i),Form("hz_mu%d",i),100,-10,10);
    hz_mu[i] = new TH1F( Form("hz_mu%d",i),Form("hz_mu%d",i),Nbins,binning);
    hz_mu[i]->GetXaxis()->SetTitle("#Delta z=z^{reco}-z^{truth} [mm]");
  }

  for(int i=0 ; i<size;i++){
    hz_mu[i]->SetLineColor(kRed+i);
    hz_mu[i]->SetLineStyle(kDashed+i);
    hz_mu[i]->Sumw2();
  }

  std::cout << "Start of the loop ..." ;
  for(int entry = 0; entry<chain->GetEntries();entry++){
    Rd.GetEntry(entry);
    hz->Fill( Rd.GetVariable("PV_z[0]-mc_vx_z[0]") );
    for( int i=0; i<(int)(mu_bins_v.size()-1);i++){
      if( Rd.GetVariable("averageIntPerXing") > mu_bins_v[i] &&
	  Rd.GetVariable("averageIntPerXing") <= mu_bins_v[i+1] ){
	hz_mu[i]->Fill( Rd.GetVariable("PV_z[0]-mc_vx_z[0]") );
      }
    }
  }
  std::cout << "... end of the loop !" << std::endl;

  hz->Scale(1./hz->Integral());

  // hz->Draw();
  for(int i=0 ; i<(int)mu_bins_v.size()-1;i++){
    hz_mu[i]->Scale(1./hz_mu[i]->Integral());
    // hz_mu[i]->Draw("same");
  }

  TFile fout("toto.root","RECREATE");
  fout.Add(hz);
  for(int i=0 ; i<(int)mu_bins_v.size()-1;i++)
    fout.Add(hz_mu[i]);
  fout.Write();
  fout.Close();

  // std::cout << "***** MEMORY CLEANING ***** "<< std::endl;
  // delete trD;
  // std::cout << "***** END OF THE RUN ***** "<< std::endl;

  return 0;

}



