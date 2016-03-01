#include <iostream>
#include <stdlib.h>
#include <map>

#include <TSystem.h>
#include <TError.h>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TBranch.h>
#include <TGraphErrors.h>

#include "ToolsCommons.h"
#include "ToolsTreeReader.h"
#include "IsolationFitter.h"
#include "IsolationGrapher.h"

using std::cerr;

int main(int argc, char ** argv){
  

  //------------- PROTECTION AGAINST WRONG NUMBER OF PARAMETERS ------------
  if( argc<13 ){
    std::cerr << "Wrong usage ! "
	      << argv[0] << " "
	      << "iso_type absetamin absetamax ptmin ptmax npvmin npvmax mumin mumax L_or_SL data" 
	      << " sample outfile \n" ;
    exit(1);
  }
  //------------------------------------------------------------------------

  //------------- INPUT PARAMETERS VALUES ------------------------------------
  TString iso_type  = argv[1]; //isolation type
  double  absetamin = atof(argv[2]); // lower |eta| bound
  double  absetamax = atof(argv[3]);// upper |eta| bound
  double  ptmin     = atof(argv[4]);// liwer pt bound
  double  ptmax     = atof(argv[5]);// upper pt bound
  double  npvmin    = atof(argv[6]);//lower NPV bound
  double  npvmax    = atof(argv[7]); //upper NPV bound
  double  mumin     = atof(argv[8]);//lower mu bound
  double  mumax     = atof(argv[9]); //upper mu bound
  TString L_or_SL   = argv[10];// Choose between leading and subleading
  bool    data      = atoi(argv[11]);// data=true mc=false
  TString sample    = argv[12];// Choose between leading and subleading
  TString outfile   = argv[13]; // output file name
  //--------------------------------------------------------------------------
  std::cerr << "Processing : " << argv[0]  << " " 
	    << iso_type << " " 
	    << absetamin << " " 
	    << absetamax << " " 
	    << ptmin << " "    
	    << ptmax << " "    
	    << npvmin << " "    
	    << npvmax << " "   
	    << mumin << " "    
	    << mumax << " "   
	    << L_or_SL << " "  
	    << data << " "     
	    << sample <<" " 
	    << outfile  
	    << "\n";

  Commons::Setup();

  //  if( iso_type !="CONE" && iso_type != "TOPO" && iso_type != "TOPO_PTDEP")
  //    Fatal("isolation_study","Wrong iso_type");

  

  TChain *treeD;
  TTree* treeMC;

  TreeReader * TR;

  if(data){
    treeD = new TChain("tree");
    treeD->Add( Commons::outdir+"runs_19_06_13/datafile_"+iso_type+"_run*.root");
    TR = new TreeReader(treeD);
  }
  
  else {
    // TFile* file = TFile::Open(Commons::outdir+"pythiagamgam_"+iso_type+".root");
    TFile* file;




    TString fname;
    if     (sample == "pythia_gg") fname = "bkg_mc12_gamgam_all_";
    else if(sample == "pythia_gg_npvrew") fname = "bkg_mc12_gamgam_all_withnpvweightbranch_";
    else if(sample == "pythia_gj") fname = "bkg_mc12_gamjet_all_";
    else if(sample == "pythia_rs") fname = "sig_mc12_Ggg_flat_";
    else if(sample == "pythia_merge") fname = "bkg_mc12_gamgam_gamjet_all_";
    else if(sample == "pythia_merge_npvrew") fname = "bkg_mc12_gamgam_gamjet_all_withnpvweightbranch_";



    std::cout<<sample<<" "<<iso_type<<std::endl;

    if( iso_type == "CONE" || iso_type == "TOPO" || iso_type == "TOPO_PTDEP") {
      //         file = TFile::Open(Commons::outdir+"prod_25_06_13/bkg_mc12_gamgam_all_"+iso_type+".root");
      
      if(sample == "pythia_gg") file = TFile::Open(Commons::outdir+"prod_25_06_13/"+fname+iso_type+".root");
      else                    	file = TFile::Open(Commons::outdir+"prod_07_10_13/"+fname+iso_type+".root");
    }
    else file = TFile::Open("/afs/in2p3.fr/home/j/jbrown/GravitonAnalysis/GravitonHunter/prod0/bkg_mc12_gamgam_all_"+iso_type+".root");
    treeMC = (TTree*)file->Get("tree");
    TR = new TreeReader(treeMC);
  }

  
  std::cout << "TOTO" << std::endl;

  if( iso_type == "TOPO_PTDEP_modifiediso" ) iso_type = "TOPO_PTDEP"  ;
  if( iso_type == "TOPO_modifiediso" )       iso_type = "TOPO"  ;


  TH1F * hpt  = new TH1F("hpt","hpt",1000,50,7000);
  TH1F * hiso = new TH1F("hiso","hiso",1000,-8,14);
//
//
//
//  double p0 =  9.74858e-03;
//  double p1 =  2.08163e-03;
//  double p2 = -1.36790e-04;
//  double p3 =  2.40864e-06;
//  double p4 = -6.65704e-02;
//  double p5 = -4.91916e+00;
//
//
//
//
//
//
        for(int i=0; i<TR->GetEntries(); i++) {
          
          TR->GetEntry(i);
      
          if(L_or_SL == "L") {
      
            if( !TR->GetVariable("IsTight_L")        ) continue;
            //      if( !TR->GetVariable("IsTight_SL")       ) continue;
      
            if( !TR->GetVariable("IsLoosePrime4_L")        ) continue;
            //      if( !TR->GetVariable("IsLoosePrime4_SL")       ) continue;
      
      
            if( TR->GetVariable("eta_L") < absetamin ) continue;
            if( TR->GetVariable("eta_L") > absetamax ) continue;
            if( TR->GetVariable("pT_L") < ptmin      ) continue;
            if( TR->GetVariable("pT_L") > ptmax      ) continue;
            if( TR->GetVariable("NPV") < npvmin      ) continue;
            if( TR->GetVariable("NPV") > npvmax      ) continue;
            if( TR->GetVariable("mu") < mumin	       ) continue;
            if( TR->GetVariable("mu") > mumax        ) continue; 
      
            //       double weight = data ? 1 : TR->GetVariable("weight");// * TR->GetVariable("gen_weight");
            double weight = 1;
            double npv = TR->GetVariable("NPV");
            
            if(!data) {
      	weight *= TR->GetVariable("weight");// * TR->GetVariable("gen_weight");
      	//	weight *= (p0+p1*npv+p2*npv*npv+p3*npv*npv*npv)*exp(p4*npv-p5);
            }
            if(iso_type == "TOPO_PTDEP")
      	hiso->Fill(TR->GetVariable("Iso_L_mod"),weight);
            else
      	hiso->Fill(TR->GetVariable("Iso_L"),weight);
            
      
            if(iso_type == "TOPO_PTDEP") {
      	if(TR->GetVariable("Iso_L_mod") < 5)
      	  hpt->Fill(TR->GetVariable("pT_L"),weight);
            }
            else {
      	if(TR->GetVariable("Iso_L") < 5) 
      	  hpt->Fill(TR->GetVariable("pT_L"),weight);
            }
      
          }
      
      
          else{
      
            //      if( !TR->GetVariable("IsTight_L")        ) continue;
            if( !TR->GetVariable("IsTight_SL")       ) continue;
      
            //      if( !TR->GetVariable("IsLoosePrime4_L")        ) continue;
            if( !TR->GetVariable("IsLoosePrime4_SL")       ) continue;
      
            if( TR->GetVariable("eta_SL") < absetamin ) continue;
            if( TR->GetVariable("eta_SL") > absetamax ) continue;
            if( TR->GetVariable("pT_SL") < ptmin      ) continue;
            if( TR->GetVariable("pT_SL") > ptmax      ) continue;
            if( TR->GetVariable("NPV") < npvmin      ) continue;
            if( TR->GetVariable("NPV") > npvmax      ) continue;
            if( TR->GetVariable("mu") < mumin	       ) continue;
            if( TR->GetVariable("mu") > mumax        ) continue; 
      
            //      double weight = data ? 1 : TR->GetVariable("weight");// * TR->GetVariable("gen_weight");
            double weight = 1;
            double npv = TR->GetVariable("NPV");
            
            if(!data) {
      	weight *= TR->GetVariable("weight");// * TR->GetVariable("gen_weight");
      	//	weight *= (p0+p1*npv+p2*npv*npv+p3*npv*npv*npv)*exp(p4*npv-p5);
            }
      
            if(iso_type == "TOPO_PTDEP") 
      	hiso->Fill(TR->GetVariable("Iso_SL_mod"),weight);
            else
      	hiso->Fill(TR->GetVariable("Iso_SL"),weight);
      
      
            if(iso_type == "TOPO_PTDEP") {
      	if(TR->GetVariable("Iso_SL_mod") < 5) 
      	  hpt->Fill(TR->GetVariable("pT_SL"),weight);
            }
            else {
      	if(TR->GetVariable("Iso_SL") < 5) 
      	  hpt->Fill(TR->GetVariable("pT_SL"),weight);
            }
      //
      //
      //      double weight = 1;
      //	if(!data) weight = TR->GetVariable("weight") * TR->GetVariable("gen_weight");
      
      
          }
         
          
        }



  IsolationFitter* IF = new IsolationFitter();
  if(data){
    IF->SetEntries(treeD);
    IF->SetTree(treeD);
  }else{
    IF->SetEntries(treeMC);
    IF->SetTree(treeMC);
  }
  IF->SetParFile("FitParameters/Iso_PtNPVeta.config");
  IF->SetStreamType(data);
//  if(data)
//    IF->SetMggBounds(100,1000);//data blinding
//  else
    IF->SetMggBounds(100,100000);
  IF->SetPtBounds(ptmin,ptmax);
  IF->SetEtaBounds(absetamin,absetamax);
  IF->SetNPVBounds(npvmin,npvmax);
  IF->SetMuBounds(mumin,mumax);

  if( iso_type == "CONE") {
    IF->SetIsoBounds(-10,25);
    IF->SetIsoNorm(10,25);
  }else {
    IF->SetIsoBounds(-8,14);
    IF->SetIsoNorm(10,14);
  }
  if(L_or_SL=="L"){
    if( iso_type == "TOPO_PTDEP" )
      IF->SetPtEtaIsoTightLoosePrimeNames("pT_L","eta_L","Iso_L_mod","IsTight_L","IsLoosePrime4_L");//leading photon
    else
      IF->SetPtEtaIsoTightLoosePrimeNames("pT_L","eta_L","Iso_L","IsTight_L","IsLoosePrime4_L");//leading photon
  }else if(L_or_SL=="SL"){
    if( iso_type == "TOPO_PTDEP" )
      IF->SetPtEtaIsoTightLoosePrimeNames("pT_SL","eta_SL","Iso_SL_mod","IsTight_SL","IsLoosePrime4_SL");//subleading photon
    else
      IF->SetPtEtaIsoTightLoosePrimeNames("pT_SL","eta_SL","Iso_SL","IsTight_SL","IsLoosePrime4_SL");//subleading photon
  }else Fatal("isolation_study", "WRONG L_or_SL flag");
  
  IF->Init_Vars();
  IF->Init_DataSets();
  IF->Fitter(true);


  std::cout<<outfile<<std::endl;

  IF->StoreToRootFile(outfile);
  TFile::Open(outfile,"UPDATE");
  hpt->Write();
  hiso->Write();

}


