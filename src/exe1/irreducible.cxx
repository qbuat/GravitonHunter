#include <iostream>
#include <map>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TChain.h>
#include <TBranch.h>
#include <TH1.h>
#include <TSystem.h>
#include <TF1.h> 
#include <TError.h>
#include "ToolsCommons.h"
#include "ToolsTreeReader.h"
#include "ToolsUtilities.h"

#include "BkgExtrapolation.h"
#include "Bkg2DFit.h"
#include "BkgEstimatorUpgrade.h"
#include "GoodRunsLists/DQHelperFunctions.h"

using std::cerr;

int main(int argc, char ** argv){



  //------------- PROTECTION AGAINST WRONG NUMBER OF PARAMETERS ------------
  if( argc<2 ){
    std::cerr << "Wrong usage ! "
	      << " data_x etacat filename \n" ;
    exit(1);
  }
  //------------------------------------------------------------------------



  //------------- INPUT PARAMETERS VALUES ------------------------------------
  TString etacat    = argv[1];// eta category : --> "NONE", "CC", "CE", "EC", "EE", "EE_S", "EE_O"
  TString filename  = argv[2];

  TString iso_type  = "TOPO_PTDEP"; //argv[1]; // iso_type --> "TOPO" or "CONE"
  double isocut_l   = 8.; //atof(argv[2]); // lead isocut value
  double isocut_sl  = 8.; //atof(argv[2]); // sublead isocut value
  double leadcut    = 50; //atoi(argv[4]); //leading photon cut
  double subleadcut = 50.; //atoi(argv[5]); //subleading photon cut
  bool doDATA       = false;
  bool doIRR        = true;
  bool doRED        = false;
  bool doPUR        = false;
  bool doSYS        = true;
  bool doFINAL      = false;
  //--------------------------------------------------------------------------
  std::cerr << "Processing : " << argv[0] << " " << iso_type << " " 
	    << Form("%1.1f",isocut_l) << " " << etacat << " " 
	    << Form("%d %d",(int)leadcut,(int)subleadcut) << " "
	    << Form("%d %d %d %d %d %d %s \n",(int)doDATA,(int)doIRR,(int)doRED,(int)doPUR,(int)doSYS,(int)doFINAL, filename.Data());


  Commons::Setup();

  // Commons::SetPurityUncert(0.02);
  //----> Normalize between 180 and 409 (bin 26 to 47) for 50-50 config
  Commons::SetNormBin(26,47);

  std::cout << "--------> RUN PROCESSING <-------" << std::endl;
  if(doDATA)  std::cout << "---> Process DATA HISTOGRAMS <---" << std::endl;
  if(doIRR)   std::cout << "---> Process IRREDUCIBLE BKG <---" << std::endl;
  if(doRED)   std::cout << "---> Process REDUCIBLE BKG <-----" << std::endl;
  if(doPUR)   std::cout << "---> Process PURITY ESTIMATE <---" << std::endl;
  if(doSYS)   std::cout << "---> Process UNCERTAINTIES <-----" << std::endl;
  if(doFINAL) std::cout << "---> Process FINAL ESTIMATE <----" << std::endl;
  std::cout << "----------------><---------------" << std::endl;

  std::cout << "-----------> RUN CONFIGURATION <----------" << std::endl;
  std::cout << Form( "--> Normalization region : [%1.2f,%1.2f]",
		     Commons::binning[Commons::norm_bin.first-1],
		     Commons::binning[Commons::norm_bin.second]) << std::endl;
  std::cout << Form( "--> Pt cuts : (%d,%d)",(int)leadcut,(int)subleadcut) << std::endl; 
  std::cout << "--> Type of isolation : " << iso_type << std::endl;
  std::cout << Form( "--> Isolation cut value = %1.1f GeV",isocut_l) << std::endl; 
  if( iso_type.Contains("PTDEP" ))
    std::cout << "--> Use pT dependent isocut" << std::endl; 
  std::cout << "--> Eta category : " << etacat << std::endl;
  std::cout << "-----------------><-----------------" << std::endl;

  
  TString iso_type_name;
  if( iso_type.Contains("TOPO") )
    iso_type_name = "TOPO";
  else if(iso_type.Contains("CONE") )
    iso_type_name = "CONE";
  else Fatal("main()","Wrong iso type !");

  //------------------ Data Handling -------------------------------------------------------
  TChain * trB = new TChain("tree");
  trB->Add("/sps/atlas/q/qbuat/DIPHOTON_2012_REF/fuerSarahEtAl/bkg_mc12_gamgam_all_TOPO_PTDEP.root");
  //  trB->Add(Commons::outdir+"prod_07_10_13/bkg_mc12_gamjet_all_"+iso_type_name+".root");

  std::cout<<trB->GetEntries()<<std::endl;


  TFile* firr_syst = new TFile("./uncerts/uncert_hists_8.0.root", "read");
  //---------------------------------------------------------------------------------------



  //------------------ Output Naming ------------------------------------------------
  TString extension = etacat + "_" + iso_type + "_" + Form("%1.1f",isocut_l);
  extension += Form("_lead%d_sublead%d", (int)leadcut, (int)subleadcut);
  extension += Form("_norm%d_%d_", Commons::norm_bin.first, Commons::norm_bin.second);
  extension += filename;

  TString irrfile      = "./rootfiles/irreducible_hist_"+extension+".root";

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //
  //                    IRREDUCIBLE IRREDUCIBLE IRREDUCIBLE
  if( doIRR ){
    //---------------- Create irreducible background mgg file --------------------------------
    std::cout << "******************  IRREDUCIBLE *******************" << std::endl;
    TreeReader Rd_B(trB);


    //    TFile * f_irrslim = new TFile("irr_pt_gt_300.root","RECREATE");
    //    TTree * irrslim = trB->CloneTree();


    int nentriesB = (int)trB->GetEntriesFast();
    TH1D* hMC = new TH1D("irreducible_shape","irreducible_shape",Commons::nBins,Commons::binning);
    hMC->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");hMC->GetYaxis()->SetTitle("Events/bin");
    hMC->Sumw2();
    std::map<TString,TH1D*> hMC_sys;
    hMC_sys["pdfset"] = new TH1D("irreducible_shape_pdfset","irreducible_shape_pdfset",Commons::nBins,Commons::binning);
    hMC_sys["pdfeig"] = new TH1D("irreducible_shape_pdfeig","irreducible_shape_pdfeig",Commons::nBins,Commons::binning);
    hMC_sys["scale"]  = new TH1D("irreducible_shape_scale","irreducible_shape_scale",Commons::nBins,Commons::binning);
    hMC_sys["iso"]    = new TH1D("irreducible_shape_iso","irreducible_shape_iso",Commons::nBins,Commons::binning);
    // new syst for irr bkg covering isolation diff between data and mc
    hMC_sys["isodatamc"]    = new TH1D("irreducible_shape_isodatamc","irreducible_shape_isodatamc",Commons::nBins,Commons::binning);

    std::map<TString,TH1D*>::iterator it;
    for ( it=hMC_sys.begin() ; it != hMC_sys.end(); it++ ) ((*it).second)->Sumw2();

    TH1D* hMC_nokfac = new TH1D("irreducible_shape_nokfac","irreducible_shape_nokfac",Commons::nBins,Commons::binning);
    hMC_nokfac->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");hMC_nokfac->GetYaxis()->SetTitle("Events/bin");
    hMC_nokfac->Sumw2();

    std::cout << "START TO LOOP OVER THE MC TREE" << std::endl;
    for(int entry=0;entry<nentriesB;entry++){
      Rd_B.GetEntry(entry);
      AnalysisTools::Processing(entry,nentriesB,(int)nentriesB/100);    
      double w_1 = Rd_B.GetVariable("weight");
      double w_g = Rd_B.GetVariable("gen_weight");

      double w_nlo;
      if(leadcut==40 && subleadcut==30) w_nlo = Commons::GetkFactor_40_30(Rd_B.GetVariable("truth_mgg"));
      else if( leadcut==50 && subleadcut==50) w_nlo = Commons::GetkFactor_50_50(Rd_B.GetVariable("truth_mgg"));
      else std::cout << "ARGH !! KFACTOR IS NOT COMPUTED !!" << std::endl;

      double mgg = Rd_B.GetVariable("mgg");
      if( !Commons::EtaCategory( Rd_B.GetVariable("eta_L"),Rd_B.GetVariable("eta_SL"),etacat) ) continue;
      if( Rd_B.GetVariable("pT_L")<leadcut) continue;
      if( Rd_B.GetVariable("pT_SL")<subleadcut) continue;
      if( (int)Rd_B.GetVariable("IsTight_L")!=1) continue;
      if( (int)Rd_B.GetVariable("IsTight_SL")!=1) continue;

      //      irrslim->Fill();

      double iso_l  = Rd_B.GetVariable("Iso_L");
      double iso_sl = Rd_B.GetVariable("Iso_SL");
      if( iso_type== "TOPO_PTDEP" ){

	iso_l  = Rd_B.GetVariable("Iso_L_mod");
	iso_sl = Rd_B.GetVariable("Iso_SL_mod");

	// iso_l  = Rd_B.GetVariable("Iso_L")  - Commons::GetTopoIsoPtcorr(Rd_B.GetVariable("pT_L"));
	// iso_sl = Rd_B.GetVariable("Iso_SL") - Commons::GetTopoIsoPtcorr(Rd_B.GetVariable("pT_SL"));
	//	iso_l += Commons::GetIsoDataMCDiff(Rd_B.GetVariable("pT_L"),"L");
	//	iso_sl += Commons::GetIsoDataMCDiff(Rd_B.GetVariable("pT_SL"),"SL");
      }
      if( iso_l>isocut_l) continue;
      if( iso_sl>isocut_sl) continue;

      hMC->Fill(mgg,w_1*w_g*w_nlo);
      hMC_nokfac->Fill(mgg,w_1*w_g);
    }
    std::cout << "END THE LOOP OVER THE MC TREE" << std::endl;
    if(doSYS){
      std::cout << "FILL SYSTEMATIC UNCERTAINTIES " << std::endl;
      for ( it=hMC_sys.begin() ; it != hMC_sys.end(); it++ ){
	TString fitname = "fit_uncert_"+(*it).first;
	std::cout << fitname;
	TF1 *fit_syst;
	fit_syst = (TF1*)firr_syst->Get(fitname);
	//if(fitname.Contains("isodatamc")) fit_syst = (TF1*)firr_syst2->Get(fitname);
	//else fit_syst = (TF1*)firr_syst->Get(fitname);
		
	for(int ibin=0;ibin<(((*it).second)->GetNbinsX()+2);ibin++){
	  double syst_int = fit_syst->Integral( ((*it).second)->GetBinLowEdge(ibin),
						((*it).second)->GetBinLowEdge(ibin+1) );
	  syst_int *= 1./(*it).second->GetBinWidth(ibin);
	  std::cout<<ibin<<" "<<syst_int<<" "<<(*it).second->GetName()<<std::endl;
	  ((*it).second)->SetBinContent(ibin,hMC->GetBinContent(ibin)+syst_int*hMC->GetBinContent(ibin));
	  //	    if(fitname.Contains("scale")) std::cout<< syst_int <<" "<<hMC->GetBinContent(ibin)+syst_int*hMC->GetBinContent(ibin)<<std::endl;
	}
	std::cout<< std::endl;
      }
    }
    TFile fout_irr(irrfile,"RECREATE");
    fout_irr.Add(hMC);
    for ( it=hMC_sys.begin() ; it != hMC_sys.end(); it++ ) fout_irr.Add(hMC_sys[(*it).first]);
    fout_irr.Add(hMC_nokfac);
    fout_irr.Write();
    fout_irr.Close();

    //    f_irrslim->cd();
    //    f_irrslim->Add(irrslim);
    //    f_irrslim->Write();
    //    f_irrslim->Close();



    std::cout << "**************  END OF IRREDUCIBLE ****************" << std::endl;
    //----------------------------------------------------------------------------------------------
  }
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "***** MEMORY CLEANING ***** "<< std::endl;
  delete trB;
  //  delete fbkg;
  delete firr_syst;
  std::cout << "***** END OF THE RUN ***** "<< std::endl;

  return 0;

}


