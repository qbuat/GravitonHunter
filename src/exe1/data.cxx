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
  bool doDATA       = true;
  bool doIRR        = false;
  bool doRED        = false;
  bool doPUR        = false;
  bool doSYS        = false;
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
  TChain *trD = new TChain("tree");
  // trD->Add( Commons::outdir+"runs_19_06_13/datafile_"+iso_type+"_run*.root");//old
  trD->Add("/sps/atlas/q/qbuat/DIPHOTON_2012_REF/fuerSarahEtAl/Look/datafile_TOPO_PTDEP_run*.root");
  // trD->Add( Commons::outdir+"runs/datafile_"+iso_type_name+"_run215027.root");
  TreeReader Rd_D(trD);
  int nentriesD = (int)trD->GetEntriesFast();

  //------------------ Output Naming ------------------------------------------------
  TString extension = etacat + "_" + iso_type + "_" + Form("%1.1f",isocut_l);
  extension += Form("_lead%d_sublead%d", (int)leadcut, (int)subleadcut);
  extension += Form("_norm%d_%d_", Commons::norm_bin.first, Commons::norm_bin.second);
  extension += filename;



  TString datafile     = "./rootfiles/data_hist_"+extension+".root";
  //-------------------------------------------------------------------------


  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //
  //                        DATA DATA DATA DATA
  if( doDATA ){
    std::cout << "*************  DATA DATA DATA DATA ****************" << std::endl;

    //    TFile * f_dataslim = new TFile("data_pt_gt_300.root","RECREATE");
    //    TTree * dataslim = trD->CloneTree();

    //---> Set the latest GRL file
    DQ::SetXMLFile("../prod0/grl/data12_8TeV.periodAllYear_DetStatus-v61-pro14-02_DQDefects-00-01-00_PHYS_StandardGRL_All_Good.xml");



    //------------- Create data mgg spectrum file -------------------------------------------
    TH1D* hD = new TH1D("mgg_data","mgg_data",Commons::nBins,Commons::binning);
    hD->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    hD->GetYaxis()->SetTitle("Events/bin");
    TH1D* hD_ptgg       = AnalysisTools::GetPtggHist("ptgg_data"); 
    TH1D* hD_dphigg     = AnalysisTools::GetAbsDeltaphiHist("dphigg_data");
    TH1D* hD_costh      = AnalysisTools::GetAbsCosthetastarHist("costhetastar_data");
    TH1D* hD_deltaetagg = AnalysisTools::GetDeltaetaggHist("deltaetagg_data");
    TH1D* hD_ygg        = AnalysisTools::GetYggHist("ygg_data");  
    TH1D* hD_deltar     = AnalysisTools::GetDeltarHist("deltar_data");
    TH1D* hD_ptl        = AnalysisTools::GetPtHist("ptl_data");	  
    TH1D* hD_etal       = AnalysisTools::GetEtaHist("etal_data");  
    TH1D* hD_phil       = AnalysisTools::GetPhiHist("phil_data");  
    TH2D* hD_etaphil    = AnalysisTools::GetEtaPhiMap("etaphil_data");
    TH1D* hD_ptsl       = AnalysisTools::GetPtHist("ptsl_data");	  
    TH1D* hD_etasl      = AnalysisTools::GetEtaHist("etasl_data");  
    TH1D* hD_phisl      = AnalysisTools::GetPhiHist("phisl_data");  
    TH2D* hD_etaphisl   = AnalysisTools::GetEtaPhiMap("etaphisl_data");
    //----------- Start the loop over the entries ----------------//
    for(int entry=0;entry<nentriesD;entry++){
      Rd_D.GetEntry(entry);
      AnalysisTools::Processing(entry,nentriesD,(int)nentriesD/100);    

      //--> GRL requirement 
      if(  !DQ::PassRunLB((int)Rd_D.GetVariable("RunNumber"),
			  (int)Rd_D.GetVariable("LumiBlock")) ) continue;



      if( !Commons::EtaCategory( Rd_D.GetVariable("eta_L"),
				 Rd_D.GetVariable("eta_SL"),
				 etacat) ) continue;
      if( Rd_D.GetVariable("pT_L")<leadcut) continue;
      if( Rd_D.GetVariable("pT_SL")<subleadcut) continue;
      if( (int)Rd_D.GetVariable("IsTight_L")!=1) continue;
      if( (int)Rd_D.GetVariable("IsTight_SL")!=1) continue;

      double iso_l  = Rd_D.GetVariable("Iso_L");
      double iso_sl = Rd_D.GetVariable("Iso_SL");
      if( iso_type== "TOPO_PTDEP" ){
	iso_l  = Rd_D.GetVariable("Iso_L_mod");
	iso_sl = Rd_D.GetVariable("Iso_SL_mod");
      }

      //      if (Rd_D.GetVariable("pT_L")>300 || Rd_D.GetVariable("pT_SL")>300 ) 
      //      dataslim->Fill();
      
      
      if( iso_l>isocut_l) continue;
      if( iso_sl>isocut_sl) continue;

      double mgg = Rd_D.GetVariable("mgg");
      // if(mgg>1000.) continue; ///---> DATA BLINDING ..............
      hD->Fill(mgg);
      if( mgg> Commons::binning[Commons::norm_bin.first] ){
	hD_ptgg       ->Fill(Rd_D.GetVariable("ptgg"));
	hD_dphigg     ->Fill( fabs(Rd_D.GetVariable("deltaphi")) );
	hD_costh      ->Fill( fabs(Rd_D.GetVariable("costhetastar")) );
	hD_deltaetagg ->Fill(Rd_D.GetVariable("eta_PV_L-eta_PV_SL"));
	hD_ygg        ->Fill(Rd_D.GetVariable("ygg"));
	hD_deltar     ->Fill(Rd_D.GetVariable("sqrt(pow(eta_L-eta_SL,2)+pow(deltaphi,2))"));
	hD_ptl        ->Fill(Rd_D.GetVariable("pT_L"));
	hD_etal       ->Fill(Rd_D.GetVariable("eta_L"));
	hD_phil       ->Fill(Rd_D.GetVariable("phi_L"));
	hD_etaphil    ->Fill(Rd_D.GetVariable("eta_L"),Rd_D.GetVariable("phi_L"));
	hD_ptsl       ->Fill(Rd_D.GetVariable("pT_SL"));
	hD_etasl      ->Fill(Rd_D.GetVariable("eta_SL"));
	hD_phisl      ->Fill(Rd_D.GetVariable("phi_SL"));
	hD_etaphisl   ->Fill(Rd_D.GetVariable("eta_SL"),Rd_D.GetVariable("phi_SL"));
      }
    }
    //--------------- End of the loop over the entries -------------
    TFile fout_data(datafile,"RECREATE");
    fout_data.Add(hD            );
    fout_data.Add(hD_ptgg       );
    fout_data.Add(hD_dphigg     );
    fout_data.Add(hD_costh      );
    fout_data.Add(hD_deltaetagg );
    fout_data.Add(hD_ygg        );
    fout_data.Add(hD_deltar     );
    fout_data.Add(hD_ptl        );
    fout_data.Add(hD_etal       );
    fout_data.Add(hD_phil       );
    fout_data.Add(hD_etaphil    );
    fout_data.Add(hD_ptsl       );
    fout_data.Add(hD_etasl      );
    fout_data.Add(hD_phisl      );
    fout_data.Add(hD_etaphisl   );
    fout_data.Write();
    fout_data.Close();
    std::cout << "*************  END OF DATA *************************" << std::endl;
    //--------------------------------------------------------------------------------------- 
    
    //    f_dataslim->cd();
    //    f_dataslim->Add(dataslim);
    //    f_dataslim->Write();
    //    f_dataslim->Close();
  }


  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  std::cout << "***** MEMORY CLEANING ***** "<< std::endl;
  delete trD;
  std::cout << "***** END OF THE RUN ***** "<< std::endl;

  return 0;

}


