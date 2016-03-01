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
  if( argc<13 ){
    std::cerr << "Wrong usage ! "
	      << " analysis_new_x iso_type isocut etacat leadcut subleadcut"
	      << " doDATA doIRR doRED doPUR doSYS doFINAL filename\n" ;
    exit(1);
  }
  //------------------------------------------------------------------------



  //------------- INPUT PARAMETERS VALUES ------------------------------------
  TString iso_type  = argv[1]; // iso_type --> "TOPO" or "CONE"
  double isocut_l   = atof(argv[2]); // lead isocut value
  double isocut_sl  = atof(argv[2]); // sublead isocut value
  TString etacat    = argv[3];// eta category : --> "NONE", "CC", "CE", "EC", "EE", "EE_S", "EE_O"
  double leadcut    = atoi(argv[4]); //leading photon cut
  double subleadcut = atoi(argv[5]); //subleading photon cut
  bool doDATA       = atoi(argv[6]);
  bool doIRR        = atoi(argv[7]);
  bool doRED        = atoi(argv[8]);
  bool doPUR        = atoi(argv[9]);
  bool doSYS        = atoi(argv[10]);
  bool doFINAL      = atoi(argv[11]);
  TString filename  = argv[12];
  //--------------------------------------------------------------------------
  std::cerr << "Processing : " << argv[0] << " " << iso_type << " " 
	    << Form("%1.1f",isocut_l) << " " << etacat << " " 
	    << Form("%d %d",(int)leadcut,(int)subleadcut) << " "
	    << Form("%d %d %d %d %d %d \n",(int)doDATA,(int)doIRR,(int)doRED,(int)doPUR,(int)doSYS,(int)doFINAL)<< " " 
	    << filename;


  Commons::Setup();

  // Commons::SetPurityUncert(0.02);
  //----> Normalize between 180 and 409 (bin 26 to 47) for 50-50 config
  if( (int)leadcut == 50 && (int)subleadcut == 50) Commons::SetNormBin(26,47);

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
  trD->Add("/sps/atlas/q/qbuat/DIPHOTON_2012_REF/fuerSarahEtAl/Look/datafile_" + iso_type + "_run*.root");
  // trD->Add( Commons::outdir+"runs/datafile_"+iso_type_name+"_run215027.root");
  TreeReader Rd_D(trD);
  int nentriesD = (int)trD->GetEntriesFast();

//  TFile* fbkg = new TFile(Commons::outdir+"prod_25_06_13/bkg_mc12_gamgam_all_"+iso_type_name+".root","read");
//  TTree * trB = (TTree*)fbkg->Get("tree");


  TChain * trB = new TChain("tree");
  trB->Add("/sps/atlas/q/qbuat/DIPHOTON_2012_REF/fuerSarahEtAl/bkg_mc12_gamgam_all_TOPO_PTDEP.root");
  //  trB->Add(Commons::outdir+"prod_07_10_13/bkg_mc12_gamjet_all_"+iso_type_name+".root");

  std::cout<<trB->GetEntries()<<std::endl;



  TFile* firr_syst;
  TFile* firr_syst2;
  if(leadcut==40 && subleadcut==30) 
    firr_syst = new TFile("./rootfiles/uncert_hists.root","read");
  else if( leadcut == 50 && subleadcut == 50) {
    //    firr_syst = new TFile("./rootfiles/uncert_hists.root","read");
    firr_syst = new TFile("./rootfiles/uncert_hists_06sep13.root","read");
    firr_syst2 = new TFile(Form("./rootfiles/uncert_hists_%1.1f.root",isocut_l),"read");
  }
  else std::cout << "ARGH !! NO IRREDUCIBLE UNCERT FILE !!" << std::endl;
  //---------------------------------------------------------------------------------------



  //------------------------- Fit Parameters Files -------------------------------
    TString fitpar_ext = etacat+"_";
    fitpar_ext += iso_type_name+"_";
    fitpar_ext += Form("%d_%d",(int)leadcut,(int)subleadcut);
    TString FitParFile_L0 = "FitParameters/FitStudy_"+fitpar_ext+"_L0.config";
    TString FitParFile_L2 = "FitParameters/FitStudy_"+fitpar_ext+"_L2.config";
    TString FitParFile_L3 = "FitParameters/FitStudy_"+fitpar_ext+"_L3.config";
    TString FitParFile_L4 = "FitParameters/FitStudy_"+fitpar_ext+"_L4.config";
    TString FitParFile_L5 = "FitParameters/FitStudy_"+fitpar_ext+"_L5.config";
  //------------------------------------------------------------------------------

  //------------------ Output Naming ------------------------------------------------
  TString extension = etacat+"_"+iso_type+"_"+Form("%1.1f",isocut_l);
  extension += Form("_lead%d_sublead%d",(int)leadcut,(int)subleadcut);
  extension += Form("_norm%d_%d",Commons::norm_bin.first,Commons::norm_bin.second);
  extension += filename;



  TString datafile     = "./rootfiles/data_hist_"+extension+".root";
  TString irrfile      = "./rootfiles/irreducible_hist_"+extension+".root";
  //------------------------------------------------------------------------
  TString redfile      = "./rootfiles/reducible_hist_"+extension+".root";
  TString redfile_L3   = "./rootfiles/reducible_hist_"+extension+"_L3.root";
  TString redfile_L4   = "./rootfiles/reducible_hist_"+extension+"_L4.root";
  TString redfile_L5   = "./rootfiles/reducible_hist_"+extension+"_L5.root";
  //------------------------------------------------------------------------
  TString yieldfile_L2 = "./rootfiles/2DFitResult_"+extension+"_L2.root";
  TString yieldfile_L3 = "./rootfiles/2DFitResult_"+extension+"_L3.root";
  TString yieldfile_L4 = "./rootfiles/2DFitResult_"+extension+"_L4.root";
  TString yieldfile_L5 = "./rootfiles/2DFitResult_"+extension+"_L5.root";
  TString fit2dfile_L2 = "./rootfiles/bkg2dfit_"+extension+"_L2.root";
  TString fit2dfile_L3 = "./rootfiles/bkg2dfit_"+extension+"_L3.root";
  TString fit2dfile_L4 = "./rootfiles/bkg2dfit_"+extension+"_L4.root";
  TString fit2dfile_L5 = "./rootfiles/bkg2dfit_"+extension+"_L5.root";
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
	fit_syst = (TF1*)firr_syst2->Get(fitname);
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
  //
  //                    REDUCIBLE REDUCIBLE REDUCIBLE

  if( doRED ){
    std::cout << "********************  REDUCIBLE *******************" << std::endl;

    //--------------- Create Reducible background mgg file --------------------------------------
    std::pair<double,double> fitrange( Commons::binning[Commons::norm_bin.first-1],1000 );
    std::pair<double,double> mggbin( Commons::binning[Commons::norm_bin.first-1],1000);

    //    pair<double,double> fitrange( Commons::binning[Commons::norm_bin.first-1],1000 );
    //    pair<double,double> mggbin( Commons::binning[Commons::norm_bin.first-1],3000);
    
    BkgExtrapolation * BE = new BkgExtrapolation();
    BE->SetTree(trD);
    BE->SetFitRange();
    BE->SetNPE(1000);
    BE->SetMggBounds(Commons::binning[0],Commons::binning[Commons::nBins]);
    BE->SetPhotonsIsoCut(isocut_l,isocut_sl);
    BE->SetPhotonsPtCut(leadcut,subleadcut);
    BE->SetEtaCategory(etacat);
    if( iso_type == "TOPO_PTDEP" )
      BE->SetPtDependentIsocut(true);

    BE->Init(mggbin,FitParFile_L0);
    BE->SetFitRange(fitrange.first,fitrange.second);
    BE->Fitter(true);
    BE->Extrapolate(true);
    BE->StoreToRootFile(redfile);
    //    BE->StoreToRootFile("./rootfiles/reducible_hist_"+extension+"data_extended_mgg_range_redfit.root");

    delete BE;
    //----------------------------------------------------------------------------------------------
    
    if( doSYS ){
      //--------------- Create Reducible systematic uncertainty background mgg file ------------------- 
      BkgExtrapolation * BE_SYS = new BkgExtrapolation();
      BE_SYS->SetTree(trD);
      BE_SYS->SetFitRange();
      BE_SYS->SetNPE(1);
      BE_SYS->SetMggBounds(Commons::binning[0],Commons::binning[Commons::nBins]);
      BE_SYS->SetPhotonsIsoCut(isocut_l,isocut_sl);
      BE_SYS->SetPhotonsPtCut(leadcut,subleadcut);
      BE_SYS->SetEtaCategory(etacat);
      if( iso_type== "TOPO_PTDEP")
	BE_SYS->SetPtDependentIsocut(true);
      BE_SYS->SetLoosePrimeType(3);
      BE_SYS->Init(mggbin,FitParFile_L3);
      BE_SYS->SetFitRange(fitrange.first,fitrange.second);
      BE_SYS->Fitter(true);
      BE_SYS->Extrapolate(true);
      BE_SYS->StoreToRootFile(redfile_L3);
      //-----------------------------------------------------------------------------
      BE_SYS = new BkgExtrapolation();
      BE_SYS->SetTree(trD);
      BE_SYS->SetFitRange();
      BE_SYS->SetNPE(1);
      BE_SYS->SetMggBounds(Commons::binning[0],Commons::binning[Commons::nBins]);
      BE_SYS->SetPhotonsIsoCut(isocut_l,isocut_sl);
      BE_SYS->SetPhotonsPtCut(leadcut,subleadcut);
      BE_SYS->SetEtaCategory(etacat);
      if( iso_type== "TOPO_PTDEP" )
	BE_SYS->SetPtDependentIsocut(true);
      BE_SYS->SetLoosePrimeType(4);
      BE_SYS->Init(mggbin,FitParFile_L4);
      BE_SYS->SetFitRange(fitrange.first,fitrange.second);
      BE_SYS->Fitter(true);
      BE_SYS->Extrapolate(true);
      BE_SYS->StoreToRootFile(redfile_L4);
      //-----------------------------------------------------------------------------
      BE_SYS = new BkgExtrapolation();
      BE_SYS->SetTree(trD);
      BE_SYS->SetFitRange();
      BE_SYS->SetNPE(1);
      BE_SYS->SetMggBounds(Commons::binning[0],Commons::binning[Commons::nBins]);
      BE_SYS->SetPhotonsIsoCut(isocut_l,isocut_sl);
      BE_SYS->SetPhotonsPtCut(leadcut,subleadcut);
      BE_SYS->SetEtaCategory(etacat);
      if( iso_type== "TOPO_PTDEP" )
	BE_SYS->SetPtDependentIsocut(true);
      BE_SYS->SetLoosePrimeType(5);
      BE_SYS->Init(mggbin,FitParFile_L5);
      BE_SYS->SetFitRange(fitrange.first,fitrange.second);
      BE_SYS->Fitter(true);
      BE_SYS->Extrapolate(true);
      BE_SYS->StoreToRootFile(redfile_L5);
      delete BE_SYS;
      //----------------------------------------------------------------------------------------------
    }
    std::cout << "********************  END OF REDUCIBLE *******************" << std::endl;
  }
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //
  //                    PURITY PURITY PURITY PURITY
  if( doPUR ){
    std::cout << "********************  SAMPLE COMPOSITION (gg/gj/jg/jj) *******************" << std::endl;
    //---------------- Determine background composition -------------------------------------------
    std::pair<double,double> mggbin_2dfit( Commons::binning[Commons::norm_bin.first-1],
					   Commons::binning[Commons::norm_bin.second] );
    Bkg2DFit * RB2DF= new Bkg2DFit();
    RB2DF->SetTree(trD);
    RB2DF->SetLoosePrimeType(4);
    RB2DF->SetEtaCategory(etacat);
    RB2DF->SetPhotonsPtCut(leadcut,subleadcut);
    RB2DF->SetPhotonsIsoCut(isocut_l,isocut_sl);
    if( iso_type.Contains("CONE") ){
      RB2DF->SetIsoBounds(-10,25);
      RB2DF->SetIsoNorm(10,25);
    } else if( iso_type.Contains("TOPO") ){
      RB2DF->SetIsoBounds(-4,14);
      RB2DF->SetIsoNorm(11,14);
    }
    if( iso_type== "TOPO_PTDEP")
	RB2DF->SetPtDependentIsocut(true);
    RB2DF->Init(mggbin_2dfit,FitParFile_L4);
    bool do1Dfits    = true;
    bool do2Dfits    = true;
    bool doSplot     = false;
    RB2DF->Fitter(do1Dfits,do2Dfits,doSplot);
    RB2DF->RandomizeYieldsResults(20000,yieldfile_L4);
    RB2DF->StoreToRootFile(fit2dfile_L4);
    //---------------------------------------------
    std::cout << "Purity (%) = "  << 100*RB2DF->GetPurity() 
	      << " +/- " << 100*RB2DF->GetPurityError() 
	      << std::endl;
    //---------------------------------------------
    delete RB2DF;

    if( doSYS ){
      Bkg2DFit * RB2DF_sys;
      //------------------------------------------------------------------------------------------------
      RB2DF_sys = new Bkg2DFit();
      RB2DF_sys->SetTree(trD);
      RB2DF_sys->SetLoosePrimeType(2);
      RB2DF_sys->SetEtaCategory(etacat);
      RB2DF_sys->SetPhotonsPtCut(leadcut,subleadcut);
      RB2DF_sys->SetPhotonsIsoCut(isocut_l,isocut_sl);
      if( iso_type.Contains("CONE") ){
	RB2DF_sys->SetIsoBounds(-10,25);
	RB2DF_sys->SetIsoNorm(10,25);
      } else if( iso_type.Contains("TOPO") ){
	RB2DF_sys->SetIsoBounds(-4,14);
	RB2DF_sys->SetIsoNorm(11,14);
      } 
      if( iso_type== "TOPO_PTDEP")
	RB2DF_sys->SetPtDependentIsocut(true);
      RB2DF_sys->Init(mggbin_2dfit,FitParFile_L2);
      RB2DF_sys->Fitter(true,true,false);
      RB2DF_sys->RandomizeYieldsResults(20000,yieldfile_L2);
      RB2DF_sys->StoreToRootFile(fit2dfile_L2);
      //------------------------------------------------------------------------------------------------
      RB2DF_sys = new Bkg2DFit();
      RB2DF_sys->SetTree(trD);
      RB2DF_sys->SetLoosePrimeType(3);
      RB2DF_sys->SetEtaCategory(etacat);
      RB2DF_sys->SetPhotonsPtCut(leadcut,subleadcut);
      RB2DF_sys->SetPhotonsIsoCut(isocut_l,isocut_sl);
      if( iso_type.Contains("CONE") ){
	RB2DF_sys->SetIsoBounds(-10,25);
	RB2DF_sys->SetIsoNorm(10,25);
      } else if( iso_type.Contains("TOPO") ){
	RB2DF_sys->SetIsoBounds(-4,14);
	RB2DF_sys->SetIsoNorm(11,14);
      }
      RB2DF_sys->Init(mggbin_2dfit,FitParFile_L3);
      RB2DF_sys->Fitter(true,true,false);
      RB2DF_sys->RandomizeYieldsResults(20000,yieldfile_L3);
      RB2DF_sys->StoreToRootFile(fit2dfile_L3);
      //------------------------------------------------------------------------------------------------
      RB2DF_sys = new Bkg2DFit();
      RB2DF_sys->SetTree(trD);
      RB2DF_sys->SetLoosePrimeType(5);
      RB2DF_sys->SetEtaCategory(etacat);
      RB2DF_sys->SetPhotonsPtCut(leadcut,subleadcut);
      RB2DF_sys->SetPhotonsIsoCut(isocut_l,isocut_sl);
      if( iso_type.Contains("CONE") ){
	RB2DF_sys->SetIsoBounds(-10,25);
	RB2DF_sys->SetIsoNorm(10,25);
      } else if( iso_type.Contains("TOPO") ){
	RB2DF_sys->SetIsoBounds(-4,14);
	RB2DF_sys->SetIsoNorm(11,14);
      }
      RB2DF_sys->Init(mggbin_2dfit,FitParFile_L5);
      RB2DF_sys->Fitter(true,true,false);
      RB2DF_sys->RandomizeYieldsResults(20000,yieldfile_L5);
      RB2DF_sys->StoreToRootFile(fit2dfile_L5);
      delete RB2DF_sys;
    }
    std::cout << "*************  END OF SAMPLE COMPOSITION (gg/gj/jg/jj) *******************" << std::endl;



  }
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //
  //                    FINAL ESTIMATE-FINAL ESTIMATE-FINAL ESTIMATE 

  

  if( doFINAL ){
    std::cout << "********************  FINAL ESTIMATE *******************" << std::endl;
    //----------- Final Background/data plot and output ---------------------------------------------
    BkgEstimatorUpgrade * BES = new BkgEstimatorUpgrade();
    BES->SetVerbose(true);
    BES->SetUncertFlag(doSYS);
    BES->SetDataFile(datafile);
    BES->SetIrrFile(irrfile) ;
    BES->SetRedFile(redfile) ;
    BES->SetYieldFile(yieldfile_L4); 
    if(doSYS){
      BES->SetRedSystFiles(redfile_L4,redfile_L5);
      BES->SetYieldSystFiles(yieldfile_L2,yieldfile_L3,yieldfile_L5);
    }
    BES->Init();
    BES->CreateBATInput("./rootfiles/Bkg_Template_"+extension+".root");
    BES->CreateBumpHunterInput("./rootfiles/Bkg_Total_"+extension+".root");  
    BES->GetFinalYieldsPerBin_NoUncert();
    delete BES;
    std::cout << "***************** END OF FINAL ESTIMATE ****************" << std::endl;
    //------------------------------------------------------------------------------------------------
  }
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "***** MEMORY CLEANING ***** "<< std::endl;
  delete trD;
  //  delete fbkg;
  delete firr_syst;
  std::cout << "***** END OF THE RUN ***** "<< std::endl;

  return 0;

}


