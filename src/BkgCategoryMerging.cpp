#include "BkgCategoryMerging.h"
#include "ToolsUtilities.h"
#include <TFile.h>
#include <TError.h>
#include <TMath.h>
#include <RooHist.h>

//////////////////////////////////////////////////////////////
BkgCategoryMerging::BkgCategoryMerging(TString type,
				       TString f0,TString f1,
				       TString f2,TString f3,
				       TString f4,TString f5,
				       TString f6,TString f7)

/////////////////////////////////////////////////////////////
{
  //Default constructor
  SetVerbose(false);
  SetFiles(f0,f1,f2,f3,f4,f5,f6,f7);
  m_type = type;
  if(m_type == "full"){
    BuildHistograms();
    BuildGraphs();
  }else if(m_type == "search")
    BuildHistograms_Search();
  else if( m_type == "2DFit")
    BuildRooPlots_2DFit();
  else if( m_type == "RedFit")
    BuildRooPlots_RedFit();
  else Fatal("BkgCategoryMerging::BkgCategoryMerging()","Wrong type !!!");
}
//////////////////////////////////////////
BkgCategoryMerging::~BkgCategoryMerging()
////////////////////////////////////////
{

}


///////////////////////////////////////////////////////
void BkgCategoryMerging::SetFiles(TString f0,TString f1,
				  TString f2,TString f3,
				  TString f4,TString f5,
				  TString f6,TString f7)
///////////////////////////////////////////////////////
{

  if( f0=="") Fatal("BkgCategoryMerging::SetFiles()",
		    "You need to specify at least one file !");

  else m_filenames.push_back(f0);

  if(f1 !="") m_filenames.push_back(f1);
  if(f2 !="") m_filenames.push_back(f2);
  if(f3 !="") m_filenames.push_back(f3);
  if(f4 !="") m_filenames.push_back(f4);
  if(f5 !="") m_filenames.push_back(f5);
  if(f6 !="") m_filenames.push_back(f6);
  if(f7 !="") m_filenames.push_back(f7);


}
/////////////////////////////////////////
bool BkgCategoryMerging::BuildHistograms()
/////////////////////////////////////////
{
  TFile *file0 = TFile::Open(m_filenames[0]);
  std::cout << "Open: " << m_filenames[0] << std::endl;
  m_hh_data         = (TH1D*)file0->Get("mgg_data");
  m_hh_bkg          = (TH1D*)file0->Get("bkg_total_gg_full");
  m_hh_red          = (TH1D*)file0->Get("bkg_reducible_gg_full");
  m_hh_gj           = (TH1D*)file0->Get("bkg_gammajet_gg_full");
  m_hh_jg           = (TH1D*)file0->Get("bkg_jetgamma_gg_full");
  m_hh_jj           = (TH1D*)file0->Get("bkg_jetjet_gg_full");
  m_hh_gg           = (TH1D*)file0->Get("bkg_irreducible_gg_full");
  // m_hh_tot_uncert   = (TH1D*)file0->Get("bkg_total_syst_gg_full");
  // m_hh_pur_uncert   = (TH1D*)file0->Get("bkg_purity_syst_gg_full");
  // m_hh_red_uncert   = (TH1D*)file0->Get("bkg_reducible_syst_gg_full");
  // m_hh_irr_uncert   = (TH1D*)file0->Get("bkg_irreducible_syst_gg_full");
  // m_hh_pdf_uncert   = (TH1D*)file0->Get("bkg_pdf_syst_gg_full");
  // m_hh_iso_uncert   = (TH1D*)file0->Get("bkg_iso_syst_gg_full");
  // m_hh_scale_uncert = (TH1D*)file0->Get("bkg_scale_syst_gg_full");
  // m_hh_bkg_syst     = (TH1D*)file0->Get("hh_bkg_syst_gg");
  // file0->Close();


  if( m_filenames.size() == 1){
    std::cout << "Only 1 file is considered for the merging !" << std::endl;
    return 0;
  }else{
    for(int i=1;i<(int)m_filenames.size();i++){
      std::cout << "Open: " << m_filenames[i] << std::endl;
      TFile *file_i = TFile::Open(m_filenames[i]);
      TH1D* temp_hh_data = (TH1D*)file_i->Get("mgg_data");
      TH1D* temp_hh_bkg  = (TH1D*)file_i->Get("bkg_total_gg_full");
      TH1D* temp_hh_red  = (TH1D*)file_i->Get("bkg_reducible_gg_full");
      TH1D* temp_hh_gj   = (TH1D*)file_i->Get("bkg_gammajet_gg_full");
      TH1D* temp_hh_jg   = (TH1D*)file_i->Get("bkg_jetgamma_gg_full");
      TH1D* temp_hh_jj   = (TH1D*)file_i->Get("bkg_jetjet_gg_full");
      TH1D* temp_hh_gg   = (TH1D*)file_i->Get("bkg_irreducible_gg_full");
      // TH1D* temp_hh_tot_uncert   = (TH1D*)file_i->Get("bkg_total_syst_gg_full");
      // TH1D* temp_hh_pur_uncert   = (TH1D*)file_i->Get("bkg_purity_syst_gg_full");
      // TH1D* temp_hh_red_uncert   = (TH1D*)file_i->Get("bkg_reducible_syst_gg_full");
      // TH1D* temp_hh_irr_uncert   = (TH1D*)file_i->Get("bkg_irreducible_syst_gg_full");
      // TH1D* temp_hh_pdf_uncert   = (TH1D*)file_i->Get("bkg_pdf_syst_gg_full");
      // TH1D* temp_hh_iso_uncert   = (TH1D*)file_i->Get("bkg_iso_syst_gg_full");
      // TH1D* temp_hh_scale_uncert = (TH1D*)file_i->Get("bkg_scale_syst_gg_full");
      // TH1D* temp_hh_bkg_syst     = (TH1D*)file_i->Get("hh_bkg_syst_gg");
      // file_i->Close();
      m_hh_data->Add(temp_hh_data);
      m_hh_bkg->Add(temp_hh_bkg);
      m_hh_red->Add(temp_hh_red);
      m_hh_gg->Add(temp_hh_gg);
      m_hh_gj->Add(temp_hh_gj);
      m_hh_jg->Add(temp_hh_jg);
      m_hh_jj->Add(temp_hh_jj);
      // m_hh_bkg_syst->Add(temp_hh_bkg_syst);
      // m_hh_pur_uncert = AnalysisTools::HistAddedUncert(*m_hh_pur_uncert,*temp_hh_pur_uncert);
      // m_hh_irr_uncert = AnalysisTools::HistAddedUncert(*m_hh_irr_uncert,*temp_hh_irr_uncert);
      // m_hh_pdf_uncert = AnalysisTools::HistAddedUncert(*m_hh_pdf_uncert,*temp_hh_pdf_uncert);
      // m_hh_iso_uncert = AnalysisTools::HistAddedUncert(*m_hh_iso_uncert,*temp_hh_iso_uncert);
      // m_hh_scale_uncert =  AnalysisTools::HistAddedUncert(*m_hh_scale_uncert,*temp_hh_scale_uncert);
      // m_hh_red_uncert =  AnalysisTools::HistAddedUncert(*m_hh_red_uncert,*temp_hh_red_uncert);
      // m_hh_tot_uncert =  AnalysisTools::HistAddedUncert(*m_hh_tot_uncert,*temp_hh_tot_uncert);

    } 
    return 1;
  }
}
/////////////////////////////////////////////////
bool BkgCategoryMerging::BuildHistograms_Search()
/////////////////////////////////////////////////
{
  TFile *file0 = TFile::Open(m_filenames[0]);
  std::cout << "Open: " << m_filenames[0] << std::endl;
  m_hh_data       = (TH1D*)file0->Get("hh_data");
  m_hh_bkg        = (TH1D*)file0->Get("bkg_total_gg");
  m_hh_red        = (TH1D*)file0->Get("bkg_reductible_gg");
  m_hh_gg         = (TH1D*)file0->Get("bkg_irreductible_gg");
  m_hh_pur_uncert = (TH1D*)file0->Get("bkg_purity_syst_gg");
  m_hh_irr_uncert = (TH1D*)file0->Get("bkg_irreducible_syst_gg");
  m_hh_red_uncert = (TH1D*)file0->Get("bkg_reducible_syst_gg");
  m_hh_tot_uncert = (TH1D*)file0->Get("bkg_total_syst_gg");
  // file0->Close();




  if( m_filenames.size() == 1){
    std::cout << "Only 1 file is considered for the merging !" << std::endl;
    return 0;
  }else{
    for(int i=1;i<(int)m_filenames.size();i++){
      std::cout << "Open: " << m_filenames[i] << std::endl;

      TFile *file_i = TFile::Open(m_filenames[i]);
      TH1D* temp_hh_data       = (TH1D*)file_i->Get("hh_data");
      TH1D* temp_hh_bkg        = (TH1D*)file_i->Get("bkg_total_gg");
      TH1D* temp_hh_red        = (TH1D*)file_i->Get("bkg_reductible_gg");
      TH1D* temp_hh_gg         = (TH1D*)file_i->Get("bkg_irreductible_gg");
      TH1D* temp_hh_pur_uncert = (TH1D*)file_i->Get("bkg_purity_syst_gg");
      TH1D* temp_hh_irr_uncert = (TH1D*)file_i->Get("bkg_irreducible_syst_gg");
      TH1D* temp_hh_red_uncert = (TH1D*)file_i->Get("bkg_reducible_syst_gg");
      TH1D* temp_hh_tot_uncert = (TH1D*)file_i->Get("bkg_total_syst_gg");
      // file_i->Close();
      m_hh_data->Add(temp_hh_data);
      m_hh_bkg->Add(temp_hh_bkg);
      m_hh_red->Add(temp_hh_red);
      m_hh_gg->Add(temp_hh_gg);
      m_hh_pur_uncert =  AnalysisTools::HistAddedUncert(*m_hh_pur_uncert,*temp_hh_pur_uncert);
      m_hh_irr_uncert =  AnalysisTools::HistAddedUncert(*m_hh_irr_uncert,*temp_hh_irr_uncert);
      m_hh_red_uncert =  AnalysisTools::HistAddedUncert(*m_hh_red_uncert,*temp_hh_red_uncert);
      m_hh_tot_uncert =  AnalysisTools::HistAddedUncert(*m_hh_tot_uncert,*temp_hh_tot_uncert);
    } 
    return 1;
  }
}

//////////////////////////////////////
bool BkgCategoryMerging::BuildGraphs()
//////////////////////////////////////
{
  m_gmgg_data = AnalysisTools::DataGraph(*m_hh_data);
  m_gmgg_bkg  = AnalysisTools::BkgGraph(*m_hh_bkg);
  m_gmgg_red  = AnalysisTools::BkgGraph(*m_hh_red);
  m_gmgg_gg   = AnalysisTools::BkgGraph(*m_hh_gg);
  m_gmgg_gj   = AnalysisTools::BkgGraph(*m_hh_gj);
  m_gmgg_jg   = AnalysisTools::BkgGraph(*m_hh_jg);
  m_gmgg_jj   = AnalysisTools::BkgGraph(*m_hh_jj);
  return 1;
}

//////////////////////////////////////////////
bool BkgCategoryMerging::BuildRooPlots_2DFit()
/////////////////////////////////////////////
{
  TFile *file0 = TFile::Open(m_filenames[0]);
  std::cout << "Open: " << m_filenames[0] << std::endl;


  if( m_filenames.size() == 1){
    frame_TiTi_L  = (RooPlot*)file0->Get("TiTiL");
    frame_TiTi_SL = (RooPlot*)file0->Get("TiTiSL");
    frame_PH_L  = (RooPlot*)file0->Get("PHL");
    frame_PH_SL = (RooPlot*)file0->Get("PHSL");
    frame_Jet_L  = (RooPlot*)file0->Get("JetL");
    frame_Jet_SL = (RooPlot*)file0->Get("JetSL");
    // file0->Close();

    std::cout << "Only 1 file is considered for the merging !" << std::endl;
    return 0;
  }else{
    //==============================================================
    RooPlot* frame0_TiTi_L  = (RooPlot*)file0->Get("TiTiL");
    RooHist*  data_L       = frame0_TiTi_L->getHist("data_L");
    RooCurve* pdf_tot_L    = frame0_TiTi_L->getCurve("pdf_tot");
    RooCurve* pdf_phph_L   = frame0_TiTi_L->getCurve("pdf_phph");
    RooCurve* pdf_phjet_L  = frame0_TiTi_L->getCurve("pdf_phjet");
    RooCurve* pdf_jetph_L  = frame0_TiTi_L->getCurve("pdf_jetph");
    RooCurve* pdf_jetjet_L = frame0_TiTi_L->getCurve("pdf_jetjet");
    frame_TiTi_L  = frame0_TiTi_L->emptyClone("TiTiL");;
    //---------------------------------------------------------------
    RooPlot* frame0_TiTi_SL  = (RooPlot*)file0->Get("TiTiSL");
    RooHist*  data_SL       = frame0_TiTi_SL->getHist("data_SL");
    RooCurve* pdf_tot_SL    = frame0_TiTi_SL->getCurve("pdf_tot");
    RooCurve* pdf_phph_SL   = frame0_TiTi_SL->getCurve("pdf_phph");
    RooCurve* pdf_phjet_SL  = frame0_TiTi_SL->getCurve("pdf_phjet");
    RooCurve* pdf_jetph_SL  = frame0_TiTi_SL->getCurve("pdf_jetph");
    RooCurve* pdf_jetjet_SL = frame0_TiTi_SL->getCurve("pdf_jetjet");
    frame_TiTi_SL  = frame0_TiTi_SL->emptyClone("TiTiSL");
    //=================================================================
    RooPlot*  frame0_PH_L  = (RooPlot*)file0->Get("PHL");
    RooHist*  data_Ti_L    = frame0_PH_L->getHist("data_TiL");
    RooCurve* pdf_PHL      = frame0_PH_L->getCurve("pdf_PHL");
    RooCurve* pdf_JL       = frame0_PH_L->getCurve("pdf_JL");
    RooCurve* pdf_TiL      = frame0_PH_L->getCurve("pdf_TiL");
    frame_PH_L = frame0_PH_L->emptyClone("PHL");
    //--------------------------------------------------------------------
    RooPlot*  frame0_PH_SL  = (RooPlot*)file0->Get("PHSL");
    RooHist*  data_Ti_SL    = frame0_PH_SL->getHist("data_TiSL");
    RooCurve* pdf_PHSL      = frame0_PH_SL->getCurve("pdf_PHSL");
    RooCurve* pdf_JSL       = frame0_PH_SL->getCurve("pdf_JSL");
    RooCurve* pdf_TiSL      = frame0_PH_SL->getCurve("pdf_TiSL");
    frame_PH_SL = frame0_PH_SL->emptyClone("PHSL");
    //=================================================================
    RooPlot* frame0_Jet_L = (RooPlot*)file0->Get("JetL");
    RooHist* data_JetL  = frame0_Jet_L->getHist("data_JL");
    RooCurve* pdf_JetL  = frame0_Jet_L->getCurve("pdf_JL");
    frame_Jet_L = frame0_Jet_L->emptyClone("JetL");
    //--------------------------------------------------------------------
    RooPlot* frame0_Jet_SL = (RooPlot*)file0->Get("JetSL");
    RooHist* data_JetSL  = frame0_Jet_SL->getHist("data_JSL");
    RooCurve* pdf_JetSL  = frame0_Jet_SL->getCurve("pdf_JSL");
    frame_Jet_SL = frame0_Jet_SL->emptyClone("JetSL");
    //=================================================================
    // file0->Close();

    for(int i=1;i<(int)m_filenames.size();i++){
      std::cout << "Open: " << m_filenames[i] << std::endl;
      TFile *file_i = TFile::Open(m_filenames[i]);
      //=================================================================
      RooPlot* frame_i_TiTi_L  = (RooPlot*)file_i->Get("TiTiL");
      RooHist*  data_L_i       = frame_i_TiTi_L->getHist("data_L");
      RooCurve* pdf_tot_L_i    = frame_i_TiTi_L->getCurve("pdf_tot");
      RooCurve* pdf_phph_L_i   = frame_i_TiTi_L->getCurve("pdf_phph");
      RooCurve* pdf_phjet_L_i  = frame_i_TiTi_L->getCurve("pdf_phjet");
      RooCurve* pdf_jetph_L_i  = frame_i_TiTi_L->getCurve("pdf_jetph");
      RooCurve* pdf_jetjet_L_i = frame_i_TiTi_L->getCurve("pdf_jetjet");
      data_L = new RooHist(*data_L,*data_L_i);
      pdf_tot_L    = new RooCurve(pdf_tot_L->GetName(),pdf_tot_L->GetTitle(),
				  *pdf_tot_L    ,*pdf_tot_L_i    );
      pdf_phph_L   = new RooCurve(pdf_phph_L->GetName(),pdf_phph_L->GetTitle(),
				    *pdf_phph_L   ,*pdf_phph_L_i   );
      pdf_phjet_L  = new RooCurve(pdf_phjet_L->GetName(),pdf_phjet_L->GetTitle(),
				  *pdf_phjet_L  ,*pdf_phjet_L_i  );
      pdf_jetph_L  = new RooCurve(pdf_jetph_L->GetName(),pdf_jetph_L->GetTitle(),
				    *pdf_jetph_L  ,*pdf_jetph_L_i  );
      pdf_jetjet_L = new RooCurve(pdf_jetjet_L->GetName(),pdf_jetjet_L->GetTitle(),
				    *pdf_jetjet_L ,*pdf_jetjet_L_i );
      //----------------------------------------------------------------------
      RooPlot* frame_i_TiTi_SL  = (RooPlot*)file_i->Get("TiTiSL");
      RooHist*  data_SL_i       = frame_i_TiTi_SL->getHist("data_SL");
      RooCurve* pdf_tot_SL_i    = frame_i_TiTi_SL->getCurve("pdf_tot");
      RooCurve* pdf_phph_SL_i   = frame_i_TiTi_SL->getCurve("pdf_phph");
      RooCurve* pdf_phjet_SL_i  = frame_i_TiTi_SL->getCurve("pdf_phjet");
      RooCurve* pdf_jetph_SL_i  = frame_i_TiTi_SL->getCurve("pdf_jetph");
      RooCurve* pdf_jetjet_SL_i = frame_i_TiTi_SL->getCurve("pdf_jetjet");
      data_SL = new RooHist(*data_SL,*data_SL_i);
      pdf_tot_SL    = new RooCurve(pdf_tot_SL->GetName(),pdf_tot_SL->GetTitle(),
				    *pdf_tot_SL    ,*pdf_tot_SL_i    );
      pdf_phph_SL   = new RooCurve(pdf_phph_SL->GetName(),pdf_phph_SL->GetTitle(),
				    *pdf_phph_SL   ,*pdf_phph_SL_i   );
      pdf_phjet_SL  = new RooCurve(pdf_phjet_SL->GetName(),pdf_phjet_SL->GetTitle(),
				    *pdf_phjet_SL  ,*pdf_phjet_SL_i  );
      pdf_jetph_SL  = new RooCurve(pdf_jetph_SL->GetName(),pdf_jetph_SL->GetTitle(),
				    *pdf_jetph_SL  ,*pdf_jetph_SL_i  );
      pdf_jetjet_SL = new RooCurve(pdf_jetjet_SL->GetName(),pdf_jetjet_SL->GetTitle(),
				    *pdf_jetjet_SL ,*pdf_jetjet_SL_i );
      //=================================================================
      RooPlot*  frame_i_PH_L   = (RooPlot*)file_i->Get("PHL");
      RooHist*  data_Ti_L_i    = frame_i_PH_L->getHist("data_TiL");
      RooCurve* pdf_PHL_i      = frame_i_PH_L->getCurve("pdf_PHL");
      RooCurve* pdf_JL_i       = frame_i_PH_L->getCurve("pdf_JL");
      RooCurve* pdf_TiL_i      = frame_i_PH_L->getCurve("pdf_TiL");
      data_Ti_L = new RooHist(*data_Ti_L,*data_Ti_L_i);
      pdf_PHL   = new RooCurve(pdf_PHL->GetName(),pdf_PHL->GetTitle(),
			       *pdf_PHL,*pdf_PHL_i);
      pdf_JL    = new RooCurve(pdf_JL->GetName(),pdf_JL->GetTitle(),
			       *pdf_JL,*pdf_JL_i);
      pdf_TiL   = new RooCurve(pdf_TiL->GetName(),pdf_TiL->GetTitle(),
			       *pdf_TiL,*pdf_TiL_i);
      //----------------------------------------------------------------
      RooPlot*  frame_i_PH_SL   = (RooPlot*)file_i->Get("PHSL");
      RooHist*  data_Ti_SL_i    = frame_i_PH_SL->getHist("data_TiSL");
      RooCurve* pdf_PHSL_i      = frame_i_PH_SL->getCurve("pdf_PHSL");
      RooCurve* pdf_JSL_i       = frame_i_PH_SL->getCurve("pdf_JSL");
      RooCurve* pdf_TiSL_i      = frame_i_PH_SL->getCurve("pdf_TiSL");
      data_Ti_SL = new RooHist(*data_Ti_SL,*data_Ti_SL_i);
      pdf_PHSL   = new RooCurve(pdf_PHSL->GetName(),pdf_PHSL->GetTitle(),
			       *pdf_PHSL,*pdf_PHSL_i);
      pdf_JSL    = new RooCurve(pdf_JSL->GetName(),pdf_JSL->GetTitle(),
			       *pdf_JSL,*pdf_JSL_i);
      pdf_TiSL   = new RooCurve(pdf_TiSL->GetName(),pdf_TiSL->GetTitle(),
			       *pdf_TiSL,*pdf_TiSL_i);
      //=====================================================================
      RooPlot* frame_i_Jet_L = (RooPlot*)file_i->Get("JetL");
      RooHist*  data_JetL_i     = frame_i_Jet_L->getHist("data_JL");
      RooCurve* pdf_JetL_i     = frame_i_Jet_L->getCurve("pdf_JL");
      data_JetL = new RooHist(*data_JetL,*data_JetL_i);
      pdf_JetL  = new RooCurve(pdf_JetL->GetName(),pdf_JetL->GetTitle(),
			       *pdf_JetL,*pdf_JetL_i); 
      //--------------------------------------------------------------------
      RooPlot* frame_i_Jet_SL = (RooPlot*)file_i->Get("JetSL");
      RooHist*  data_JetSL_i     = frame_i_Jet_SL->getHist("data_JSL");
      RooCurve* pdf_JetSL_i     = frame_i_Jet_SL->getCurve("pdf_JSL");
      data_JetSL = new RooHist(*data_JetSL,*data_JetSL_i);
      pdf_JetSL  = new RooCurve(pdf_JetSL->GetName(),pdf_JetSL->GetTitle(),
				*pdf_JetSL,*pdf_JetSL_i); 
      //=================================================================
      // file_i->Close();
    
    }

    data_L->SetMarkerSize(1);
    pdf_tot_L->SetLineColor(kBlack);
    pdf_phph_L->SetLineColor(kRed);
    pdf_phjet_L->SetLineColor(kRed);
    pdf_jetjet_L->SetLineColor(kGreen);
    pdf_jetph_L->SetLineColor(kGreen);
    pdf_phjet_L->SetLineStyle(kDashed);
    pdf_jetph_L->SetLineStyle(kDashed);
    pdf_tot_L->SetLineColor(kBlack);
    frame_TiTi_L->addPlotable(data_L       ,"P");
    frame_TiTi_L->addPlotable(pdf_tot_L    );
    frame_TiTi_L->addPlotable(pdf_phph_L   );
    frame_TiTi_L->addPlotable(pdf_phjet_L  );
    frame_TiTi_L->addPlotable(pdf_jetph_L  );
    frame_TiTi_L->addPlotable(pdf_jetjet_L );

    data_SL->SetMarkerSize(1);
    pdf_tot_SL->SetLineColor(kBlack);
    pdf_phph_SL->SetLineColor(kRed);
    pdf_jetph_SL->SetLineColor(kRed);
    pdf_phjet_SL->SetLineColor(kGreen);
    pdf_jetjet_SL->SetLineColor(kGreen);
    pdf_phjet_SL->SetLineStyle(kDashed);
    pdf_jetph_SL->SetLineStyle(kDashed);
    frame_TiTi_SL->addPlotable(data_SL       ,"P");
    frame_TiTi_SL->addPlotable(pdf_tot_SL    );
    frame_TiTi_SL->addPlotable(pdf_phph_SL   );
    frame_TiTi_SL->addPlotable(pdf_phjet_SL  );
    frame_TiTi_SL->addPlotable(pdf_jetph_SL  );
    frame_TiTi_SL->addPlotable(pdf_jetjet_SL );


    data_Ti_L    ->SetMarkerSize(1);
    pdf_PHL      ->SetLineColor(kRed);
    pdf_JL       ->SetLineColor(kGreen);
    pdf_TiL      ->SetLineColor(kBlue);
    frame_PH_L   ->addPlotable(data_Ti_L,"P");
    frame_PH_L   ->addPlotable(pdf_PHL  );
    frame_PH_L   ->addPlotable(pdf_JL   );
    frame_PH_L   ->addPlotable(pdf_TiL   );

    data_Ti_SL    ->SetMarkerSize(1);
    pdf_PHSL      ->SetLineColor(kRed);
    pdf_JSL       ->SetLineColor(kGreen);
    pdf_TiSL      ->SetLineColor(kBlue);
    frame_PH_SL   ->addPlotable(data_Ti_SL,"P");
    frame_PH_SL   ->addPlotable(pdf_PHSL  );
    frame_PH_SL   ->addPlotable(pdf_JSL   );
    frame_PH_SL   ->addPlotable(pdf_TiSL   );


    data_JetL->SetMarkerSize(1);
    data_JetL->SetMarkerColor(kGreen);
    pdf_JetL->SetLineColor(kGreen);
    frame_Jet_L->addPlotable(data_JetL,"P");
    frame_Jet_L->addPlotable(pdf_JetL);

    data_JetSL->SetMarkerSize(1);
    data_JetSL->SetMarkerColor(kGreen);
    pdf_JetSL->SetLineColor(kGreen);
    frame_Jet_SL->addPlotable(data_JetSL,"P");
    frame_Jet_SL->addPlotable(pdf_JetSL);


    return 1;
  }
}

//////////////////////////////////////////////
bool BkgCategoryMerging::BuildRooPlots_RedFit()
/////////////////////////////////////////////
{
  TFile *file0 = TFile::Open(m_filenames[0]);
  std::cout << "Open: " << m_filenames[0] << std::endl;


  if( m_filenames.size() == 1){
    frame_lead_fake    = (RooPlot*)file0->Get("lead_fake");
    frame_sublead_fake = (RooPlot*)file0->Get("sublead_fake");
    frame_double_fake  = (RooPlot*)file0->Get("double_fake");
    std::cout << "Only 1 file is considered for the merging !" << std::endl;
    // file0->Close();
    return 0;
  }else{

    //==============================================================
    RooPlot* frame0_lead_fake  = (RooPlot*)file0->Get("lead_fake");
    RooHist*  LoTi          = frame0_lead_fake->getHist("LoTi");
    RooCurve* lead_fake_pdf = frame0_lead_fake->getCurve("lead_fake_pdf");
    //    RooCurve* lead_fake_pdf = frame0_lead_fake->getCurve("LoTi");
    frame_lead_fake  = frame0_lead_fake->emptyClone("lead_fake");
    //=================================================================
    RooPlot* frame0_sublead_fake  = (RooPlot*)file0->Get("sublead_fake");
    RooHist*  TiLo          = frame0_sublead_fake->getHist("TiLo");
    RooCurve* sublead_fake_pdf = frame0_sublead_fake->getCurve("sublead_fake_pdf");
    frame_sublead_fake  = frame0_sublead_fake->emptyClone("sublead_fake");
    //=================================================================
    RooPlot* frame0_double_fake  = (RooPlot*)file0->Get("double_fake");
    RooHist*  LoLo          = frame0_double_fake->getHist("LoLo");
    RooCurve* double_fake_pdf = frame0_double_fake->getCurve("double_fake_pdf");
    frame_double_fake  = frame0_double_fake->emptyClone("double_fake");
    //=================================================================
    // file0->Close();

//
//    //==============================================================
//    RooPlot* frame0_lead_fake  = (RooPlot*)file0->Get("lead_fake");
//    RooHist*  LoTi          = frame0_lead_fake->getHist("LoTi");
//    RooCurve* lead_fake_pdf = frame0_lead_fake->getCurve("lead_fake_pdf");
//    frame_lead_fake  = frame0_lead_fake->emptyClone("lead_fake");
//    //=================================================================
//    RooPlot* frame0_sublead_fake  = (RooPlot*)file0->Get("sublead_fake");
//    RooHist*  TiLo          = frame0_sublead_fake->getHist("TiLo");
//    RooCurve* sublead_fake_pdf = frame0_sublead_fake->getCurve("sublead_fake_pdf");
//    frame_sublead_fake  = frame0_sublead_fake->emptyClone("sublead_fake");
//    //=================================================================
//    RooPlot* frame0_double_fake  = (RooPlot*)file0->Get("double_fake");
//    RooHist*  LoLo          = frame0_double_fake->getHist("LoLo");
//    RooCurve* double_fake_pdf = frame0_double_fake->getCurve("double_fake_pdf");
//    frame_double_fake  = frame0_double_fake->emptyClone("double_fake");
//    //=================================================================
//    // file0->Close();


    for(int i=1;i<(int)m_filenames.size();i++){
      std::cout << "Open: " << m_filenames[i] << std::endl;
      TFile *file_i = TFile::Open(m_filenames[i]);

      //==============================================================
      RooPlot* frame_i_lead_fake  = (RooPlot*)file_i->Get("lead_fake");
      RooHist*  LoTi_i            = frame_i_lead_fake->getHist("LoTi");
      RooCurve* lead_fake_pdf_i   = frame_i_lead_fake->getCurve("lead_fake_pdf");
      LoTi = new RooHist(*LoTi,*LoTi_i);
      lead_fake_pdf = new RooCurve(lead_fake_pdf->GetName(),lead_fake_pdf->GetTitle(),
				   *lead_fake_pdf,*lead_fake_pdf_i);


      //=================================================================
      RooPlot* frame_i_sublead_fake  = (RooPlot*)file_i->Get("sublead_fake");
      RooHist*  TiLo_i               = frame_i_sublead_fake->getHist("TiLo");
      RooCurve* sublead_fake_pdf_i   = frame_i_sublead_fake->getCurve("sublead_fake_pdf");
      TiLo = new RooHist(*TiLo,*TiLo_i);
      sublead_fake_pdf = new RooCurve(sublead_fake_pdf->GetName(),sublead_fake_pdf->GetTitle(),
				      *sublead_fake_pdf,*sublead_fake_pdf_i);
      //=================================================================
      RooPlot* frame_i_double_fake  = (RooPlot*)file_i->Get("double_fake");
      RooHist*  LoLo_i              = frame_i_double_fake->getHist("LoLo");
      RooCurve* double_fake_pdf_i   = frame_i_double_fake->getCurve("double_fake_pdf");
      LoLo = new RooHist(*LoLo,*LoLo_i);
      double_fake_pdf = new RooCurve(double_fake_pdf->GetName(),double_fake_pdf->GetTitle(),
				     *double_fake_pdf,*double_fake_pdf_i);
      //=================================================================
      // file_i->Close();
      
    }






//
//
//
//    //    y_arr = h_data.GetY()
//    for(int i = 0; i<450; i++) {
//      if( (LoTi->GetY()[i]) < 1e-9 ) 	LoTi->SetPointError(i,0.,0.,0.,0.);
//      if( (TiLo->GetY()[i]) < 1e-9 )    TiLo->SetPointError(i,0.,0.,0.,0.);
//      if( (LoLo->GetY()[i]) < 1e-9 )    LoLo->SetPointError(i,0.,0.,0.,0.);
//    }
//

    
    LoTi->SetMarkerSize(1);
    lead_fake_pdf->SetLineColor(kGreen);
    frame_lead_fake->addPlotable(LoTi,"P");
    frame_lead_fake->addPlotable(lead_fake_pdf);
    //    frame_lead_fake->SetMinimum(1e-9);

    TiLo->SetMarkerSize(1);
    sublead_fake_pdf->SetLineColor(kRed);
    frame_sublead_fake->addPlotable(TiLo,"P");
    frame_sublead_fake->addPlotable(sublead_fake_pdf);
    //    frame_sublead_fake->SetMinimum(1e-9);

    LoLo->SetMarkerSize(1);
    double_fake_pdf->SetLineColor(kBlue);
    frame_double_fake->addPlotable(LoLo,"P");
    frame_double_fake->addPlotable(double_fake_pdf);
    //    frame_double_fake->SetMinimum(1e-9);

    return 1;
  }
  

}


////////////////////////////////////////////////////
void BkgCategoryMerging::StoreToRootFile(TString st)
////////////////////////////////////////////////////
{
  TFile fout(st,"RECREATE");


  //----------------------------
  std::cout << "Add histos to output" << std::endl;
  fout.Add(m_hh_data         );
  fout.Add(m_hh_bkg          );
  fout.Add(m_hh_red          );
  fout.Add(m_hh_gj           );
  fout.Add(m_hh_jg           );
  fout.Add(m_hh_jj           );
  fout.Add(m_hh_gg           );
  // fout.Add(m_hh_tot_uncert   );
  // fout.Add(m_hh_pur_uncert   );
  // fout.Add(m_hh_red_uncert   );
  // fout.Add(m_hh_irr_uncert   );
  // fout.Add(m_hh_pdf_uncert   );
  // fout.Add(m_hh_iso_uncert   );
  // fout.Add(m_hh_scale_uncert );
  // fout.Add(m_hh_bkg_syst     );
  //----------------------------
  std::cout << "Add graphs to output" << std::endl;
  m_gmgg_data->SetName("gmgg_data");
  fout.Add(m_gmgg_data);
  m_gmgg_bkg->SetName("gmgg_bkg");
  fout.Add(m_gmgg_bkg);
  m_gmgg_red->SetName("gmgg_red");
  fout.Add(m_gmgg_red);
  fout.Add(m_gmgg_gj);
  fout.Add(m_gmgg_jg);
  fout.Add(m_gmgg_jj);
  fout.Add(m_gmgg_gg);
  //----------------------------
  fout.Write();
  fout.Close();
}




///////////////////////////////////////////////////////////
void BkgCategoryMerging::StoreToRootFile_Search(TString st)
/////////////////////////////////////////////////////////
{
  TFile fout(st,"RECREATE");

  fout.Add(m_hh_data      );
  fout.Add(m_hh_bkg       );
  fout.Add(m_hh_red       );
  fout.Add(m_hh_gg        );
  // fout.Add(m_hh_pur_uncert);
  // fout.Add(m_hh_irr_uncert);
  // fout.Add(m_hh_red_uncert);
  // fout.Add(m_hh_tot_uncert);

  //----------------------------
  std::cout << "Add histos to output" << std::endl;
  //----------------------------
  fout.Write();
  fout.Close();
}

///////////////////////////////////////////////////////////
void BkgCategoryMerging::StoreToRootFile_2DFit(TString st)
/////////////////////////////////////////////////////////
{
  TFile fout(st,"RECREATE");

  //----------------------------
  std::cout << "Add Rooplots to output" << std::endl;
  //----------------------------
  fout.Add(frame_TiTi_L);
  fout.Add(frame_TiTi_SL);
  fout.Add(frame_PH_L);
  fout.Add(frame_PH_SL);
  fout.Add(frame_Jet_L);
  fout.Add(frame_Jet_SL);
  fout.Write();
  fout.Close();



}
///////////////////////////////////////////////////////////
void BkgCategoryMerging::StoreToRootFile_RedFit(TString st)
/////////////////////////////////////////////////////////
{
  TFile fout(st,"RECREATE");

  //----------------------------
  std::cout << "Add Rooplots to output" << std::endl;
  //----------------------------
  fout.Add(frame_lead_fake);
  fout.Add(frame_sublead_fake);
  fout.Add(frame_double_fake);
  fout.Write();
  fout.Close();

 
}

