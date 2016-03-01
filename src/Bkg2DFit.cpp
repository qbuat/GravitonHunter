#include "Bkg2DFit.h"
#include "ToolsCommons.h"
#include "ToolsUtilities.h"
#include "ToolsSignificanceHist.h"
#include <TError.h>
#include <TAxis.h>
#include <RooMsgService.h>
#include <TIterator.h>
//////////////////////////////////////////////////////////////////////////////
Bkg2DFit::Bkg2DFit(TTree* tree,
		   std::pair<double,double > mggbin,
		     TString FitParFile) : BkgFitFramework(tree,mggbin,FitParFile)
//////////////////////////////////////////////////////////////////////////////
{
  // Constructor: Set the m_tree variable (See BkgFitFramework)
  // and the range of the isolation fit + call the Init method
  SetIsoNorm(10,25);
  Init(mggbin,FitParFile);
}
//////////////////////////////////////////
Bkg2DFit::Bkg2DFit() : BkgFitFramework()
//////////////////////////////////////////
{
  m_mggbin         = std::make_pair(0.,10000.);
  m_do1Dfits        = false;
  m_do2Dfits        = false;
  m_doSplot         = false;

}
//////////////////////////////////////////////////
void Bkg2DFit::SetIsoNorm(double Xmin,double Xmax)
/////////////////////////////////////////////////
{
  m_isonormrange.first  = Xmin;
  m_isonormrange.second = Xmax;
}

///////////////////////////////////////////////////////////////////
void Bkg2DFit::Init(std::pair<double,double> mggbin, TString FitParFile)
///////////////////////////////////////////////////////////////////
{
  // Call the BkgFitFramework to create the DataSet (m_DataSet_mggcut)
  // with all the criteria applied
  BkgFitFramework::Init(mggbin,FitParFile);
  m_do1Dfits   = false;
  m_do2Dfits   = false;
  m_doSplot    = false;

  // On top of the criteria applied in BkgFitFramework
  // apply the isolation criteria and built the different control regions
  std::cout << "Parameters from : " << m_fitparfile << std::endl;
  //---------------------------------------------------------------------------------------------
  m_DataSet_TiTi        = (RooDataSet*)m_DataSet_mggcut->reduce("IsTight_L==1 && IsTight_SL==1");
  m_DataSet_LoLo        = (RooDataSet*)m_DataSet_mggcut->reduce("IsTight_L==0 && IsTight_SL==0");
  m_DataSet_LoTi        = (RooDataSet*)m_DataSet_mggcut->reduce("IsTight_L==0 && IsTight_SL==1");
  m_DataSet_TiLo        = (RooDataSet*)m_DataSet_mggcut->reduce("IsTight_L==1 && IsTight_SL==0");
  m_DataSet_Ti_L        = (RooDataSet*)m_DataSet_mggcut->reduce("IsTight_L==1");
  m_DataSet_Ti_SL       = (RooDataSet*)m_DataSet_mggcut->reduce("IsTight_SL==1");
  std::cout << m_isocut_st << std::endl;
  m_DataSet_TiIsoTiIso  = (RooDataSet*)m_DataSet_TiTi->reduce(m_isocut_st);
  //---------------------------------------------------------------------------------------------
}
///////////////////////////////////////////
void Bkg2DFit::LeadingJetFit(bool doFit)
///////////////////////////////////////////
{
  // Perform the fit of the isolation on the AntiTight-Tight
  // sample. This corresponds to the extraction of the jet
  // template for the leading candidate

  int NbinsFit = 2*abs((int)m_Isomax-(int)m_Isomin);

  RooRealVar* peak_L  = new RooRealVar("L_Jet_peak","#mu",0,"[GeV]");
  RooRealVar* width_L = new RooRealVar("L_Jet_width","#sigma",0,"[GeV]");
  RooRealVar* tail_L  = new RooRealVar("L_Jet_tail","tail",0);

  m_JetL  = new RooArgSet(*peak_L,*width_L,*tail_L,"JetL_Set");
  m_JetL->readFromFile(m_fitparfile,"READ","Leading Jet" );
  m_LeadJetPdf = new RooNovosibirsk( "LeadJetPdf","LeadJetPdf",
				     *m_Iso_L,*peak_L,*width_L,*tail_L );
  if(doFit){
    m_FitRes_Jet_L = m_LeadJetPdf->fitTo( *m_DataSet_LoTi,
					  RooFit::Range(m_Isomin,m_Isomax),
					  RooFit::Save(kTRUE) );
    m_FitRes_Jet_L->Print("v");
  }
  frame_JetL = m_Iso_L->frame(m_Isomin,m_Isomax,NbinsFit);
  frame_JetL->SetName("JetL");
  m_DataSet_LoTi->plotOn( frame_JetL,
			  RooFit::LineColor(3),
			  RooFit::MarkerColor(3),
			  RooFit::Name("data_JL"));
  m_LeadJetPdf->plotOn( frame_JetL,
			RooFit::LineColor(3),
			RooFit::Name("pdf_JL"));
}
//////////////////////////////////////////////
void Bkg2DFit::SubLeadingJetFit(bool doFit)
//////////////////////////////////////////////
{
  // Perform the fit of the isolation on the Tight-AntiTight
  // sample. This corresponds to the extraction of the jet
  // template for the subleading candidate
  int NbinsFit = 2*abs((int)m_Isomax-(int)m_Isomin);
  RooRealVar* peak_SL  = new RooRealVar("SL_Jet_peak","#mu",0,"[GeV]");
  RooRealVar* width_SL = new RooRealVar("SL_Jet_width","#sigma",0,"[GeV]");
  RooRealVar* tail_SL  = new RooRealVar("SL_Jet_tail","tail",0);
  m_JetSL  = new RooArgSet(*peak_SL,*width_SL,*tail_SL,"JetSL_Set"); 
  m_JetSL->readFromFile(m_fitparfile,"READ","SubLeading Jet" );
  m_SubLeadJetPdf = new RooNovosibirsk( "SubLeadJetPdf","SubLeadJetPdf",
					*m_Iso_SL,*peak_SL,*width_SL,*tail_SL );
  if(doFit){
    m_FitRes_Jet_SL = m_SubLeadJetPdf->fitTo( *m_DataSet_TiLo,
					      RooFit::Range(m_Isomin,m_Isomax),
					      RooFit::Save(kTRUE) );
    m_FitRes_Jet_SL->Print("v");
  }
  frame_JetSL= m_Iso_SL->frame(m_Isomin,m_Isomax,NbinsFit);
  frame_JetSL->SetName("JetSL");
  m_DataSet_TiLo->plotOn( frame_JetSL,
			  RooFit::LineColor(3),
			  RooFit::MarkerColor(3),
			  RooFit::Name("data_JSL") );
  m_SubLeadJetPdf->plotOn( frame_JetSL,
			   RooFit::LineColor(3),
			   RooFit::Name("pdf_JSL") );

}
///////////////////////////////////////
void Bkg2DFit::JetJetFit(bool doFit)
//////////////////////////////////////
{
  // Compute the 2D template of the jetjet component
  // The template is taken from the AntiTight-AntiTight
  // sample and smoothed using the RooNDKeysPdf class
  int Nbins = 2*abs((int)m_Isomax-(int)m_Isomin);
  m_hJetJet    = (TH2F*)m_DataSet_LoLo->createHistogram( *m_Iso_L,*m_Iso_SL,
							 Nbins,Nbins,"","hJetJet" );
  m_hJetJet->GetXaxis()->SetTitle(m_Iso_L->GetTitle());
  m_hJetJet->GetYaxis()->SetTitle(m_Iso_SL->GetTitle());
  m_hJetJet->SetName("JetJetHist");
  m_hJetJet->SetTitle("JetJetHist");
  if(doFit){
    m_JetJetPdf  = new RooNDKeysPdf( "JetJetPdf","JetJetPdf",RooArgList(*m_Iso_L,*m_Iso_SL),
				     *m_DataSet_LoLo,RooKeysPdf::NoMirror,2 );
    TString histname = m_Iso_L->GetName();
    histname        += ":";
    histname        += m_Iso_SL->GetName();

    m_hJetJetFit = (TH2F*)m_JetJetPdf->createHistogram(histname,Nbins,Nbins);
    m_hJetJetFit->SetName("JetJetPdf");
    m_hJetJetFit->SetTitle("JetJetPdf");
  }
}
///////////////////////////////////////////
void Bkg2DFit::LeadPhotonFit(bool doFit)
///////////////////////////////////////////
{
  // Perform a fit of the leading photon candidate
  // in the Tight-Tight sample. The leading jet template 
  // is normalized to data in the high isolation values.
  // Then a fit of a Crystal-Ball is performed to create
  // the leading photon template

  int NbinsFit = 2*abs((int)m_Isomax-(int)m_Isomin);

  TString cut_L = m_Iso_L->GetName();
  cut_L        += Form(">%f",m_isonormrange.first);
  cut_L        += "&&";
  cut_L        += m_Iso_L->GetName();
  cut_L        += Form("<%f",m_isonormrange.second); 
  std::cout << cut_L << std::endl;
  double num_L     = m_DataSet_TiTi->sumEntries(cut_L);
  double denom_L   = m_DataSet_LoTi->sumEntries(cut_L);
  double scale_L   = num_L/denom_L;
  double Ntot_LoTi = m_DataSet_LoTi->sumEntries();
  RooRealVar N_Jet_L("N_Jet_L","N_Jet_L",Ntot_LoTi*scale_L);
  RooRealVar N_PH_L("L_PH_N","L_PH_N",1000.);

  RooRealVar* meanCB_L  = new RooRealVar("L_PH_meanCB","meanCB",4);
  RooRealVar* sigmaCB_L = new RooRealVar("L_PH_sigmaCB","sigmaCB",4);
  RooRealVar* alphaCB_L = new RooRealVar("L_PH_alphaCB","alphaCB",4);
  RooRealVar* nCB_L     = new RooRealVar("L_PH_nCB","nCB",4);
  RooRealVar* fracCB_L  = new RooRealVar("L_PH_fracCB","fracCB",0.5);
  RooRealVar* sigmaG_L  = new RooRealVar("L_PH_sigmaG","sigmaG",4);

  m_PHL     = new RooArgSet(*meanCB_L,*sigmaCB_L,*alphaCB_L,*nCB_L,*fracCB_L,
			    *sigmaG_L,"PHL_Set");
  RooArgSet TiL(*m_PHL);
  TiL.add(N_PH_L);
  TiL.readFromFile(m_fitparfile,"READ","Leading Tight" );
  TiL.Print();
  std::cout << N_PH_L  .getVal() << std::endl;
  std::cout << Ntot_LoTi << std::endl;
  std::cout << scale_L << std::endl;
  std::cout << N_Jet_L  .getVal() << std::endl;
  std::cout << meanCB_L  ->getVal() << std::endl;
  std::cout << sigmaCB_L ->getVal() << std::endl;
  std::cout << alphaCB_L ->getVal() << std::endl;
  std::cout << nCB_L     ->getVal() << std::endl;
  std::cout << fracCB_L  ->getVal() << std::endl;
  std::cout << sigmaG_L  ->getVal() << std::endl;


  // m_LeadPHPdf= new RooCBShape( "LeadPHPdf","LeadPHPdf",
  // 			       *m_Iso_L,*meanCB_L,*sigmaCB_L,*alphaCB_L,*nCB_L );
  RooCBShape* LeadPH_CB = new RooCBShape( "LeadPH_CB","LeadPH_CB",
				*m_Iso_L,*meanCB_L,*sigmaCB_L,*alphaCB_L,*nCB_L );
  RooGaussian* LeadPH_G = new RooGaussian( "LeadPH_G","LeadPH_G",*m_Iso_L,*meanCB_L,*sigmaG_L);
  m_LeadPHPdf = new RooAddPdf( "LeadPHPdf","LeadPHPdf",
			       *LeadPH_CB,*LeadPH_G, *fracCB_L);

  RooAddPdf* LeadTightPdf = new RooAddPdf( "LeadTightPdf","LeadTightPdf",
					   RooArgList(*m_LeadPHPdf,*m_LeadJetPdf),
					   RooArgList(N_PH_L,N_Jet_L) );
  if(doFit){
    m_FitRes_Ti_L = LeadTightPdf->fitTo( *m_DataSet_TiTi,
					 RooFit::Range(m_Isomin,m_Isomax),
					 RooFit::Extended(kTRUE),
					 RooFit::Save(kTRUE) );
    m_FitRes_Ti_L->Print("v");
  }
  
  frame_PHL= m_Iso_L->frame(m_Isomin,m_Isomax,NbinsFit);
  frame_PHL->SetName("PHL");
  m_DataSet_TiTi->plotOn( frame_PHL,
			  RooFit::LineColor(1),
			  RooFit::MarkerColor(1),
			  RooFit::Name("data_TiL") );
  LeadTightPdf->plotOn( frame_PHL,
			RooFit::LineColor(2),
			RooFit::Components(*m_LeadPHPdf),
			RooFit::Name("pdf_PHL") );
  LeadTightPdf->plotOn( frame_PHL,
			RooFit::LineColor(3),
			RooFit::Components(*m_LeadJetPdf),
			RooFit::Name("pdf_JL") );
  
  LeadTightPdf->plotOn( frame_PHL,
			RooFit::LineColor(4),
			RooFit::Name("pdf_TiL") );
}
/////////////////////////////////////////////
void Bkg2DFit::SubLeadPhotonFit(bool doFit)
/////////////////////////////////////////////
{
  // Perform a fit of the subleading photon candidate
  // in the Tight-Tight sample. The leading jet template 
  // is normalized to data in the high isolation values.
  // Then a fit of a Crystal-Ball is performed to create
  // the subleading photon template

  int NbinsFit = 2*abs((int)m_Isomax-(int)m_Isomin);

  TString cut_SL = m_Iso_SL->GetName();
  cut_SL        += Form(">%f",m_isonormrange.first);
  cut_SL        += "&&";
  cut_SL        += m_Iso_SL->GetName();
  cut_SL        += Form("<%f",m_isonormrange.second); 
  std::cout << cut_SL << std::endl;


  double num_SL    = m_DataSet_TiTi->sumEntries(cut_SL);
  double denom_SL  = m_DataSet_TiLo->sumEntries(cut_SL);
  double scale_SL  = num_SL/denom_SL;
  double Ntot_TiLo = m_DataSet_TiLo->sumEntries();
  RooRealVar N_Jet_SL("N_Jet_SL","N_Jet_SL",Ntot_TiLo*scale_SL);
  RooRealVar N_PH_SL("SL_PH_N","SL_PH_N",1000.);

  RooRealVar* meanCB_SL  = new RooRealVar("SL_PH_meanCB","meanCB",4);
  RooRealVar* sigmaCB_SL = new RooRealVar("SL_PH_sigmaCB","sigmaCB",4);
  RooRealVar* alphaCB_SL = new RooRealVar("SL_PH_alphaCB","alphaCB",4);
  RooRealVar* nCB_SL     = new RooRealVar("SL_PH_nCB","nCB",4);
  RooRealVar* fracCB_SL  = new RooRealVar("SL_PH_fracCB","fracCB",0.5);
  RooRealVar* sigmaG_SL  = new RooRealVar("SL_PH_sigmaG","sigmaG",4);

  m_PHSL     = new RooArgSet(*meanCB_SL,*sigmaCB_SL,*alphaCB_SL,*nCB_SL,*fracCB_SL,
			     *sigmaG_SL,"PHSL_Set");
  RooArgSet TiSL(*m_PHSL);
  TiSL.add(N_PH_SL);
  TiSL.readFromFile(m_fitparfile,"READ","SubLeading Tight" );
  // m_SubLeadPHPdf = new RooCBShape( "SubLeadPHPdf","SubLeadPHPdf",
  // 				   *m_Iso_SL,*meanCB_SL,*sigmaCB_SL,*alphaCB_SL,*nCB_SL);
  RooCBShape* SubLeadPH_CB = new RooCBShape( "SubLeadPH_CB","SubLeadPH_CB",
					     *m_Iso_SL,*meanCB_SL,*sigmaCB_SL,*alphaCB_SL,*nCB_SL);
  RooGaussian* SubLeadPH_G = new RooGaussian( "SubLeadPH_G","SubLeadPH_G",
					      *m_Iso_SL,*meanCB_SL,*sigmaG_SL);
  m_SubLeadPHPdf = new RooAddPdf( "SubLeadPHPdf", "SubLeadPHPdf" ,
				  *SubLeadPH_CB,*SubLeadPH_G,
				  *fracCB_SL );



  //---------------------------------------------------------------------------------
  RooAddPdf* SubLeadTightPdf = new RooAddPdf( "SubLeadTightPdf","SubLeadTightPdf",
					      RooArgList(*m_SubLeadPHPdf,*m_SubLeadJetPdf),
					      RooArgList(N_PH_SL,N_Jet_SL));

  if(doFit){
    m_FitRes_Ti_SL = SubLeadTightPdf->fitTo( *m_DataSet_TiTi,
					     RooFit::Range(m_Isomin,m_Isomax),
					     RooFit::Extended(kTRUE),
					     RooFit::Save(kTRUE) );
    m_FitRes_Ti_SL->Print("v");
  }
  frame_PHSL = m_Iso_SL->frame(m_Isomin,m_Isomax,NbinsFit);
  frame_PHSL->SetName("PHSL");
  m_DataSet_TiTi->plotOn( frame_PHSL,
			  RooFit::LineColor(1),
			  RooFit::MarkerColor(1),
			  RooFit::Name("data_TiSL") );
  SubLeadTightPdf->plotOn( frame_PHSL,
			   RooFit::LineColor(2),
			   RooFit::Components(*m_SubLeadPHPdf),
			   RooFit::Name("pdf_PHSL") );
  SubLeadTightPdf->plotOn( frame_PHSL,
			   RooFit::LineColor(3),
			   RooFit::Components(*m_SubLeadJetPdf),
			   RooFit::Name("pdf_JSL") );
  SubLeadTightPdf->plotOn( frame_PHSL,
			   RooFit::LineColor(4),
			   RooFit::Name("pdf_TiSL") );
}
////////////////////////////////////////////////////////
void Bkg2DFit::TightTightFit(bool doFit,bool doFrame)
///////////////////////////////////////////////////////
{
  // Perform the 2D fit of the isolations variables
  // on the TightTight sample
  // For gg, gj and jg the 2D template is a 1DX1D template
  // obtained with RooProdPdf

  int NbinsFit = 2*abs((int)m_Isomax-(int)m_Isomin);

  //-------------------------------------------------------------------------------------------
  m_PHPHPdf  = new RooProdPdf("PHPHPdf" ,"PHPHPdf" ,RooArgList(*m_LeadPHPdf,*m_SubLeadPHPdf) );  
  m_PHJetPdf = new RooProdPdf("PHJetPdf","PHJetPdf",RooArgList(*m_LeadPHPdf,*m_SubLeadJetPdf) );  
  m_JetPHPdf = new RooProdPdf("JetPHPdf","JetPHPdf",RooArgList(*m_LeadJetPdf,*m_SubLeadPHPdf) );
  RooArgList PdfList(*m_PHPHPdf,*m_PHJetPdf,*m_JetPHPdf,*m_JetJetPdf);
  //----------------------------------------------------------------
  Ngamgam = new RooRealVar("Ngamgam","Ngamgam",10000,0,1000000000.);
  Ngamjet = new RooRealVar("Ngamjet","Ngamjet",10000,0,1000000000.);
  Njetgam = new RooRealVar("Njetgam","Njetgam",10000,0,1000000000.);
  Njetjet = new RooRealVar("Njetjet","Njetjet",10000,0,1000000000.);
  RooArgList CoeList(*Ngamgam,*Ngamjet,*Njetgam,*Njetjet);
  RooArgSet CoeSet(*Ngamgam,*Ngamjet,*Njetgam, *Njetjet);
  //----------------------------------------------------------
  m_TiTiPdf = new RooAddPdf("TiTiPdf","TiTiPdf",PdfList,CoeList );
  RooArgSet TiTiSet(*m_PHL,*m_PHSL,"TiTi_Set");
  TiTiSet.add(*m_JetL); 
  TiTiSet.add(*m_JetSL);
  TiTiSet.add(CoeSet);
  TiTiSet.readFromFile(m_fitparfile,"READ","Tight Tight");
  std::cout << "Tight Tight Fit starts here" << std::endl;
  if(doFit){
    m_FitRes_TiTi= m_TiTiPdf->fitTo(*m_DataSet_TiTi,RooFit::Extended(kTRUE),RooFit::Save(kTRUE) );
    m_FitRes_TiTi->Print("v");
  }
  //------------------------------------------------------------------------------//
  if(doFrame){
    frame_TiTi_L = m_Iso_L->frame(m_Isomin,m_Isomax,NbinsFit);
    frame_TiTi_L->SetName("TiTiL");
    m_DataSet_TiTi->plotOn( frame_TiTi_L,RooFit::Name("data_L"),
			    RooFit::LineColor(1),RooFit::MarkerColor(1) );
    m_TiTiPdf->plotOn( frame_TiTi_L,RooFit::LineColor(1),RooFit::Name("pdf_tot") );
    m_TiTiPdf->plotOn( frame_TiTi_L,RooFit::LineColor(2),
		       RooFit::Components(*m_PHPHPdf),RooFit::Name("pdf_phph") );
    m_TiTiPdf->plotOn( frame_TiTi_L,RooFit::LineColor(2),RooFit::LineStyle(kDashed),
		       RooFit::Components(*m_PHJetPdf),RooFit::Name("pdf_phjet") );
    m_TiTiPdf->plotOn( frame_TiTi_L,RooFit::LineColor(3),
		       RooFit::Components(*m_JetJetPdf),RooFit::Name("pdf_jetjet") );
    m_TiTiPdf->plotOn( frame_TiTi_L,RooFit::LineColor(3),RooFit::LineStyle(kDashed),
		       RooFit::Components(*m_JetPHPdf),RooFit::Name("pdf_jetph") );
    m_TiTiPdf->plotOn( frame_TiTi_L,RooFit::LineColor(3),RooFit::Invisible(),
		       RooFit::Components(RooArgSet(*m_PHJetPdf,*m_PHPHPdf)),
		       RooFit::Name("pdf_phph_phjet") );
    m_TiTiPdf->plotOn( frame_TiTi_L,RooFit::LineColor(3),RooFit::Invisible(),
		       RooFit::Components(RooArgSet(*m_JetPHPdf,*m_JetJetPdf)),
		       RooFit::Name("pdf_jetph_jetjet") );
    //------------------------------------------------------------------------------//
    frame_TiTi_SL= m_Iso_SL->frame(m_Isomin,m_Isomax,NbinsFit);
    frame_TiTi_SL->SetName("TiTiSL");
    m_DataSet_TiTi->plotOn( frame_TiTi_SL,RooFit::Name("data_SL"),
			    RooFit::LineColor(1),RooFit::MarkerColor(1) );
    m_TiTiPdf->plotOn( frame_TiTi_SL,RooFit::LineColor(1),RooFit::Name("pdf_tot") );
    m_TiTiPdf->plotOn( frame_TiTi_SL,RooFit::LineColor(2),
		       RooFit::Components(*m_PHPHPdf),RooFit::Name("pdf_phph") );
    m_TiTiPdf->plotOn( frame_TiTi_SL,RooFit::LineColor(2),RooFit::LineStyle(kDashed),
		       RooFit::Components(*m_JetPHPdf),RooFit::Name("pdf_jetph") );
    m_TiTiPdf->plotOn( frame_TiTi_SL,RooFit::LineColor(3),
		       RooFit::Components(*m_JetJetPdf),RooFit::Name("pdf_jetjet") );
    m_TiTiPdf->plotOn( frame_TiTi_SL,RooFit::LineColor(3),RooFit::LineStyle(kDashed),
		       RooFit::Components(*m_PHJetPdf),RooFit::Name("pdf_phjet") );
    m_TiTiPdf->plotOn( frame_TiTi_SL,RooFit::LineColor(3),RooFit::Invisible(),
		       RooFit::Components(RooArgSet(*m_JetPHPdf,*m_PHPHPdf)),
		       RooFit::Name("pdf_phph_jetph") );
    m_TiTiPdf->plotOn( frame_TiTi_SL,RooFit::LineColor(3),RooFit::Invisible(),
		       RooFit::Components(RooArgSet(*m_PHJetPdf,*m_JetJetPdf)),
		       RooFit::Name("pdf_phjet_jetjet") );
  }
}
///////////////////////////////////////////////////////////////////
void Bkg2DFit::Fitter(bool do1Dfits,bool do2Dfits,bool doSplot)
///////////////////////////////////////////////////////////////////
{
  //--> Main method to run the whole procedure
  // Run all the fit successively

  m_do1Dfits        = do1Dfits;
  m_do2Dfits        = do2Dfits;
  m_doSplot         = doSplot;

  //----------------------------
  LeadingJetFit(m_do1Dfits);
  SetConstantParameters(*m_JetL);
  //-----------------------------
  SubLeadingJetFit(m_do1Dfits);
  SetConstantParameters(*m_JetSL);
  //-----------------------------
  LeadPhotonFit(m_do1Dfits);
  SetConstantParameters(*m_PHL);
  //-----------------------------
  SubLeadPhotonFit(m_do1Dfits);
  SetConstantParameters(*m_PHSL);
  //----------------------------
  JetJetFit(m_do1Dfits);
  std::cout << "JetJet Fit performed" << std::endl;
  //-----------------------------------------------
  TightTightFit(m_do2Dfits);
  //------------------------
  if(m_do2Dfits)
    YieldsExtraction();
  //--------------------
  if(m_doSplot)
    ExtractShapesFrom_sPlot();
}
/////////////////////////////////
void Bkg2DFit::YieldsExtraction()
/////////////////////////////////
{
  //--> Compute the 4 differents yields from
  // the 2D pdf TiTiPdf(iso_l,iso_sl)
  // By applying the isolation criteria

  m_Iso_L->setRange("TiIso",m_Isomin,m_IsoCut_L);
  m_Iso_SL->setRange("TiIso",m_Isomin,m_IsoCut_SL);
  RooArgSet Var(*m_Iso_L,*m_Iso_SL);

  iPHPH_TiIso   = new RooRealVar("iPHPH_TiIso","iPHPH_TiIso",0);
  iPHJET_TiIso  = new RooRealVar("iPHJET_TiIso","iPHJET_TiIso",0);
  iJETPH_TiIso  = new RooRealVar("iJETPH_TiIso","iJETPH_TiIso",0);
  iJETJET_TiIso = new RooRealVar("iJETJET_TiIso","iJETJET_TiIso",0);
  iPHPH_TiIso   = (RooRealVar*)m_PHPHPdf->createIntegral( Var,RooFit::NormSet(Var),
							  RooFit::Range("TiIso") );
  iPHJET_TiIso  = (RooRealVar*)m_PHJetPdf->createIntegral( Var,RooFit::NormSet(Var),
							   RooFit::Range("TiIso") );
  iJETPH_TiIso  = (RooRealVar*)m_JetPHPdf->createIntegral( Var,RooFit::NormSet(Var),
							   RooFit::Range("TiIso") );
  iJETJET_TiIso = (RooRealVar*)m_JetJetPdf->createIntegral( Var,RooFit::NormSet(Var),
							    RooFit::Range("TiIso") );

  double iPHPH_TiIso_val   = iPHPH_TiIso->getVal();
  double iPHJET_TiIso_val  = iPHJET_TiIso->getVal();
  double iJETPH_TiIso_val  = iJETPH_TiIso->getVal();
  double iJETJET_TiIso_val = iJETJET_TiIso->getVal();

  m_NgamgamYield      = iPHPH_TiIso_val*Ngamgam->getVal();
  m_NgamjetYield      = iPHJET_TiIso_val*Ngamjet->getVal();
  m_NjetgamYield      = iJETPH_TiIso_val*Njetgam->getVal();
  m_NjetjetYield      = iJETJET_TiIso_val*Njetjet->getVal();
  m_NgamgamYieldError = iPHPH_TiIso_val*Ngamgam->getError();
  m_NgamjetYieldError = iPHJET_TiIso_val*Ngamjet->getError();
  m_NjetgamYieldError = iJETPH_TiIso_val*Njetgam->getError();
  m_NjetjetYieldError = iJETJET_TiIso_val*Njetjet->getError();

  std::cout << "**************************************************" << std::endl;
  std::cout << "*************** TIGHT REGION *********************" << std::endl;
  std::cout << "Ngamgam = " << Ngamgam->getVal() << " +/- "<< Ngamgam->getError() << std::endl;
  std::cout << "Ngamjet = " << Ngamjet->getVal() << " +/- "<< Ngamjet->getError() << std::endl;
  std::cout << "Njetgam = " << Njetgam->getVal() << " +/- "<< Njetgam->getError() << std::endl;
  std::cout << "Njetjet = " << Njetjet->getVal() << " +/- "<< Njetjet->getError() << std::endl;
  std::cout << "**************************************************" << std::endl;
  std::cout << "*************** TIGHT+ISOLATED REGION ************" << std::endl;
  std::cout << "Ngamgam = " <<iPHPH_TiIso_val*Ngamgam->getVal()  << " +/- "
	    << iPHPH_TiIso_val*Ngamgam->getError()<< std::endl;
  std::cout << "Ngamjet = " <<iPHJET_TiIso_val*Ngamjet->getVal() << " +/- "
	    << iPHJET_TiIso_val*Ngamjet->getError()<< std::endl;
  std::cout << "Njetgam = " <<iJETPH_TiIso_val*Njetgam->getVal() << " +/- "
	    << iJETPH_TiIso_val*Njetgam->getError()<< std::endl;
  std::cout << "Njetjet = " <<iJETJET_TiIso_val*Njetjet->getVal()<< " +/- "
	    << iJETJET_TiIso_val*Njetjet->getError()<< std::endl;
}
///////////////////////////////
void Bkg2DFit::PlotResults()
///////////////////////////////
{
  TString mggbin_name= Form("_mgg%d%d",(int)m_mggbin.first,(int)m_mggbin.second);
  if(m_do1Dfits){
    TCanvas *cJL=new TCanvas("cJL","Lead Jet",800,600);
    cJL->cd();
    frame_JetL->Draw();
    cJL->SaveAs("./plots/LeadJet"+mggbin_name+".eps");
    TCanvas* cJSL=new TCanvas("cJSL","SubLead Jet",800,600);
    cJSL->cd();
    frame_JetSL->Draw();
    cJSL->SaveAs("./plots/SubLeadJet"+mggbin_name+".eps");
    TCanvas* cPHL=new TCanvas("cPHL","Lead Tight",800,600);
    cPHL->cd();
    frame_PHL->Draw();
    cPHL->SaveAs("./plots/LeadPhoton"+mggbin_name+".eps");
    TCanvas* cPHSL=new TCanvas("cPHSL","SubLead Tight",800,600);
    cPHSL->cd();
    frame_PHSL->Draw();
    cPHSL->SaveAs("./plots/SubLeadPhoton"+mggbin_name+".eps");
    TCanvas* cJetJet = new TCanvas("cJetJet","cJetJet",800,600);
    cJetJet->cd();
    m_hJetJet->Draw("lego");
    cJetJet->SaveAs("./plots/JetJet"+mggbin_name+".eps");
    TCanvas* cJetJetFit = new TCanvas("cJetJetFit","cJetJetFit",800,600);
    cJetJetFit->cd();
    m_hJetJetFit->Draw("lego");
    cJetJetFit->SaveAs("./plots/JetJetFit"+mggbin_name+".eps");
  }
  //================================================================//
  if(m_do2Dfits){
    TCanvas* cTiTiL=new TCanvas("cTiTiL","Lead Tight final ",800,600);
    cTiTiL->cd();
    frame_TiTi_L->Draw();
    TLegend* leg_TiTi_L = new TLegend(0.5,0.5,0.8,0.7);
    TLatex * mytex = new TLatex();
    mytex->SetTextSize(0.04);
    leg_TiTi_L->SetFillColor(0);
    TString names_TiTi_L[] = { "pdf_tot","pdf_phph",
			       "pdf_phjet","pdf_jetph",
			       "pdf_jetjet", "" };
    TString signs_TiTi_L[] = { "#gamma#gamma+#gammaj+j#gamma+jj",
			       "#gamma#gamma","#gammaj",
			       "j#gamma","jj", };
    Int_t i_TiTi_L=-1;
    while ( names_TiTi_L[++i_TiTi_L] != "" ) {
      TObject *obj = frame_TiTi_L->findObject(names_TiTi_L[i_TiTi_L].Data());
      if (!obj) {
	Warning("fitBi4",Form("Can't find item='%s' in the frame2!\n",names_TiTi_L[i_TiTi_L].Data()));
	continue;
      }
      leg_TiTi_L->AddEntry(obj,signs_TiTi_L[i_TiTi_L],"l");
    }
    leg_TiTi_L->Draw("same");
    mytex->DrawLatex(5 ,1600,Form("#int L dt = %1.2f fb^{-1}",Commons::int_lumi_fb) );
    mytex->DrawLatex(17,1600,"#sqrt{s} = 8 TeV");
    cTiTiL->SaveAs("./plots/LeadFinal"+mggbin_name+".eps");
    //====================================================================//
    TCanvas* cTiTiSL=new TCanvas("cTiTiSL","SubLead Tight final ",800,600);
    cTiTiSL->cd();
    frame_TiTi_SL->Draw();
    TLegend* leg_TiTi_SL = new TLegend(0.5,0.5,0.8,0.7);
    leg_TiTi_SL->SetFillColor(0);
    TString names_TiTi_SL[] = { "pdf_tot","pdf_phph",
				"pdf_jetph","pdf_phjet",
				"pdf_jetjet","" };
    TString signs_TiTi_SL[] = { "#gamma#gamma+#gammaj+j#gamma+jj",
				"#gamma#gamma","j#gamma",
				"#gammaj","jj", };
    Int_t i_TiTi_SL=-1;
    while ( names_TiTi_SL[++i_TiTi_SL] != "" ) {
      TObject *obj = frame_TiTi_SL->findObject(names_TiTi_SL[i_TiTi_SL].Data());
      if (!obj) {
	Warning("fitBi4",Form("Can't find item='%s' in the frame!\n",names_TiTi_SL[i_TiTi_SL].Data()));
	continue;
      }
      leg_TiTi_SL->AddEntry(obj,signs_TiTi_SL[i_TiTi_SL],"l");
    }
    leg_TiTi_SL->Draw("same");
    mytex->DrawLatex(5 ,1300,Form("#int L dt = %1.2f fb^{-1}",Commons::int_lumi_fb) );
    mytex->DrawLatex(17,1300,"#sqrt{s} = 8 TeV");
    cTiTiSL->SaveAs("./plots/SubLeadFinal"+mggbin_name+".eps");
  }
  //=============================================================================//
  if(m_doSplot){
    //----------------------------------------------------------//
    hmgg_Ngg->SetLineColor(2);
    hmgg_Ngg->SetFillColor(2);
    hmgg_Ngg->SetMarkerColor(2);
    
    hmgg_Ngj->SetLineColor(3);
    hmgg_Ngj->SetFillColor(3);
    hmgg_Ngj->SetMarkerColor(3);
    
    hmgg_Njg->SetLineColor(4);
    hmgg_Njg->SetFillColor(4);
    hmgg_Njg->SetMarkerColor(4);

    hmgg_Njj->SetLineColor(5);
    hmgg_Njj->SetFillColor(5);
    hmgg_Njj->SetMarkerColor(5);
    //----------------------------------------------------------//
    THStack* hs=new THStack("hs","hs");
    hs->SetMinimum(1);
    hs->Add(hmgg_Njj);
    hs->Add(hmgg_Njg);
    hs->Add(hmgg_Ngj);
    hs->Add(hmgg_Ngg);
    TCanvas* c_Splot2=new TCanvas("c_Splot2","c_Splot2",800,600);
    c_Splot2->cd();
    c_Splot2->SetLogy();
    hs->Draw("");
    hs->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    hs->Draw("sameHIST");
    hmgg->Draw("samePE");
    c_Splot2->SaveAs("./plots/Splot2.eps" );
    //----------------------------------------------------------//
    TCanvas* c_Splot3=new TCanvas("c_Splot3","c_Splot3",800,600);
    c_Splot3->cd();
    THStack* hs1=new THStack("hs1","hs1");
    hs1->Add(hmgg_red);
    hs1->Add(hmgg_Ngg);
    hmgg->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    hmgg->Draw("PE");
    hs1->Draw("sameHIST");
    hmgg->Draw("samePE");
    c_Splot3->SaveAs("./plots/Splot3.eps" );
    //----------------------------------------------------------//
    TCanvas* c_Splot4=new TCanvas("c_Splot4","c_Splot4",800,600);
    c_Splot4->Divide(2,2);
    c_Splot4->cd(1);
    hmgg_Ngg->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    hmgg_Ngg->Draw("PE");
    c_Splot4->cd(2);
    hmgg_Ngj->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    hmgg_Ngj->Draw("PE");
    hmgg_Njg->Draw("samePE");
    hmgg_Njj->Draw("samePE");
    c_Splot4->cd(3);
    hmgg_Ngjjg->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    hmgg_Ngjjg->Draw("PE");
    hmgg_Njj->Draw("samePE");
    c_Splot4->cd(4);
    hmgg_Njj->Draw("PE");
    c_Splot4->SaveAs("./plots/Splot4.eps" );
    //----------------------------------------------------------//
  }
  //=============================================================================//
}
////////////////////////////////////////////
void Bkg2DFit::ExtractShapesFrom_sPlot()
///////////////////////////////////////////
{
  RooArgList yieldsList(*Ngamgam,*Ngamjet,*Njetgam,*Njetjet);
  RooStats::SPlot mySplot("mySPlot","mySPlot",*m_DataSet_TiTi,m_TiTiPdf,yieldsList);

  m_DataSet_TiTi->Print("v");
  m_DataSet_TiTi_Ngamgamw = new RooDataSet( m_DataSet_TiTi->GetName(),
					    m_DataSet_TiTi->GetTitle(),
					    m_DataSet_TiTi,
					    *m_DataSet_TiTi->get(),
					    0,"Ngamgam_sw" );
  m_DataSet_TiTi_Ngamjetw = new RooDataSet( m_DataSet_TiTi->GetName(),
					    m_DataSet_TiTi->GetTitle(),
					    m_DataSet_TiTi,
					    *m_DataSet_TiTi->get(),
					    0,"Ngamjet_sw" );
  m_DataSet_TiTi_Njetgamw = new RooDataSet( m_DataSet_TiTi->GetName(),
					    m_DataSet_TiTi->GetTitle(),
					    m_DataSet_TiTi,
					    *m_DataSet_TiTi->get(),
					    0,"Njetgam_sw" );
  m_DataSet_TiTi_Njetjetw = new RooDataSet( m_DataSet_TiTi->GetName(),
					    m_DataSet_TiTi->GetTitle(),
					    m_DataSet_TiTi,
					    *m_DataSet_TiTi->get(),
					    0,"Njetjet_sw" );

  m_DataSet_TiIsoTiIso_Ngamgamw =(RooDataSet*)m_DataSet_TiTi_Ngamgamw->reduce(m_isocut_st);
  m_DataSet_TiIsoTiIso_Ngamjetw =(RooDataSet*)m_DataSet_TiTi_Ngamjetw->reduce(m_isocut_st);
  m_DataSet_TiIsoTiIso_Njetgamw =(RooDataSet*)m_DataSet_TiTi_Njetgamw->reduce(m_isocut_st);
  m_DataSet_TiIsoTiIso_Njetjetw =(RooDataSet*)m_DataSet_TiTi_Njetjetw->reduce(m_isocut_st);

  m_DataSet_TiTi_Ngamgamw->Print("v");
  m_DataSet_TiIsoTiIso_Ngamgamw->Print("v");
  m_DataSet_TiTi_Ngamjetw->Print("v");
  m_DataSet_TiIsoTiIso_Ngamjetw->Print("v");
  m_DataSet_TiTi_Njetgamw->Print("v");
  m_DataSet_TiIsoTiIso_Njetgamw->Print("v");
  m_DataSet_TiTi_Njetjetw->Print("v");
  m_DataSet_TiIsoTiIso_Njetjetw->Print("v");

  int    Nbins = 13;
  double Xmin  = Commons::norm_val.first;
  double Xmax  = Commons::norm_val.second;
  hmgg       = new TH1F("hmgg","mgg data",Nbins,Xmin,Xmax);
  hmgg_Ngg   = new TH1F("hmgg_gg","mgg gg",Nbins,Xmin,Xmax);
  hmgg_Ngj   = new TH1F("hmgg_gj","mgg gj",Nbins,Xmin,Xmax);
  hmgg_Njg   = new TH1F("hmgg_jg","mgg jg",Nbins,Xmin,Xmax);
  hmgg_Njj   = new TH1F("hmgg_jj","mgg jj",Nbins,Xmin,Xmax);
  hmgg_Ngjjg = new TH1F("hmgg_gjjg","mgg gj&jg",Nbins,Xmin,Xmax);
  hmgg_red   = new TH1F("hmgg_red","mgg red bkg",Nbins,Xmin,Xmax);

  m_DataSet_TiTi->fillHistogram(hmgg,RooArgList(*m_mgg));
  m_DataSet_TiTi_Ngamgamw->fillHistogram(hmgg_Ngg,RooArgList(*m_mgg));
  m_DataSet_TiTi_Ngamjetw->fillHistogram(hmgg_Ngj,RooArgList(*m_mgg));
  m_DataSet_TiTi_Njetgamw->fillHistogram(hmgg_Njg,RooArgList(*m_mgg));
  m_DataSet_TiTi_Njetjetw->fillHistogram(hmgg_Njj,RooArgList(*m_mgg));

  hmgg_iso       = new TH1F("hmgg_iso","mgg iso data",Nbins,Xmin,Xmax);
  hmgg_Ngg_iso   = new TH1F("hmgg_gg_iso","mgg iso gg",Nbins,Xmin,Xmax);
  hmgg_Ngj_iso   = new TH1F("hmgg_gj_iso","mgg iso gj",Nbins,Xmin,Xmax);
  hmgg_Njg_iso   = new TH1F("hmgg_jg_iso","mgg iso jg",Nbins,Xmin,Xmax);
  hmgg_Njj_iso   = new TH1F("hmgg_jj_iso","mgg iso jj",Nbins,Xmin,Xmax);
  hmgg_Ngjjg_iso = new TH1F("hmgg_gjjg_iso","mgg iso gj&jg",Nbins,Xmin,Xmax);
  hmgg_red_iso   = new TH1F("hmgg_red_iso","mgg iso red bkg",Nbins,Xmin,Xmax);

  m_DataSet_TiIsoTiIso->fillHistogram(hmgg_iso,RooArgList(*m_mgg));
  m_DataSet_TiIsoTiIso_Ngamgamw->fillHistogram(hmgg_Ngg_iso,RooArgList(*m_mgg));
  m_DataSet_TiIsoTiIso_Ngamjetw->fillHistogram(hmgg_Ngj_iso,RooArgList(*m_mgg));
  m_DataSet_TiIsoTiIso_Njetgamw->fillHistogram(hmgg_Njg_iso,RooArgList(*m_mgg));
  m_DataSet_TiIsoTiIso_Njetjetw->fillHistogram(hmgg_Njj_iso,RooArgList(*m_mgg));

  hmgg->Sumw2(); 
  hmgg_Ngg->Sumw2(); 
  hmgg_Ngj->Sumw2();
  hmgg_Njg->Sumw2(); 
  hmgg_Njj->Sumw2();
  hmgg_iso->Sumw2();
  hmgg_Ngg_iso->Sumw2();
  hmgg_Ngj_iso->Sumw2();
  hmgg_Njg_iso->Sumw2(); 
  hmgg_Njj_iso->Sumw2();

  hmgg_Ngjjg=(TH1F*)hmgg_Ngj->Clone();
  hmgg_Ngjjg->Add(hmgg_Njg);

  hmgg_Ngjjg_iso=(TH1F*)hmgg_Ngj_iso->Clone();
  hmgg_Ngjjg_iso->Add(hmgg_Njg);

  hmgg_red=(TH1F*)hmgg_Njj->Clone();
  hmgg_red->Add(hmgg_Ngj,1);
  hmgg_red->Add(hmgg_Njg,1);
  hmgg_red->SetFillColor(5);

  hmgg_red_iso=(TH1F*)hmgg_Njj_iso->Clone();
  hmgg_red_iso->Add(hmgg_Ngj_iso,1);
  hmgg_red_iso->Add(hmgg_Njg_iso,1);
  hmgg_red_iso->SetFillColor(5);

}
///////////////////////////////////////////
void Bkg2DFit::StoreToRootFile(TString out)
///////////////////////////////////////////
{
  //--> Store results to a rootfile
  TString ext = Form("_mgg%d%d",(int)m_mggbin.first,(int)m_mggbin.second);
  TString filename; 
  if(out=="none")
    filename = "rootfiles/redbkg2dfit"+ext+".root";
  else
    filename = out ;
  TFile fout(filename,"RECREATE");
  fout.cd();
  if(m_do1Dfits){
    fout.Add( frame_JetL );
    fout.Add( frame_JetSL );
    fout.Add( frame_PHL );
    fout.Add( frame_PHSL );
    fout.Add( m_hJetJet );
    fout.Add( m_hJetJetFit );
  }
  if(m_do2Dfits){
    fout.Add( frame_TiTi_L );
    fout.Add( frame_TiTi_SL );
  }
  if(m_doSplot){
    fout.Add( hmgg );
    fout.Add( hmgg_Ngg );
    fout.Add( hmgg_Ngj );
    fout.Add( hmgg_Njg );
    fout.Add( hmgg_Ngjjg );
    fout.Add( hmgg_Njj );
    fout.Add( hmgg_red );
    fout.Add( hmgg_iso );
    fout.Add( hmgg_Ngg_iso );
    fout.Add( hmgg_Ngj_iso );
    fout.Add( hmgg_Njg_iso );
    fout.Add( hmgg_Njj_iso );
    fout.Add( hmgg_red_iso );
  }
  fout.Write();
  fout.Close();
}
/////////////////////////////
double Bkg2DFit::GetPurity()
////////////////////////////
{
  double Ntot    = m_NgamgamYield+m_NgamjetYield+m_NjetgamYield+m_NjetjetYield;
  double Purity  = m_NgamgamYield/(Ntot) ;
  return Purity;
}
/////////////////////////////////
double Bkg2DFit::GetPurityError()
/////////////////////////////////
{

  //Compute the stat uncertainty on the purity by 
  // propagation of the error matrix to the Purity 
  // formula
  RooMsgService::instance().saveState();
  RooMsgService::instance().deleteStream(1);
  //Calculation using the getPropagatedError method from RooFormulaVar class
  RooArgList List( *iPHPH_TiIso,*iPHJET_TiIso,*iJETPH_TiIso,*iJETJET_TiIso,
		   *Ngamgam,*Ngamjet,*Njetgam,*Njetjet );
  RooFormulaVar Purity("RFV","RFV","@0*@4/(@0*@4+@1*@5+@2*@6+@3*@7)",List);
  double error = Purity.getPropagatedError(*m_FitRes_TiTi);
  RooMsgService::instance().restoreState();
  return error;
}
/////////////////////////////////
double Bkg2DFit::GetSumofYields()
/////////////////////////////////
{
  double sum = m_NgamgamYield + m_NgamjetYield +
    m_NjetgamYield + m_NjetjetYield;
  return sum;
}
///////////////////////////////////////////////
double Bkg2DFit::GetSumofYields_NoIso()
///////////////////////////////////////////////
{
  double sum = Ngamgam->getVal() + Ngamjet->getVal() + 
    Njetgam->getVal() + Njetjet->getVal();
  return sum;
}
//////////////////////////////////////
double Bkg2DFit::GetSumofYieldsError()
//////////////////////////////////////
{
  RooMsgService::instance().saveState();
  RooMsgService::instance().deleteStream(1);
  //Calculation using the getPropagatedError method from RooFormulaVar class
  RooArgList List( *iPHPH_TiIso,*iPHJET_TiIso,*iJETPH_TiIso,*iJETJET_TiIso,
		   *Ngamgam,*Ngamjet,*Njetgam,*Njetjet );
  RooFormulaVar sum("RFV","RFV","(@0*@4+@1*@5+@2*@6+@3*@7)",List);
  double error = sum.getPropagatedError(*m_FitRes_TiTi);
  RooMsgService::instance().restoreState();
  return error;
}
////////////////////////////////////////////
double Bkg2DFit::GetSumofYieldsError_NoIso()
////////////////////////////////////////////
{
  //Calculation using the getPropagatedError method from RooFormulaVar class
  RooArgList List(*Ngamgam,*Ngamjet,*Njetgam,*Njetjet);
  RooFormulaVar sum("RFV","RFV","@0+@1+@2+@3",List);
  double error = sum.getPropagatedError(*m_FitRes_TiTi);
  return error;
}
//////////////////////////////////////////////
void Bkg2DFit::RandomizeFitResults(int NPE)
/////////////////////////////////////////////
{

}
//////////////////////////////////////////////////////////
void Bkg2DFit::RandomizeYieldsResults(int NPE,TString out)
/////////////////////////////////////////////////////////
{
  // Apply a randomization of the yield result 
  // and store it in a ttree.
  // The result will be fitted by a gaussian 
  // to extract the final yield in BkgEstimatorUpgrade

  double iPHPH_TiIso_val = iPHPH_TiIso->getVal();
  double iPHJET_TiIso_val = iPHJET_TiIso->getVal();
  double iJETPH_TiIso_val = iJETPH_TiIso->getVal();
  double iJETJET_TiIso_val = iJETJET_TiIso->getVal();

  TString filename;
  if(out== "none")
    filename = "rootfiles/2DFitResult.root";
  else
    filename = out;
  
  TFile fout(filename,"RECREATE");
  TTree *t_out = new TTree("yieldstree","yields");
  double Ngg = -999.;
  double Ngj = -999.;
  double Njg = -999.;
  double Njj = -999.;
  t_out->Branch("Ngg",&Ngg,"Ngg/D");
  t_out->Branch("Ngj",&Ngj,"Ngj/D");
  t_out->Branch("Njg",&Njg,"Njg/D");
  t_out->Branch("Njj",&Njj,"Njj/D");
  for( int ipe=0 ; ipe<NPE;ipe++){
    RooArgList Yields_List( *Ngamgam,*Ngamjet,*Njetgam,*Njetjet );
    Yields_List  = m_FitRes_TiTi->randomizePars();
    Ngg = iPHPH_TiIso_val*Ngamgam->getVal();
    Ngj = iPHJET_TiIso_val*Ngamjet->getVal();
    Njg = iJETPH_TiIso_val*Njetgam->getVal();
    Njj = iJETJET_TiIso_val*Njetjet->getVal();
    t_out->Fill();
  }
  t_out->Write();
  fout.Close();

}
///////////////////////////////////////
void Bkg2DFit::GetChiSquares(int Npar)
///////////////////////////////////////
{
  double JET_L_chi2  = frame_JetL->chiSquare("pdf_JL","data_JL",3);
  double JET_SL_chi2 = frame_JetSL->chiSquare("pdf_JSL","data_JSL",3);
  std::cout << "------ Jet fit ----------" << std::endl;
  std::cout << "Lead     ---> Chi2/ndf = " << JET_L_chi2 << std::endl;
  std::cout << "Sublead  ---> Chi2/ndf = " << JET_SL_chi2 << std::endl;

  double PH_L_chi2  = frame_PHL->chiSquare("pdf_TiL","data_TiL",6);
  double PH_SL_chi2 = frame_PHSL->chiSquare("pdf_TiSL","data_TiSL",6);
  std::cout << "------ Photon fit --------" << std::endl;
  std::cout << "Lead     ---> Chi2/ndf = " << PH_L_chi2 << std::endl;
  std::cout << "Sublead  ---> Chi2/ndf = " << PH_SL_chi2 << std::endl;

  double TiTi_L_chi2    = frame_TiTi_L->chiSquare("pdf_tot","data_L",Npar);
  double TiTi_SL_chi2   = frame_TiTi_SL->chiSquare("pdf_tot","data_SL",Npar);
  std::cout << "------ Final fit --------" << std::endl;
  std::cout << "Lead     ---> Chi2/ndf = " << TiTi_L_chi2 << std::endl;
  std::cout << "Sublead  ---> Chi2/ndf = " << TiTi_SL_chi2 << std::endl;
}

/////////////////////////////////////////////////////////////
double Bkg2DFit::GetNgamgamYield() { return m_NgamgamYield; }
/////////////////////////////////////////////////////////////
double Bkg2DFit::GetNgamjetYield() { return m_NgamjetYield; }
//////////////////////////////////////////////////////////////
double Bkg2DFit::GetNjetgamYield() { return m_NjetgamYield; }
//////////////////////////////////////////////////////////////
double Bkg2DFit::GetNjetjetYield() { return m_NjetjetYield; }
///////////////////////////////////////////////////////////////
double Bkg2DFit::GetNgamgamYieldError() { return m_NgamgamYieldError; }
////////////////////////////////////////////////////////////////////////
double Bkg2DFit::GetNgamjetYieldError() { return m_NgamjetYieldError; }
////////////////////////////////////////////////////////////////////////
double Bkg2DFit::GetNjetgamYieldError() { return m_NjetgamYieldError; }
////////////////////////////////////////////////////////////////////////
double Bkg2DFit::GetNjetjetYieldError() { return m_NjetjetYieldError; }
////////////////////////////////////////////////////////////////////////
double Bkg2DFit::GetNgamgamYield_NoIso() { return Ngamgam->getVal(); }
///////////////////////////////////////////////////////////////////////
double Bkg2DFit::GetNgamjetYield_NoIso() { return Ngamjet->getVal(); }
///////////////////////////////////////////////////////////////////////
double Bkg2DFit::GetNjetgamYield_NoIso() { return Njetgam->getVal(); }
///////////////////////////////////////////////////////////////////////
double Bkg2DFit::GetNjetjetYield_NoIso() { return Njetjet->getVal(); }
/////////////////////////////////////////////////////////////////////////////
double Bkg2DFit::GetNgamgamYieldError_NoIso() { return Ngamgam->getError(); }
/////////////////////////////////////////////////////////////////////////////
double Bkg2DFit::GetNgamjetYieldError_NoIso() { return Ngamjet->getError(); }
/////////////////////////////////////////////////////////////////////////////
double Bkg2DFit::GetNjetgamYieldError_NoIso() { return Njetgam->getError(); }
/////////////////////////////////////////////////////////////////////////////
double Bkg2DFit::GetNjetjetYieldError_NoIso() { return Njetjet->getError(); }
/////////////////////////////////////////////////////////////////////////////
TH1F* Bkg2DFit::GetDataHist() { return hmgg; }
//////////////////////////////////////////////////////
TH1F* Bkg2DFit::GetIsoDataHist() { return hmgg_iso; }
/////////////////////////////////////////////////////
RooFitResult* Bkg2DFit::GetTiTiFitResult() { return m_FitRes_TiTi; }
/////////////////////////////////////////////////////////////////////

///////////////////////////////
void Bkg2DFit::Test2DIntegral()
///////////////////////////////
{
  //--> Implement a quick test
  // of a 2D integral computation
  RooRealVar x("x","x",2,-10,10);
  RooRealVar y("y","y",2,-10,10);
  RooRealVar sigmax("sigmax","sigmax",3);
  RooRealVar meanx("meanx","meanx",0);
  RooRealVar sigmay("sigmay","sigmay",2);
  RooRealVar meany("meany","meany",1);
  //--------------------
  x.setRange("sig",-5,5);
  y.setRange("sig",-5,5);
  //------------------------------------------------
  RooGaussian gx("gx","gaussian PDF",x,meanx,sigmax);
  RooGaussian gy("gy","gaussian PDF",y,meany,sigmay);
  RooAbsReal* Gx = gx.createIntegral(x,RooFit::NormSet(x),RooFit::Range("sig"));
  RooAbsReal* Gy = gy.createIntegral(y,RooFit::NormSet(y),RooFit::Range("sig"));
  std::cout << "Gx = " << Gx->getVal() << std::endl;
  std::cout << "Gy = " << Gy->getVal() << std::endl;
  //---------------------------------------------------
  RooProdPdf gaussxy("gxy","gx*gy",RooArgList(gx,gy));
  RooArgList Var(x,y);
  RooAbsReal* Gxy = gaussxy.createIntegral(Var,RooFit::NormSet(Var),RooFit::Range("sig") );
  std::cout << "Gxy = " << Gxy->getVal() << std::endl;
}

/////////////////////
Bkg2DFit::~Bkg2DFit()
/////////////////////
{
  //destructor
  // delete m_tree;
  delete Ngamgam;
  delete Ngamjet;
  delete Njetgam;
  delete Njetjet;
  delete iPHPH_TiIso;
  delete iPHJET_TiIso;
  delete iJETPH_TiIso;
  delete iJETJET_TiIso;
  delete m_DataSet_TiTi; 
  delete m_DataSet_LoLo; 
  delete m_DataSet_Ti_L; 
  delete m_DataSet_Ti_SL;
  delete m_DataSet_TiLo; 
  delete m_DataSet_LoTi;
  delete m_DataSet_TiIsoTiIso; 
  delete m_LeadJetPdf;
  delete m_JetL;
  delete m_SubLeadJetPdf;
  delete m_JetSL;
  delete m_JetJetPdf;
  delete m_LeadPHPdf;
  delete m_PHL;
  delete m_SubLeadPHPdf;
  delete m_PHSL;
  delete m_PHPHPdf;
  delete m_PHJetPdf;
  delete m_JetPHPdf;
  delete m_TiTiPdf;
  delete m_FitRes_Jet_L;
  delete m_FitRes_Jet_SL;
  delete m_FitRes_Ti_L;
  delete m_FitRes_Ti_SL;
  delete m_FitRes_TiTi;

 }
