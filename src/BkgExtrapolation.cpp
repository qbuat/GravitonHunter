#include "BkgExtrapolation.h"

#include <RooMsgService.h>
#include <RooProdPdf.h>
#include <RooWorkspace.h>
#include <TError.h>
#include <TLine.h>
#include <TFormula.h>
#include <TF1.h>

#include "ToolsCommons.h"
#include "ToolsHPDInterval.h"
#include "ToolsUtilities.h"
#include "ToolsSignificanceHist.h"

#include "BkgRooDijetFunction.h"
#include "BkgRooErfDijetFunction.h"
#include "BkgRooDoubleErfDijetFunction.h"

#include "BkgRooErfFunction.h"
// #include "BkgRooReducibleFunction.h"

////////////////////////////////////////////////////////////////////////////////////////////////
BkgExtrapolation::BkgExtrapolation(TTree* tree,
				   std::pair<double,double > mggbin,
				     TString FitParFile) : BkgFitFramework(tree,mggbin,FitParFile)
////////////////////////////////////////////////////////////////////////////////////////////////
{
  // Constructor
  // Take a TTree with all the events (used to construct the dataset)
  // and initialize the m_tree internal variable (See BkgFitFramework)
  // mggbin is used to define the range where the invariant mass is 
  // considered for the input

  BkgFitFramework::Init(mggbin,FitParFile);
  SetFitRange();
  m_dofits         = false;
  m_doextrapol     = false;
}
//////////////////////////////////////////////////////////
BkgExtrapolation::BkgExtrapolation() : BkgFitFramework()
//////////////////////////////////////////////////////////
{
  //Default constructor
  m_dofits         = false;
  m_doextrapol     = false;

}

///////////////////////////////////////////////////////////
void BkgExtrapolation::SetFitRange(double Xmin,double Xmax)
///////////////////////////////////////////////////////////
{ m_Xmin = Xmin ; m_Xmax = Xmax ; }
///////////////////////////////////////////////////////////////////////////
void BkgExtrapolation::Init(std::pair<double,double> mggbin, TString FitParFile)
///////////////////////////////////////////////////////////////////////////
{
  // Call the BkgFitFramework to create the DataSet (m_DataSet_mggcut)
  // with all the criteria applied 
  BkgFitFramework::Init(mggbin,FitParFile);

  // On top of the criteria applied in BkgFitFramework
  // apply the isolation criteria and built the different control regions
  // Tight-AntiTight, AntiTight-Tight, AntiTight-AntiTight and the sum of the 3.
  std::cout << m_isocut_st << std::endl;
  std::cout << "Parameters from : " << m_fitparfile << std::endl;
  m_DataSet_IsoIso      = (RooDataSet*)m_DataSet_mggcut->reduce(m_isocut_st);
  m_DataSet_LoLo        = (RooDataSet*)m_DataSet_IsoIso->reduce("IsTight_L==0 && IsTight_SL==0");
  m_DataSet_LoTi        = (RooDataSet*)m_DataSet_IsoIso->reduce("IsTight_L==0 && IsTight_SL==1");
  m_DataSet_TiLo        = (RooDataSet*)m_DataSet_IsoIso->reduce("IsTight_L==1 && IsTight_SL==0");
  m_DataSet_n_fake      = (RooDataSet*)m_DataSet_IsoIso->reduce("IsTight_L==0 || IsTight_SL==0");
}
//////////////////////////////////////////////////
void BkgExtrapolation::Lead_fake_fit(bool dofit)
/////////////////////////////////////////////////
{
  // Perform the fit on the AntiTight-Tight sample
  // which is the modeling of the jet-gamma component
  // The definition of the function to be used depend on the eta 
  // category considered.
  //-------------------------------------------------------------------------------------------
  std::cout << "DEFINE THE FUNCTION " << std::endl;
  // RooArgSet *ParamSet;//("lead_fake_parameters");
  if( m_etacat == "NONE"|| m_etacat == "CC"){
    RooRealVar* lead_fake_k1 = new RooRealVar("lead_fake_k1","k1",1);
    RooRealVar* lead_fake_k2 = new RooRealVar("lead_fake_k2","k2",1);
    m_lead_fake_pdf = new RooDijetFunction( "lead_fake","lead_fake",
					    *m_mgg,*lead_fake_k1,*lead_fake_k2);
    m_lead_fake_params = new RooArgSet(*lead_fake_k1,*lead_fake_k2);
  }else if( m_etacat == "EE_O"){
    RooRealVar* lead_fake_k1     = new RooRealVar("lead_fake_k1","k1",1);
    RooRealVar* lead_fake_k2     = new RooRealVar("lead_fake_k2","k2",1);
    RooRealVar* lead_fake_mean1  = new RooRealVar("lead_fake_mean1","lead_fake_mean1",1);
    RooRealVar* lead_fake_mean2  = new RooRealVar("lead_fake_mean2","lead_fake_mean1",1);
    RooRealVar* lead_fake_width1 = new RooRealVar("lead_fake_width1","width1",1);
    RooRealVar* lead_fake_width2 = new RooRealVar("lead_fake_width2","width2",1);
    m_lead_fake_pdf = new RooDoubleErfDijetFunction("lead_fake","lead_fake",
						    *m_mgg,*lead_fake_k1,*lead_fake_k2,
						    *lead_fake_mean1,*lead_fake_width1,
						    *lead_fake_mean2,*lead_fake_width2 );

    m_lead_fake_params = new RooArgSet(*lead_fake_k1,*lead_fake_k2,
				       *lead_fake_mean1,*lead_fake_mean2,
				       *lead_fake_width1,*lead_fake_width2);
  }else{
    RooRealVar* lead_fake_k1    = new RooRealVar("lead_fake_k1","k1",1);
    RooRealVar* lead_fake_k2    = new RooRealVar("lead_fake_k2","k2",1);
    RooRealVar* lead_fake_mean  = new RooRealVar("lead_fake_mean","k3",1);
    RooRealVar* lead_fake_width = new RooRealVar("lead_fake_width","width",1);
    m_lead_fake_pdf = new RooErfDijetFunction("lead_fake","lead_fake",
					      *m_mgg,
					      *lead_fake_k1,*lead_fake_k2,
					      *lead_fake_mean,*lead_fake_width);
    m_lead_fake_params = new RooArgSet(*lead_fake_k1,*lead_fake_k2,
				       *lead_fake_mean,*lead_fake_width);
  }
  //-------------------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------------------
  

  //-------------------------------------------------------------------------------
  std::cout << "READ PARAMETERS " << std::endl;
  m_lead_fake_params->readFromFile(m_fitparfile,"READ","ANTITIGHT-TIGHT");
  if(dofit){
    std::cout << "FIT START " << std::endl;
    m_FitRes_lead_fake = m_lead_fake_pdf->fitTo( *m_DataSet_LoTi,
						 RooFit::Range(m_Xmin,m_Xmax),
						 RooFit::Save(kTRUE) );
    m_FitRes_lead_fake->Print("v");
  }


  // Create the RooPlot where the graphical representation of the fit
  // is stored
  int NbinsFit = 150;
  //  frame_lead_fake = m_mgg->frame(m_Xmin,m_Xmax,NbinsFit);
  frame_lead_fake = m_mgg->frame(m_mggbin.first,m_mggbin.second,NbinsFit);
  frame_lead_fake->SetName("lead_fake");
  m_DataSet_LoTi->plotOn( frame_lead_fake,
			  RooFit::LineColor(1),
			  RooFit::MarkerColor(1),
			  RooFit::Name("LoTi") ,
			  RooFit::Range(m_mggbin.first,m_mggbin.second));
  m_lead_fake_pdf->plotOn( frame_lead_fake,
			   RooFit::LineColor(3),
			   RooFit::Name("lead_fake_pdf"),
			   RooFit::Range(m_mggbin.first,m_mggbin.second));
  //  			   RooFit::Range(180,3000));

  //---------------------------------------------------------------------------------
}
////////////////////////////////////////////////////
void BkgExtrapolation::Sublead_fake_fit(bool dofit)
////////////////////////////////////////////////////
{
  // Perform the fit on the Tight-AntiTight sample
  // which is the modeling of the gamma-jet component
  // The definition of the function to be used depend on the eta 
  // category considered.
  std::cout << "DEFINE THE FUNCTION " << std::endl;

  // RooArgSet *ParamSet;
  if( m_etacat == "NONE"|| m_etacat == "CC"){
    RooRealVar* sublead_fake_k1 = new RooRealVar("sublead_fake_k1","k1",1);
    RooRealVar* sublead_fake_k2 = new RooRealVar("sublead_fake_k2","k2",1);
    m_sublead_fake_pdf = new RooDijetFunction( "sublead_fake","sublead_fake",
					    *m_mgg,*sublead_fake_k1,*sublead_fake_k2);
    m_sublead_fake_params = new RooArgSet(*sublead_fake_k1,*sublead_fake_k2);
  }else if( m_etacat == "EE_O"){
    RooRealVar* sublead_fake_k1     = new RooRealVar("sublead_fake_k1","k1",1);
    RooRealVar* sublead_fake_k2     = new RooRealVar("sublead_fake_k2","k2",1);
    RooRealVar* sublead_fake_mean1  = new RooRealVar("sublead_fake_mean1","mean1",1);
    RooRealVar* sublead_fake_mean2  = new RooRealVar("sublead_fake_mean2","mean2",1);
    RooRealVar* sublead_fake_width1 = new RooRealVar("sublead_fake_width1","width1",1);
    RooRealVar* sublead_fake_width2 = new RooRealVar("sublead_fake_width2","width2",1);
 
   m_sublead_fake_pdf = new RooDoubleErfDijetFunction("sublead_fake","sublead_fake",
						      *m_mgg,
						      *sublead_fake_k1,*sublead_fake_k2,
						      *sublead_fake_mean1,*sublead_fake_width1,
						      *sublead_fake_mean2,*sublead_fake_width2 ); 

    m_sublead_fake_params = new RooArgSet( *sublead_fake_k1,*sublead_fake_k2,
					   *sublead_fake_mean1,*sublead_fake_mean2,
					   *sublead_fake_width1,*sublead_fake_width2 );
  }else{
    RooRealVar* sublead_fake_k1    = new RooRealVar("sublead_fake_k1","k1",1);
    RooRealVar* sublead_fake_k2    = new RooRealVar("sublead_fake_k2","k2",1);
    RooRealVar* sublead_fake_mean  = new RooRealVar("sublead_fake_mean","k3",1);
    RooRealVar* sublead_fake_width = new RooRealVar("sublead_fake_width","width",1);
    m_sublead_fake_pdf = new RooErfDijetFunction("sublead_fake","sublead_fake",
						 *m_mgg,
						 *sublead_fake_k1,*sublead_fake_k2,
						 *sublead_fake_mean,*sublead_fake_width);
    m_sublead_fake_params = new RooArgSet( *sublead_fake_k1,*sublead_fake_k2,
					   *sublead_fake_mean,*sublead_fake_width );
  }
    




  //------------------------------------------------------------------------------------
  std::cout << "READ PARAMETERS " << std::endl;
  m_sublead_fake_params->readFromFile(m_fitparfile,"READ","TIGHT-ANTITIGHT");
  if(dofit){
    std::cout << "FIT START " << std::endl;
    m_FitRes_sublead_fake = m_sublead_fake_pdf->fitTo( *m_DataSet_TiLo,
						       RooFit::Range(m_Xmin,m_Xmax),
						       RooFit::Save(kTRUE) );
    m_FitRes_sublead_fake->Print("v");
  }

  // Create the RooPlot where the graphical representation of the fit
  // is stored
  int NbinsFit = 150;
  //  frame_sublead_fake =  m_mgg->frame(m_Xmin,m_Xmax,NbinsFit);
  frame_sublead_fake = m_mgg->frame(m_mggbin.first,m_mggbin.second,NbinsFit*3);
  frame_sublead_fake->SetName("sublead_fake");
  m_DataSet_TiLo->plotOn( frame_sublead_fake,
			  RooFit::LineColor(1),
			  RooFit::MarkerColor(1),
			  RooFit::Name("TiLo") ,
			  //			  RooFit::Range(180,3000));
			  RooFit::Range(m_mggbin.first,m_mggbin.second));

  m_sublead_fake_pdf->plotOn( frame_sublead_fake,
			      RooFit::LineColor(2),
			      RooFit::Name("sublead_fake_pdf") ,
			      //			      RooFit::Range(180,3000));
			  RooFit::Range(m_mggbin.first,m_mggbin.second));

  //------------------------------------------------------------------------------------
}
//////////////////////////////////////////////////
void BkgExtrapolation::Double_fake_fit(bool dofit)
/////////////////////////////////////////////////
{
  // Perform the fit on the AntiTight-AntiTight sample
  // which is the modeling of the jet-jet component
  // The definition of the function to be used depend on the eta 
  // category considered.
  std::cout << "DEFINE THE FUNCTION " << std::endl;

  // RooArgSet* ParamSet;
  if( m_etacat == "NONE"|| m_etacat == "CC"){
    RooRealVar* double_fake_k1 = new RooRealVar("double_fake_k1","k1",1);
    RooRealVar* double_fake_k2 = new RooRealVar("double_fake_k2","k2",1);
    m_double_fake_pdf = new RooDijetFunction( "double_fake","double_fake",
					      *m_mgg,*double_fake_k1,*double_fake_k2);
    m_double_fake_params = new RooArgSet(*double_fake_k1,*double_fake_k2);
  }else if( m_etacat == "EE_O"){
    RooRealVar* double_fake_k1     = new RooRealVar("double_fake_k1","k1",1);
    RooRealVar* double_fake_k2     = new RooRealVar("double_fake_k2","k2",1);
    RooRealVar* double_fake_mean1  = new RooRealVar("double_fake_mean1","mean1",1);
    RooRealVar* double_fake_mean2  = new RooRealVar("double_fake_mean2","mean2",1);
    RooRealVar* double_fake_width1 = new RooRealVar("double_fake_width1","width1",1);
    RooRealVar* double_fake_width2 = new RooRealVar("double_fake_width2","width2",1);
    m_double_fake_pdf = new RooDoubleErfDijetFunction( "double_fake","double_fake",
						       *m_mgg,
						       *double_fake_k1,*double_fake_k2,
						       *double_fake_mean1,*double_fake_width1,
						       *double_fake_mean2,*double_fake_width2 ); 
    m_double_fake_params = new RooArgSet( *double_fake_k1,*double_fake_k2,
					  *double_fake_mean1,*double_fake_mean2,
					  *double_fake_width1,*double_fake_width2 );
  }else{
    RooRealVar* double_fake_k1    = new RooRealVar("double_fake_k1","k1",1);
    RooRealVar* double_fake_k2    = new RooRealVar("double_fake_k2","k2",1);
    RooRealVar* double_fake_mean  = new RooRealVar("double_fake_mean","k3",1);
    RooRealVar* double_fake_width = new RooRealVar("double_fake_width","width",1);
    m_double_fake_pdf = new RooErfDijetFunction("double_fake","double_fake",
						 *m_mgg,
						 *double_fake_k1,*double_fake_k2,
						 *double_fake_mean,*double_fake_width);
    m_double_fake_params = new RooArgSet(*double_fake_k1,*double_fake_k2,*double_fake_mean,*double_fake_width);
  }


 
  //------------------------------------------------------------------------------------
  std::cout << "READ PARAMETERS " << std::endl;
  m_double_fake_params->readFromFile(m_fitparfile,"READ","ANTITIGHT-ANTITIGHT");
 if(dofit){
    std::cout << "FIT START " << std::endl;
    m_FitRes_double_fake = m_double_fake_pdf->fitTo( *m_DataSet_LoLo,
						     RooFit::Range(m_Xmin,m_Xmax),
						     RooFit::Save(kTRUE) );
    m_FitRes_double_fake->Print("v");
  }

  // Create the RooPlot where the graphical representation of the fit
  // is stored
  int NbinsFit = 150;
  //  frame_double_fake =  m_mgg->frame(m_Xmin,m_Xmax,NbinsFit);
  frame_double_fake = m_mgg->frame(m_mggbin.first,m_mggbin.second,NbinsFit*3);
  frame_double_fake->SetName("double_fake");
  m_DataSet_LoLo->plotOn( frame_double_fake,
			  RooFit::LineColor(1),
			  RooFit::MarkerColor(1),
			  RooFit::Name("LoLo"),
			  //			  RooFit::Range(180,300));
			  RooFit::Range(m_mggbin.first,m_mggbin.second));

  m_double_fake_pdf->plotOn( frame_double_fake,
			     RooFit::LineColor(4),
			     RooFit::Name("double_fake_pdf"),
			     //			     RooFit::Range(180,300) );
			  RooFit::Range(m_mggbin.first,m_mggbin.second));

  //------------------------------------------------------------------------------------
}
/////////////////////////////////////////////
void BkgExtrapolation::n_fake_fit(bool dofit)
/////////////////////////////////////////////
{

  // Perform the fit on the sum of the 3 control samples
  // which was the modeling of the reducible bkg in the 
  // first 7TeV analysis
  // The definition of the function to be used depend on the eta 
  // category considered.
  std::cout << "DEFINE THE FUNCTION " << std::endl;

  // RooArgSet* ParamSet;
  if( m_etacat == "NONE"|| m_etacat == "CC"){
    RooRealVar* n_fake_k1 = new RooRealVar("n_fake_k1","k1",1);
    RooRealVar* n_fake_k2 = new RooRealVar("n_fake_k2","k2",1);
    m_n_fake_pdf = new RooDijetFunction( "n_fake","n_fake",
					    *m_mgg,*n_fake_k1,*n_fake_k2);
    m_n_fake_params = new RooArgSet(*n_fake_k1,*n_fake_k2);
  }else if( m_etacat == "EE_O"){
    RooRealVar* n_fake_k1     = new RooRealVar("n_fake_k1","k1",1);
    RooRealVar* n_fake_k2     = new RooRealVar("n_fake_k2","k2",1);
    RooRealVar* n_fake_mean1  = new RooRealVar("n_fake_mean1","mean1",1);
    RooRealVar* n_fake_mean2  = new RooRealVar("n_fake_mean2","mean2",1);
    RooRealVar* n_fake_width1 = new RooRealVar("n_fake_width1","width1",1);
    RooRealVar* n_fake_width2 = new RooRealVar("n_fake_width2","width2",1);
    m_n_fake_pdf = new RooDoubleErfDijetFunction( "n_fake","n_fake",
						  *m_mgg,
						  *n_fake_k1,*n_fake_k2,
						  *n_fake_mean1,*n_fake_width1,
						  *n_fake_mean2,*n_fake_width2 ); 
    m_n_fake_params = new RooArgSet( *n_fake_k1,*n_fake_k2,
				     *n_fake_mean1,*n_fake_mean2,
				     *n_fake_width1,*n_fake_width2 );

  }else{
    RooRealVar* n_fake_k1       = new RooRealVar("n_fake_k1","k1",1);
    RooRealVar* n_fake_k2       = new RooRealVar("n_fake_k2","k2",1);
    RooRealVar* n_fake_mean     = new RooRealVar("n_fake_mean","k3",1);
    RooRealVar* n_fake_width    = new RooRealVar("n_fake_width","width",1);
    m_n_fake_pdf = new RooErfDijetFunction("n_fake","n_fake",
						 *m_mgg,
						 *n_fake_k1,*n_fake_k2,
						 *n_fake_mean,*n_fake_width);
    m_n_fake_params = new RooArgSet(*n_fake_k1,*n_fake_k2,*n_fake_mean,*n_fake_width);

  }


  //------------------------------------------------------------------------------------
  std::cout << "READ PARAMETERS " << std::endl;
  m_n_fake_params->readFromFile(m_fitparfile,"READ","n_fake");
  if(dofit){
    std::cout << "FIT START " << std::endl;
    m_FitRes_n_fake = m_n_fake_pdf->fitTo( *m_DataSet_n_fake,
					   RooFit::Range(m_Xmin,m_Xmax),
					   RooFit::Save(kTRUE) );
    m_FitRes_n_fake->Print("v");
  }

  // Create the RooPlot where the graphical representation of the fit
  // is stored
  int NbinsFit = 150;
  //  frame_n_fake = m_mgg->frame(m_Xmin,m_Xmax,NbinsFit);
  frame_n_fake = m_mgg->frame(m_mggbin.first,m_mggbin.second,NbinsFit*3);
  frame_n_fake->SetName("n_fake");
  m_DataSet_n_fake->plotOn( frame_n_fake,
			    RooFit::LineColor(1),
			    RooFit::MarkerColor(1),
			    RooFit::Name("n_fake_hist")
			    );
  m_n_fake_pdf->plotOn( frame_n_fake,
			RooFit::LineColor(1),
			RooFit::Name("n_fake_pdf") );
  //------------------------------------------------------------------------------------

}
//////////////////////////////////////////
void BkgExtrapolation::Fitter(bool dofits)
//////////////////////////////////////////
{
  //This method calls the four methods
  // above. It is meant to be called outside the class
  // but the individual fit can also be called separately

  m_dofits = dofits;
  //----------------------------
  std::cout << " LEAD FAKE FIT " << std::endl;
  Lead_fake_fit(dofits); 
  std::cout << " SUBLEAD FAKE FIT " << std::endl;
  Sublead_fake_fit(dofits);
  std::cout << " DOUBLE FAKE FIT " << std::endl;
  Double_fake_fit(dofits);
  std::cout << " N FAKE FIT " << std::endl;
  n_fake_fit(dofits);  
  m_FitRes_lead_fake->Print();
  m_FitRes_sublead_fake->Print();
  m_FitRes_double_fake->Print();
  m_FitRes_n_fake->Print();
}
///////////////////////////////////////////////////
void BkgExtrapolation::Extrapolate(bool doextrapol)
///////////////////////////////////////////////////
{
  // Perfom the extrapolation of the fit up to the
  // end of the mgg spectrum (defined by 
  // Commons::binning[Commons::nBins])
  // The method used is BuildHistogram and is private to 
  // this class.
  m_doextrapol   = doextrapol;

  std::cout << "HISTOGRAM BUILDING" << std::endl;
  //----------------------------------------------------------------------
  h_lead_fake    = BuildHistogram(*m_lead_fake_pdf, "h_leading_fake");
  h_sublead_fake = BuildHistogram(*m_sublead_fake_pdf, "h_subleading_fake");
  h_double_fake  = BuildHistogram(*m_double_fake_pdf, "h_double_fake");
  h_n_fake       = BuildHistogram(*m_n_fake_pdf, "h_nfake");
  //----------------------------------------------------------------------
  std::cout << "RANDOMIZATION" << std::endl;
  //-----------------------------
  Randomizer(false);

  //------------------------------
  h_n_fake->SetLineColor(1);
  h_sublead_fake->SetLineColor(2);
  h_lead_fake->SetLineColor(3);
  h_double_fake->SetLineColor(4);
  //------------------------------


 
}
///////////////////////////////////////////////
void BkgExtrapolation::PlotResults(bool dosave)
//////////////////////////////////////////////
{
}
//////////////////////////////////////////////
void BkgExtrapolation::GetChiSquares(int Npar)
///////////////////////////////////////////////
{

}
////////////////////////////////////////////////////////
void BkgExtrapolation::StoreToRootFile(TString filename)
////////////////////////////////////////////////////////
{
  // Dump the results in a rootfile as histograms 
  // to be used by other classes (BkgEstimatorUpgrade)
  // Dump also graphs whith the stat uncertainties of the fit
  // are propagated.
  // Store RooPlot to be drawn.
  // Store also the PDFs in a RooWorspace

  TFile f(filename,"RECREATE");
  //--------------------
  f.Add(h_n_fake);
  f.Add(h_lead_fake);
  f.Add(h_sublead_fake);
  f.Add(h_double_fake);
  //------------------------
  f.Add(graph_n_fake);
  f.Add(graph_lead_fake);
  f.Add(graph_sublead_fake);
  f.Add(graph_double_fake);
  //------------------------
  f.Add(frame_lead_fake);
  f.Add(frame_sublead_fake);
  f.Add(frame_double_fake);
  f.Write();

  RooWorkspace w("workspace");
  w.import(*m_lead_fake_pdf);
  w.import(*m_sublead_fake_pdf);
  w.import(*m_double_fake_pdf);
  w.Write();

  f.Close();
}

////////////////////////////////////////////////
void BkgExtrapolation::Randomizer(bool VERBOSE)
////////////////////////////////////////////////
{
  // Randomize the fit result using the covariance matrix
  // to estimate the stat error on the shape computation
  RooMsgService::instance().saveState();
  RooMsgService::instance().deleteStream(1);
  std::cout << "n_fake error computation" << std::endl;
  graph_n_fake = FitErrorPropagator( *h_n_fake,*m_n_fake_pdf,
				     *m_DataSet_n_fake,*m_FitRes_n_fake,
				     VERBOSE );
  std::cout << "lead_fake error computation" << std::endl;
  graph_lead_fake = FitErrorPropagator( *h_lead_fake,*m_lead_fake_pdf,
					*m_DataSet_LoTi,*m_FitRes_lead_fake,
					VERBOSE );
  std::cout << "sublead_fake error computation" << std::endl;
  graph_sublead_fake = FitErrorPropagator( *h_sublead_fake,*m_sublead_fake_pdf,
					   *m_DataSet_TiLo,*m_FitRes_sublead_fake,
					   VERBOSE );
  std::cout << "double_fake error computation" << std::endl;
  graph_double_fake = FitErrorPropagator( *h_double_fake,*m_double_fake_pdf,
					  *m_DataSet_LoLo,*m_FitRes_double_fake,
					  VERBOSE );
  graph_n_fake->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  graph_lead_fake->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  graph_double_fake->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  graph_sublead_fake->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  RooMsgService::instance().restoreState();

}
////////////////////////////////////////////////////////////////////////////////////
TGraphAsymmErrors* BkgExtrapolation::FitErrorPropagator( const TH1D& hist, 
							  const RooAbsPdf& pdf,
							  const RooAbsData& data,
							  const RooFitResult& fitres,
							  bool VERBOSE )
/////////////////////////////////////////////////////////////////////////////////////
{
  //This method randomize the pdf and perform the extrapolation
  // The number of PE is controled by m_NPE.
  // Each mgg bin is populated and the 
  // Highest Probability Interval(@68%) is taken as the uncertainty 
  // on the histogram hist.

  TString histname = hist.GetName();
  //----------------------------------------------------
  TH1D up_dev   = *(TH1D*)hist.Clone("up_deviation");
  TH1D down_dev = *(TH1D*)hist.Clone("down_deviation");
  //----------------------------------------------------
  TH1D*  dev[Commons::nBins] ;
  TLine* line[Commons::nBins] ;
  TLine* line_xdown[Commons::nBins];
  TLine* line_xup[Commons::nBins];
  TH1D*  hrand[m_NPE];
  //------------------------------------------
  RooArgSet* argset = pdf.getParameters(data);
  RooArgList par(*argset);

  //====================== Loop over all the Pseudo Experiments ====================//
  for( int ipe = 0 ; ipe <m_NPE;ipe++){
    if( ipe%1000==0 && ipe!=0 )
      std::cout << "PE NUMBER " << ipe << std::endl;
    par  = fitres.randomizePars();
    hrand[ipe] = BuildHistogram(pdf,histname+Form("_h%d",ipe));
  }//---------------- End of the loop over all Pseudo Experiments --------------------

  //====================== Loop over all the bins ================================//
  for ( int ibin=0;ibin<=up_dev.GetNbinsX();ibin++){
    double center = hist.GetBinContent(ibin);
    std::vector<double> integrals;
    double Xmin = hrand[0]->GetBinContent(ibin);
    double Xmax = hrand[0]->GetBinContent(ibin);
    for( int ipe = 0 ; ipe < m_NPE; ipe++){
      double I = hrand[ipe]->GetBinContent(ibin);
      if(I>Xmax) Xmax=I;
      if(I<Xmin) Xmin=I;
      integrals.push_back(I);
    }
    //------------------------------------------------------------------
    TString title = Form("%d",(int)up_dev.GetBinCenter(ibin));
    dev[ibin] = new TH1D(histname+Form("dev%d",ibin),"deviation",200,Xmin,Xmax);
    dev[ibin]->GetXaxis()->SetTitle(title);
    //------------------------------------
    for( int ipe = 0 ; ipe < m_NPE; ipe++)
      dev[ibin]->Fill(integrals[ipe]);
    //------------------------------------------------------------------
    line[ibin] = new TLine(center,0,center,
			   dev[ibin]->GetBinContent(dev[ibin]->GetMaximumBin()));
    if( VERBOSE ){
      std::cout << "bin n0 : " << ibin << " Xmin = " << Xmin << " Xmax = " << Xmax 
		<< " Delta X/ref(%) = " << 100*(Xmax-Xmin)/center 
		<< " Underflow : " << dev[ibin]->GetBinContent(0)
		<< " (Mean-ref)/ref(%) : " << 100*(dev[ibin]->GetMean()-center)/center 
		<< std::endl;
    }
    HPDInterval HPDI(*dev[ibin]);//See ToolsHPDInterval
    std::pair<double,double> interval = HPDI.GetHPDInterval();
    line_xdown[ibin] = new TLine( interval.first,0,interval.first,
				  dev[ibin]->GetBinContent(dev[ibin]->GetMaximumBin()) );
    line_xup[ibin] = new TLine( interval.second,0,interval.second,
				dev[ibin]->GetBinContent(dev[ibin]->GetMaximumBin()) );
    up_dev.SetBinContent(ibin,(interval.second) );
    down_dev.SetBinContent(ibin,(interval.first) );
  }//-------------------- End of the loop over all the bins ---------------------------------
  for( int ipe = 0 ; ipe < m_NPE; ipe++) delete hrand[ipe];


  //---------- Build output TGraphAsymmErrors() ---------
  TGraphAsymmErrors* gr_temp = new TGraphAsymmErrors();
  gr_temp->SetName( "graph_"+histname );
  for(int ibin=1;ibin<=hist.GetNbinsX();ibin++){
    double center = hist.GetBinCenter(ibin);
    double value  = hist.GetBinContent(ibin);
    double ex_up  = hist.GetBinLowEdge(ibin+1)-center;
    double ex_do  = center-hist.GetBinLowEdge(ibin);
    double ey_up  = up_dev.GetBinContent(ibin) - value;
    double ey_do  = value - down_dev.GetBinContent(ibin);
    gr_temp->SetPoint(ibin-1, center, value);
    gr_temp->SetPointError(ibin-1, ex_do, ex_up, ey_do, ey_up);
  }//------------------------------------------------------------

  if(VERBOSE){
    for(int ibin = 1 ; ibin<=up_dev.GetNbinsX();ibin++){
      TString cname = histname+Form("_c%d",ibin);
      TCanvas c(cname,cname,800,600);
      dev[ibin]->Draw();
      line[ibin]->SetLineColor(2);
      line_xdown[ibin]->SetLineColor(4);
      line_xup[ibin]->SetLineColor(4);
      line[ibin]->Draw("same");
      line_xdown[ibin]->Draw("same");
      line_xup[ibin]->Draw("same");
      c.SaveAs("./plots/randomized/"+cname+".eps");
      c.Close();
    }
  }
  for ( int ibin=0;ibin<=up_dev.GetNbinsX();ibin++){
    delete dev[ibin] ;
    delete line[ibin] ;
    delete line_xdown[ibin];
    delete line_xup[ibin];
  }
  return gr_temp;
}
/////////////////////////////////////////////////////////////////////////
TH1D* BkgExtrapolation::BuildHistogram(const RooAbsPdf& pdf,TString name)
/////////////////////////////////////////////////////////////////////////
{
  // Build the histogram from a TF1
  // build using the asTF method of RooAbsPDf
  // The stat error is set to zero.
  TH1D *h = new TH1D(name,name,Commons::nBins,Commons::binning);
  h->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  h->GetYaxis()->SetTitle("Arbitrary unit");


  TF1* f = pdf.asTF(RooArgList(*m_mgg));
  for(int ibin=0;ibin<Commons::nBins;ibin++)
    h->SetBinContent(ibin+1,f->Integral(Commons::binning[ibin],Commons::binning[ibin+1]));

  // pdf.fillHistogram(h,RooArgList(*m_mgg));//--> WRONG BEHAVIOUR !!!!
  // double norm = h->Integral( Commons::norm_bin.first,
  // 			     Commons::norm_bin.second );
  // h->Scale(1./norm);
  for(int ibin=0;ibin<=Commons::nBins+1;ibin++)
    h->SetBinError(ibin,0);
  return h;
}
/////////////////////////////////////
BkgExtrapolation::~BkgExtrapolation()
////////////////////////////////////
{
  delete m_DataSet_IsoIso; 
  delete m_double_fake_pdf;
  delete m_lead_fake_pdf;
  delete m_sublead_fake_pdf;
  delete m_n_fake_pdf;
  delete m_FitRes_double_fake;
  delete m_FitRes_lead_fake;
  delete m_FitRes_sublead_fake;
  delete m_FitRes_n_fake;
}


