
#include "SignalTemplateFitter.h"
#include "SignalTemplateCreator.h"
#include <RooCategory.h>

/////////////////////////////////////////////
SignalTemplateFitter::SignalTemplateFitter()
/////////////////////////////////////////////
{
  //Default constructor
}
//////////////////////////////////////////////
SignalTemplateFitter::~SignalTemplateFitter()
//////////////////////////////////////////////
{
  //Default destructor
  delete m_DataSet;
  delete m_mass;
  delete m_BW_pdf;
  delete m_var_width;
  delete m_FitRes_BW;
  delete m_frame_BW;

}
///////////////////////////////////////////////
void SignalTemplateFitter::SetTree(TTree * tree)
///////////////////////////////////////////////
{
  m_tree = tree;
  m_mass  = new RooRealVar("mc_m","True mass [GeV]",0,0,10000);

  RooRealVar m_pileupweight("weight","pileup weight",1);
  RooRealVar m_mgg("mgg","Reconstructed mass [GeV]",1000);
  RooRealVar m_NPV("NPV","Number of primary vertices",5);


  RooRealVar polemass("polemass","polemass",m_polemass);
  RooRealVar coupling("coupling","coupling",m_coupling);
  RooAbsReal* m_genweight = new RooFormulaVar("genweight","SignalTemplateCreator::GetGravitonWeight(@0,@1,@2)",RooArgList(polemass,*m_mass,coupling));




  RooRealVar m_Iso_L("Iso_L","Leading photon isolation [GeV]",0);
  RooRealVar m_Iso_SL("Iso_SL","Subleading photon isolation [GeV]",0);

  RooCategory m_IsTight_L("IsTight_L","IsTight_L");
  m_IsTight_L.defineType("Tight",1);
  m_IsTight_L.defineType("NotTight",0);
  RooCategory m_IsTight_SL("IsTight_SL","IsTight_SL");
  m_IsTight_SL.defineType("Tight",1);
  m_IsTight_SL.defineType("NotTight",0);

  RooRealVar m_DeltaMass("mgg-mc_m","Reconstructed mass - True mass [GeV]",0);

  //----------------------------------------------------------------
  m_DataSet = new RooDataSet("DataSet","DataSet",
			     RooArgSet(m_Iso_L,m_IsTight_L,
				       m_Iso_SL,m_IsTight_SL,
				       *m_mass,m_mgg,m_NPV,
				       m_pileupweight,*m_genweight),
			     RooFit::Import(*m_tree),
			     RooFit::WeightVar(m_pileupweight),
			     RooFit::WeightVar(*(RooRealVar*)m_genweight));

}
///////////////////////////////////////////////////////////////
void SignalTemplateFitter::SetFitRange(double Xmin,double Xmax)
//////////////////////////////////////////////////////////////
{
  m_Xmin = Xmin;
  m_Xmax = Xmax;
}
/////////////////////////////////////////////
void SignalTemplateFitter::BW_Fitter(bool doFit)
/////////////////////////////////////////////
{
  RooRealVar m_var_polemass("polemass","polemass",m_polemass,"[GeV]");
  m_var_width    = new RooRealVar("width","width",3,0,100);
  m_BW_pdf = new RooBreitWigner("BW_pdf","BW_pdf",*m_mass,m_var_polemass,*m_var_width);
  if(doFit){
    m_FitRes_BW = m_BW_pdf->fitTo(*m_DataSet,
				  RooFit::Range(m_Xmin,m_Xmax),
				  RooFit::Save(kTRUE) );
    m_FitRes_BW->Print("v");
  }
  int NbinsFit = 200;
  m_frame_BW = m_mass->frame(m_Xmin,m_Xmax,NbinsFit);
  m_frame_BW->SetName("BW");
  m_DataSet->plotOn( m_frame_BW,RooFit::LineColor(1),RooFit::MarkerColor(1));
  m_BW_pdf->plotOn(m_frame_BW,RooFit::LineColor(2));
}

//////////////////////////////////////////
TCanvas* SignalTemplateFitter::PlotBWFit()
//////////////////////////////////////////
{
  TCanvas * c1 = new TCanvas("c1","c1",800,600);
  c1->cd();
  m_frame_BW->Draw();
  return c1;
}

//////////////////////////////////////////////////////////////////
void SignalTemplateFitter::PlotBWFit_old(bool doSave,TString name)
//////////////////////////////////////////////////////////////////
{
  TCanvas * c1 = new TCanvas("c1","c1",800,600);
  c1->cd();
  m_frame_BW->Draw();
  if(doSave)
    c1->SaveAs(Form("./plots/BWFit%d_"+name+".eps",(int)m_polemass));
}
