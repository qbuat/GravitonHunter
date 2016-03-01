#include <TIterator.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TFile.h>

#include <RooNovosibirsk.h>
#include <RooAddPdf.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>

#include "SingleParticleIsolationFitter.h"
#include "ToolsUtilities.h"
#include "ToolsSignificanceHist.h"

/////////////////////////////////////////////////////////////
SingleParticleIsolationFitter::SingleParticleIsolationFitter(TTree* tree,TString ParFile)
/////////////////////////////////////////////////////////////
{

  if (tree == 0)
    Fatal("SingleParticleIsolationFitter::SingleParticleIsolationFitter", "No input tree specified !!");
  SetEntries(tree);
  SetTree(tree);

  SetParFile();
  SetPtEtaIsoTightConvNames("_pt","_eta","_iso40","_istight","_isconv");
  
  SetIsoBounds();
  SetPtBounds();
  SetEtaBounds();
  SetIsConv();

  Init_Vars();
  Init_DataSets();

}
/////////////////////////////////
SingleParticleIsolationFitter::SingleParticleIsolationFitter()
//////////////////////////////////
{
  //Default constructor
  m_tree           = 0;
  m_nentries       = 0;
  m_nentriesforfit = 0;
  m_FitParFile     = "";

}
///////////////////////////////////////////////////////////////
void SingleParticleIsolationFitter::SetEntries(TTree* tree,int entriesforfit)
///////////////////////////////////////////////////////////////
{
  m_nentries = (int)tree->GetEntriesFast();
  if(entriesforfit==-1)
    m_nentriesforfit = m_nentries;
  else
    m_nentriesforfit = entriesforfit;
}
///////////////////////////////////////////
void SingleParticleIsolationFitter::SetTree(TTree* tree)
/////////////////////////////////////////
{
  if( m_nentriesforfit == m_nentries) m_tree = tree;
  else                                m_tree = tree->CloneTree(m_nentriesforfit);
  if (tree == 0)  Fatal("SingleParticleIsolationFitter::SetTree", "No dataset tree specified !!");
}

///////////////////////////////////////////////////////////
void SingleParticleIsolationFitter::SetPtBounds(double Xmin,double Xmax)
//////////////////////////////////////////////////////////
{
  m_ptrange.first  = Xmin;
  m_ptrange.second = Xmax;
}
///////////////////////////////////////////////////////////
void SingleParticleIsolationFitter::SetEtaBounds(double Xmin,double Xmax)
///////////////////////////////////////////////////////////
{
  m_etarange.first  = Xmin;
  m_etarange.second = Xmax;
}
/////////////////////////////////////////////////////
void SingleParticleIsolationFitter::SetIsConv(int conv)
/////////////////////////////////////////////////////
{
  m_conv = conv;
}
///////////////////////////////////////////////////////////
void SingleParticleIsolationFitter::SetIsoBounds(double Xmin,double Xmax)
///////////////////////////////////////////////////////////
{
  m_isorange.first  = Xmin;
  m_isorange.second = Xmax;
}
/////////////////////////////////////////////////
void SingleParticleIsolationFitter::SetParFile(TString ParFile)
////////////////////////////////////////////////
{m_FitParFile = ParFile; }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SingleParticleIsolationFitter::SetPtEtaIsoTightConvNames(TString pt,TString eta,TString iso,TString ti,TString conv)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
{ 
  m_pt_name         = pt;
  m_eta_name        = eta;
  m_iso_name        = iso;
  m_istight_name    = ti;
  m_isconv_name     = conv;
}

/////////////////////////////////
void SingleParticleIsolationFitter::Init_Vars()
/////////////////////////////////
{

  m_dofit = false;
  m_pt      = new RooRealVar(m_pt_name,"E_{T} [GeV]",0,0,1000000);
  m_eta     = new RooRealVar(m_eta_name,"#eta",0,-3,3);
  m_iso     = new RooRealVar(m_iso_name,"E^{iso}_{T} [GeV]",0,m_isorange.first,m_isorange.second);
  m_istight = new RooRealVar(m_istight_name,"IsTight",0,0,2);
  m_isconv  = new RooRealVar(m_isconv_name,"IsConverted",0,0,2);
  //-------------------------------------------------------
  std::cout << "READ THE DATASET"  << std::endl;
  RooArgSet Variables_MC(*m_iso,*m_pt,*m_eta,*m_isconv,"mc");
  Variables_MC.add(*m_istight);

  m_DataSet = new RooDataSet( "DataSet","DataSet",
			      Variables_MC,
			      RooFit::Import(*m_tree) );
  
  std::cout << "DATASET SIZE = " << m_DataSet->sumEntries() << std::endl;
}
/////////////////////////////////////
void SingleParticleIsolationFitter::Init_DataSets()
/////////////////////////////////////
{

  TString ptcut  = Form(m_pt_name+">=%f && "+m_pt_name+"<%f",m_ptrange.first,m_ptrange.second);
  TString etacut = Form("abs("+m_eta_name+") >=%f && abs("+m_eta_name+")<%f",m_etarange.first,m_etarange.second);
  TString isconvcut = Form(m_isconv_name+"==%i",m_conv);

  TString mastercut = etacut;
  mastercut += "&&";
  mastercut += ptcut;
  mastercut += "&&";
  mastercut += isconvcut;

  //--------------------------------------------------------
  std::cout << "******* CUT APPLIED ********" << std::endl;
  std::cout << etacut << std::endl;
  std::cout << ptcut  << std::endl;
  std::cout << isconvcut  << std::endl;
  std::cout << "****************************" << std::endl;
  //-------------------------------------------------------------
  m_DataSet_cut = 0;
  m_DataSet_Ti  = 0;
  m_DataSet_cut = (RooDataSet*)m_DataSet->reduce(mastercut);
  std::cout << "REDUCED DATASET SIZE = " << m_DataSet_cut->sumEntries() << std::endl;
  // delete m_DataSet;
  // TString tightcut = m_tight_name+"==1"; 
  TString tightcut = "_istight==1";
  std::cout << tightcut << std::endl;
  m_DataSet_Ti  = (RooDataSet*)m_DataSet_cut->reduce(tightcut);
  std::cout << "TIGHT DATASET SIZE = " << m_DataSet_Ti->sumEntries() << std::endl;
}
///////////////////////////////////////////
void SingleParticleIsolationFitter::PhotonFit(bool dofit)
///////////////////////////////////////////
{
  //  int NbinsFit = 2*abs((int)m_isorange.second-(int)m_isorange.first);
  int NbinsFit = 120;
  TString cut =  Form( m_iso_name+">%f &&"+m_iso_name+"<%f" );


  RooRealVar N_PH("N_PH","N_PH",1000.);
  meanCB  = new RooRealVar("meanCB","meanCB",4);
  sigmaCB = new RooRealVar("sigmaCB","sigmaCB",4);
  alphaCB = new RooRealVar("alphaCB","alphaCB",4);
  nCB     = new RooRealVar("nCB","nCB",4);
  m_PHSet = new RooArgSet(*meanCB,*sigmaCB,*alphaCB,*nCB,"PHSet");
  RooArgSet TiSet(*m_PHSet,"TiSet");
  TiSet.add(N_PH);
  TiSet.readFromFile(m_FitParFile,"READ","Tight" );
  m_PHPdf= new RooCBShape( "PHPdf","PHPdf",*m_iso,
			   *meanCB,*sigmaCB,*alphaCB,*nCB );

  std::cout << "Parameters values :" << N_PH.getVal() << " "
	    << meanCB->getVal() << " " << sigmaCB->getVal() << " "
	    << nCB->getVal() << " " << alphaCB->getVal() << std::endl;


  if(dofit){
    m_FitRes_Ti = m_PHPdf->fitTo( *m_DataSet_Ti,
				  RooFit::Range(m_isorange.first,m_isorange.second),
				  RooFit::Save(kTRUE),
				  RooFit::SumW2Error(kTRUE) );
    m_FitRes_Ti->Print("v");
  }
  
  frame_PH = m_iso->frame(m_isorange.first,m_isorange.second,NbinsFit);
  frame_PH->SetName("PH");
  m_DataSet_Ti->plotOn( frame_PH,
			RooFit::LineColor(1),
			RooFit::MarkerColor(1),
			RooFit::Name("data_tight") );
    m_PHPdf->plotOn( frame_PH,
		     RooFit::LineColor(4),
		     RooFit::Name("pdf_tot") );

  RooHist* RH = frame_PH->getHist("data_tight");
  m_h_iso = RooHistToHist(*RH);


}
////////////////////////////////////////
void SingleParticleIsolationFitter::Fitter(bool dofit)
///////////////////////////////////////
{
  PhotonFit(dofit);
}

/////////////////////////////////////////////////////////////////
void SingleParticleIsolationFitter::SetConstantParameters(const RooArgSet& set)
/////////////////////////////////////////////////////////////////
{
  TIterator *iterat= set.createIterator();
  RooRealVar *next = 0;
  while( (0 != (next= (RooRealVar*)iterat->Next())) )
    next->setConstant();

  delete iterat;  
}

/////////////////////////////////////////////////
void SingleParticleIsolationFitter::StoreToRootFile(TString st)
/////////////////////////////////////////////////
{
  TFile fout(st,"RECREATE");
  fout.Add(m_h_iso);
  fout.Add(frame_PH);
  fout.Write();
  RooWorkspace w("workspace");
    w.import(*m_PHPdf);
 
  w.Write();
  fout.Close();

}


///////////////////////////////////////////////////////
TH1F* SingleParticleIsolationFitter::RooHistToHist(const RooHist& RH)
///////////////////////////////////////////////////////
{
  double x=-999; double y= -999;

  std::vector<double> low_edge;
  std::vector<double> high_edge;

  for(int ip=0;ip<RH.GetN();ip++){
    int point = RH.GetPoint(ip,x,y);
    if( point == -1 ) Fatal("SingleParticleIsolationFitter::RooHistToHist()","Wrong input point 1 !");
    low_edge.push_back( x-RH.GetErrorXlow(ip) );
    high_edge.push_back( x+RH.GetErrorXhigh(ip) );
  }    

  // std::cout << " low edge = " ;
  // for(int i=0; i<(int)low_edge.size();i++){
  //   std::cout << low_edge[i] << ", ";
  // }
  // std::cout << std::endl;

  // std::cout << " high edge = " ;
  // for(int i=0; i<low_edge.size();i++){
  //   std::cout << high_edge[i] << ", ";
  // }
  // std::cout << std::endl;

  std::vector<double> binning = low_edge;
  binning.push_back( high_edge[high_edge.size()-1] );
  // for(int i=0; i<binning.size();i++){
  //   std::cout << binning[i] << ", ";
  // }

  TString name = "h_";
  name+= RH.GetName();
  TH1F* h = new TH1F( name, RH.GetTitle(),binning.size()-1, &binning[0]);
  for(int ibin=1 ; ibin<=h->GetNbinsX() ; ibin++){
    int p = RH.GetPoint(ibin-1,x,y);
    // std::cout << " p = " << p 
    // 	      << ", ibin-1 = " << ibin-1 
    // 	      << ", Nbins = " << h->GetNbinsX()
    // 	      << ", bin low edge = "<< h->GetBinLowEdge(ibin)
    // 	      << ", x = " << x 
    // 	      << ", y = " << y << std::endl;
    if( p == -1) Fatal("SingleParticleIsolationFitter::RooHistToHist()","Wrong input point !");
    h->SetBinContent(ibin,y);
  }

  //--> Set underflow and overflow to 0 
  //--> This information is lost in a TGraph.
  h->SetBinContent(0,0);
  h->SetBinContent(h->GetNbinsX()+1,0);

  return h;
}

///////////////////////////////////
SingleParticleIsolationFitter::~SingleParticleIsolationFitter()
///////////////////////////////////
{
  //destructor
  delete m_pt;
  delete m_eta;
  delete m_iso;
  delete m_istight;
  delete m_DataSet_cut;
  delete m_DataSet_Ti;
  delete m_DataSet;
  delete m_PHSet;
  delete m_TightPdf;
  delete meanCB;
  delete alphaCB;
  delete sigmaCB;
  delete nCB;
  delete m_FitRes_Ti;
  delete frame_PH;
}
