#include <TIterator.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TFile.h>

#include <RooNovosibirsk.h>
#include <RooAddPdf.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>

#include "IsolationFitter.h"
#include "ToolsUtilities.h"
#include "ToolsSignificanceHist.h"

/////////////////////////////////////////////////////////////
IsolationFitter::IsolationFitter(TTree* tree,TString ParFile)
/////////////////////////////////////////////////////////////
{

  if (tree == 0)
    Fatal("IsolationFitter::IsolationFitter", "No input tree specified !!");
  SetEntries(tree);
  SetTree(tree);
  SetStreamType();
  SetParFile();
  SetPtEtaIsoTightLoosePrimeNames("pT_L","eta_L","Iso_L","Tight_L","IsLoosePrime4_L");
  
  SetMggBounds();
  SetIsoBounds();
  SetPtBounds();
  SetEtaBounds();
  SetNPVBounds();
  Init_Vars();
  Init_DataSets();

}
/////////////////////////////////
IsolationFitter::IsolationFitter()
//////////////////////////////////
{
  //Default constructor
  m_tree           = 0;
  m_nentries       = 0;
  m_nentriesforfit = 0;
  m_FitParFile     = "";

}
///////////////////////////////////////////////////////////////
void IsolationFitter::SetEntries(TTree* tree,int entriesforfit)
///////////////////////////////////////////////////////////////
{
  m_nentries = (int)tree->GetEntriesFast();
  if(entriesforfit==-1)
    m_nentriesforfit = m_nentries;
  else
    m_nentriesforfit = entriesforfit;
}
///////////////////////////////////////////
void IsolationFitter::SetTree(TTree* tree)
/////////////////////////////////////////
{
  if( m_nentriesforfit == m_nentries)
    m_tree = tree;
  else
    m_tree = tree->CloneTree(m_nentriesforfit);
  if (tree == 0)
    Fatal("IsolationFitter::SetTree", "No dataset tree specified !!");
}

///////////////////////////////////////////////////////////
void IsolationFitter::SetMggBounds(double Xmin,double Xmax)
//////////////////////////////////////////////////////////
{
  m_mggrange.first  = Xmin;
  m_mggrange.second = Xmax;
}
///////////////////////////////////////////////////////////
void IsolationFitter::SetPtBounds(double Xmin,double Xmax)
//////////////////////////////////////////////////////////
{
  m_ptrange.first  = Xmin;
  m_ptrange.second = Xmax;
}
///////////////////////////////////////////////////////////
void IsolationFitter::SetEtaBounds(double Xmin,double Xmax)
///////////////////////////////////////////////////////////
{
  m_etarange.first  = Xmin;
  m_etarange.second = Xmax;
}
/////////////////////////////////////////////////////
void IsolationFitter::SetNPVBounds(int Xmin,int Xmax)
/////////////////////////////////////////////////////
{
  m_npvrange.first  = Xmin;
  m_npvrange.second = Xmax;
}
/////////////////////////////////////////////////////
void IsolationFitter::SetMuBounds(int Xmin,int Xmax)
/////////////////////////////////////////////////////
{
  m_murange.first  = Xmin;
  m_murange.second = Xmax;
}
///////////////////////////////////////////////////////////
void IsolationFitter::SetIsoBounds(double Xmin,double Xmax)
///////////////////////////////////////////////////////////
{
  m_isorange.first  = Xmin;
  m_isorange.second = Xmax;
}
/////////////////////////////////////////////////////////
void IsolationFitter::SetIsoNorm(double Xmin,double Xmax)
/////////////////////////////////////////////////////////
{
  m_isonormrange.first  = Xmin;
  m_isonormrange.second = Xmax;
}
/////////////////////////////////////////////////
void IsolationFitter::SetParFile(TString ParFile)
////////////////////////////////////////////////
{m_FitParFile = ParFile; }
/////////////////////////////////////////////////
void IsolationFitter::SetStreamType(bool data)
////////////////////////////////////////////////
{ m_data = data; }

//////////////////////////////////////////////////////////////////////////////////////////
void IsolationFitter::SetPtEtaIsoTightNames(TString pt,TString eta,TString iso,TString ti)
/////////////////////////////////////////////////////////////////////////////////////////
{ 
  m_pt_name    = pt;
  m_eta_name   = eta;
  m_iso_name   = iso;
  m_tight_name = ti;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void IsolationFitter::SetPtEtaIsoTightLoosePrimeNames(TString pt,TString eta,TString iso,TString ti,TString lp)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
{ 
  m_pt_name         = pt;
  m_eta_name        = eta;
  m_iso_name        = iso;
  m_tight_name      = ti;
  m_looseprime_name = lp;
}

/////////////////////////////////
void IsolationFitter::Init_Vars()
/////////////////////////////////
{

  m_dofit = false;
  m_pt     = new RooRealVar(m_pt_name,"E_{T} [GeV]",10,0,10000);
  m_eta    = new RooRealVar(m_eta_name,"#eta",0,-3,3);
  m_Iso    = new RooRealVar(m_iso_name,"E^{iso}_{T} [GeV]",0,m_isorange.first,m_isorange.second);
  m_mgg    = new RooRealVar("mgg","m_{#gamma#gamma} [GeV]",100,0,4000);
  m_npv    = new RooRealVar("NPV","NPV",100,0,4000);
  m_mu     = new RooRealVar("mu","mu",100,0,4000);
  m_weight = new RooRealVar("weight","weight",1,0,1000);
  //-------------------------------------------------------
  m_IsTight_L = new RooCategory("IsTight_L","IsTight_L");
  m_IsTight_L->defineType("Tight",1);
  m_IsTight_L->defineType("NotTight",0);
  m_IsTight_SL = new RooCategory("IsTight_SL","IsTight_SL");
  m_IsTight_SL->defineType("Tight",1);
  m_IsTight_SL->defineType("NotTight",0);
  //--------------------------------------------------------
  m_IsLoosePrime_L = new RooCategory("IsLoosePrime4_L","IsLoosePrime4_L");
  m_IsLoosePrime_L->defineType("LoosePrime",1);
  m_IsLoosePrime_L->defineType("NotLoosePrime",0);
  m_IsLoosePrime_SL = new RooCategory("IsLoosePrime4_SL","IsLoosePrime4_SL");
  m_IsLoosePrime_SL->defineType("LoosePrime",1);
  m_IsLoosePrime_SL->defineType("NotLoosePrime",0);
  //--------------------------------------------------------
  std::cout << "READ THE DATASET"  << std::endl;
  RooArgSet Variables_D(*m_Iso,*m_pt,*m_eta,*m_mgg,*m_npv,*m_mu,"data");
  RooArgSet Variables_MC(*m_Iso,*m_pt,*m_eta,*m_mgg,*m_npv,*m_mu,*m_weight,"mc");
  Variables_D.add(*m_IsLoosePrime_L);
  Variables_D.add(*m_IsLoosePrime_SL);
  Variables_D.add(*m_IsTight_L);
  Variables_D.add(*m_IsTight_SL);
  Variables_MC.add(*m_IsLoosePrime_L);
  Variables_MC.add(*m_IsLoosePrime_SL);
  Variables_MC.add(*m_IsTight_L);
  Variables_MC.add(*m_IsTight_SL);

  if( m_data ) {
    m_DataSet = new RooDataSet( "DataSet","DataSet",
				Variables_D,
				RooFit::Import(*m_tree) );
  }else{ 
    m_DataSet = new RooDataSet( "DataSet","DataSet",
				Variables_MC,
				RooFit::Import(*m_tree),
				RooFit::WeightVar(*m_weight) );
  }
  std::cout << "DATASET SIZE = " << m_DataSet->sumEntries() << std::endl;
}
/////////////////////////////////////
void IsolationFitter::Init_DataSets()
/////////////////////////////////////
{
  TString mggcut = Form("mgg>=%f && mgg<%f",m_mggrange.first,m_mggrange.second);
  TString npvcut = Form("NPV>=%d && NPV<%d",m_npvrange.first,m_npvrange.second);
  TString mucut  = Form("mu>=%f && mu<%f"  ,m_murange.first ,m_murange.second);
  TString ptcut  = Form(m_pt_name+">=%f && "+m_pt_name+"<%f",m_ptrange.first,m_ptrange.second);
  TString etacut = Form("abs("+m_eta_name+") >=%f && abs("+m_eta_name+")<%f",m_etarange.first,m_etarange.second);
  TString lprimecut = "IsLoosePrime4_L==1 && IsLoosePrime4_SL==1";

  TString mastercut = mggcut ;
  mastercut += "&&";
  mastercut += npvcut;
  mastercut += "&&";
  mastercut += mucut;
  mastercut += "&&";
  mastercut += etacut;
  mastercut += "&&";
  mastercut += ptcut;
  mastercut += "&&";
  mastercut += lprimecut;

  //--------------------------------------------------------
  std::cout << "******* CUT APPLIED ********" << std::endl;
  std::cout << mggcut << std::endl;
  std::cout << npvcut << std::endl;
  std::cout << mucut << std::endl;
  std::cout << etacut << std::endl;
  std::cout << ptcut  << std::endl;
  std::cout << "****************************" << std::endl;
  //-------------------------------------------------------------
  m_DataSet_cut = 0;
  m_DataSet_Ti  = 0;
  m_DataSet_Lo  = 0;
  m_DataSet_cut = (RooDataSet*)m_DataSet->reduce(mastercut);
  std::cout << "REDUCED DATASET SIZE = " << m_DataSet_cut->sumEntries() << std::endl;
  // delete m_DataSet;
  // TString tightcut = m_tight_name+"==1"; 
  TString tightcut = "IsTight_L==1 && IsTight_SL==1";//Need to require tight for leading and subleading for comparison with gg sm mc.
  std::cout << tightcut << std::endl;
  m_DataSet_Ti  = (RooDataSet*)m_DataSet_cut->reduce(tightcut);
  std::cout << "TIGHT DATASET SIZE = " << m_DataSet_Ti->sumEntries() << std::endl;
  // TString antitightcut  = m_tight_name+"==0";
  TString antitightcut ;
  if( m_pt_name.Contains("_L") )
    antitightcut = "IsTight_L==0";
  else
    antitightcut = "IsTight_SL==0";

  std::cout << antitightcut << std::endl;
  m_DataSet_Lo  = (RooDataSet*)m_DataSet_cut->reduce(antitightcut);
  std::cout << "ANTI-TIGHT DATASET SIZE = " << m_DataSet_Lo->sumEntries() << std::endl;
  //------------------------------------------------------------------
}

////////////////////////////////////////
void IsolationFitter::JetFit(bool dofit)
////////////////////////////////////////
{
  int NbinsFit = 2*abs((int)m_isorange.second-(int)m_isorange.first);
  //  int NbinsFit = 10;

  RooRealVar* peak  = new RooRealVar("peak","#mu",0,"[GeV]");
  RooRealVar* width = new RooRealVar("width","#sigma",0,"[GeV]");
  RooRealVar* tail  = new RooRealVar("tail","tail",0);
  m_JetSet = new RooArgSet(*peak,*width,*tail,"JetSet");
  m_JetSet->readFromFile(m_FitParFile,"READ","Jet" );
  m_JetPdf = new RooNovosibirsk( "JetPdf","JetPdf",*m_Iso,*peak,*width,*tail );

  if(dofit){
    if(m_data)
      m_FitRes_Jet = m_JetPdf->fitTo( *m_DataSet_Lo,
				      RooFit::Range(m_isorange.first,m_isorange.second),
				      RooFit::Save(kTRUE) );
    else 
      m_FitRes_Jet = m_JetPdf->fitTo( *m_DataSet_Lo,
				      RooFit::Range(m_isorange.first,m_isorange.second),
				      RooFit::SumW2Error(kTRUE),
				      RooFit::Save(kTRUE) );

    m_FitRes_Jet->Print("v");
  }
  frame_Jet = m_Iso->frame(m_isorange.first,m_isorange.second,NbinsFit);
  frame_Jet->SetName("Jet");
  m_DataSet_Lo->plotOn( frame_Jet,
			RooFit::LineColor(3),
			RooFit::MarkerColor(3),
			RooFit::Name("data_jet"));
  m_JetPdf->plotOn( frame_Jet,
		    RooFit::LineColor(3),
		    RooFit::Name("pdf_jet") );
}
///////////////////////////////////////////
void IsolationFitter::PhotonFit(bool dofit)
///////////////////////////////////////////
{
  int NbinsFit = 2*abs((int)m_isorange.second-(int)m_isorange.first);
  //int NbinsFit = 10;
  TString cut =  Form( m_iso_name+">%f &&"+m_iso_name+"<%f",
		       m_isonormrange.first,
		       m_isonormrange.second );
  double num     = m_DataSet_Ti->sumEntries(cut);
  double denom   = m_DataSet_Lo->sumEntries(cut);
  std::cout << "num / denom = " << num << "/" << denom << std::endl;
  double scale   = num/denom;
  double Ntot_Lo = m_DataSet_Lo->sumEntries();
  RooRealVar N_Jet("N_Jet","N_Jet",Ntot_Lo*scale );
  RooRealVar N_PH("N_PH","N_PH",1000.);
  meanCB  = new RooRealVar("meanCB","meanCB",4);
  sigmaCB = new RooRealVar("sigmaCB","sigmaCB",4);
  alphaCB = new RooRealVar("alphaCB","alphaCB",4);
  nCB     = new RooRealVar("nCB","nCB",4);
  m_PHSet = new RooArgSet(*meanCB,*sigmaCB,*alphaCB,*nCB,"PHSet");
  RooArgSet TiSet(*m_PHSet,"TiSet");
  TiSet.add(N_PH);
  TiSet.readFromFile(m_FitParFile,"READ","Tight" );
  m_PHPdf= new RooCBShape( "PHPdf","PHPdf",*m_Iso,
			   *meanCB,*sigmaCB,*alphaCB,*nCB );
  
  if(m_data)
    m_TightPdf = new RooAddPdf( "TightPdf","TightPdf",
				RooArgList(*m_PHPdf,*m_JetPdf),
				RooArgList(N_PH,N_Jet) );

  std::cout << "Parameters values :" << N_PH.getVal() 
	    << " " << N_Jet.getVal() <<  " "
	    << meanCB->getVal() << " " << sigmaCB->getVal() << " "
	    << nCB->getVal() << " " << alphaCB->getVal() << std::endl;


  if(dofit){
    if(m_data) 
      m_FitRes_Ti = m_TightPdf->fitTo( *m_DataSet_Ti,
				       RooFit::Range(m_isorange.first,m_isorange.second),
				       RooFit::Extended(kTRUE),
				       RooFit::Save(kTRUE) );
    else
      m_FitRes_Ti = m_PHPdf->fitTo( *m_DataSet_Ti,
				    RooFit::Range(m_isorange.first,m_isorange.second),
				    RooFit::Save(kTRUE),
				    RooFit::SumW2Error(kTRUE) );
    m_FitRes_Ti->Print("v");
  }
  
  frame_PH = m_Iso->frame(m_isorange.first,m_isorange.second,NbinsFit);
  frame_PH->SetName("PH");
  m_DataSet_Ti->plotOn( frame_PH,
			RooFit::LineColor(1),
			RooFit::MarkerColor(1),
			RooFit::Name("data_tight") );
  if(m_data){
    m_TightPdf->plotOn( frame_PH,
			RooFit::LineColor(2),
			RooFit::Components(*m_PHPdf),
			RooFit::Name("pdf_ph") );
    m_TightPdf->plotOn( frame_PH,
			RooFit::LineColor(3),
			RooFit::Components(*m_JetPdf),
			RooFit::Name("pdf_jet") );
    m_TightPdf->plotOn( frame_PH,
			RooFit::LineColor(4),
			RooFit::Name("pdf_tot") );
  }else{  
    m_PHPdf->plotOn( frame_PH,
		     RooFit::LineColor(4),
		     RooFit::Name("pdf_tot") );
  }
  RooHist* RH = frame_PH->getHist("data_tight");
  m_h_iso = RooHistToHist(*RH);


}
////////////////////////////////////////
void IsolationFitter::Fitter(bool dofit)
///////////////////////////////////////
{
  if(m_data){
    JetFit(dofit);
    SetConstantParameters(*m_JetSet);
  }
  PhotonFit(dofit);
}

//////////////////////////////////////
TCanvas* IsolationFitter::GetJetPlot()
//////////////////////////////////////
{

  // TString st_iso = Form( "iso_%d-%d_",(int)m_isorange.first,(int)m_isorange.second);
  // TString st_pt  = Form( "pt_%d-%d_",(int)m_ptrange.first,(int)m_ptrange.second);
  // TString st_eta = Form( "eta_%1.1f-%1.1f_",m_etarange.first,m_etarange.second);
  // TString st_npv = Form( "npv_%d-%d_",(int)m_npvrange.first,(int)m_npvrange.second);
  TString name = "Plot_Jet_";
  // name += st_iso;
  // name += st_pt;
  // name += st_eta;
  // name += st_npv;
  // if( m_iso_name.Contains("_L")) name+= "leading";
  // else if( m_iso_name.Contains("_SL")) name += "subleading";
  // else Fatal("IsolationFitter::GetJetPlot()","Wrong m_iso_name !!");

  // TLatex * lat = new TLatex();
  // lat->SetNDC(true);
  // lat->SetTextSize(0.025);
  TCanvas* cJ = new TCanvas(name,name,800,600);
  // cJ->cd();
  // frame_Jet->Draw();
  // lat->DrawLatex(0.18,0.90,Form("p_{T} [GeV] #in [%d,%d]",(int)m_ptrange.first,(int)m_ptrange.second));
  // lat->DrawLatex(0.18,0.86,Form("|#eta| #in [%1.1f,%1.1f]",m_etarange.first,m_etarange.second));
  // lat->DrawLatex(0.18,0.82,Form("NPV #in [%d,%d]",(int)m_npvrange.first,(int)m_npvrange.second));
  // lat->DrawLatex(0.18,0.78,Form("#bf{#chi^{2}/ndf = %1.1f}",JetChiSquare()) );
  return cJ;
}
/////////////////////////////////////////
TCanvas* IsolationFitter::GetPhotonPlot()
/////////////////////////////////////////
{
  // TString st_iso = Form( "iso_%d-%d_",(int)m_isorange.first,(int)m_isorange.second);
  // TString st_pt  = Form( "pt_%d-%d_",(int)m_ptrange.first,(int)m_ptrange.second);
  // TString st_eta = Form( "eta_%1.1f-%1.1f_",m_etarange.first,m_etarange.second);
  // TString st_npv = Form( "npv_%d-%d_",(int)m_npvrange.first,(int)m_npvrange.second);
  TString name = "Plot_PH_";
  // name += st_iso;
  // name += st_pt;
  // name += st_eta;
  // name += st_npv;
  // if( m_iso_name.Contains("_L")) name+= "leading";
  // else if( m_iso_name.Contains("_SL")) name += "subleading";
  // else Fatal("IsolationFitter::GetPhotonPlot()","Wrong m_iso_name !!");


  // TLatex * lat = new TLatex();
  // lat->SetNDC(true);
  // lat->SetTextSize(0.025);
  TCanvas* cPH = new TCanvas(name,name,800,600);
  // cPH->cd();
  // frame_PH->Draw();
  // lat->DrawLatex(0.80,0.90,Form("p_{T} [GeV] #in [%d,%d]",(int)m_ptrange.first,(int)m_ptrange.second));
  // lat->DrawLatex(0.83,0.86,Form("|#eta| #in [%1.1f,%1.1f]",m_etarange.first,m_etarange.second));
  // lat->DrawLatex(0.83,0.82,Form("NPV #in [%d,%d]",(int)m_npvrange.first,(int)m_npvrange.second));
  // lat->DrawLatex(0.83,0.78,Form("#bf{#chi^{2}/ndf = %1.1f}",TotChiSquare()) );

  return cPH;

}
///////////////////////////////////////
double IsolationFitter::TotChiSquare()
//////////////////////////////////////
{

  double chi2    = frame_PH->chiSquare("pdf_tot","data_tight", 
				       m_FitRes_Ti->floatParsFinal().getSize()  );

  return chi2;
}
///////////////////////////////////////
double IsolationFitter::JetChiSquare()
//////////////////////////////////////
{

  double chi2    = frame_Jet->chiSquare("pdf_jet","data_jet", 
					m_FitRes_Jet->floatParsFinal().getSize()  );

  return chi2;
}


/////////////////////////////////////////////////////////////////
void IsolationFitter::SetConstantParameters(const RooArgSet& set)
/////////////////////////////////////////////////////////////////
{
  TIterator *iterat= set.createIterator();
  RooRealVar *next = 0;
  while( (0 != (next= (RooRealVar*)iterat->Next())) )
    next->setConstant();

  delete iterat;  
}

/////////////////////////////////////////////////
void IsolationFitter::StoreToRootFile(TString st)
/////////////////////////////////////////////////
{
  TFile fout(st,"RECREATE");
  fout.Add(m_h_iso);
  fout.Add(frame_PH);
  if(m_data)
    fout.Add(frame_Jet);
  fout.Write();
  RooWorkspace w("workspace");
    w.import(*m_PHPdf);
 
  if(m_data){
    w.import(*m_JetPdf);
  }//else
  w.Write();
  fout.Close();

}


///////////////////////////////////////////////////////
TH1F* IsolationFitter::RooHistToHist(const RooHist& RH)
///////////////////////////////////////////////////////
{
  double x=-999; double y= -999;

  std::vector<double> low_edge;
  std::vector<double> high_edge;

  for(int ip=0;ip<RH.GetN();ip++){
    int point = RH.GetPoint(ip,x,y);
    if( point == -1 ) Fatal("IsolationFitter::RooHistToHist()","Wrong input point 1 !");
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
    if( p == -1) Fatal("IsolationFitter::RooHistToHist()","Wrong input point !");
    h->SetBinContent(ibin,y);
  }

  //--> Set underflow and overflow to 0 
  //--> This information is lost in a TGraph.
  h->SetBinContent(0,0);
  h->SetBinContent(h->GetNbinsX()+1,0);

  return h;
}

///////////////////////////////////
IsolationFitter::~IsolationFitter()
///////////////////////////////////
{
  //destructor
  delete m_mgg;
  delete m_npv;
  delete m_mu;
  delete m_pt;
  delete m_eta;
  delete m_Iso;
  delete m_IsTight_L;
  delete m_IsTight_SL;
  delete m_weight;
  delete m_DataSet_cut;
  delete m_DataSet_Ti;
  delete m_DataSet_Lo;
  delete m_DataSet;
  delete m_JetSet;
  delete m_JetPdf;
  delete m_PHSet;
  delete m_TightPdf;
  delete meanCB;
  delete alphaCB;
  delete sigmaCB;
  delete nCB;
  delete m_FitRes_Jet;
  delete m_FitRes_Ti;
  delete frame_Jet;
  delete frame_PH;
}
