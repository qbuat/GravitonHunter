#include "IsolationGrapher.h"
#include "ToolsUtilities.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
IsolationGrapher::IsolationGrapher(IsolationFitter* IF, int Nbins, const double* bins,TString type)
///////////////////////////////////////////////////////////////////////////////////////////////////
{
  //Default constructor
  m_IF = IF;
  m_IF->SetPtBounds();
  m_IF->SetNPVBounds();
  // m_IF->SetEtaBounds();
  SetNbins(Nbins);
  SetType(type);
  for ( int i=0;i<m_Nbins;i++){
    double width = (bins[i+1]-bins[i])/2.;
    bin_width.push_back(width);
    bin_center.push_back(bins[i]+width);
  }

  m_canvas = new TObjArray();

}

///////////////////////////////////
IsolationGrapher::~IsolationGrapher()
///////////////////////////////////
{
  //destructor
  // delete m_IF;
  delete m_canvas;
}


////////////////////////////////////////////////////////
void IsolationGrapher::SetBounds(double min, double max)
///////////////////////////////////////////////////////
{
  if( m_type == "pT" )
    m_IF->SetPtBounds(min,max);
  else if( m_type == "NPV" )
    m_IF->SetNPVBounds(min,max);
  else if( m_type == "eta")
    m_IF->SetEtaBounds(min,max);
  else 
    Fatal( "IsolationGrapher::SetBounds", 
	   "Wrong boundary type (pT/NPV/eta)!!");
 
}

///////////////////////////////////
void IsolationGrapher::FillGraphs()
///////////////////////////////////
{
  for( int ibin = 0; ibin< m_Nbins;ibin++ ){
    SetBounds( bin_center[ibin]-bin_width[ibin],
	       bin_center[ibin]+bin_width[ibin] );
    m_IF->Init_DataSets();
    m_IF->Fitter(true);
    iso_fitmean.push_back(m_IF->GetPhotonFitMean());
    iso_fitmean_error.push_back(m_IF->GetPhotonFitMeanError());
    iso_fitwidth.push_back(m_IF->GetPhotonFitWidth());
    iso_fitwidth_error.push_back(m_IF->GetPhotonFitWidthError());
    if( m_IF->GetStreamType()==true)
      m_canvas->Add( m_IF->GetJetPlot() );
    m_canvas->Add( m_IF->GetPhotonPlot() );
    TH1F* h_iso = m_IF->GetPhotonHist();
    iso_mean.push_back( h_iso->GetMean() );
    iso_mean_error.push_back( h_iso->GetMeanError() );
    iso_rms.push_back( h_iso->GetRMS() );
    iso_rms_error.push_back( h_iso->GetRMSError() );
  }
}

//////////////////////////////////////////////////
TGraphErrors* IsolationGrapher::GetFitMeanGraph()
/////////////////////////////////////////////////
{
  TGraphErrors* gr = new TGraphErrors( m_Nbins,
				       &bin_center[0],
				       &iso_fitmean[0],
				       &bin_width[0],
				       &iso_fitmean_error[0] );

  gr->GetYaxis()->SetTitle("meanCB [GeV]");
  TString xtitle;
  if( m_type == "pT")
    xtitle = "p_{T} [GeV]";
  else if( m_type == "NPV" )
    xtitle = "Number of Primary Vertices";
  else if( m_type == "eta" )
    xtitle = "#eta";
  gr->GetXaxis()->SetTitle(xtitle);
  return gr;
}
//////////////////////////////////////////////
TGraphErrors* IsolationGrapher::GetFitWidthGraph()
//////////////////////////////////////////////
{
  TGraphErrors* gr = new TGraphErrors( m_Nbins,
				       &bin_center[0],
				       &iso_fitwidth[0],
				       &bin_width[0],
				       &iso_fitwidth_error[0] );

  gr->GetYaxis()->SetTitle("widthCB [GeV]");
  TString xtitle;
  if( m_type == "pT")
    xtitle = "p_{T} [GeV]";
  else if( m_type == "NPV" )
    xtitle = "Number of Primary Vertices";
  else if( m_type == "eta" )
    xtitle = "#eta";
  gr->GetXaxis()->SetTitle(xtitle);
  return gr;
}

//////////////////////////////////////////////////
TGraphErrors* IsolationGrapher::GetHistMeanGraph()
/////////////////////////////////////////////////
{
  TGraphErrors* gr = new TGraphErrors( m_Nbins,
				       &bin_center[0],
				       &iso_mean[0],
				       &bin_width[0],
				       &iso_mean_error[0] );

  gr->GetYaxis()->SetTitle("mean [GeV]");
  TString xtitle;
  if( m_type == "pT")
    xtitle = "p_{T} [GeV]";
  else if( m_type == "NPV" )
    xtitle = "Number of Primary Vertices";
  else if( m_type == "eta" )
    xtitle = "#eta";
  gr->GetXaxis()->SetTitle(xtitle);
  return gr;
}
//////////////////////////////////////////////
TGraphErrors* IsolationGrapher::GetHistRMSGraph()
//////////////////////////////////////////////
{
  TGraphErrors* gr = new TGraphErrors( m_Nbins,
				       &bin_center[0],
				       &iso_rms[0],
				       &bin_width[0],
				       &iso_rms_error[0] );

  gr->GetYaxis()->SetTitle("RMS [GeV]");
  TString xtitle;
  if( m_type == "pT")
    xtitle = "p_{T} [GeV]";
  else if( m_type == "NPV" )
    xtitle = "Number of Primary Vertices";
  else if( m_type == "eta" )
    xtitle = "#eta";
  gr->GetXaxis()->SetTitle(xtitle);
  return gr;
}

