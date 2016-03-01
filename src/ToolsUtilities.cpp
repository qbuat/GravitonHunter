#define ANALYSISTOOLS_cpp
#include "ToolsUtilities.h"
#include <iostream>
#include <TMath.h>
#include <TPad.h>
#include <TError.h>
//////////////////////////////
AnalysisTools::AnalysisTools()
//////////////////////////////
{

}

///////////////////////////////////////////
AnalysisTools::~AnalysisTools()
///////////////////////////////////////////
{
  //Destructor
}

//////////////////////////////////////////////////////
double AnalysisTools::efficiency(double k,double n)
//////////////////////////////////////////////////////
{
  //--> Returns the efficiency computed following
  //--> http://arxiv.org/abs/0908.0130
  double temp=(k+1)/(n+2);
  return temp;
}

////////////////////////////////////////////////////////////
double AnalysisTools::sigma_efficiency(double k,double n)
////////////////////////////////////////////////////////////
{
  //--> Returns the erro on the efficiency computed following
  //--> http://arxiv.org/abs/0908.0130
  double r1=(k+1)/(n+2);
  double r2=(k+2)/(n+3);
  double r3=(k+1)/(n+2);
  double temp=sqrt(r1*r2-r3*r3);
  return temp;
}

/////////////////////////////////////////////////////
double AnalysisTools::RoundUpDouble(double x,int pre)
/////////////////////////////////////////////////////
{
  //--> Rounding method. Choose the precision with pre
  double temp=TMath::Ceil(pre*x)/(double)pre;
  return temp;
}
//////////////////////////////////////////////////////////////////////////////////
void AnalysisTools::ShowNumberOfProcessedEvents(int jentry,int nentries,int period)
//////////////////////////////////////////////////////////////////////////////////
{
  //------------------------------------------------------
  if(jentry%period==0 && jentry!=0){
    int r= (int)( (double)jentry/(double)nentries*100 );
    std::cout << "\E[31;1m" 
	      << jentry 
	      << " (" << r << ")%"
	      << " events processed so far\E[m" 
	      << std::endl;
  }
  //-------------------------------------------------------
}

////////////////////////////////////////////////////////////////////
TH1D* AnalysisTools::HistAddedUncert(const TH1D& h1,const TH1D& h2)
///////////////////////////////////////////////////////////////////
{
  //--> Add two histograms in quadrature and compute the uncertainty on the 
  //--> resulting histogram
  //--> Returns the quadratic summed histogram
  TH1D* htot = (TH1D*)h1.Clone(h1.GetName());
  for(int ibin=0;ibin<=htot->GetNbinsX()+1;ibin++){

    double bin1_square = h1.GetBinContent(ibin)*h1.GetBinContent(ibin);
    double bin2_square = h2.GetBinContent(ibin)*h2.GetBinContent(ibin);
    double error1_square = h1.GetBinError(ibin)*h1.GetBinError(ibin);
    double error2_square = h2.GetBinError(ibin)*h2.GetBinError(ibin);
    
    double bin_content = sqrt(bin1_square+bin2_square);
    double bin_error   = 0;
    if( bin_content!=0 )
      bin_error = sqrt( (bin1_square*error1_square+bin2_square*error2_square) / (bin1_square+bin2_square) );
    htot->SetBinContent(ibin, bin_content);
    htot->SetBinError(ibin,bin_error);
  }
  return htot;
}
/////////////////////////////////////////////////////////////////
TGraphAsymmErrors*  AnalysisTools::DataGraph(const TH1D& hd)
/////////////////////////////////////////////////////////////////
{
  //--> Transform an histogram into an asymmetric graph
  //--> The x errors are set to 0
  //--> The y errors are the 1sigma band of the poisson distribution
  //--> with parameter the bin content.


  TString name = "graph_";
  name += hd.GetName();
  TGraphAsymmErrors* graph = new TGraphAsymmErrors();
  graph->SetName(name);
  graph->SetTitle(name);
  graph->GetXaxis()->SetTitle( hd.GetXaxis()->GetTitle() );
  graph->GetYaxis()->SetTitle( hd.GetYaxis()->GetTitle() );

  for (int ibin = 1 ; ibin <= hd.GetNbinsX()+1;ibin++){
    double value = hd.GetBinContent(ibin);
    if( value!=0){
      double y1 = value + 1.0;
      double d1 = 1.0 - 1.0/(9.0*y1) + 1.0/(3*sqrt(y1));
      double error_poisson_up = y1*d1*d1*d1 - value;
      double y2 = value;
      double d2 = 1.0 - 1.0/(9.0*y2) - 1.0/(3.0*sqrt(y2));
      double error_poisson_down = value - y2*d2*d2*d2;
      graph->SetPoint(ibin-1, hd.GetBinCenter(ibin), value);
      graph->SetPointError(ibin-1, 0, 0, error_poisson_down, error_poisson_up);
    }      
  }
  return graph;

}
////////////////////////////////////////////////////////////////
TGraphAsymmErrors*  AnalysisTools::BkgGraph(const TH1D& hb)
////////////////////////////////////////////////////////////////
{
  //--> Transform a TH1D into an asymmetric graph 
  //--> The bin width is used as x error and the bin error is used as y error

  TString name = "graph_";
  name += hb.GetName();
  TGraphAsymmErrors* graph = new TGraphAsymmErrors();
  graph->SetName(name);
  graph->SetTitle(name);
  graph->GetXaxis()->SetTitle( hb.GetXaxis()->GetTitle() );
  graph->GetYaxis()->SetTitle( hb.GetYaxis()->GetTitle() );

  for (int ibin = 1 ; ibin <= hb.GetNbinsX()+1;ibin++){
    double value  = hb.GetBinContent(ibin);
    double erry   = hb.GetBinError(ibin);
    double errx   = hb.GetBinWidth(ibin)/2;
    if(value !=0){
      graph->SetPoint(ibin-1, hb.GetBinCenter(ibin), value);
      graph->SetPointError(ibin-1, errx, errx, erry, erry);
    }      
  }
  return graph;

}

//////////////////////////////////////////////////////////////////////////////////////
TH1D* AnalysisTools::GetTruncatedHist(const TH1D& hfull, int nbins,const double* bins)
//////////////////////////////////////////////////////////////////////////////////////
{
  //--> Truncate an histogram by returning only the bins given by
  //--> the user in the bins vector
  //--> Returns the truncated hist
  TString histname = hfull.GetName();
  histname += "_truncated";
  TH1D* htemp = new TH1D(histname,histname,nbins,bins);

  if( nbins > hfull.GetNbinsX() ) 
    Fatal("GetTruncatedHist()","Number of bins is too large !" ); 

  int bin_start = -99;
  for( int ibin=1;ibin<(hfull.GetNbinsX()+1);ibin++ ){
    if( hfull.GetBinLowEdge(ibin) == bins[0] ) 
      bin_start = ibin ;
  }      
  if ( bin_start == -99 ) 
    Fatal("GetTruncatedHist()","Didn't find the first matching bin" ); 
  for( int i=0;i<(nbins+1);i++){
    // std::cout <<  hfull.GetBinLowEdge(bin_start+i) << "=/=" << bins[i] << std::endl;
    // if( hfull.GetBinLowEdge(bin_start+i) != bins[i] ) 
    //   Fatal("GetTruncatedHist()","Error in the bin matching !" ); 
    htemp->SetBinContent(i+1,hfull.GetBinContent(bin_start+i));
  }
  return htemp;

}






//////////////////////////////////////////////////////////////////
void AnalysisTools::Processing(int jentry,int nentries,int period)
//////////////////////////////////////////////////////////////////
{
  //------------------------------------------------------
  if(jentry%period==0 && jentry!=0){
    int r= (int)( (double)jentry/(double)nentries*100 );
    std::cout << "\E[31;1m" 
	      << jentry 
	      << " events processed so far" 
	      << " ---> (" << r << ")%"
	      << "\E[m"
	      << std::endl;
  }
  //-------------------------------------------------------
}

/////////////////////////////////////////////////////////////////
TCanvas* AnalysisTools::Get2PadPlot(TString name, TString title)
////////////////////////////////////////////////////////////////
{
  TCanvas* canv = new TCanvas(name,title,800,800);
  canv->Divide(1,2);
  TPad *p1 = (TPad*)canv->GetPad(1);
  p1->SetPad(0.05,0.30,0.97,0.97);
  p1->SetBottomMargin(0.00);
  TPad *p2 = (TPad*)canv->GetPad(2);
  p2->SetPad(0.05,0.0,0.97,0.30);
  p2->SetTopMargin(0.0);
  p2->SetBottomMargin(0.40);
  p1->SetTicks(1,1);
  p2->SetTicks(1,1);
  canv->cd();
  p1->Update();
  canv->cd();
  p2->cd();
  // p2->SetGridx();
  p2->SetGridy();
  p2->Update();

  return canv;
}


/////////////////////////////////////////////////////
TH1D* AnalysisTools::GetPtggHist(TString name)
/////////////////////////////////////////////////////
{
  TH1D* h = new TH1D(name,name,100,0,800);
  h->GetXaxis()->SetTitle("p_{T,#gamma#gamma} [GeV]");
  h->GetYaxis()->SetTitle("Events/bin");
  return h;
}
/////////////////////////////////////////////////////
TH1D* AnalysisTools::GetAbsDeltaphiHist(TString name)
/////////////////////////////////////////////////////
{  
  const int Nbins = 18;
  double bins[] = 
    { 0, 0.5, 1.0, 1.5, 1.75, 2.0, 2.25, 2.35, 2.45, 2.55, 2.65,
      2.70, 2.75, 2.80, 2.85, 2.90, 2.95, 3.1, 3.1416 };
  TH1D* h = new TH1D(name,name,Nbins,bins);
  h->GetXaxis()->SetTitle("|#Delta#phi_{#gamma#gamma}| [rad]");
  h->GetYaxis()->SetTitle("Events/bin");
 return h;
}
//////////////////////////////////////////////////////////
TH1D* AnalysisTools::GetAbsCosthetastarHist(TString name)
/////////////////////////////////////////////////////////
{
  const int Nbins = 10;
  double bins[] = { 0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
  TH1D* h = new TH1D(name,name,Nbins,bins);
  h->GetXaxis()->SetTitle("|cos(#theta*)|");
  h->GetYaxis()->SetTitle("Events/bin");
  return h;
}
/////////////////////////////////////////////////////
TH1D* AnalysisTools::GetDeltaetaggHist(TString name)
/////////////////////////////////////////////////////
{
  TH1D* h = new TH1D(name,name,10,-5,5);
  h->GetXaxis()->SetTitle("#Delta#eta");
  h->GetYaxis()->SetTitle("Events/bin");
  return h;
}
/////////////////////////////////////////////////////
TH1D* AnalysisTools::GetYggHist(TString name)
/////////////////////////////////////////////////////
{
  TH1D* h = new TH1D(name,name,60,-3,3);
  h->GetXaxis()->SetTitle("Y_{#gamma#gamma}");
  h->GetYaxis()->SetTitle("Events/bin");
  return h;
}
/////////////////////////////////////////////////////
TH1D* AnalysisTools::GetDeltarHist(TString name)
/////////////////////////////////////////////////////
{
  TH1D* h = new TH1D(name,name,60,0,6);
  h->GetXaxis()->SetTitle("p_{T,#gamma#gamma} [GeV]");
  h->GetYaxis()->SetTitle("Events/bin");
  return h;
}
/////////////////////////////////////////////////////
TH1D* AnalysisTools::GetPtHist(TString name)
/////////////////////////////////////////////////////
{
  TH1D* h = new TH1D(name,name,100,0,800);
  h->GetXaxis()->SetTitle("p_{T,#gamma#gamma} [GeV]");
  h->GetYaxis()->SetTitle("Events/bin");
  return h;
}
/////////////////////////////////////////////////////
TH1D* AnalysisTools::GetEtaHist(TString name)
/////////////////////////////////////////////////////
{ 
  const int Nbins = 34;
  double bins[] = 
    {-2.4,-2.37,-2.27,-2.17,-2.07,-1.97,-1.87,-1.77,-1.67,-1.57,-1.52,-1.37,-1.27,-1.17,-1.07,-0.97,-0.6,0,
     0.6, 0.97, 1.07, 1.17, 1.27, 1.37, 1.52, 1.57, 1.67, 1.77, 1.87, 1.97, 2.07, 2.17, 2.27, 2.37, 2.4 };

  TH1D* h = new TH1D(name,name,Nbins,bins);
  h->GetXaxis()->SetTitle("#eta_{#gamma}");
  h->GetYaxis()->SetTitle("Events/bin");
  return h;
}
/////////////////////////////////////////////////////
TH1D* AnalysisTools::GetPhiHist(TString name)
/////////////////////////////////////////////////////
{  
  const int Nbins = 66;
  double bins[] = 
    { -3.2,-3.1416,-3.1,-3,-2.9,-2.8,-2.7,-2.6,-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,
      -1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,
      0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,
      2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.1416,3.2};

  TH1D* h = new TH1D(name,name,Nbins,bins);
  h->GetXaxis()->SetTitle("#eta_{#gamma}");
  h->GetYaxis()->SetTitle("Events/bin");
  return h;

}
/////////////////////////////////////////////////////
TH2D* AnalysisTools::GetEtaPhiMap(TString name)
/////////////////////////////////////////////////////
{  


  TH2D* h = new TH2D(name,name,480,-2.40,2.40,630,-3.15,3.15);
  h->GetXaxis()->SetTitle("#eta_{#gamma}");
  h->GetYaxis()->SetTitle("#phi_{#gamma}");
  return h;
}
