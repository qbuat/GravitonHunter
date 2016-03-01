#include "PlottingEfficiencies.h"

#include "ToolsUtilities.h"
#include "ToolsSignificanceHist.h"
#include "ToolsCommons.h"

///////////////////////////////////////////////////////////
PlottingEfficiencies::PlottingEfficiencies(TString filemaps)
///////////////////////////////////////////////////////////
{
  InitMaps(filemaps);
  m_lat = new TLatex();

}
//////////////////////////////////////////////
PlottingEfficiencies::~PlottingEfficiencies()
/////////////////////////////////////////////
{
  delete m_filemaps;
  delete m_lat;
}


/////////////////////////////////////////////////////
void PlottingEfficiencies::InitMaps(TString filemaps)
/////////////////////////////////////////////////////
{
  if(filemaps=="NONE") m_filemaps = 0;
  else  m_filemaps = new TFile(filemaps,"read");
}
///////////////////////////////////////////////////////////////////
void PlottingEfficiencies::SetPtBinning(int Nbins, double *edges)
//////////////////////////////////////////////////////////////////
{
  for(int i=0;i<Nbins;i++)
    m_pt_binedges.push_back( edges[i] );
}

//////////////////////////////////////////////////////////////////////////////////////////////
std::vector<TGraphErrors*> PlottingEfficiencies::BuildPtEfficiencyGraphs(const TH2 & m_den,
									 const TH2 & m_num )
//////////////////////////////////////////////////////////////////////////////////////////////
{
  AnalysisTools AT;
  std::vector<TGraphErrors*> vec_out;
  for(int ietabin = 0;ietabin<m_den.GetNbinsY();ietabin++){

    std::vector<double> pt_values;
    std::vector<double> pt_errors;
    std::vector<double> eff_values;
    std::vector<double> eff_errors;

    for(int ipt=0;ipt<(int)m_pt_binedges.size()-1;ipt++){
      std::cout << "[ "<< m_pt_binedges[ipt] << ", "<< m_pt_binedges[ipt+1] << "]" << std::endl;
      double err_pt = (m_pt_binedges[ipt+1]-m_pt_binedges[ipt])/2;
      pt_values.push_back( m_pt_binedges[ipt]+ err_pt);
      pt_errors.push_back( err_pt );
      double num_y =0; double den_y=0;
      for(int iptbin = 0;iptbin<m_den.GetNbinsX();iptbin++){
	if( m_den.GetXaxis()->GetBinLowEdge(iptbin)< m_pt_binedges[ipt] ) continue;
	if( m_den.GetXaxis()->GetBinLowEdge(iptbin)>= m_pt_binedges[ipt+1] ) continue;
	num_y += m_num.GetBinContent(iptbin,ietabin);
	den_y += m_den.GetBinContent(iptbin,ietabin);
      }
      std::cout << num_y << " / "<< den_y << " = " << 100*AT.efficiency(num_y,den_y) << std::endl;
      eff_values.push_back( 100*AT.efficiency(num_y,den_y) ) ;
      eff_errors.push_back( 100*AT.sigma_efficiency(num_y,den_y) );
    }
    TGraphErrors* gr = new TGraphErrors((int)pt_values.size(),
					&pt_values[0],&eff_values[0],
					&pt_errors[0],&eff_errors[0]);
    gr->SetName( Form("graph%d",ietabin) );
    gr->SetTitle( Form("graph%d",ietabin) );
    vec_out.push_back(gr);
  }
  return vec_out;
}
//////////////////////////////////////////////////////////////////
TCanvas* PlottingEfficiencies::RecoTightEfficiencyPlot(int etabin)
//////////////////////////////////////////////////////////////////
{
  //  TH2F* m_num = (TH2F*)m_filemaps->Get("map_lead_recotightiso");
  TH2F* m_num = (TH2F*)m_filemaps->Get("map_lead_reco");
  TH2F* m_den = (TH2F*)m_filemaps->Get("map_lead_truth");
  std::vector<TGraphErrors*> grs = BuildPtEfficiencyGraphs(*m_den,*m_num);
  TCanvas * c = new TCanvas();
  grs[etabin]->GetYaxis()->SetTitle("Efficiency (%)");
  grs[etabin]->GetXaxis()->SetTitle("p_{T}^{truth}[GeV]");
  grs[etabin]->GetYaxis()->SetRangeUser(0,110);
  grs[etabin]->Draw("AP");

  TFile * savefile;
  TString samplename = (TString(m_filemaps->GetName()).Contains("RSGG")) ? "RSGG" : "SMGG";
  TString isoname = (TString(m_filemaps->GetName()).Contains("TOPO_PTDEP")) ? "TOPO_PTDEP" : "CONE";
  
  TString fullname = "ph_eff_"+TString(samplename)+"_"+TString(isoname)+"_etabin";
  savefile = new TFile(Form("%s%i.root",fullname.Data(),etabin),"RECREATE");

  grs[etabin]->Write();
  c->Write();

  return c;

}
