#include "ToolsExtendedCanvas.h"
#include <iostream>


ExtendedCanvas::ExtendedCanvas(const char *name, const char *title,
			       Int_t ww, Int_t wh,int npads) : TCanvas(name,title,ww,wh)
{
  m_haslumilabel = false;
  m_hascdmlabel = false;
  m_hasatlaslabel = false;
  m_lat = new TLatex();
  m_lat->SetTextSize(0.042);
  m_lat->SetNDC(true);
  if(npads !=1)
    SetMultiplePads(npads);
}

ExtendedCanvas::ExtendedCanvas(Bool_t build,int npads) : TCanvas(build)
{
  m_haslumilabel = false;
  m_hascdmlabel = false;
  m_hasatlaslabel = false;
  m_lat = new TLatex();
  m_lat->SetNDC(true);
  m_lat->SetTextSize(0.042);
  if(npads !=1)
    SetMultiplePads(npads);
}
ExtendedCanvas::~ExtendedCanvas()
{
  //Default destructor
  delete m_lat;
}


void ExtendedCanvas::SetLumiLabel(double X, double Y,double lumi)
{
  m_haslumilabel = true;
  m_lat->DrawLatex(X,Y,Form("#int Ldt = %1.1f fb^{-1}",lumi) );
}
void ExtendedCanvas::SetCdmLabel(double X, double Y,TString status)
{
  m_hascdmlabel = true;
  if (status == "2010" || status == "2011" )
    m_lat->DrawLatex(X,Y,"#sqrt{s} = 7 TeV");
  else if( status == "2012" )
    m_lat->DrawLatex(X,Y,"#sqrt{s} = 8 TeV");
  else Fatal("ExtendedCanvas::SetCdmLabel()","Wrong status !");
}
void ExtendedCanvas::SetAtlasLabel(double X, double Y,TString status)
{
  m_hasatlaslabel = true;
  if (status == "internal")
    m_lat->DrawLatex(X,Y,"#font[72]{ATLAS Internal}");
  else if( status == "progress" )
    m_lat->DrawLatex(X,Y,"#font[72]{ATLAS Work In Progress}");
  else if( status == "approval")
    m_lat->DrawLatex(X,Y,"#font[72]{ATLAS For Approval}");
  else if( status == "conference" )
    m_lat->DrawLatex(X,Y,"#font[72]{ATLAS Preliminary}");
  else if( status == "paper" )
    m_lat->DrawLatex(X,Y,"#font[72]{ATLAS}");
  else Fatal("ExtendedCanvas::SetAtlasLabel()","Wrong status !");
}

void ExtendedCanvas::SetGenericLabel(double X, double Y,TString label)
{  m_lat->DrawLatex(X,Y,label); }

void ExtendedCanvas::SetEtaCategoryLabel(double X, double Y,TString cat)
{  
  m_lat->SetTextSize(0.037);
  if(cat=="CC") m_lat->DrawLatex(X,Y,"#it{Central-Central}");
  else if(cat=="CE") m_lat->DrawLatex(X,Y,"#it{Central-Endcap}");
  else if(cat=="EC") m_lat->DrawLatex(X,Y,"#it{Endcap-Central}");
  else if(cat=="EE") m_lat->DrawLatex(X,Y,"#it{Endcap-Endcap}");
  else if(cat=="EE_S") m_lat->DrawLatex(X,Y,"#it{Endcap-Endcap (same sign)}");
  else if(cat=="EE_O") m_lat->DrawLatex(X,Y,"#it{Endcap-Endcap (opposite sign)}");
  else if(cat=="NONE") m_lat->DrawLatex(X,Y,"#it{No categorization}");
  else if(cat=="allcat") m_lat->DrawLatex(X,Y,"All categories");
  else Fatal("ExtendedCanvas::SetEtaCategoryLabel","Wrong Eta category !!");
  m_lat->SetTextSize(0.042);

 }

void ExtendedCanvas::SetMultiplePads(int npads)
{
  if( npads !=2) Fatal("ExtendedCanvas::SetMultiplePads()","Number of pads not supported" );
  m_lat->SetTextSize(0.062);
  this->Divide(1,2);
  TPad *p1 = (TPad*)this->GetPad(1);
  p1->SetPad(0.05,0.30,0.97,0.97);
  p1->SetBottomMargin(0.00);
  TPad *p2 = (TPad*)this->GetPad(2);
  p2->SetPad(0.05,0.0,0.97,0.30);
  p2->SetTopMargin(0.0);
  p2->SetBottomMargin(0.40);
  p1->SetTicks(1,1);
  p2->SetTicks(1,1);
  this->cd();
  p1->Update();
  this->cd();
  p2->cd();
  // p2->SetGridx();
  p2->SetGridy();
  p2->Update();

}

