#ifndef ExtendedCanvas_h
#define ExtendedCanvas_h

#include <TCanvas.h>
#include <TLatex.h>
#include <TString.h>

class ExtendedCanvas : public TCanvas{

  private:

    bool m_haslumilabel;
    bool m_hascdmlabel;
    bool m_hasatlaslabel;
    TLatex* m_lat;

    void SetMultiplePads(int npads = 2);

  public:

    //--> Overloaded only those two constructors
    // from TCanvas constructors
    ExtendedCanvas(Bool_t build=kTRUE,int npads =1);
    ExtendedCanvas(const char *name, const char *title, 
		   Int_t ww, Int_t wh,int npads=1);
    virtual ~ExtendedCanvas();

    void SetLumiLabel(double X, double Y,double lumi);
    void SetCdmLabel(double X,double Y,TString status= "2012");
    void SetAtlasLabel(double X,double Y,TString status="progress");
    void SetGenericLabel(double X,double Y,TString label);
    void SetEtaCategoryLabel(double X,double Y,TString cat);
    ClassDef(ExtendedCanvas,0);


};

#endif
