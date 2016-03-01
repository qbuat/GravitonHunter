import uuid
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)



if __name__ == '__main__':

    NTUPLES = '/sps/atlas/e/escalier/ExchangeInfo/Publication/Graviton/Graviton_0.0.2/TreesData_2012/Total.root'

    ch = ROOT.TChain('tree_selected')
    ch.Add(NTUPLES)


    h2 = ROOT.TH1F("h2", "h2", 20, -0.1, 0.1)
    h2.GetXaxis().SetTitle("#Delta(m_{#gamma#gamma})/m_{#gamma#gamma}")
    h2.SetLineWidth(2)
    ch.Draw("(invariant_mass-invariant_mass_HighestSumPt2)/invariant_mass>>h2", "invariant_mass > 200 && invariant_mass < 300")
    h2.Scale(1. / h2.Integral())
    h2.SetLineColor(ROOT.kBlack)
    h3 = ROOT.TH1F("h3", "h3", 20, -0.1, 0.1)
    h3.GetXaxis().SetTitle("m_{#gamma#gamma} (MVA vertex) [GeV]")
    h3.GetYaxis().SetTitle("#Delta(m_{#gamma#gamma}) [GeV]")
    h3.SetLineColor(ROOT.kRed)
    h3.SetLineWidth(2)
    ch.Draw("(invariant_mass-invariant_mass_HighestSumPt2)/invariant_mass>>h3", "invariant_mass > 300 && invariant_mass < 500")
    h3.Scale(1. / h3.Integral())


    h4 = ROOT.TH1F("h4", "h4", 20, -0.1, 0.1)
    h4.GetXaxis().SetTitle("m_{#gamma#gamma} (MVA vertex) [GeV]")
    h4.GetYaxis().SetTitle("#Delta(m_{#gamma#gamma}) [GeV]")
    h4.SetLineColor(ROOT.kBlue)
    h4.SetLineWidth(2)
    ch.Draw("(invariant_mass-invariant_mass_HighestSumPt2)/invariant_mass>>h4", "invariant_mass > 500 && invariant_mass < 1000")
    h4.Scale(1. / h4.Integral())

    canv = ROOT.TCanvas()
    canv.SetLogy(True)
    h2.Draw("HIST")
    h3.Draw("SAMEHIST")
    h4.Draw("SAMEHIST")
    h2.Draw("SAMEHIST")


    leg = ROOT.TLegend(
        0.7, 0.7,
        1 - canv.GetRightMargin(),
        1 - canv.GetTopMargin())
    leg.SetTextSize(0.9*leg.GetTextSize())
    leg.SetFillColor(0)
    leg.AddEntry(h2, '200 < m_{#gamma#gamma} < 300', 'f')
    leg.AddEntry(h3, '300 < m_{#gamma#gamma} < 500', 'f')
    leg.AddEntry(h4, '500 < m_{#gamma#gamma} < 1000', 'f')
    leg.Draw('same')
    canv.SaveAs('./plots/dummy_mgg_res.pdf')
