import uuid
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)




def ratio_hist(tree, cut_den, cut_nom, name=None):
    """
    Make a ratio between two selection, return the histogram
    """
    h_den = ROOT.TH1F(uuid.uuid1().hex, 'mgg', 50, 0, 1000)
    tree.Draw('mgg>>{0}'.format(h_den.GetName()), cut_den)
    h_den.Sumw2()

    h_nom = ROOT.TH1F(uuid.uuid1().hex, 'mgg', 50, 0, 1000)
    tree.Draw('mgg>>{0}'.format(h_nom.GetName()), cut_nom)
    h_nom.Sumw2()

    h_ratio = ROOT.TH1F(uuid.uuid1().hex, 'mgg', 50, 0, 1000)
    h_ratio.Divide(h_nom, h_den)
    if name is not None:
        h_ratio.SetName(name)

    return h_ratio

def ratio_plot(
    r1, r2, xtitle='X', 
    ytitle='Cut Efficiency',
    yrange=(0, 1.1),
    saving_name='./plots/dummy.pdf'):
    """
    Make the fancy (lol) plot
    """

    canv = ROOT.TCanvas()
    canv.SetName(uuid.uuid1().hex)
    r1.GetXaxis().SetTitle(xtitle)
    r1.GetYaxis().SetTitle(ytitle)
    r1.GetYaxis().SetRangeUser(yrange[0], yrange[1])
    r2.SetLineColor(ROOT.kRed)
    r1.SetLineWidth(2)
    r2.SetLineWidth(2)
    r1.Draw('HISTE')
    r2.Draw('sameHISTE')
    leg = ROOT.TLegend(
        0.7, 0.7,
        1 - canv.GetRightMargin(),
        1 - canv.GetTopMargin())
    leg.SetTextSize(0.9*leg.GetTextSize())
    leg.SetFillColor(0)
    leg.AddEntry(r1, 'published', 'f')
    leg.AddEntry(r2, 'latest', 'f')
    leg.Draw('same')
    canv.SaveAs(saving_name)
    # return canv

def compare_shape(
    h1, h2, xtitle='X',
    normalize=True,
    logy=False,
    saving_name='./plots/dummy.pdf'):
    h1.Sumw2()
    h2.Sumw2()

    h2.SetLineColor(ROOT.kRed)
    if normalize:
        h1.Scale(1. / h1.Integral())
        h2.Scale(1. / h2.Integral())

    h1.GetXaxis().SetTitle(xtitle)
    h2.GetXaxis().SetTitle(xtitle)
    h1.GetYaxis().SetTitle('Arbitrary Units')
    h2.GetYaxis().SetTitle('Arbitrary Units')

    # h1.GetYaxis().SetRangeUser(0.1, 1.1 * max(h1.GetBinContent(h1.GetMaximumBin()), h2.GetBinContent(h2.GetMaximumBin())))
    h1.SetLineWidth(2)
    h2.SetLineWidth(2)

    canv = ROOT.TCanvas()
    canv.SetName(uuid.uuid1().hex)
    h1.Draw('HIST')
    h2.Draw('sameHIST')
    leg = ROOT.TLegend(
        0.7, 0.7,
        1 - canv.GetRightMargin(),
        1 - canv.GetTopMargin())
    leg.SetTextSize(0.9*leg.GetTextSize())
    leg.SetFillColor(0)
    leg.AddEntry(h1, 'published', 'f')
    leg.AddEntry(h2, 'latest', 'f')
    leg.Draw('same')
    if logy:
        canv.SetLogy(True)
    canv.SaveAs(saving_name)

if __name__ == '__main__':

    RUN1_NTUPLES = '/sps/atlas/q/qbuat/DIPHOTON_2012_REF/fuerSarahEtAl/Look/datafile_TOPO_PTDEP_run*.root'
    # RUN1_NTUPLES = '/sps/atlas/q/qbuat/BUMPY_BUMP_OLD/datafile_TOPO_PTDEP_run*.root'
    # RUN1_NTUPLES = '/sps/atlas/q/qbuat/BUMPY_BUMP_OLD_NOCALIBEFORISO/datafile_TOPO_PTDEP_run*.root'
    RUN2_NTUPLES = '/sps/atlas/q/qbuat/BUMPY_BUMP/datafile_TOPO_PTDEP_run*.root'


    cut_pt = ROOT.TCut("pT_L >= 50 && pT_SL >= 50")
    cut_id = ROOT.TCut("IsTight_L == 1 && IsTight_SL == 1")
    cut_iso = ROOT.TCut("Iso_L_mod <= 8 && Iso_SL_mod <= 8")
    cut_mass = ROOT.TCut("mgg >= 150")
    # cut_mass = ROOT.TCut("mgg >= 179.1034055")

    cut_noid = cut_pt + cut_mass
    cut_noiso = cut_pt + cut_id + cut_mass
    cut_all = cut_pt + cut_id + cut_iso + cut_mass
    cut_all.Print()
    ch_r1 = ROOT.TChain('tree')
    ch_r1.Add(RUN1_NTUPLES)

    ch_r1_geo21 = ROOT.TChain('tree')
    ch_r1_geo21.Add(RUN2_NTUPLES)
    
    ratio_r1_id = ratio_hist(ch_r1, cut_noid, cut_noiso, name='run1')
    ratio_r1_geo21_id = ratio_hist(ch_r1_geo21, cut_noid, cut_noiso, name='run1_geo21')


    
    canv = ratio_plot(ratio_r1_id, ratio_r1_geo21_id, xtitle='m_{#gamma#gamma} [GeV]', yrange=(0, 0.2), saving_name='./plots/dummy_id.pdf')


    ratio_r1 = ratio_hist(ch_r1, cut_noiso, cut_all, name='run1')
    ratio_r1_geo21 = ratio_hist(ch_r1_geo21, cut_noiso, cut_all, name='run1_geo21')
 
    canv = ratio_plot(ratio_r1, ratio_r1_geo21, xtitle='m_{#gamma#gamma} [GeV]', saving_name='./plots/dummy_iso.pdf')



    # PT LEADING
    h1 = ROOT.TH1F("h_pt_l_1", "pt_l", 50, 50, 500)
    h2 = ROOT.TH1F("h_pt_l_2", "pt_l", 50, 50, 500)
    ch_r1.Draw('pT_L>>h_pt_l_1', cut_all)
    ch_r1_geo21.Draw('pT_L>>h_pt_l_2', cut_all)
    canv = compare_shape(h1, h2, xtitle='Leading Photon p_{T} [GeV]', saving_name='./plots/dummy_pt_l_shape.pdf')

    # PT SUBLEADING
    h1 = ROOT.TH1F("h_pt_sl_1", "pt_sl", 50, 50, 500)
    h2 = ROOT.TH1F("h_pt_sl_2", "pt_sl", 50, 50, 500)
    ch_r1.Draw('pT_SL>>h_pt_sl_1', cut_all)
    ch_r1_geo21.Draw('pT_SL>>h_pt_sl_2', cut_all)
    canv = compare_shape(h1, h2, xtitle='Subleading Photon p_{T} [GeV]', saving_name='./plots/dummy_pt_sl_shape.pdf')

    # CONV LEADING
    h1 = ROOT.TH1F("h_conv_l_1", "conv_l", 3, 0, 3)
    h2 = ROOT.TH1F("h_conv_l_2", "conv_l", 3, 0, 3)
    ch_r1.Draw('IsConv_L>>h_conv_l_1', cut_noiso)
    ch_r1_geo21.Draw('IsConv_L>>h_conv_l_2', cut_noiso)
    compare_shape(h1, h2, xtitle='Leading Photon Conversion Status', saving_name='./plots/dummy_conv_l_shape.pdf')

    # CONV SUBLEADING
    h1 = ROOT.TH1F("h_conv_sl_1", "conv_sl", 3, 0, 3)
    h2 = ROOT.TH1F("h_conv_sl_2", "conv_sl", 3, 0, 3)
    ch_r1.Draw('IsConv_SL>>h_conv_sl_1', cut_noiso)
    ch_r1_geo21.Draw('IsConv_SL>>h_conv_sl_2', cut_noiso)
    compare_shape(h1, h2, xtitle='Subleading Photon Conversion Status', saving_name='./plots/dummy_conv_sl_shape.pdf')

    # # ISOL LEADING
    # h1 = ROOT.TH1F("h_pt_corr_l_1", "isol_l", 50, -5, 20)
    # h2 = ROOT.TH1F("h_pt_corr_l_2", "isol_l", 50, -5, 20)
    # ch_r1.Draw('Iso_ptcorr_L>>h_pt_corr_l_1', cut_noiso)
    # ch_r1_geo21.Draw('Iso_ptcorr_L>>h_pt_corr_l_2', cut_noiso)
    # print h1.Integral(), h2.Integral()
    # canv = compare_shape(h1, h2, xtitle='Leading Photon Isolation [GeV] (pt leakage correction)', saving_name='./plots/dummy_pt_corr_l_shape.pdf')

    # # ISOL SUBLEADING
    # h1 = ROOT.TH1F("h_pt_corr_sl_1", "isol_sl", 50, -5, 20)
    # h2 = ROOT.TH1F("h_pt_corr_sl_2", "isol_sl", 50, -5, 20)
    # ch_r1.Draw('Iso_ptcorr_SL>>h_pt_corr_sl_1', cut_noiso)
    # ch_r1_geo21.Draw('Iso_ptcorr_SL>>h_pt_corr_sl_2', cut_noiso)
    # print h1.Integral(), h2.Integral()
    # canv = compare_shape(h1, h2, xtitle='Subleading Photon Isolation [GeV] (pt leakage correction)', saving_name='./plots/dummy_pt_corr_sl_shape.pdf')

    # # ISOL LEADING
    # h1 = ROOT.TH1F("h_ed_med_l_1", "isol_l", 50, -5, 20)
    # h2 = ROOT.TH1F("h_ed_med_l_2", "isol_l", 50, -5, 20)
    # ch_r1.Draw('ED_med_L>>h_ed_med_l_1', cut_noiso)
    # ch_r1_geo21.Draw('ED_med_L>>h_ed_med_l_2', cut_noiso)
    # print h1.Integral(), h2.Integral()
    # canv = compare_shape(h1, h2, xtitle='Leading Photon Isolation [GeV] (energy density)', saving_name='./plots/dummy_ed_med_l_shape.pdf')

    # # ISOL SUBLEADING
    # h1 = ROOT.TH1F("h_ed_med_sl_1", "isol_sl", 50, -5, 20)
    # h2 = ROOT.TH1F("h_ed_med_sl_2", "isol_sl", 50, -5, 20)
    # ch_r1.Draw('ED_med_SL>>h_ed_med_sl_1', cut_noiso)
    # ch_r1_geo21.Draw('ED_med_SL>>h_ed_med_sl_2', cut_noiso)
    # print h1.Integral(), h2.Integral()
    # canv = compare_shape(h1, h2, xtitle='Subleading Photon Isolation [GeV] (energy density)', saving_name='./plots/dummy_ed_med_sl_shape.pdf')

    # # ISOL LEADING
    # h1 = ROOT.TH1F("h_isol_uncor_l_1", "isol_l", 50, -5, 20)
    # h2 = ROOT.TH1F("h_isol_uncor_l_2", "isol_l", 50, -5, 20)
    # ch_r1.Draw('Iso_uncor_L>>h_isol_uncor_l_1', cut_noiso)
    # ch_r1_geo21.Draw('Iso_uncor_L>>h_isol_uncor_l_2', cut_noiso)
    # print h1.Integral(), h2.Integral()
    # canv = compare_shape(h1, h2, xtitle='Leading Photon Isolation [GeV] (uncorrected)', saving_name='./plots/dummy_iso_uncor_l_shape.pdf')

    # # ISOL SUBLEADING
    # h1 = ROOT.TH1F("h_isol_uncor_sl_1", "isol_sl", 50, -5, 20)
    # h2 = ROOT.TH1F("h_isol_uncor_sl_2", "isol_sl", 50, -5, 20)
    # ch_r1.Draw('Iso_uncor_SL>>h_isol_uncor_sl_1', cut_noiso)
    # ch_r1_geo21.Draw('Iso_uncor_SL>>h_isol_uncor_sl_2', cut_noiso)
    # print h1.Integral(), h2.Integral()
    # canv = compare_shape(h1, h2, xtitle='Leading Photon Isolation [GeV] (uncorrected)', saving_name='./plots/dummy_iso_uncor_sl_shape.pdf')


    # ISOL LEADING
    h1 = ROOT.TH1F("h_isol_l_1", "isol_l", 50, -5, 20)
    h2 = ROOT.TH1F("h_isol_l_2", "isol_l", 50, -5, 20)
    # Quentin fix this MeV -> GeV conversion, please
    ch_r1.Draw('Iso_L>>h_isol_l_1', cut_noiso)
    ch_r1_geo21.Draw('Iso_L>>h_isol_l_2', cut_noiso)
    print h1.Integral(), h2.Integral()
    canv = compare_shape(h1, h2, xtitle='Leading Photon Isolation [GeV] (corrected)', saving_name='./plots/dummy_iso_corr_l_shape.pdf')


    # ISO SUBLEADING
    h1 = ROOT.TH1F("h_isol_sl_1", "isol_sl", 50, -5, 20)
    h2 = ROOT.TH1F("h_isol_sl_2", "isol_sl", 50, -5, 20)
    ch_r1.Draw('Iso_SL>>h_isol_sl_1', cut_noiso)
    ch_r1_geo21.Draw('Iso_SL>>h_isol_sl_2', cut_noiso)
    print h1.Integral(), h2.Integral()
    print ch_r1.GetEntries(cut_all.GetTitle())/ float(ch_r1.GetEntries(cut_noiso.GetTitle())) * 100
    print ch_r1_geo21.GetEntries(cut_all.GetTitle()) / float(ch_r1_geo21.GetEntries(cut_noiso.GetTitle())) * 100

    canv = compare_shape(h1, h2, xtitle='Subleading Photon Isolation [GeV] (corrected)', saving_name='./plots/dummy_iso_corr_sl_shape.pdf')

    # ISOL LEADING
    h1 = ROOT.TH1F("h_isol_l_1", "isol_l", 50, -5, 20)
    h2 = ROOT.TH1F("h_isol_l_2", "isol_l", 50, -5, 20)
    # Quentin fix this MeV -> GeV conversion, please
    ch_r1.Draw('Iso_L_mod>>h_isol_l_1', cut_noiso)
    ch_r1_geo21.Draw('Iso_L_mod>>h_isol_l_2', cut_noiso)
    print h1.Integral(), h2.Integral()
    canv = compare_shape(h1, h2, xtitle='Leading Photon Isolation [GeV] (run1 graviton style)', saving_name='./plots/dummy_iso_l_shape.pdf')

    # ISO SUBLEADING
    h1 = ROOT.TH1F("h_isol_sl_1", "isol_sl", 50, -5, 20)
    h2 = ROOT.TH1F("h_isol_sl_2", "isol_sl", 50, -5, 20)
    ch_r1.Draw('Iso_SL_mod>>h_isol_sl_1', cut_noiso)
    ch_r1_geo21.Draw('Iso_SL_mod>>h_isol_sl_2', cut_noiso)
    print h1.Integral(), h2.Integral()
    print ch_r1.GetEntries(cut_all.GetTitle())/ float(ch_r1.GetEntries(cut_noiso.GetTitle())) * 100
    print ch_r1_geo21.GetEntries(cut_all.GetTitle()) / float(ch_r1_geo21.GetEntries(cut_noiso.GetTitle())) * 100
    canv = compare_shape(h1, h2, xtitle='Subleading Photon Isolation [GeV] (run1 graviton style)', saving_name='./plots/dummy_iso_sl_shape.pdf')

    
    # MGG 
    h1 = ROOT.TH1F("mgg_old", "mgg", 100, 0, 1000)
    h2 = ROOT.TH1F("mgg_new", "mgg", 100, 0, 1000)
    ch_r1.Draw("mgg>>mgg_old", cut_all)
    ch_r1_geo21.Draw("mgg>>mgg_new", cut_all)
    canv = compare_shape(
        h1, h2, 
        xtitle='m_{#gamma#gamma} [GeV]', 
        normalize=False,
        logy=True,
        saving_name='./plots/dummy_mgg_comparison.pdf')
