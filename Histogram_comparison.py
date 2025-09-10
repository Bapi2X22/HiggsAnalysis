import ROOT
import os


def compare_hists(
    histname, 
    file1, dir1, xsec1, lumi1,
    file2, dir2, xsec2, lumi2,
    file3=None, dir3=None, xsec3=1.0, lumi3=1.0,  # optional
    rebin=1,
    leg1="Sample1", leg2="Sample2", leg3="Sample3",
    xlabel="X-axis", png_name="comparison.png",
    outprefix="hist1_vs_hist2",
    outroot="Histo_comparison.root",
    do_scale=True
):

    # --- Load files and hists ---
    f1 = ROOT.TFile.Open(file1)
    f2 = ROOT.TFile.Open(file2)
    h1 = f1.Get(f"{dir1}/{histname}")
    h2 = f2.Get(f"{dir2}/{histname}")

    if not h1 or not h2:
        raise RuntimeError("Could not load histograms from file1/file2!")

    h1 = h1.Clone("h1_clone")
    h2 = h2.Clone("h2_clone")

    h3 = None
    if file3 and dir3:
        f3 = ROOT.TFile.Open(file3)
        h3 = f3.Get(f"{dir3}/{histname}")
        if not h3:
            raise RuntimeError("Could not load histogram from file3!")
        h3 = h3.Clone("h3_clone")

    # --- Rebin ---
    if rebin > 1:
        h1.Rebin(rebin)
        h2.Rebin(rebin)
        if h3: h3.Rebin(rebin)

    # --- Normalize ---
    if h1.Integral() > 0: h1.Scale(1.0 / h1.Integral())
    if h2.Integral() > 0: h2.Scale(1.0 / h2.Integral())
    if h3 and h3.Integral() > 0: h3.Scale(1.0 / h3.Integral())

    # --- Cross section × lumi scaling ---
    if do_scale:
        h1.Scale(xsec1 * lumi1)
        h2.Scale(xsec2 * lumi2)
        if h3: h3.Scale(xsec3 * lumi3)

    tag = "scaled" if do_scale else "raw"

    # --- Styling ---
    h1.SetLineColor(ROOT.kRed);   h1.SetLineWidth(2)
    h2.SetLineColor(ROOT.kBlue);  h2.SetLineWidth(2)
    if h3:
        h3.SetLineColor(ROOT.kGreen+2); h3.SetLineWidth(2)

    ROOT.gStyle.SetOptStat(0)

    h1.GetXaxis().SetTitle(xlabel)
    h1.GetYaxis().SetTitle("Expected Events")

    max_y = max(h1.GetMaximum(), h2.GetMaximum(), (h3.GetMaximum() if h3 else 0))
    h1.SetMaximum(1.5 * max_y)

    # --- Linear canvas ---
    c_lin = ROOT.TCanvas("c_lin", "Linear scale", 800, 600)
    h1.Draw("HIST")
    h2.Draw("HIST SAME")
    if h3: h3.Draw("HIST SAME")

    leg = ROOT.TLegend(0.55, 0.7, 0.85, 0.85)
    leg.SetTextSize(0.03); leg.SetBorderSize(0); leg.SetFillStyle(0)
    leg.AddEntry(h1, leg1, "l")
    leg.AddEntry(h2, leg2, "l")
    if h3: leg.AddEntry(h3, leg3, "l")
    leg.Draw()

    outdir = "/eos/user/b/bbapi/www/Shape_comparison/"
    c_lin.SaveAs(os.path.join(outdir, f"{histname}_{leg1}_vs_{leg2}{'_vs_'+leg3 if h3 else ''}_{tag}_lin.png"))
    c_lin.SaveAs(os.path.join(outdir, f"{histname}_{leg1}_vs_{leg2}{'_vs_'+leg3 if h3 else ''}_{tag}_lin.pdf"))

    # --- Log canvas ---
    c_log = ROOT.TCanvas("c_log", "Log scale", 800, 600)
    c_log.SetLogy()
    h1.Draw("HIST")
    h2.Draw("HIST SAME")
    if h3: h3.Draw("HIST SAME")
    leg.Draw()
    c_log.SaveAs(os.path.join(outdir, f"{histname}_{leg1}_vs_{leg2}{'_vs_'+leg3 if h3 else ''}_{tag}_log.png"))
    c_log.SaveAs(os.path.join(outdir, f"{histname}_{leg1}_vs_{leg2}{'_vs_'+leg3 if h3 else ''}_{tag}_log.pdf"))

    # --- Save ROOT outputs ---
    mode = "UPDATE" if os.path.exists(outroot) else "RECREATE"
    fout = ROOT.TFile(outroot, mode)
    fout.cd()
    h1.Write(f"{histname}_{leg1}_hist", ROOT.TObject.kOverwrite)
    h2.Write(f"{histname}_{leg2}_hist", ROOT.TObject.kOverwrite)
    if h3: h3.Write(f"{histname}_{leg3}_hist", ROOT.TObject.kOverwrite)
    c_lin.Write(f"{histname}_{leg1}_vs_{leg2}{'_vs_'+leg3 if h3 else ''}_{tag}_lin", ROOT.TObject.kOverwrite)
    c_log.Write(f"{histname}_{leg1}_vs_{leg2}{'_vs_'+leg3 if h3 else ''}_{tag}_log", ROOT.TObject.kOverwrite)
    fout.Close()

compare_hists(
    histname="h_higgs_pt",
    file1="hist_output_WH.root",
    dir1="M20_RunIISummer20UL17NanoAODv9",
    xsec1=1.373,  # WH (pb)
    lumi1=41500,  # Run2016 UL (pb^-1)

    file2="hist_output_ggH.root",
    dir2="M20_Run3Summer22NanoAODv13",
    xsec2=48.58,  # ggH (pb)
    lumi2=34700,  # Run3 2022 (pb^-1)

    file3="hist_output_ZH.root",
    dir3="M20_RunIISummer20UL17NanoAODv9",
    xsec3=0.8839,  # ZH (pb)
    lumi3=41500,   # Run2017 UL (pb^-1)

    leg1="WH (2017 UL)",
    leg2="ggH (2022 Run3)",
    leg3="ZH (2017 UL)",

    xlabel="Higgs p_{T} [GeV]",
    rebin=2,
    outprefix="Leading_Photon_pt_WH_vs_ggH_vs_ZH",
    outroot="Histo_comparison.root",
    do_scale=False
)

with open("commands.txt") as f:
    code = f.read()   # read whole file
    exec(code)        # run all compare_hists(...) calls