import numpy as np
import awkward as ak
import uproot
import ROOT

# --- Config ---
parquet_file = "AllDatasets.parquet"
data_pu_files = {
    "2016": "/eos/user/b/bbapi/CMSSW_14_0_15/src/pileup_histos/pileup_2016_Golden.root",
    "2017": "/eos/user/b/bbapi/CMSSW_14_0_15/src/pileup_histos/pileup_2017_Golden.root",
    "2018": "/eos/user/b/bbapi/CMSSW_14_0_15/src/pileup_histos/pileup_2018_Golden.root",
}
output_root = "PUPlots_AllYears.root"

ROOT.TH1.SetDefaultSumw2(True)  # Proper error handling for histograms

# Load events once
Events = ak.from_parquet(parquet_file)

# Open ROOT output file once
fout = ROOT.TFile(output_root, "RECREATE")

for year, data_pu_file in data_pu_files.items():
    print(f"\n=== Processing {year} ===")

    # Filter datasets for this year
   # mask = np.char.find(ak.to_numpy(Events["dataset"]), f"UL{year}") >= 0
    mask = np.char.find(ak.to_numpy(Events["dataset"]), f"UL{str(year)[2:]}") >= 0
    Events_year = Events[mask]
    mass_points = sorted(set(ak.to_numpy(Events_year["dataset"])))

    # Load Data pileup histogram
    with uproot.open(data_pu_file) as f:
        data_vals, bin_edges = f["pileup"].to_numpy()
    data_pdf = data_vals / np.sum(data_vals)

    # Create directory for this year in the output ROOT file
    fout.mkdir(year)
    fout.cd(year)

    for mass in mass_points:
        print(f"  → {mass}")

        # Select events for this dataset
        M_MC = Events_year[ak.to_numpy(Events_year["dataset"]) == mass]
        mc_nTrueInt = ak.to_numpy(ak.flatten(M_MC["pu_true"]))

        # MC pileup histogram (normalized)
        mc_vals, _ = np.histogram(mc_nTrueInt, bins=bin_edges)
        mc_pdf = mc_vals / np.sum(mc_vals)

        # PU weights
        eps = 1e-8
        pu_weights_per_bin = np.where(mc_pdf > eps, data_pdf / mc_pdf, 0.0)

        # MC after reweighting
        bin_idx = np.digitize(mc_nTrueInt, bin_edges) - 1
        bin_idx = np.clip(bin_idx, 0, len(pu_weights_per_bin) - 1)
        pu_weight = pu_weights_per_bin[bin_idx]
        mc_vals_rw, _ = np.histogram(mc_nTrueInt, bins=bin_edges, weights=pu_weight)
        mc_pdf_rw = mc_vals_rw / np.sum(mc_vals_rw)

        # Ratio after reweighting
        ratio_rw = np.where(data_pdf > eps, mc_pdf_rw / data_pdf, 0)

        # Unique histogram names — detached from ROOT directories
        h_mc = ROOT.TH1F(f"h_mc_{mass}", "", len(bin_edges)-1, bin_edges)
        h_mc.SetDirectory(0)
        h_data = ROOT.TH1F(f"h_data_{mass}", "", len(bin_edges)-1, bin_edges)
        h_data.SetDirectory(0)
        h_weights = ROOT.TH1F(f"h_weights_{mass}", "", len(bin_edges)-1, bin_edges)
        h_weights.SetDirectory(0)
        h_mc_rw = ROOT.TH1F(f"h_mc_rw_{mass}", "", len(bin_edges)-1, bin_edges)
        h_mc_rw.SetDirectory(0)
        h_ratio = ROOT.TH1F(f"h_ratio_{mass}", "", len(bin_edges)-1, bin_edges)
        h_ratio.SetDirectory(0)

        for i in range(len(mc_pdf)):
            h_mc.SetBinContent(i+1, mc_pdf[i])
            h_data.SetBinContent(i+1, data_pdf[i])
            h_weights.SetBinContent(i+1, pu_weights_per_bin[i])
            h_mc_rw.SetBinContent(i+1, mc_pdf_rw[i])
            h_ratio.SetBinContent(i+1, ratio_rw[i])

        # Write histograms
        fout.cd(year)
        h_mc.Write()
        h_data.Write()
        h_weights.Write()
        h_mc_rw.Write()
        h_ratio.Write()

        # --- Create Canvases ---
        # 1. MC vs Data pileup
        c1 = ROOT.TCanvas(f"c_mc_data_{mass}", f"MC vs Data Pileup - {mass}", 800, 600)
        h_mc.SetLineColor(ROOT.kRed)
        h_mc.SetLineWidth(2)
        h_mc.Draw("HIST")
        h_data.SetMarkerStyle(20)
        h_data.SetMarkerSize(1.0)
        h_data.SetMarkerColor(ROOT.kBlue)
        h_data.SetLineColor(ROOT.kBlue)
        h_data.Draw("EP SAME")
        leg1 = ROOT.TLegend(0.6, 0.7, 0.88, 0.88)
        leg1.AddEntry(h_mc, "MC", "l")
        leg1.AddEntry(h_data, "Data", "p")
        leg1.Draw()
        c1.Update()
        c1.Write()

        # 2. PU Weights
        c2 = ROOT.TCanvas(f"c_weights_{mass}", f"PU Weights - {mass}", 800, 600)
        h_weights.SetLineColor(ROOT.kBlack)
        h_weights.SetLineWidth(2)
        h_weights.Draw("HIST")
        c2.Update()
        c2.Write()

        # 3. MC after reweighting vs Data
        c3 = ROOT.TCanvas(f"c_mc_rw_data_{mass}", f"MC after PU Reweighting - {mass}", 800, 600)
        h_mc_rw.SetLineColor(ROOT.kGreen+2)
        h_mc_rw.SetLineWidth(2)
        h_mc_rw.Draw("HIST")
        h_data.Draw("EP SAME")
        leg3 = ROOT.TLegend(0.6, 0.7, 0.88, 0.88)
        leg3.AddEntry(h_mc_rw, "MC reweighted", "l")
        leg3.AddEntry(h_data, "Data", "p")
        leg3.Draw()
        c3.Update()
        c3.Write()

        # 4. Ratio after reweighting
        c4 = ROOT.TCanvas(f"c_ratio_{mass}", f"Ratio MC_rw/Data - {mass}", 800, 600)
        h_ratio.SetLineColor(ROOT.kMagenta)
        h_ratio.SetLineWidth(2)
        h_ratio.Draw("HIST")
        c4.Update()
        c4.Write()

fout.Close()
print(f"\nSaved PU plots and histograms for 2016, 2017, and 2018 in {output_root}")
