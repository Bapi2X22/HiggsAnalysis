import ROOT
import numpy as np
import os
import awkward as ak


Events_ak = ak.from_parquet("Mass_20_2018_WH.parquet")

# higgs_pt = ak.to_numpy(Events_ak["higgs_pt"])

higgs_pt = ak.to_numpy(Events_ak["higgs_pt"])
higgs_eta = ak.to_numpy(Events_ak["higgs_eta"])
higgs_phi = ak.to_numpy(Events_ak["higgs_phi"])
A_pt_1 = ak.to_numpy(Events_ak["A_pt_1"])
A_pt_2 = ak.to_numpy(Events_ak["A_pt_2"])
A_eta_1 = ak.to_numpy(Events_ak["A_eta_1"])
A_eta_2 = ak.to_numpy(Events_ak["A_eta_2"])
A_phi_1 = ak.to_numpy(Events_ak["A_phi_1"])
A_phi_2 = ak.to_numpy(Events_ak["A_phi_2"])
leading_A_pt = ak.to_numpy(Events_ak["leading_A_pt"])
subleading_A_pt = ak.to_numpy(Events_ak["subleading_A_pt"])
leading_A_eta = ak.to_numpy(Events_ak["leading_A_eta"])
subleading_A_eta = ak.to_numpy(Events_ak["subleading_A_eta"])
leading_A_phi = ak.to_numpy(Events_ak["leading_A_phi"])
subleading_A_phi = ak.to_numpy(Events_ak["subleading_A_phi"])
pho_from_a_pt_1 = ak.to_numpy(Events_ak["pho_from_a_pt_1"])
pho_from_a_pt_2 = ak.to_numpy(Events_ak["pho_from_a_pt_2"])
pho_from_a_eta_1 = ak.to_numpy(Events_ak["pho_from_a_eta_1"])
pho_from_a_eta_2 = ak.to_numpy(Events_ak["pho_from_a_eta_2"])
pho_from_a_phi_1 = ak.to_numpy(Events_ak["pho_from_a_phi_1"])
pho_from_a_phi_2 = ak.to_numpy(Events_ak["pho_from_a_phi_2"])
lead_pt_pho_gen = ak.to_numpy(Events_ak["lead_pt_pho_gen"])
sublead_pt_pho_gen = ak.to_numpy(Events_ak["sublead_pt_pho_gen"])
lead_eta_pho_gen = ak.to_numpy(Events_ak["lead_eta_pho_gen"])
sublead_eta_pho_gen = ak.to_numpy(Events_ak["sublead_eta_pho_gen"])
lead_phi_pho_gen = ak.to_numpy(Events_ak["lead_phi_pho_gen"])
sublead_phi_pho_gen = ak.to_numpy(Events_ak["sublead_phi_pho_gen"])
bquark_from_a_pt_1 = ak.to_numpy(Events_ak["bquark_from_a_pt_1"])
bquark_from_a_pt_2 = ak.to_numpy(Events_ak["bquark_from_a_pt_2"])
bquark_from_a_eta_1 = ak.to_numpy(Events_ak["bquark_from_a_eta_1"])
bquark_from_a_eta_2 = ak.to_numpy(Events_ak["bquark_from_a_eta_2"])
bquark_from_a_phi_1 = ak.to_numpy(Events_ak["bquark_from_a_phi_1"])
bquark_from_a_phi_2 = ak.to_numpy(Events_ak["bquark_from_a_phi_2"])
lead_pt_bquark_gen = ak.to_numpy(Events_ak["lead_pt_bquark_gen"])
sublead_pt_bquark_gen = ak.to_numpy(Events_ak["sublead_pt_bquark_gen"])
lead_eta_bquark_gen = ak.to_numpy(Events_ak["lead_eta_bquark_gen"])
sublead_eta_bquark_gen = ak.to_numpy(Events_ak["sublead_eta_bquark_gen"])
lead_phi_bquark_gen = ak.to_numpy(Events_ak["lead_phi_bquark_gen"])
sublead_phi_bquark_gen = ak.to_numpy(Events_ak["sublead_phi_bquark_gen"])
Reco_pho_pt = ak.to_numpy(Events_ak["Reco_pho_pt"])
Reco_pho_eta = ak.to_numpy(Events_ak["Reco_pho_eta"])
Reco_pho_phi = ak.to_numpy(Events_ak["Reco_pho_phi"])
Reco_lead_pho_pt = ak.to_numpy(Events_ak["Reco_lead_pho_pt"])
Reco_sublead_pho_pt = ak.to_numpy(Events_ak["Reco_sublead_pho_pt"])
Reco_lead_pho_eta = ak.to_numpy(Events_ak["Reco_lead_pho_eta"])
Reco_sublead_pho_eta = ak.to_numpy(Events_ak["Reco_sublead_pho_eta"])
Reco_lead_pho_phi = ak.to_numpy(Events_ak["Reco_lead_pho_phi"])
Reco_sublead_pho_phi = ak.to_numpy(Events_ak["Reco_sublead_pho_phi"])
gen_photon_from_a_1_pt = ak.to_numpy(Events_ak["gen_photon_from_a_1_pt"])
gen_photon_from_a_2_pt = ak.to_numpy(Events_ak["gen_photon_from_a_2_pt"])
Gen_photon_pt_1 = ak.to_numpy(Events_ak["Gen_photon_pt_1"])
Gen_photon_pt_2 = ak.to_numpy(Events_ak["Gen_photon_pt_2"])
Genmatched_pho_1_pt = ak.to_numpy(Events_ak["Genmatched_pho_1_pt"])
Genmatched_pho_2_pt = ak.to_numpy(Events_ak["Genmatched_pho_2_pt"])
Genmatched_pho_1_eta = ak.to_numpy(Events_ak["Genmatched_pho_1_eta"])
Genmatched_pho_2_eta = ak.to_numpy(Events_ak["Genmatched_pho_2_eta"])
Genmatched_pho_1_phi = ak.to_numpy(Events_ak["Genmatched_pho_1_phi"])
Genmatched_pho_2_phi = ak.to_numpy(Events_ak["Genmatched_pho_2_phi"])
gen_lead_b_pt = ak.to_numpy(Events_ak["gen_lead_b_pt"])
gen_lead_b_eta = ak.to_numpy(Events_ak["gen_lead_b_eta"])
gen_lead_b_phi = ak.to_numpy(Events_ak["gen_lead_b_phi"])
gen_sublead_b_pt = ak.to_numpy(Events_ak["gen_sublead_b_pt"])
gen_sublead_b_eta = ak.to_numpy(Events_ak["gen_sublead_b_eta"])
gen_sublead_b_phi = ak.to_numpy(Events_ak["gen_sublead_b_phi"])
gen_b1_pt = ak.to_numpy(Events_ak["gen_b1_pt"])
gen_b1_eta = ak.to_numpy(Events_ak["gen_b1_eta"])
gen_b1_phi = ak.to_numpy(Events_ak["gen_b1_phi"])
gen_b2_pt = ak.to_numpy(Events_ak["gen_b2_pt"])
gen_b2_eta = ak.to_numpy(Events_ak["gen_b2_eta"])
gen_b2_phi = ak.to_numpy(Events_ak["gen_b2_phi"])
reco_lead_bjet_pt = ak.to_numpy(Events_ak["reco_lead_bjet_pt"])
reco_lead_bjet_eta = ak.to_numpy(Events_ak["reco_lead_bjet_eta"])
reco_lead_bjet_phi = ak.to_numpy(Events_ak["reco_lead_bjet_phi"])
reco_sublead_bjet_pt = ak.to_numpy(Events_ak["reco_sublead_bjet_pt"])
reco_sublead_bjet_eta = ak.to_numpy(Events_ak["reco_sublead_bjet_eta"])
reco_sublead_bjet_phi = ak.to_numpy(Events_ak["reco_sublead_bjet_phi"])
matched_bjet1_pt = ak.to_numpy(Events_ak["matched_bjet1_pt"])
matched_bjet1_eta = ak.to_numpy(Events_ak["matched_bjet1_eta"])
matched_bjet1_phi = ak.to_numpy(Events_ak["matched_bjet1_phi"])
matched_bjet2_pt = ak.to_numpy(Events_ak["matched_bjet2_pt"])
matched_bjet2_eta = ak.to_numpy(Events_ak["matched_bjet2_eta"])
matched_bjet2_phi = ak.to_numpy(Events_ak["matched_bjet2_phi"])
gen_invmasses_diphoton = ak.to_numpy(Events_ak["gen_invmasses_diphoton"])
gen_invmasses_bb = ak.to_numpy(Events_ak["gen_invmasses_bb"])
invmasses_diphoton = ak.to_numpy(Events_ak["invmasses_diphoton"])
invmasses_bb = ak.to_numpy(Events_ak["invmasses_bb"])

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\------Creating Hitogram Part---\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


branch_arrays = {
    "higgs_pt": higgs_pt,
    "higgs_eta": higgs_eta,
    "higgs_phi": higgs_phi,
    "A_pt_1": A_pt_1,
    "A_pt_2": A_pt_2,
    "A_eta_1": A_eta_1,
    "A_eta_2": A_eta_2,
    "A_phi_1": A_phi_1,
    "A_phi_2": A_phi_2,
    "leading_A_pt": leading_A_pt,
    "subleading_A_pt": subleading_A_pt,
    "leading_A_eta": leading_A_eta,
    "subleading_A_eta": subleading_A_eta,
    "leading_A_phi": leading_A_phi,
    "subleading_A_phi": subleading_A_phi,
    "pho_from_a_pt_1": pho_from_a_pt_1,
    "pho_from_a_pt_2": pho_from_a_pt_2,
    "pho_from_a_eta_1": pho_from_a_eta_1,
    "pho_from_a_eta_2": pho_from_a_eta_2,
    "pho_from_a_phi_1": pho_from_a_phi_1,
    "pho_from_a_phi_2": pho_from_a_phi_2,
    "lead_pt_pho_gen": lead_pt_pho_gen,
    "sublead_pt_pho_gen": sublead_pt_pho_gen,
    "lead_eta_pho_gen": lead_eta_pho_gen,
    "sublead_eta_pho_gen": sublead_eta_pho_gen,
    "lead_phi_pho_gen": lead_phi_pho_gen,
    "sublead_phi_pho_gen": sublead_phi_pho_gen,
    "bquark_from_a_pt_1": bquark_from_a_pt_1,
    "bquark_from_a_pt_2": bquark_from_a_pt_2,
    "bquark_from_a_eta_1": bquark_from_a_eta_1,
    "bquark_from_a_eta_2": bquark_from_a_eta_2,
    "bquark_from_a_phi_1": bquark_from_a_phi_1,
    "bquark_from_a_phi_2": bquark_from_a_phi_2,
    "lead_pt_bquark_gen": lead_pt_bquark_gen,
    "sublead_pt_bquark_gen": sublead_pt_bquark_gen,
    "lead_eta_bquark_gen": lead_eta_bquark_gen,
    "sublead_eta_bquark_gen": sublead_eta_bquark_gen,
    "lead_phi_bquark_gen": lead_phi_bquark_gen,
    "sublead_phi_bquark_gen": sublead_phi_bquark_gen,
    "Reco_pho_pt": Reco_pho_pt,
    "Reco_pho_eta": Reco_pho_eta,
    "Reco_pho_phi": Reco_pho_phi,
    "Reco_lead_pho_pt": Reco_lead_pho_pt,
    "Reco_sublead_pho_pt": Reco_sublead_pho_pt,
    "Reco_lead_pho_eta": Reco_lead_pho_eta,
    "Reco_sublead_pho_eta": Reco_sublead_pho_eta,
    "Reco_lead_pho_phi": Reco_lead_pho_phi,
    "Reco_sublead_pho_phi": Reco_sublead_pho_phi,
    "gen_photon_from_a_1_pt": gen_photon_from_a_1_pt,
    "gen_photon_from_a_2_pt": gen_photon_from_a_2_pt,
    "Gen_photon_pt_1": Gen_photon_pt_1,
    "Gen_photon_pt_2": Gen_photon_pt_2,
    "Genmatched_pho_1_pt": Genmatched_pho_1_pt,
    "Genmatched_pho_2_pt": Genmatched_pho_2_pt,
    "Genmatched_pho_1_eta": Genmatched_pho_1_eta,
    "Genmatched_pho_2_eta": Genmatched_pho_2_eta,
    "Genmatched_pho_1_phi": Genmatched_pho_1_phi,
    "Genmatched_pho_2_phi": Genmatched_pho_2_phi,
    "gen_lead_b_pt": gen_lead_b_pt,
    "gen_lead_b_eta": gen_lead_b_eta,
    "gen_lead_b_phi": gen_lead_b_phi,
    "gen_sublead_b_pt": gen_sublead_b_pt,
    "gen_sublead_b_eta": gen_sublead_b_eta,
    "gen_sublead_b_phi": gen_sublead_b_phi,
    "gen_b1_pt": gen_b1_pt,
    "gen_b1_eta": gen_b1_eta,
    "gen_b1_phi": gen_b1_phi,
    "gen_b2_pt": gen_b2_pt,
    "gen_b2_eta": gen_b2_eta,
    "gen_b2_phi": gen_b2_phi,
    "reco_lead_bjet_pt": reco_lead_bjet_pt,
    "reco_lead_bjet_eta": reco_lead_bjet_eta,
    "reco_lead_bjet_phi": reco_lead_bjet_phi,
    "reco_sublead_bjet_pt": reco_sublead_bjet_pt,
    "reco_sublead_bjet_eta": reco_sublead_bjet_eta,
    "reco_sublead_bjet_phi": reco_sublead_bjet_phi,
    "matched_bjet1_pt": matched_bjet1_pt,
    "matched_bjet1_eta": matched_bjet1_eta,
    "matched_bjet1_phi": matched_bjet1_phi,
    "matched_bjet2_pt": matched_bjet2_pt,
    "matched_bjet2_eta": matched_bjet2_eta,
    "matched_bjet2_phi": matched_bjet2_phi,
    "gen_invmasses_diphoton": gen_invmasses_diphoton,
    "gen_invmasses_bb": gen_invmasses_bb,
    "invmasses_diphoton": invmasses_diphoton,
    "invmasses_bb": invmasses_bb,
}


bin_settings = {
    # "nPhoton": (100, 0, 10),

    # Higgs
    "higgs_pt":  (500, 0, 500),
    "higgs_eta": (100, 1, -1),
    "higgs_phi": (64, 1, -1),

    "A_pt_1": (500, 0, 500),
    "A_pt_2": (500, 0, 500),
    "A_eta_1": (100, 1, -1),
    "A_eta_2": (100, 1, -1),
    "A_phi_1": (64, 1, -1),
    "A_phi_2": (64, 1, -1),
    "leading_A_pt": (500, 0, 500),
    "subleading_A_pt": (500, 0, 500),
    "leading_A_eta": (100, 1, -1),
    "subleading_A_eta": (100, 1, -1),
    "leading_A_phi": (64, 1, -1),
    "subleading_A_phi": (64, 1, -1),

    # Unsorted photons from a
    "pho_from_a_pt_1": (500, 0, 500),
    "pho_from_a_pt_2": (500, 0, 500),
    "pho_from_a_eta_1": (100, 1, -1),
    "pho_from_a_eta_2": (100, 1, -1),
    "pho_from_a_phi_1": (64, 1, -1),
    "pho_from_a_phi_2": (64, 1, -1),

    # Sorted photons from a
    "lead_pt_pho_gen":     (500, 0, 500),
    "sublead_pt_pho_gen":  (500, 0, 500),
    "lead_eta_pho_gen":    (100, 1, -1),
    "sublead_eta_pho_gen": (100, 1, -1),
    "lead_phi_pho_gen":    (64, 1, -1),
    "sublead_phi_pho_gen": (64, 1, -1),

    # Unsorted b-quarks from a
    "bquark_from_a_pt_1": (500, 0, 500),
    "bquark_from_a_pt_2": (500, 0, 500),
    "bquark_from_a_eta_1": (100, 1, -1),
    "bquark_from_a_eta_2": (100, 1, -1),
    "bquark_from_a_phi_1": (64, 1,-1),
    "bquark_from_a_phi_2": (64, 1,-1),

    # Sorted b-quarks from a
    "lead_pt_bquark_gen":     (500, 0, 500),
    "sublead_pt_bquark_gen":  (500, 0, 500),
    "lead_eta_bquark_gen":    (100, 1, -1),
    "sublead_eta_bquark_gen": (100, 1, -1),
    "lead_phi_bquark_gen":    (64, 1,-1),
    "sublead_phi_bquark_gen": (64, 1,-1),

    "Reco_pho_pt":  (500, 0, 500),
    "Reco_pho_eta": (200, -10, 10), 
    "Reco_pho_phi": (64, 1, -1),   


    "Reco_lead_pho_pt":     (500, 0, 500),
    "Reco_sublead_pho_pt":  (500, 0, 500),
    "Reco_lead_pho_eta":    (100, 1, -1),
    "Reco_sublead_pho_eta": (100, 1, -1),
    "Reco_lead_pho_phi":    (64, 1, -1),
    "Reco_sublead_pho_phi": (64, 1, -1),

    "gen_photon_from_a_1_pt": (500, 0, 500),
    "gen_photon_from_a_2_pt": (500, 0, 500),

    "Gen_photon_pt_1": (500, 0, 500),
    "Gen_photon_pt_2": (500, 0, 500),

    "Genmatched_pho_1_pt":  (500, 0, 500),
    "Genmatched_pho_2_pt":  (500, 0, 500),
    "Genmatched_pho_1_eta": (100, 1, -1),
    "Genmatched_pho_2_eta": (100, 1, -1),
    "Genmatched_pho_1_phi": (64, 1, -1),
    "Genmatched_pho_2_phi": (64, 1, -1),

        # Gen-level b-quarks sorted
    "gen_lead_b_pt":     (500, 0, 500),
    "gen_lead_b_eta":    (100, 1, -1),
    "gen_lead_b_phi":    (64, 1, -1),
    "gen_sublead_b_pt":  (500, 0, 500),
    "gen_sublead_b_eta": (100, 1, -1),
    "gen_sublead_b_phi": (64, 1, -1),

    # Gen-level b-quarks unsorted
    "gen_b1_pt":  (500, 0, 500),
    "gen_b1_eta": (100, 1, -1),
    "gen_b1_phi": (64, 1, -1),
    "gen_b2_pt":  (500, 0, 500),
    "gen_b2_eta": (100, 1, -1),
    "gen_b2_phi": (64, 1, -1),

    # Reco-level b-jets sorted
    "reco_lead_bjet_pt":     (500, 0, 500),
    "reco_lead_bjet_eta":    (100, 1, -1),
    "reco_lead_bjet_phi":    (64, 1, -1),
    "reco_sublead_bjet_pt":  (500, 0, 500),
    "reco_sublead_bjet_eta": (100, 1, -1),
    "reco_sublead_bjet_phi": (64, 1, -1),

    # Gen-matched reco b-jets
    "matched_bjet1_pt":  (500, 0, 500),
    "matched_bjet1_eta": (100, 1, -1),
    "matched_bjet1_phi": (64, 1, -1),
    "matched_bjet2_pt":  (500, 0, 500),
    "matched_bjet2_eta": (100, 1, -1),
    "matched_bjet2_phi": (64, 1, -1),

    "gen_invmasses_diphoton": (200,19,21),
    "gen_invmasses_bb": (200, 17, 21),
    "invmasses_diphoton": (500, 15, 25),
    "invmasses_bb": (500, 5, 50),
}


# ----------------------------------------------------------------------Plotting Part -------------------------------------------------------------------------------------------


# Mode options: "create", "update", "selected"
mode = "selected"

# If using 'selected' mode, specify which ones to process
selected_branches = [("gen_photon_from_a_1_pt", "Genmatched_pho_1_pt"), ("gen_photon_from_a_2_pt", "Genmatched_pho_2_pt"), ("Gen_photon_pt_1", "Genmatched_pho_1_pt"), ("Gen_photon_pt_2", "Genmatched_pho_2_pt"), ("gen_b1_pt", "matched_bjet1_pt"), ("gen_b2_pt", "matched_bjet2_pt"), "matched_bjet2_eta", "invmasses_bb",     ("leading_A_pt", "subleading_A_pt"),
    ("lead_pt_pho_gen", "sublead_pt_pho_gen"),
    ("lead_pt_bquark_gen", "sublead_pt_bquark_gen"),
    ("Reco_lead_pho_pt", "Reco_sublead_pho_pt"),
    ("reco_lead_bjet_pt", "reco_sublead_bjet_pt"),
    ("A_pt_1", "A_pt_2")]  # example

# Select which branches to process
if mode == "selected":
    branches_to_process = selected_branches
else:
    branches_to_process = list(branch_arrays.keys())

# Open ROOT file
file_mode = "UPDATE" if os.path.exists("hist_output_WH_20_2018_new.root") else "RECREATE"
f = ROOT.TFile("hist_output_WH_20_2018_new.root", file_mode)

for branch in branches_to_process:
    is_2d = isinstance(branch, tuple) and len(branch) == 2
    if is_2d:
        b1, b2 = branch
        if b1 not in branch_arrays or b2 not in branch_arrays:
            print(f"Skipping 2D histogram: '{branch}' not found.")
            continue
        hist_name = f"h_{b1}_vs_{b2}"
    else:
        if branch not in branch_arrays:
            print(f"Skipping: Branch '{branch}' not provided.")
            continue
        hist_name = f"h_{branch}"

    hist_exists = bool(f.Get(hist_name))

    if mode == "create":
        if hist_exists:
            print(f"Already exists: {branch} skipping (create mode)")
            continue
        print(f"Creating new histogram: {branch}")
    elif mode == "update":
        if not hist_exists:
            print(f"Missing histogram: {branch} skipping (update mode)")
            continue
        print(f"Updating histogram: {branch}")
        f.Delete(f"{hist_name};*")
    elif mode == "selected":
        if hist_exists:
            print(f"Overwriting histogram: {branch}")
            f.Delete(f"{hist_name};*")
        else:
            print(f"Creating histogram: {branch}")
    else:
        raise ValueError("Invalid mode: must be 'create', 'update', or 'selected'")

    # # Save histogram in ROOT file
    # hist.Write()
    outdir = "/eos/user/b/bbapi/www/WHToAA/2018/M_20/"
    os.makedirs(outdir, exist_ok=True)
    canvas = ROOT.TCanvas("canvas", "", 800, 600)

    if is_2d:
        data_x = np.asarray(branch_arrays[b1])
        data_y = np.asarray(branch_arrays[b2])
        nbinsx, xmin, xmax = bin_settings.get(b1, (100, float(np.min(data_x)), float(np.max(data_x))))
        nbinsy, ymin, ymax = bin_settings.get(b2, (100, float(np.min(data_y)), float(np.max(data_y))))
        hist = ROOT.TH2F(hist_name, f"{b1} vs {b2};{b1};{b2}", nbinsx, xmin, xmax, nbinsy, ymin, ymax)
        hist.Sumw2()
        for x, y in zip(data_x, data_y):
            hist.Fill(x, y)
        hist.Draw("COLZ")

        # Save 2D plots
        png_path = os.path.join(outdir, f"{hist_name}_2D.png")
        pdf_path = os.path.join(outdir, f"{hist_name}_2D.pdf")
        canvas.SaveAs(png_path)
        canvas.SaveAs(pdf_path)

    else:
        data = np.asarray(branch_arrays[branch])
        nbins, xmin, xmax = bin_settings.get(branch, (100, float(np.min(data)), float(np.max(data))))
        hist = ROOT.TH1F(hist_name, f"{branch};{branch};Events", nbins, xmin, xmax)
        hist.Sumw2()
        for val in data:
            hist.Fill(val)
        hist.Draw()

        # Save normal version
        png_path = os.path.join(outdir, f"{hist_name}.png")
        pdf_path = os.path.join(outdir, f"{hist_name}.pdf")
        canvas.SaveAs(png_path)
        canvas.SaveAs(pdf_path)

        # Save log y-scale version if "pt" is in branch name
        if "pt" in branch.lower():
            canvas.SetLogy()
            hist.Draw()
            log_png_path = os.path.join(outdir, f"{hist_name}_logy.png")
            log_pdf_path = os.path.join(outdir, f"{hist_name}_logy.pdf")
            canvas.SaveAs(log_png_path)
            canvas.SaveAs(log_pdf_path)
            canvas.SetLogy(0)

    # Save histogram in ROOT file
    hist.Write()


f.Close()
print(f"\nDone writing histograms in '{mode}' mode.")





