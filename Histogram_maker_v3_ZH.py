import ROOT
import numpy as np
import os
import awkward as ak

# Read the full parquet file
# all_events = ak.from_parquet("WH_all_datasets.parquet")
all_events = ak.from_parquet("AllDatasets_ZH.parquet")

# Extract dataset names
unique_datasets = ak.to_list(all_events.dataset)

# Define branch plotting configuration
bin_settings = {
    "PileUp": (100, 0, 100),
    "nPhoton": (10, 0, 10), 

    "higgs_pt":  (500, 0, 500),
    "higgs_eta": (100, 1, -1),
    "higgs_phi": (64, -3.2, 3.2),

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

    "Reco_photon_lead_pt_all": (500, 0, 500),
    "Reco_photon_sublead_pt_all": (500, 0, 500),
    "Reco_photon_lead_eta_all": (100, 1, -1),
    "Reco_photon_sublead_eta_all": (100, 1, -1),
    "Reco_photon_lead_phi_all": (64, 1, -1),
    "Reco_photon_sublead_phi_all": (64, 1, -1),

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

    "gen_invmasses_diphoton": (10000,0,100),
    "gen_invmasses_bb": (10000, 0, 100),
    "invmasses_diphoton": (500, 0, 100),
    "invmasses_bb": (500, 0, 100),
}

# Mode: "create", "update", "selected"
mode = "selected"

# # If selected, specify# Selected branches
# selected_branches = [("gen_photon_from_a_1_pt", "Genmatched_pho_1_pt"), ("gen_photon_from_a_2_pt", "Genmatched_pho_2_pt"), ("Gen_photon_pt_1", "Genmatched_pho_1_pt"), ("Gen_photon_pt_2", "Genmatched_pho_2_pt"), ("gen_b1_pt", "matched_bjet1_pt"), ("gen_b2_pt", "matched_bjet2_pt"), ("leading_A_pt", "subleading_A_pt"),
#     ("lead_pt_pho_gen", "sublead_pt_pho_gen"),
#     ("lead_pt_bquark_gen", "sublead_pt_bquark_gen"),
#     ("Reco_lead_pho_pt", "Reco_sublead_pho_pt"),
#     ("reco_lead_bjet_pt", "reco_sublead_bjet_pt"),
#     ("A_pt_1", "A_pt_2")]  # example

selected_branches = [("Genmatched_pho_1_pt", "Genmatched_pho_2_pt")]

# selected_branches = ["Reco_photon_lead_pt_all", "Reco_photon_sublead_pt_all", "Reco_photon_lead_eta_all", "Reco_photon_sublead_eta_all", "Reco_photon_lead_phi_all", "Reco_photon_sublead_phi_all"]

# selected_branches = ["gen_invmasses_diphoton", "gen_invmasses_bb"] 
# Open ROOT file
root_filename = "hist_output_ZH.root"
file_mode = "UPDATE" if os.path.exists(root_filename) else "RECREATE"
f = ROOT.TFile(root_filename, file_mode)

# Loop over datasets
for dataset_name in unique_datasets:
    print(f"\n=== Processing dataset: {dataset_name} ===")

    # Filter events for this dataset
    # dataset_events = all_events[all_events["dataset"] == dataset_name].data
    dataset_events = all_events[all_events["dataset"] == dataset_name]

    # Build branch_arrays for this dataset
    branch_arrays = {}
    for branch in dataset_events.fields[:-1]:
        branch_arrays[branch] = ak.to_numpy(dataset_events[branch])[0]

    # Create output directory for plots
    outdir = f"/eos/user/b/bbapi/www/ZHToAA/{dataset_name}/"
    os.makedirs(outdir, exist_ok=True)

    # Create a subdirectory in ROOT file
    f.mkdir(dataset_name) if not f.GetDirectory(dataset_name) else None
    f.cd(dataset_name)

    # Select branches to process
    if mode == "selected":
        branches_to_process = selected_branches
    else:
        branches_to_process = list(branch_arrays.keys())

    # Loop over branches
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

        hist_exists = bool(f.Get(f"{dataset_name}/{hist_name}"))

        # Handle modes
        if mode == "create" and hist_exists:
            print(f"Already exists: {branch} skipping")
            continue
        elif mode == "update" and not hist_exists:
            print(f"Missing histogram: {branch} skipping")
            continue
        elif mode in ("selected", "update") and hist_exists:
            print(f"Overwriting histogram: {branch}")
            f.Delete(f"{dataset_name}/{hist_name};*")
        else:
            print(f"Creating histogram: {branch}")

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

            canvas.SaveAs(os.path.join(outdir, f"{hist_name}_2D.png"))
            canvas.SaveAs(os.path.join(outdir, f"{hist_name}_2D.pdf"))

        else:
            data = np.asarray(branch_arrays[branch])
            data = data[~np.isnan(data)]
            nbins, xmin, xmax = bin_settings.get(branch, (100, float(np.min(data)), float(np.max(data))))
            hist = ROOT.TH1F(hist_name, f"{branch};{branch};Events", nbins, xmin, xmax)
            hist.Sumw2()
            for val in data:
                hist.Fill(val)
                        # Style settings to match PNG/PDF
            hist.SetLineColor(ROOT.kBlue)
            hist.SetLineWidth(2)
            hist.SetMarkerStyle(20)
            hist.SetMarkerSize(0.8)
            hist.SetDrawOption("HIST E1")
            hist.Draw("HIST E1")

            # Save normal version
            canvas.SaveAs(os.path.join(outdir, f"{hist_name}.png"))
            canvas.SaveAs(os.path.join(outdir, f"{hist_name}.pdf"))

            # Log scale if needed
            if "pt" in branch.lower():
                canvas.SetLogy()
                hist.Draw("HIST E1")
                canvas.SaveAs(os.path.join(outdir, f"{hist_name}_logy.png"))
                canvas.SaveAs(os.path.join(outdir, f"{hist_name}_logy.pdf"))
                canvas.SetLogy(0)

        hist.Write()
        canvas.Close()

f.Close()
print("\nAll datasets processed successfully.")

