import ROOT
import numpy as np
import os
import awkward as ak

# Example input: provide preloaded awkward arrays (you'll define these outside)
# Example:
Events = ak.from_parquet("Events.parquet")
nPhoton = Events["nPhoton"]

gen = Events.GenPart
higgs_mask = (gen.pdgId == 25) & (gen.status == 22)
higgs = gen[higgs_mask]

higgs_pt = ak.flatten(higgs.pt)
higgs_eta = ak.flatten(higgs.eta)
higgs_phi = ak.flatten(higgs.phi)

# Step 1: Select all prompt photons (pdgId == 22)
# photons = gen[gen.pdgId == 22]
photons = gen[(gen.pdgId == 22) & (gen.status == 1)]

# Step 2: Get mother indices for each photon
mother_idx = photons.genPartIdxMother  # shape: (n_events, variable)

# Step 3: Check if mother's pdgId == 35
from_a_mask = gen[mother_idx].pdgId == 35  # shape: (n_events, variable)

# Step 4: Apply mask to get only photons from a
photons_from_a = photons[from_a_mask]

mask_two_photons = ak.num(photons_from_a.eta) >= 2
photons_from_a = photons_from_a[mask_two_photons]

pho_from_a_pt = photons_from_a.pt
pho_from_a_eta = photons_from_a.eta 
pho_from_a_phi = photons_from_a.phi

pho_from_a_pt_1 = pho_from_a_pt[:,0]
pho_from_a_pt_2 = pho_from_a_pt[:,1]
pho_from_a_eta_1 = pho_from_a_eta[:,0]  
pho_from_a_eta_2 = pho_from_a_eta[:,1]
pho_from_a_phi_1 = pho_from_a_phi[:,0]
pho_from_a_phi_2 = pho_from_a_phi[:,1]

# Sort the selected photons by pt in descending order (per event)
sorted_photons = photons_from_a[ak.argsort(photons_from_a.pt, axis=1, ascending=False)]

# Extract pt, eta, phi for leading and subleading photons
lead_pt_pho_gen  = sorted_photons.pt[:, 0]
sublead_pt_pho_gen = sorted_photons.pt[:, 1]

lead_eta_pho_gen = sorted_photons.eta[:, 0]
sublead_eta_pho_gen = sorted_photons.eta[:, 1]

lead_phi_pho_gen = sorted_photons.phi[:, 0]
sublead_phi_pho_gen = sorted_photons.phi[:, 1]


# Step 1: Select all prompt photons (pdgId == 22)
bquarks = gen[abs(gen.pdgId) == 5]

# Step 2: Get mother indices for each photon
mother_idx = bquarks.genPartIdxMother  # shape: (n_events, variable)

# Step 3: Check if mother's pdgId == 35
from_a_mask = gen[mother_idx].pdgId == 35  # shape: (n_events, variable)

# Step 4: Apply mask to get only photons from a
bquarks_from_a = bquarks[from_a_mask]

bquark_from_a_pt = bquarks_from_a.pt
bquark_from_a_eta = bquarks_from_a.eta 
bquark_from_a_phi = bquarks_from_a.phi

bquark_from_a_pt_1 = bquark_from_a_pt[:,0]
bquark_from_a_pt_2 = bquark_from_a_pt[:,1]
bquark_from_a_eta_1 = bquark_from_a_eta[:,0]  
bquark_from_a_eta_2 = bquark_from_a_eta[:,1]
bquark_from_a_phi_1 = bquark_from_a_phi[:,0]
bquark_from_a_phi_2 = bquark_from_a_phi[:,1]

# Sort the selected bquarks by pt in descending order (per event)
sorted_bquarks = bquarks_from_a[ak.argsort(bquarks_from_a.pt, axis=1, ascending=False)]

# Extract pt, eta, phi for leading and subleading photons
lead_pt_bquark_gen  = sorted_bquarks.pt[:, 0]
sublead_pt_bquark_gen = sorted_bquarks.pt[:, 1]

lead_eta_bquark_gen = sorted_bquarks.eta[:, 0]
sublead_eta_bquark_gen = sorted_bquarks.eta[:, 1]

lead_phi_bquark_gen = sorted_bquarks.phi[:, 0]
sublead_phi_bquark_gen = sorted_bquarks.phi[:, 1]


reco_photons = Events.Photon[mask_two_photons]

Reco_pho_pt = ak.flatten(reco_photons.pt)
Reco_pho_eta = ak.flatten(reco_photons.eta)
Reco_pho_phi = ak.flatten(reco_photons.phi)

Reco_bquark_pt = ak.flatten(reco_photons.pt)
Reco_bquark_eta = ak.flatten(reco_photons.eta)
Reco_bquark_phi = ak.flatten(reco_photons.phi)

Reco_pho_pt_uf = reco_photons.pt
Reco_pho_eta_uf = reco_photons.eta
Reco_pho_phi_uf = reco_photons.phi

Reco_bquark_pt_uf = reco_photons.pt  # <-- Update if needed
Reco_bquark_eta_uf = reco_photons.eta
Reco_bquark_phi_uf = reco_photons.phi

# Sort by pt in descending order per event
sorted_reco_photons = reco_photons[ak.argsort(reco_photons.pt, axis=1, ascending=False)]

has_lead = ak.num(sorted_reco_photons.pt) > 0
Reco_lead_pho_pt  = sorted_reco_photons.pt[has_lead][:, 0]
Reco_lead_pho_eta = sorted_reco_photons.eta[has_lead][:, 0]
Reco_lead_pho_phi = sorted_reco_photons.phi[has_lead][:, 0]

has_sublead = ak.num(sorted_reco_photons.pt) > 1
Reco_sublead_pho_pt  = sorted_reco_photons.pt[has_sublead][:, 1]
Reco_sublead_pho_eta = sorted_reco_photons.eta[has_sublead][:, 1]
Reco_sublead_pho_phi = sorted_reco_photons.phi[has_sublead][:, 1]


def dR(Pho_gen_eta, Pho_gen_phi, Pho_reco_eta, Pho_reco_phi):
    d_eta = Pho_gen_eta - Pho_reco_eta
    d_phi = Pho_gen_phi - Pho_reco_phi
    return np.sqrt(d_eta**2 + d_phi**2)

dR_pho_1 = dR(pho_from_a_eta_1, pho_from_a_phi_1, Reco_pho_eta_uf, Reco_pho_phi_uf)
dR_pho_2 = dR(pho_from_a_eta_2, pho_from_a_phi_2, Reco_pho_eta_uf, Reco_pho_phi_uf)

# Step 1: Index of photon with minimum dR per event
min_idx_1 = ak.argmin(dR_pho_1, axis=1)
min_idx_2 = ak.argmin(dR_pho_2, axis=1)

# Step 2: Build photon index array
photon_idx = ak.local_index(reco_photons)

# Step 3: Mask where photon index matches minimum dR index
mask_idx_1 = photon_idx == min_idx_1[:, None]   
mask_idx_2 = photon_idx == min_idx_2[:, None]

# Step 4: Also require dR < 0.1
mask_dR_1 = dR_pho_1 < 0.1
mask_dR_2 = dR_pho_2 < 0.1

# Step 5: Combine masks
mask_sel_1 = mask_idx_1 & mask_dR_1
mask_sel_2 = mask_idx_2 & mask_dR_2

# Step 6: Select photons satisfying both conditions
selected_photon_1 = ak.firsts(reco_photons[mask_sel_1])
selected_photon_2 = ak.firsts(reco_photons[mask_sel_2])

gen_photon_sync_1 = ak.to_numpy(pho_from_a_pt_1[~ak.is_none(selected_photon_1.pt)])
gen_photon_sync_2 = ak.to_numpy(pho_from_a_pt_2[~ak.is_none(selected_photon_2.pt)])


Genmatched_photon_1_clean = selected_photon_1[~ak.is_none(selected_photon_1)]
Genmatched_photon_2_clean = selected_photon_2[~ak.is_none(selected_photon_2)]

Genmatched_pho_1_pt  = ak.to_numpy(Genmatched_photon_1_clean.pt)
Genmatched_pho_2_pt  = ak.to_numpy(Genmatched_photon_2_clean.pt)
Genmatched_pho_1_eta = ak.to_numpy(Genmatched_photon_1_clean.eta)
Genmatched_pho_2_eta = ak.to_numpy(Genmatched_photon_2_clean.eta)
Genmatched_pho_1_phi = ak.to_numpy(Genmatched_photon_1_clean.phi)
Genmatched_pho_2_phi = ak.to_numpy(Genmatched_photon_2_clean.phi)

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\------Creating Hitogram Part---\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


branch_arrays = {
    "nPhoton": nPhoton,
    "higgs_pt": higgs_pt,
    "higgs_eta": higgs_eta,
    "higgs_phi": higgs_phi,

    # Prompt photons from a (unsorted)
    "pho_from_a_pt_1": pho_from_a_pt_1,
    "pho_from_a_pt_2": pho_from_a_pt_2,
    "pho_from_a_eta_1": pho_from_a_eta_1,
    "pho_from_a_eta_2": pho_from_a_eta_2,
    "pho_from_a_phi_1": pho_from_a_phi_1,
    "pho_from_a_phi_2": pho_from_a_phi_2,

    # Sorted (leading/subleading) photons from a
    "lead_pt_pho_gen": lead_pt_pho_gen,
    "sublead_pt_pho_gen": sublead_pt_pho_gen,
    "lead_eta_pho_gen": lead_eta_pho_gen,
    "sublead_eta_pho_gen": sublead_eta_pho_gen,
    "lead_phi_pho_gen": lead_phi_pho_gen,
    "sublead_phi_pho_gen": sublead_phi_pho_gen,

    # Prompt b-quarks from a (unsorted)
    "bquark_from_a_pt_1": bquark_from_a_pt_1,
    "bquark_from_a_pt_2": bquark_from_a_pt_2,
    "bquark_from_a_eta_1": bquark_from_a_eta_1,
    "bquark_from_a_eta_2": bquark_from_a_eta_2,
    "bquark_from_a_phi_1": bquark_from_a_phi_1,
    "bquark_from_a_phi_2": bquark_from_a_phi_2,

    # Sorted (leading/subleading) b-quarks from a
    "lead_pt_bquark_gen": lead_pt_bquark_gen,
    "sublead_pt_bquark_gen": sublead_pt_bquark_gen,
    "lead_eta_bquark_gen": lead_eta_bquark_gen,
    "sublead_eta_bquark_gen": sublead_eta_bquark_gen,
    "lead_phi_bquark_gen": lead_phi_bquark_gen,
    "sublead_phi_bquark_gen": sublead_phi_bquark_gen,


    # Reco photons
    "Reco_pho_pt":  Reco_pho_pt,
    "Reco_pho_eta": Reco_pho_eta,
    "Reco_pho_phi": Reco_pho_phi,

    # Reco b-quarks (temporarily using Photon info; replace if needed)
    "Reco_bquark_pt":  Reco_bquark_pt,
    "Reco_bquark_eta": Reco_bquark_eta,
    "Reco_bquark_phi": Reco_bquark_phi,

    "Reco_lead_pho_pt":  Reco_lead_pho_pt,
    "Reco_sublead_pho_pt": Reco_sublead_pho_pt,
    "Reco_lead_pho_eta": Reco_lead_pho_eta,
    "Reco_sublead_pho_eta": Reco_sublead_pho_eta,
    "Reco_lead_pho_phi": Reco_lead_pho_phi,
    "Reco_sublead_pho_phi": Reco_sublead_pho_phi,

    "gen_photon_from_a_1_pt": gen_photon_sync_1,
    "gen_photon_from_a_2_pt": gen_photon_sync_2,

    "Genmatched_pho_1_pt":  Genmatched_pho_1_pt,
    "Genmatched_pho_2_pt":  Genmatched_pho_2_pt,
    "Genmatched_pho_1_eta": Genmatched_pho_1_eta,
    "Genmatched_pho_2_eta": Genmatched_pho_2_eta,
    "Genmatched_pho_1_phi": Genmatched_pho_1_phi,
    "Genmatched_pho_2_phi": Genmatched_pho_2_phi


}


# bin_settings = {
#     "nPhoton": (100, 0, 10),

#     # Higgs
#     "higgs_pt":  (500, 0, 400),
#     "higgs_eta": (100, -10, 10),
#     "higgs_phi": (64, 1, -1),

#     # Unsorted photons from a
#     "pho_from_a_pt_1": (100, 0, 200),
#     "pho_from_a_pt_2": (100, 0, 200),
#     "pho_from_a_eta_1": (100, -10, 10),
#     "pho_from_a_eta_2": (100, -10, 10),
#     "pho_from_a_phi_1": (64, 1, -1),
#     "pho_from_a_phi_2": (64, 1, -1),

#     # Sorted photons from a
#     "lead_pt_pho_gen":     (100, 0, 200),
#     "sublead_pt_pho_gen":  (100, 0, 200),
#     "lead_eta_pho_gen":    (100, -10, 10),
#     "sublead_eta_pho_gen": (100, -10, 10),
#     "lead_phi_pho_gen":    (64, 1, -1),
#     "sublead_phi_pho_gen": (64, 1, -1),

#     # Unsorted b-quarks from a
#     "bquark_from_a_pt_1": (100, 0, 200),
#     "bquark_from_a_pt_2": (100, 0, 200),
#     "bquark_from_a_eta_1": (100, -10, 10),
#     "bquark_from_a_eta_2": (100, -10, 10),
#     "bquark_from_a_phi_1": (64, 1,-1),
#     "bquark_from_a_phi_2": (64, 1,-1),

#     # Sorted b-quarks from a
#     "lead_pt_bquark_gen":     (100, 0, 200),
#     "sublead_pt_bquark_gen":  (100, 0, 200),
#     "lead_eta_bquark_gen":    (100, -10, 10),
#     "sublead_eta_bquark_gen": (100, -10, 10),
#     "lead_phi_bquark_gen":    (64, 1,-1),
#     "sublead_phi_bquark_gen": (64, 1,-1),

#     "Reco_pho_pt":  (100, 0, 200),
#     "Reco_pho_eta": (100, 1, -1), 
#     "Reco_pho_phi": (64, 1, -1),   


#     "Reco_bquark_pt":  (100, 0, 200),
#     "Reco_bquark_eta": (100, 1, -1),  # reversed axis
#     "Reco_bquark_phi": (64, 1, -1),

#     "Reco_lead_pho_pt":     (100, 0, 200),
#     "Reco_sublead_pho_pt":  (100, 0, 200),
#     "Reco_lead_pho_eta":    (100, 1, -1),
#     "Reco_sublead_pho_eta": (100, 1, -1),
#     "Reco_lead_pho_phi":    (64, 1, -1),
#     "Reco_sublead_pho_phi": (64, 1, -1),

#     "gen_photon_from_a_1_pt": (100, 0, 100),
#     "gen_photon_from_a_2_pt": (100, 0, 100),

#     "Genmatched_pho_1_pt":  (100, 0, 100),
#     "Genmatched_pho_2_pt":  (100, 0, 100),
#     "Genmatched_pho_1_eta": (50, -3.0, 3.0),
#     "Genmatched_pho_2_eta": (50, -3.0, 3.0),
#     "Genmatched_pho_1_phi": (64, 1, -1),
#     "Genmatched_pho_2_phi": (64, 1, -1),
# }

bin_settings = {
    # "nPhoton": (100, 0, 10),

    # Higgs
    "higgs_pt":  (500, 0, 500),
    "higgs_eta": (500, -10, 10),
    "higgs_phi": (64, 1, -1),

    # Unsorted photons from a
    "pho_from_a_pt_1": (500, 0, 400),
    "pho_from_a_pt_2": (500, 0, 400),
    "pho_from_a_eta_1": (500, -10, 10),
    "pho_from_a_eta_2": (500, -10, 10),
    "pho_from_a_phi_1": (64, 1, -1),
    "pho_from_a_phi_2": (64, 1, -1),

    # Sorted photons from a
    "lead_pt_pho_gen":     (500, 0, 400),
    "sublead_pt_pho_gen":  (500, 0, 400),
    "lead_eta_pho_gen":    (500, -10, 10),
    "sublead_eta_pho_gen": (500, -10, 10),
    "lead_phi_pho_gen":    (64, 1, -1),
    "sublead_phi_pho_gen": (64, 1, -1),

    # Unsorted b-quarks from a
    "bquark_from_a_pt_1": (500, 0, 400),
    "bquark_from_a_pt_2": (500, 0, 400),
    "bquark_from_a_eta_1": (500, -10, 10),
    "bquark_from_a_eta_2": (500, -10, 10),
    "bquark_from_a_phi_1": (64, 1,-1),
    "bquark_from_a_phi_2": (64, 1,-1),

    # Sorted b-quarks from a
    "lead_pt_bquark_gen":     (500, 0, 400),
    "sublead_pt_bquark_gen":  (500, 0, 400),
    "lead_eta_bquark_gen":    (500, -10, 10),
    "sublead_eta_bquark_gen": (500, -10, 10),
    "lead_phi_bquark_gen":    (64, 1,-1),
    "sublead_phi_bquark_gen": (64, 1,-1),

    "Reco_pho_pt":  (500, 0, 400),
    "Reco_pho_eta": (500, 1, -1), 
    "Reco_pho_phi": (64, 1, -1),   


    "Reco_bquark_pt":  (500, 0, 400),
    "Reco_bquark_eta": (500, 1, -1),  # reversed axis
    "Reco_bquark_phi": (64, 1, -1),

    "Reco_lead_pho_pt":     (500, 0, 400),
    "Reco_sublead_pho_pt":  (500, 0, 400),
    "Reco_lead_pho_eta":    (500, 1, -1),
    "Reco_sublead_pho_eta": (500, 1, -1),
    "Reco_lead_pho_phi":    (64, 1, -1),
    "Reco_sublead_pho_phi": (64, 1, -1),

    "gen_photon_from_a_1_pt": (500, 0, 400),
    "gen_photon_from_a_2_pt": (500, 0, 400),

    "Genmatched_pho_1_pt":  (500, 0, 400),
    "Genmatched_pho_2_pt":  (500, 0, 400),
    "Genmatched_pho_1_eta": (500, -3.0, 3.0),
    "Genmatched_pho_2_eta": (500, -3.0, 3.0),
    "Genmatched_pho_1_phi": (64, 1, -1),
    "Genmatched_pho_2_phi": (64, 1, -1),
}


# ----------------------------------------------------------------------Plotting Part -------------------------------------------------------------------------------------------


# Mode options: "create", "update", "selected"
mode = "selected"

# If using 'selected' mode, specify which ones to process
selected_branches = [("gen_photon_from_a_1_pt", "Genmatched_pho_1_pt"), ("gen_photon_from_a_2_pt", "Genmatched_pho_2_pt"), "higgs_pt"]  # example

# Select which branches to process
if mode == "selected":
    branches_to_process = selected_branches
else:
    branches_to_process = list(branch_arrays.keys())

# Open ROOT file
file_mode = "UPDATE" if os.path.exists("hist_output.root") else "RECREATE"
f = ROOT.TFile("hist_output.root", file_mode)

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

    if is_2d:
        data_x = np.asarray(branch_arrays[b1])
        data_y = np.asarray(branch_arrays[b2])
        nbinsx, xmin, xmax = bin_settings.get(b1, (100, float(np.min(data_x)), float(np.max(data_x))))
        nbinsy, ymin, ymax = bin_settings.get(b2, (100, float(np.min(data_y)), float(np.max(data_y))))
        hist = ROOT.TH2F(hist_name, f"{b1} vs {b2};{b1};{b2}", nbinsx, xmin, xmax, nbinsy, ymin, ymax)
        for x, y in zip(data_x, data_y):
            hist.Fill(x, y)
    else:
        data = np.asarray(branch_arrays[branch])
        nbins, xmin, xmax = bin_settings.get(branch, (100, float(np.min(data)), float(np.max(data))))
        hist = ROOT.TH1F(hist_name, f"{branch};{branch};Events", nbins, xmin, xmax)
        for val in data:
            hist.Fill(val)

    hist.Write()

f.Close()
print(f"\nDone writing histograms in '{mode}' mode.")





