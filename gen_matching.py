import awkward as ak 
import numpy as np
import matplotlib.pyplot as plt
import ROOT

Events = ak.from_parquet("Events.parquet")

gen = Events.GenPart

# Step 1: Select all prompt photons (pdgId == 22)
photons = gen[gen.pdgId == 22]

print(photons.pt[0])

# Step 2: Get mother indices for each photon
mother_idx = photons.genPartIdxMother  # shape: (n_events, variable)

print(mother_idx[0])

# Step 3: Check if mother's pdgId == 35
from_a_mask = gen[mother_idx].pdgId == 35  # shape: (n_events, variable)

print(from_a_mask[0])

# Step 4: Apply mask to get only photons from a
photons_from_a = photons[from_a_mask]

pho_from_a_pt = photons_from_a.pt
pho_from_a_eta = photons_from_a.eta 
pho_from_a_phi = photons_from_a.phi


pho_from_a_pt_1 = pho_from_a_pt[:,0]
pho_from_a_pt_2 = pho_from_a_pt[:,1]
pho_from_a_eta_1 = pho_from_a_eta[:,0]  
pho_from_a_eta_2 = pho_from_a_eta[:,1]
pho_from_a_phi_1 = pho_from_a_phi[:,0]
pho_from_a_phi_2 = pho_from_a_phi[:,1]

Reco_pho_pt = Events.Photon.pt
Reco_pho_eta = Events.Photon.eta
Reco_pho_phi = Events.Photon.phi

print("Reco_pho_pt:\n", Reco_pho_pt)
print("Reco_pho_eta:\n", Reco_pho_eta)
print("Reco_pho_phi:\n", Reco_pho_phi)



def dR(Pho_gen_eta, Pho_gen_phi, Pho_reco_eta, Pho_reco_phi):
    """
    Calculate the delta R between photons and b quarks.
    """
    d_eta = Pho_gen_eta - Pho_reco_eta
    d_phi = Pho_gen_phi - Pho_reco_phi
    
    return np.sqrt(d_eta**2 + d_phi**2)

dR_pho_1 = dR(pho_from_a_eta_1, pho_from_a_phi_1, Reco_pho_eta, Reco_pho_phi)
dR_pho_2 = dR(pho_from_a_eta_2, pho_from_a_phi_2, Reco_pho_eta, Reco_pho_phi)

print("dR_pho_1",dR_pho_1)
print("dR_pho_2",dR_pho_2)

# Step 1: Index of photon with minimum dR per event
min_idx_1 = ak.argmin(dR_pho_1, axis=1)
min_idx_2 = ak.argmin(dR_pho_2, axis=1)

print("min_idx_1" ,min_idx_1)
print("min_idx_2" ,min_idx_2)

# Step 2: Build photon index array
photon_idx = ak.local_index(Events.Photon)

print("photon_idx", photon_idx)

# # Step 3: Mask where photon index matches minimum dR index
# mask_idx_1 = photon_idx == min_idx_1[:, None]
# mask_idx_2 = photon_idx == min_idx_2[:, None]


# # Step 4: Also require dR < 0.1
# mask_dR_1 = dR_pho_1 < 0.1
# mask_dR_2 = dR_pho_2 < 0.1

# # Step 5: Combine masks
# mask_sel_1 = mask_idx_1 & mask_dR_1
# mask_sel_2 = mask_idx_2 & mask_dR_2

# # Step 6: Select photons satisfying both conditions
# selected_photon_1 = ak.firsts(Events.Photon[mask_sel_1])
# selected_photon_2 = ak.firsts(Events.Photon[mask_sel_2])

# Step 3: Mask where photon index matches minimum dR index
mask_idx_1 = photon_idx == min_idx_1[:, None]   
mask_idx_2 = photon_idx == min_idx_2[:, None]
print("mask_idx_1:\n", mask_idx_1)
print("mask_idx_2:\n", mask_idx_2)

# Step 4: Also require dR < 0.1
mask_dR_1 = dR_pho_1 < 0.1
mask_dR_2 = dR_pho_2 < 0.1
print("mask_dR_1:\n", mask_dR_1)
print("mask_dR_2:\n", mask_dR_2)

# Step 5: Combine masks
mask_sel_1 = mask_idx_1 & mask_dR_1
mask_sel_2 = mask_idx_2 & mask_dR_2
print("mask_sel_1:\n", mask_sel_1)
print("mask_sel_2:\n", mask_sel_2)

# Step 6: Select photons satisfying both conditions
selected_photon_1 = ak.firsts(Events.Photon[mask_sel_1])
selected_photon_2 = ak.firsts(Events.Photon[mask_sel_2])
print("selected_photon_1:\n", selected_photon_1)
print("selected_photon_2:\n", selected_photon_2)

print("selected_photon_1.pt:\n", selected_photon_1.pt)
print("selected_photon_2.pt:\n", selected_photon_2.pt)


selected_pho_pt_1 = selected_photon_1.pt
selected_pho_pt_2 = selected_photon_2.pt
selected_pho_eta_1 = selected_photon_1.eta
selected_pho_eta_2 = selected_photon_2.eta
selected_pho_phi_1 = selected_photon_1.phi
selected_pho_phi_1 = selected_photon_1.phi

# Convert awkward arrays to flat numpy arrays, removing invalid entries (e.g., missing matches)
gen1 = ak.to_numpy(pho_from_a_pt_1[~ak.is_none(selected_pho_pt_1)])
sel1 = ak.to_numpy(selected_pho_pt_1[~ak.is_none(selected_pho_pt_1)])

gen2 = ak.to_numpy(pho_from_a_pt_2[~ak.is_none(selected_pho_pt_2)])
sel2 = ak.to_numpy(selected_pho_pt_2[~ak.is_none(selected_pho_pt_2)])

# Create histograms
h2_gen1_vs_sel1 = ROOT.TH2F("h2_gen1_vs_sel1", "Gen1 vs Sel1; Gen pT [GeV]; Reco pT [GeV]", 
                            100, 0, 100, 100, 0, 100)

h2_gen2_vs_sel2 = ROOT.TH2F("h2_gen2_vs_sel2", "Gen2 vs Sel2; Gen pT [GeV]; Reco pT [GeV]", 
                            100, 0, 100, 100, 0, 100)

# Fill histograms
for g, s in zip(gen1, sel1):
    h2_gen1_vs_sel1.Fill(g, s)

for g, s in zip(gen2, sel2):
    h2_gen2_vs_sel2.Fill(g, s)

# Draw histograms
canvas1 = ROOT.TCanvas("c1", "Gen1 vs Sel1", 800, 600)
h2_gen1_vs_sel1.Draw("COLZ")
canvas1.SaveAs("gen1_vs_sel1.png")

canvas2 = ROOT.TCanvas("c2", "Gen2 vs Sel2", 800, 600)
h2_gen2_vs_sel2.Draw("COLZ")
canvas2.SaveAs("gen2_vs_sel2.png")