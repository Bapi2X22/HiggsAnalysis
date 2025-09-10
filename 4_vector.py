import uproot
import awkward as ak
import numpy as np
import vector

# Your dataset dictionary
# datasets = {
#     "M20_RunIISummer20UL18NanoAODv9": [
#         "root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/WHToAA_AToBB_AToGG_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2530000/FB6A4557-3B3E-5E4A-900D-45A77C107EA2.root",
#         "root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/WHToAA_AToBB_AToGG_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2530000/52967538-6671-C748-9CEC-C21D276D640B.root",
#     ],
#     "M40_RunIISummer20UL18NanoAODv9": [
#         "root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/WHToAA_AToBB_AToGG_M-40_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2520000/668672C9-F3A9-5647-836A-6A6CA9B5B66C.root",
#         "root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/WHToAA_AToBB_AToGG_M-40_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2520000/86D3D6A3-EE02-5849-8FC0-C89B9BE7D930.root",
#     ],
#     "M60_RunIISummer20UL18NanoAODv9": [
#         "root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/WHToAA_AToBB_AToGG_M-60_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2520000/5E5C0417-1992-D84B-9109-62ABCEA530BF.root",
#         "root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/WHToAA_AToBB_AToGG_M-60_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2520000/68E0B14A-B99F-F84C-B660-19688E15967B.root",
#     ],
# }

datasets = {
    "M20_Run3Summer22NanoAODv13": [
        "root://cmsxrootd.fnal.gov//store/mc/Run3Summer22NanoAODv13/HAHMHto2A_Ato2G_MA-20_TuneCP5_13p6TeV_madgraphMLM-pythia8/NANOAODSIM/133X_mcRun3_2022_realistic_ForNanov13_v1-v2/70000/4322332b-469d-40ae-babd-1a8db37b06e6.root",
        "root://cmsxrootd.fnal.gov//store/mc/Run3Summer22NanoAODv13/HAHMHto2A_Ato2G_MA-20_TuneCP5_13p6TeV_madgraphMLM-pythia8/NANOAODSIM/133X_mcRun3_2022_realistic_ForNanov13_v1-v2/70000/551dfa5b-97e5-461d-972f-2b49e7366c59.root",
        "root://cmsxrootd.fnal.gov//store/mc/Run3Summer22NanoAODv13/HAHMHto2A_Ato2G_MA-20_TuneCP5_13p6TeV_madgraphMLM-pythia8/NANOAODSIM/133X_mcRun3_2022_realistic_ForNanov13_v1-v2/70000/05ed9656-16b4-4897-a9bd-bd3508a11683.root",
        "root://cmsxrootd.fnal.gov//store/mc/Run3Summer22NanoAODv13/HAHMHto2A_Ato2G_MA-20_TuneCP5_13p6TeV_madgraphMLM-pythia8/NANOAODSIM/133X_mcRun3_2022_realistic_ForNanov13_v1-v2/70000/a816aaaa-fdd7-49e3-bc74-bf7c82e8e456.root",
        "root://cmsxrootd.fnal.gov//store/mc/Run3Summer22NanoAODv13/HAHMHto2A_Ato2G_MA-20_TuneCP5_13p6TeV_madgraphMLM-pythia8/NANOAODSIM/133X_mcRun3_2022_realistic_ForNanov13_v1-v2/70000/841a4f27-cf66-4c1a-bc7d-59b0832053ff.root",
        "root://cmsxrootd.fnal.gov//store/mc/Run3Summer22NanoAODv13/HAHMHto2A_Ato2G_MA-20_TuneCP5_13p6TeV_madgraphMLM-pythia8/NANOAODSIM/133X_mcRun3_2022_realistic_ForNanov13_v1-v2/70000/b0961b19-2f84-416e-842c-cd0c3f54b639.root",
    ],
    "M40_Run3Summer22NanoAODv13": [
        "root://cmsxrootd.fnal.gov//store/mc/Run3Summer22NanoAODv13/HAHMHto2A_Ato2G_MA-40_TuneCP5_13p6TeV_madgraphMLM-pythia8/NANOAODSIM/133X_mcRun3_2022_realistic_ForNanov13_v1-v2/70000/818051f1-9f9e-4146-b14e-360f250baea7.root",
        "root://cmsxrootd.fnal.gov//store/mc/Run3Summer22NanoAODv13/HAHMHto2A_Ato2G_MA-40_TuneCP5_13p6TeV_madgraphMLM-pythia8/NANOAODSIM/133X_mcRun3_2022_realistic_ForNanov13_v1-v2/70000/e15e0d2d-d8d8-46fe-a8e5-6a10b9b4b538.root",
        "root://cmsxrootd.fnal.gov//store/mc/Run3Summer22NanoAODv13/HAHMHto2A_Ato2G_MA-40_TuneCP5_13p6TeV_madgraphMLM-pythia8/NANOAODSIM/133X_mcRun3_2022_realistic_ForNanov13_v1-v2/70000/8f08287e-0cb6-4280-97ca-49565a6f44f0.root",
        "root://cmsxrootd.fnal.gov//store/mc/Run3Summer22NanoAODv13/HAHMHto2A_Ato2G_MA-40_TuneCP5_13p6TeV_madgraphMLM-pythia8/NANOAODSIM/133X_mcRun3_2022_realistic_ForNanov13_v1-v2/70000/6aed8e12-1118-43e1-9eb5-5e49e237a05f.root",
        "root://cmsxrootd.fnal.gov//store/mc/Run3Summer22NanoAODv13/HAHMHto2A_Ato2G_MA-40_TuneCP5_13p6TeV_madgraphMLM-pythia8/NANOAODSIM/133X_mcRun3_2022_realistic_ForNanov13_v1-v2/70000/56b8b46e-525f-4ede-8012-933087ae85c0.root",
        "root://cmsxrootd.fnal.gov//store/mc/Run3Summer22NanoAODv13/HAHMHto2A_Ato2G_MA-40_TuneCP5_13p6TeV_madgraphMLM-pythia8/NANOAODSIM/133X_mcRun3_2022_realistic_ForNanov13_v1-v2/70000/55eacf13-1c6d-4342-a74f-68c94ee5fe8f.root",
    ],
    "M60_Run3Summer22NanoAODv13": [
        "root://cmsxrootd.fnal.gov//store/mc/Run3Summer22NanoAODv13/HAHMHto2A_Ato2G_MA-60_TuneCP5_13p6TeV_madgraphMLM-pythia8/NANOAODSIM/133X_mcRun3_2022_realistic_ForNanov13_v1-v2/2820000/1487e377-20cd-4a73-97de-b5e019014dfd.root",
        "root://cmsxrootd.fnal.gov//store/mc/Run3Summer22NanoAODv13/HAHMHto2A_Ato2G_MA-60_TuneCP5_13p6TeV_madgraphMLM-pythia8/NANOAODSIM/133X_mcRun3_2022_realistic_ForNanov13_v1-v2/2820000/dc1ca030-1c14-4125-855c-ada98de8dd70.root",
        "root://cmsxrootd.fnal.gov//store/mc/Run3Summer22NanoAODv13/HAHMHto2A_Ato2G_MA-60_TuneCP5_13p6TeV_madgraphMLM-pythia8/NANOAODSIM/133X_mcRun3_2022_realistic_ForNanov13_v1-v2/2820000/ad1d6ece-874d-4733-b63d-b0f24a0c8afe.root",
        "root://cmsxrootd.fnal.gov//store/mc/Run3Summer22NanoAODv13/HAHMHto2A_Ato2G_MA-60_TuneCP5_13p6TeV_madgraphMLM-pythia8/NANOAODSIM/133X_mcRun3_2022_realistic_ForNanov13_v1-v2/2820000/71fa2532-d7cf-4008-9fb7-b2ef60102f2b.root",
        "root://cmsxrootd.fnal.gov//store/mc/Run3Summer22NanoAODv13/HAHMHto2A_Ato2G_MA-60_TuneCP5_13p6TeV_madgraphMLM-pythia8/NANOAODSIM/133X_mcRun3_2022_realistic_ForNanov13_v1-v2/2820000/67051b5b-3156-48a3-a0fd-b8a12b12fe31.root",
        "root://cmsxrootd.fnal.gov//store/mc/Run3Summer22NanoAODv13/HAHMHto2A_Ato2G_MA-60_TuneCP5_13p6TeV_madgraphMLM-pythia8/NANOAODSIM/133X_mcRun3_2022_realistic_ForNanov13_v1-v2/2820000/b77ced61-14a1-42aa-953a-b2e61d078fd7.root",
        "root://cmsxrootd.fnal.gov//store/mc/Run3Summer22NanoAODv13/HAHMHto2A_Ato2G_MA-60_TuneCP5_13p6TeV_madgraphMLM-pythia8/NANOAODSIM/133X_mcRun3_2022_realistic_ForNanov13_v1-v2/2820000/195a8b07-7523-4e4e-b74f-a5a91c4f097d.root",
        "root://cmsxrootd.fnal.gov//store/mc/Run3Summer22NanoAODv13/HAHMHto2A_Ato2G_MA-60_TuneCP5_13p6TeV_madgraphMLM-pythia8/NANOAODSIM/133X_mcRun3_2022_realistic_ForNanov13_v1-v2/2820000/8a4935df-0a75-41bb-a80a-a8156b0d4554.root",
        "root://cmsxrootd.fnal.gov//store/mc/Run3Summer22NanoAODv13/HAHMHto2A_Ato2G_MA-60_TuneCP5_13p6TeV_madgraphMLM-pythia8/NANOAODSIM/133X_mcRun3_2022_realistic_ForNanov13_v1-v2/2820000/fa0ce2d3-d469-4f1d-a29f-655c366cec33.root",
        "root://cmsxrootd.fnal.gov//store/mc/Run3Summer22NanoAODv13/HAHMHto2A_Ato2G_MA-60_TuneCP5_13p6TeV_madgraphMLM-pythia8/NANOAODSIM/133X_mcRun3_2022_realistic_ForNanov13_v1-v2/2820000/146ea23b-ccc5-448d-ae1f-4e3a28781cbd.root",
        "root://cmsxrootd.fnal.gov//store/mc/Run3Summer22NanoAODv13/HAHMHto2A_Ato2G_MA-60_TuneCP5_13p6TeV_madgraphMLM-pythia8/NANOAODSIM/133X_mcRun3_2022_realistic_ForNanov13_v1-v2/2820000/d08a086d-3718-4fd6-b3eb-76bda22edf48.root",
        "root://cmsxrootd.fnal.gov//store/mc/Run3Summer22NanoAODv13/HAHMHto2A_Ato2G_MA-60_TuneCP5_13p6TeV_madgraphMLM-pythia8/NANOAODSIM/133X_mcRun3_2022_realistic_ForNanov13_v1-v2/2820000/890e13a3-ef43-4247-a9c6-2aad3ebae3db.root",
        "root://cmsxrootd.fnal.gov//store/mc/Run3Summer22NanoAODv13/HAHMHto2A_Ato2G_MA-60_TuneCP5_13p6TeV_madgraphMLM-pythia8/NANOAODSIM/133X_mcRun3_2022_realistic_ForNanov13_v1-v2/2820000/6dfbfe63-9ec7-4ff3-b4b2-1817c114e19f.root",
    ],
}


#////////////////////////////////////////////////////////////////////////ggH Sample\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


def process_dataset(name, files):
    print(f"Processing {name}...")

    def to_numpy(arr):
        return ak.to_numpy(ak.fill_none(arr, np.nan))
    # --- load all Events ---
    events_list = []
    for fname in files:
        with uproot.open(fname, timeout=120) as Hfile:
            Tree = Hfile["Events"]
            events = Tree.arrays(library="ak", how="zip")
            events_list.append(events)
    Events = ak.concatenate(events_list)

    gen = Events.GenPart

    # --- Higgs ---
    higgs_mask = (gen.pdgId == 25) & (gen.status == 62)
    higgs = gen[higgs_mask]
    H_pt   = to_numpy(ak.pad_none(higgs.pt, 1)[:, 0])
    H_eta  = to_numpy(ak.pad_none(higgs.eta, 1)[:, 0])
    H_phi  = to_numpy(ak.pad_none(higgs.phi, 1)[:, 0])
    H_mass = to_numpy(ak.pad_none(higgs.mass, 1)[:, 0])

    # --- A bosons (2 per event) ---
    A = gen[abs(gen.pdgId) == 35]
    A_pt   = ak.pad_none(A.pt, 2, axis=1)
    A_eta  = ak.pad_none(A.eta, 2, axis=1)
    A_phi  = ak.pad_none(A.phi, 2, axis=1)
    A_mass = ak.pad_none(A.mass, 2, axis=1)

    A1_pt   = to_numpy(A_pt[:, 0])
    A2_pt   = to_numpy(A_pt[:, 1])
    A1_eta  = to_numpy(A_eta[:, 0])
    A2_eta  = to_numpy(A_eta[:, 1])
    A1_phi  = to_numpy(A_phi[:, 0])
    A2_phi  = to_numpy(A_phi[:, 1])
    A1_mass = to_numpy(A_mass[:, 0])
    A2_mass = to_numpy(A_mass[:, 1])

    # --- photons from A (4 per event) ---
    photons = gen[(gen.pdgId == 22) & (gen.status == 1)]
    mother_idx = photons.genPartIdxMother
    from_a_mask = gen[mother_idx].pdgId == 35
    photons_from_a = photons[from_a_mask]

    pho_pt  = ak.pad_none(photons_from_a.pt, 4, axis=1)
    pho_eta = ak.pad_none(photons_from_a.eta, 4, axis=1)
    pho_phi = ak.pad_none(photons_from_a.phi, 4, axis=1)

    pho1_pt  = to_numpy(pho_pt[:, 0])
    pho2_pt  = to_numpy(pho_pt[:, 1])
    pho3_pt  = to_numpy(pho_pt[:, 2])
    pho4_pt  = to_numpy(pho_pt[:, 3])

    pho1_eta = to_numpy(pho_eta[:, 0])
    pho2_eta = to_numpy(pho_eta[:, 1])
    pho3_eta = to_numpy(pho_eta[:, 2])
    pho4_eta = to_numpy(pho_eta[:, 3])

    pho1_phi = to_numpy(pho_phi[:, 0])
    pho2_phi = to_numpy(pho_phi[:, 1])
    pho3_phi = to_numpy(pho_phi[:, 2])
    pho4_phi = to_numpy(pho_phi[:, 3])

    pho1_mass = np.zeros_like(pho1_pt)
    pho2_mass = np.zeros_like(pho2_pt)
    pho3_mass = np.zeros_like(pho3_pt)
    pho4_mass = np.zeros_like(pho4_pt)

    branches = {
    "H_pt": H_pt, "H_eta": H_eta, "H_phi": H_phi, "H_mass": H_mass,
    "A1_pt": A1_pt, "A1_eta": A1_eta, "A1_phi": A1_phi, "A1_mass": A1_mass,
    "A2_pt": A2_pt, "A2_eta": A2_eta, "A2_phi": A2_phi, "A2_mass": A2_mass,
    "pho1_pt": pho1_pt, "pho1_eta": pho1_eta, "pho1_phi": pho1_phi, "pho1_mass": pho1_mass,
    "pho2_pt": pho2_pt, "pho2_eta": pho2_eta, "pho2_phi": pho2_phi, "pho2_mass": pho2_mass,
    "pho3_pt": pho3_pt, "pho3_eta": pho3_eta, "pho3_phi": pho3_phi, "pho3_mass": pho3_mass,
    "pho4_pt": pho4_pt, "pho4_eta": pho4_eta, "pho4_phi": pho4_phi, "pho4_mass": pho4_mass,
}

    outname = f"ggH_{name}_genlevel.root"
    with uproot.recreate(outname) as fout:
        fout["Events"] = branches

    print(f"Saved {outname}")


# Loop over datasets
for name, files in datasets.items():
    process_dataset(name, files)


#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\WH sample\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

# # Fixed PDG b-quark mass
# m_b_pdg = 4.18

# def process_dataset(name, files):
#     print(f"Processing {name}...")

#     # --- load all Events from all files ---
#     events_list = []
#     for fname in files:
#         with uproot.open(fname, timeout=120) as Hfile:
#             Tree = Hfile["Events"]
#             events = Tree.arrays(library="ak", how="zip")
#             events_list.append(events)
#     Events = ak.concatenate(events_list)

#     # --- same as before, now using Events ---
#     gen = Events.GenPart

#     higgs_mask = (gen.pdgId == 25) & (gen.status == 62)
#     higgs = gen[higgs_mask]

#     higgs_pt = ak.flatten(higgs.pt)
#     higgs_eta = ak.flatten(higgs.eta)
#     higgs_phi = ak.flatten(higgs.phi)

#     is_A = (abs(gen.pdgId) == 35)
#     A = gen[is_A]

#     A_pt = A.pt
#     A_eta = A.eta
#     A_phi = A.phi

#     A1_pt = A_pt[:, 0]
#     A2_pt = A_pt[:, 1]
#     A1_eta = A_eta[:, 0]
#     A2_eta = A_eta[:, 1]
#     A1_phi = A_phi[:, 0]
#     A2_phi = A_phi[:, 1]

#     photons = gen[(gen.pdgId == 22) & (gen.status == 1)]
#     mother_idx = photons.genPartIdxMother
#     from_a_mask = gen[mother_idx].pdgId == 35
#     photons_from_a = photons[from_a_mask]

#     # Pad to at least 2 photons per event
#     photons_from_a_padded = ak.pad_none(photons_from_a, 2, axis=1, clip=True)

#     pho_from_a_pt  = photons_from_a_padded.pt
#     pho_from_a_eta = photons_from_a_padded.eta
#     pho_from_a_phi = photons_from_a_padded.phi

#     pho_from_a_pt_1  = ak.fill_none(pho_from_a_pt[:, 0], np.nan)
#     pho_from_a_pt_2  = ak.fill_none(pho_from_a_pt[:, 1], np.nan)
#     pho_from_a_eta_1 = ak.fill_none(pho_from_a_eta[:, 0], np.nan)
#     pho_from_a_eta_2 = ak.fill_none(pho_from_a_eta[:, 1], np.nan)
#     pho_from_a_phi_1 = ak.fill_none(pho_from_a_phi[:, 0], np.nan)
#     pho_from_a_phi_2 = ak.fill_none(pho_from_a_phi[:, 1], np.nan)

#     bquarks = gen[abs(gen.pdgId) == 5]
#     mother_idx = bquarks.genPartIdxMother
#     from_a_mask = gen[mother_idx].pdgId == 35
#     bquarks_from_a = bquarks[from_a_mask]

#     bquark_from_a_pt = bquarks_from_a.pt
#     bquark_from_a_eta = bquarks_from_a.eta
#     bquark_from_a_phi = bquarks_from_a.phi

#     bquark_from_a_pt_1 = bquark_from_a_pt[:, 0]
#     bquark_from_a_pt_2 = bquark_from_a_pt[:, 1]
#     bquark_from_a_eta_1 = bquark_from_a_eta[:, 0]
#     bquark_from_a_eta_2 = bquark_from_a_eta[:, 1]
#     bquark_from_a_phi_1 = bquark_from_a_phi[:, 0]
#     bquark_from_a_phi_2 = bquark_from_a_phi[:, 1]

#     def to_numpy(arr): return ak.to_numpy(ak.fill_none(arr, np.nan))

#     branches = {
#         "H_pt": to_numpy(higgs_pt),
#         "H_eta": to_numpy(higgs_eta),
#         "H_phi": to_numpy(higgs_phi),
#         "H_mass": to_numpy(ak.flatten(higgs.mass)),

#         "A1_pt": to_numpy(A1_pt),
#         "A1_eta": to_numpy(A1_eta),
#         "A1_phi": to_numpy(A1_phi),
#         "A1_mass": to_numpy(A.mass[:,0]),

#         "A2_pt": to_numpy(A2_pt),
#         "A2_eta": to_numpy(A2_eta),
#         "A2_phi": to_numpy(A2_phi),
#         "A2_mass": to_numpy(A.mass[:,1]),

#         "pho1_pt": to_numpy(pho_from_a_pt_1),
#         "pho1_eta": to_numpy(pho_from_a_eta_1),
#         "pho1_phi": to_numpy(pho_from_a_phi_1),
#         "pho1_mass": np.zeros_like(to_numpy(pho_from_a_pt_1)),

#         "pho2_pt": to_numpy(pho_from_a_pt_2),
#         "pho2_eta": to_numpy(pho_from_a_eta_2),
#         "pho2_phi": to_numpy(pho_from_a_phi_2),
#         "pho2_mass": np.zeros_like(to_numpy(pho_from_a_pt_2)),

#         "b1_pt": to_numpy(bquark_from_a_pt_1),
#         "b1_eta": to_numpy(bquark_from_a_eta_1),
#         "b1_phi": to_numpy(bquark_from_a_phi_1),
#         "b1_mass": np.full_like(to_numpy(bquark_from_a_pt_1), m_b_pdg),

#         "b2_pt": to_numpy(bquark_from_a_pt_2),
#         "b2_eta": to_numpy(bquark_from_a_eta_2),
#         "b2_phi": to_numpy(bquark_from_a_phi_2),
#         "b2_mass": np.full_like(to_numpy(bquark_from_a_pt_2), m_b_pdg),
#     }

#     outname = f"{name}_genlevel.root"
#     with uproot.recreate(outname) as fout:
#         fout["Events"] = branches

#     print(f"Saved {outname}")


# # Loop over datasets
# for name, files in datasets.items():
#     process_dataset(name, files)