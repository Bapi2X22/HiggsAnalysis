import awkward as ak
import numpy as np
from coffea.processor import ProcessorABC, dict_accumulator
# from coffea.processor.accumulator import column_accumulator
from coffea.processor.accumulator import column_accumulator, dict_accumulator
from ROOT import TLorentzVector


class HiggsAnalysisProcessor(ProcessorABC):
    def __init__(self):
        # Wrap your accumulators inside dict_accumulator
        self._accumulator = dict_accumulator({
            # "nPhoton": column_accumulator(np.array([], dtype=np.float32)),
            "pu_true": column_accumulator(np.array([], dtype=np.float32)),
            "nPhoton": column_accumulator(np.array([], dtype=np.float32)),
            "higgs_pt": column_accumulator(np.array([], dtype=np.float32)),
            "higgs_eta": column_accumulator(np.array([], dtype=np.float32)),
            "higgs_phi": column_accumulator(np.array([], dtype=np.float32)),

            "A_pt_1": column_accumulator(np.array([], dtype=np.float32)),
            "A_pt_2": column_accumulator(np.array([], dtype=np.float32)),
            "A_eta_1": column_accumulator(np.array([], dtype=np.float32)),
            "A_eta_2": column_accumulator(np.array([], dtype=np.float32)),
            "A_phi_1": column_accumulator(np.array([], dtype=np.float32)),
            "A_phi_2": column_accumulator(np.array([], dtype=np.float32)),

            # Sorted A's
            "leading_A_pt": column_accumulator(np.array([], dtype=np.float32)),
            "subleading_A_pt": column_accumulator(np.array([], dtype=np.float32)),
            "leading_A_eta": column_accumulator(np.array([], dtype=np.float32)),
            "subleading_A_eta": column_accumulator(np.array([], dtype=np.float32)),
            "leading_A_phi": column_accumulator(np.array([], dtype=np.float32)),
            "subleading_A_phi": column_accumulator(np.array([], dtype=np.float32)),

            # Prompt photons from a (unsorted)
            "pho_from_a_pt_1": column_accumulator(np.array([], dtype=np.float32)),
            "pho_from_a_pt_2": column_accumulator(np.array([], dtype=np.float32)),
            "pho_from_a_eta_1": column_accumulator(np.array([], dtype=np.float32)),
            "pho_from_a_eta_2": column_accumulator(np.array([], dtype=np.float32)),
            "pho_from_a_phi_1": column_accumulator(np.array([], dtype=np.float32)),
            "pho_from_a_phi_2": column_accumulator(np.array([], dtype=np.float32)),

            # Sorted (leading/subleading) photons from a
            "lead_pt_pho_gen": column_accumulator(np.array([], dtype=np.float32)),
            "sublead_pt_pho_gen": column_accumulator(np.array([], dtype=np.float32)),
            "lead_eta_pho_gen": column_accumulator(np.array([], dtype=np.float32)),
            "sublead_eta_pho_gen": column_accumulator(np.array([], dtype=np.float32)),
            "lead_phi_pho_gen": column_accumulator(np.array([], dtype=np.float32)),
            "sublead_phi_pho_gen": column_accumulator(np.array([], dtype=np.float32)),

            # Prompt b-quarks from a (unsorted)
            "bquark_from_a_pt_1": column_accumulator(np.array([], dtype=np.float32)),
            "bquark_from_a_pt_2": column_accumulator(np.array([], dtype=np.float32)),
            "bquark_from_a_eta_1": column_accumulator(np.array([], dtype=np.float32)),
            "bquark_from_a_eta_2": column_accumulator(np.array([], dtype=np.float32)),
            "bquark_from_a_phi_1": column_accumulator(np.array([], dtype=np.float32)),
            "bquark_from_a_phi_2": column_accumulator(np.array([], dtype=np.float32)),

            # Sorted (leading/subleading) b-quarks from a
            "lead_pt_bquark_gen": column_accumulator(np.array([], dtype=np.float32)),
            "sublead_pt_bquark_gen": column_accumulator(np.array([], dtype=np.float32)),
            "lead_eta_bquark_gen": column_accumulator(np.array([], dtype=np.float32)),
            "sublead_eta_bquark_gen": column_accumulator(np.array([], dtype=np.float32)),
            "lead_phi_bquark_gen": column_accumulator(np.array([], dtype=np.float32)),
            "sublead_phi_bquark_gen": column_accumulator(np.array([], dtype=np.float32)),

            # Reco photons
            "Reco_pho_pt": column_accumulator(np.array([], dtype=np.float32)),
            "Reco_pho_eta": column_accumulator(np.array([], dtype=np.float32)),
            "Reco_pho_phi": column_accumulator(np.array([], dtype=np.float32)),

            "Reco_photon_lead_pt_all": column_accumulator(np.array([], dtype=np.float32)),
            "Reco_photon_sublead_pt_all": column_accumulator(np.array([], dtype=np.float32)),
            "Reco_photon_lead_eta_all": column_accumulator(np.array([], dtype=np.float32)),
            "Reco_photon_sublead_eta_all": column_accumulator(np.array([], dtype=np.float32)),
            "Reco_photon_lead_phi_all": column_accumulator(np.array([], dtype=np.float32)),
            "Reco_photon_sublead_phi_all": column_accumulator(np.array([], dtype=np.float32)),


            "Reco_lead_pho_pt": column_accumulator(np.array([], dtype=np.float32)),
            "Reco_sublead_pho_pt": column_accumulator(np.array([], dtype=np.float32)),
            "Reco_lead_pho_eta": column_accumulator(np.array([], dtype=np.float32)),
            "Reco_sublead_pho_eta": column_accumulator(np.array([], dtype=np.float32)),
            "Reco_lead_pho_phi": column_accumulator(np.array([], dtype=np.float32)),
            "Reco_sublead_pho_phi": column_accumulator(np.array([], dtype=np.float32)),

            "gen_photon_from_a_1_pt": column_accumulator(np.array([], dtype=np.float32)),
            "gen_photon_from_a_2_pt": column_accumulator(np.array([], dtype=np.float32)),

            "Gen_photon_pt_1": column_accumulator(np.array([], dtype=np.float32)),
            "Gen_photon_pt_2": column_accumulator(np.array([], dtype=np.float32)),

            "Genmatched_pho_1_pt": column_accumulator(np.array([], dtype=np.float32)),
            "Genmatched_pho_2_pt": column_accumulator(np.array([], dtype=np.float32)),
            "Genmatched_pho_1_eta": column_accumulator(np.array([], dtype=np.float32)),
            "Genmatched_pho_2_eta": column_accumulator(np.array([], dtype=np.float32)),
            "Genmatched_pho_1_phi": column_accumulator(np.array([], dtype=np.float32)),
            "Genmatched_pho_2_phi": column_accumulator(np.array([], dtype=np.float32)),

                # Gen-level b-quarks from A (pt-sorted)
            "gen_lead_b_pt": column_accumulator(np.array([], dtype=np.float32)),
            "gen_lead_b_eta": column_accumulator(np.array([], dtype=np.float32)),
            "gen_lead_b_phi": column_accumulator(np.array([], dtype=np.float32)),

            "gen_sublead_b_pt": column_accumulator(np.array([], dtype=np.float32)),
            "gen_sublead_b_eta": column_accumulator(np.array([], dtype=np.float32)),
            "gen_sublead_b_phi": column_accumulator(np.array([], dtype=np.float32)),

            # Gen b-quarks unordered
            "gen_b1_pt": column_accumulator(np.array([], dtype=np.float32)),
            "gen_b1_eta": column_accumulator(np.array([], dtype=np.float32)),
            "gen_b1_phi": column_accumulator(np.array([], dtype=np.float32)),

            "gen_b2_pt": column_accumulator(np.array([], dtype=np.float32)),
            "gen_b2_eta": column_accumulator(np.array([], dtype=np.float32)),
            "gen_b2_phi": column_accumulator(np.array([], dtype=np.float32)),

            # Reco b-jets (pt-sorted)
            "reco_lead_bjet_pt": column_accumulator(np.array([], dtype=np.float32)),
            "reco_lead_bjet_eta": column_accumulator(np.array([], dtype=np.float32)),
            "reco_lead_bjet_phi": column_accumulator(np.array([], dtype=np.float32)),

            "reco_sublead_bjet_pt": column_accumulator(np.array([], dtype=np.float32)),
            "reco_sublead_bjet_eta": column_accumulator(np.array([], dtype=np.float32)),
            "reco_sublead_bjet_phi": column_accumulator(np.array([], dtype=np.float32)),

            # Reco b-jets matched to gen b-quarks
            "matched_bjet1_pt": column_accumulator(np.array([], dtype=np.float32)),
            "matched_bjet1_eta": column_accumulator(np.array([], dtype=np.float32)),
            "matched_bjet1_phi": column_accumulator(np.array([], dtype=np.float32)),

            "matched_bjet2_pt": column_accumulator(np.array([], dtype=np.float32)),
            "matched_bjet2_eta": column_accumulator(np.array([], dtype=np.float32)),
            "matched_bjet2_phi": column_accumulator(np.array([], dtype=np.float32)),
            "gen_invmasses_diphoton": column_accumulator(np.array([], dtype=np.float32)),
            "gen_invmasses_bb": column_accumulator(np.array([], dtype=np.float32)),
            "invmasses_diphoton": column_accumulator(np.array([], dtype=np.float32)),
            "invmasses_bb": column_accumulator(np.array([], dtype=np.float32)),
        })

    def accumulator(self):
        # Return the accumulator prototype for coffea
        return self._accumulator

    def process(self, events):

        Events = events

        dataset = Events.metadata["dataset"]  # get dataset name dynamically

        output = self.accumulator()
        # nPhoton = Events["nPhoton"]
        
        pu_true = Events.Pileup.nTrueInt
        # nPhoton = Events.Photon.nPhoton
        nPhoton = ak.num(events.Photon)

        # nPhoton = Events["nPhoton"]

        gen = Events.GenPart
        higgs_mask = (gen.pdgId == 25) & (gen.status == 62)
        higgs = gen[higgs_mask]

        higgs_pt = ak.flatten(higgs.pt)
        higgs_eta = ak.flatten(higgs.eta)
        higgs_phi = ak.flatten(higgs.phi)

        is_A = (abs(gen.pdgId) == 35)

        A = gen[is_A]

        A_pt = A.pt
        A_eta = A.eta
        A_phi = A.phi

        A_pt_1 = A_pt[:,0]
        A_pt_2 = A_pt[:,1]
        A_eta_1 = A_eta[:,0]
        A_eta_2 = A_eta[:,1]
        A_phi_1 = A_phi[:,0]
        A_phi_2 = A_phi[:,1]

        sorted_As = A[ak.argsort(A_pt, axis=1, ascending=False)]

        leading_A_pt = sorted_As.pt[:,0]
        subleading_A_pt = sorted_As.pt[:,1]
        leading_A_eta = sorted_As.eta[:,0]
        subleading_A_eta = sorted_As.eta[:,1]
        leading_A_phi = sorted_As.phi[:,0]
        subleading_A_phi = sorted_As.phi[:,1]

        photons = gen[(gen.pdgId == 22) & (gen.status == 1)]
        mother_idx = photons.genPartIdxMother
        from_a_mask = gen[mother_idx].pdgId == 35
        photons_from_a = photons[from_a_mask]

        mask_two_photons = ak.num(photons_from_a.eta) >= 2
        photons_from_a = photons_from_a[mask_two_photons]

        pho_from_a_pt = photons_from_a.pt
        pho_from_a_eta = photons_from_a.eta
        pho_from_a_phi = photons_from_a.phi

        pho_from_a_pt_1 = pho_from_a_pt[:, 0]
        pho_from_a_pt_2 = pho_from_a_pt[:, 1]
        pho_from_a_eta_1 = pho_from_a_eta[:, 0]
        pho_from_a_eta_2 = pho_from_a_eta[:, 1]
        pho_from_a_phi_1 = pho_from_a_phi[:, 0]
        pho_from_a_phi_2 = pho_from_a_phi[:, 1]

        sorted_photons = photons_from_a[ak.argsort(photons_from_a.pt, axis=1, ascending=False)]
        lead_pt_pho_gen = sorted_photons.pt[:, 0]
        sublead_pt_pho_gen = sorted_photons.pt[:, 1]
        lead_eta_pho_gen = sorted_photons.eta[:, 0]
        sublead_eta_pho_gen = sorted_photons.eta[:, 1]
        lead_phi_pho_gen = sorted_photons.phi[:, 0]
        sublead_phi_pho_gen = sorted_photons.phi[:, 1]

        bquarks = gen[abs(gen.pdgId) == 5]
        mother_idx = bquarks.genPartIdxMother
        from_a_mask = gen[mother_idx].pdgId == 35
        bquarks_from_a = bquarks[from_a_mask]

        bquark_from_a_pt = bquarks_from_a.pt
        bquark_from_a_eta = bquarks_from_a.eta
        bquark_from_a_phi = bquarks_from_a.phi

        bquark_from_a_pt_1 = bquark_from_a_pt[:, 0]
        bquark_from_a_pt_2 = bquark_from_a_pt[:, 1]
        bquark_from_a_eta_1 = bquark_from_a_eta[:, 0]
        bquark_from_a_eta_2 = bquark_from_a_eta[:, 1]
        bquark_from_a_phi_1 = bquark_from_a_phi[:, 0]
        bquark_from_a_phi_2 = bquark_from_a_phi[:, 1]

        sorted_bquarks = bquarks_from_a[ak.argsort(bquarks_from_a.pt, axis=1, ascending=False)]
        lead_pt_bquark_gen = sorted_bquarks.pt[:, 0]
        sublead_pt_bquark_gen = sorted_bquarks.pt[:, 1]
        lead_eta_bquark_gen = sorted_bquarks.eta[:, 0]
        sublead_eta_bquark_gen = sorted_bquarks.eta[:, 1]
        lead_phi_bquark_gen = sorted_bquarks.phi[:, 0]
        sublead_phi_bquark_gen = sorted_bquarks.phi[:, 1]

        reco_photons = Events.Photon[mask_two_photons]
        has_photon = ak.num(reco_photons) >= 1
        Reco_photons_all = reco_photons[has_photon]
        Reco_pho_pt = ak.flatten(reco_photons.pt)
        Reco_pho_eta = ak.flatten(reco_photons.eta)
        Reco_pho_phi = ak.flatten(reco_photons.phi)

        Reco_pho_pt_uf = reco_photons.pt
        Reco_pho_eta_uf = reco_photons.eta
        Reco_pho_phi_uf = reco_photons.phi

        sorted_reco_photons = reco_photons[ak.argsort(reco_photons.pt, axis=1, ascending=False)]
        Sorted_reco_photons_all = Reco_photons_all[ak.argsort(Reco_photons_all.pt, axis=1, ascending=False)]

        Reco_photon_lead_pt_all = Sorted_reco_photons_all.pt[:,0]
        Reco_photon_sublead_pt_all = ak.pad_none(Sorted_reco_photons_all.pt, 2)[:, 1]
        Reco_photon_lead_eta_all = Sorted_reco_photons_all.eta[:,0]
        Reco_photon_sublead_eta_all = ak.pad_none(Sorted_reco_photons_all.eta, 2)[:, 1]
        Reco_photon_lead_phi_all = Sorted_reco_photons_all.phi[:,0]
        Reco_photon_sublead_phi_all = ak.pad_none(Sorted_reco_photons_all.phi, 2)[:, 1]
        # has_lead = ak.num(sorted_reco_photons.pt) > 0
        # Reco_lead_pho_pt = sorted_reco_photons.pt[has_lead][:, 0]
        # Reco_lead_pho_eta = sorted_reco_photons.eta[has_lead][:, 0]
        # Reco_lead_pho_phi = sorted_reco_photons.phi[has_lead][:, 0]

        # has_sublead = ak.num(sorted_reco_photons.pt) > 1
        # Reco_sublead_pho_pt = sorted_reco_photons.pt[has_sublead][:, 1]
        # Reco_sublead_pho_eta = sorted_reco_photons.eta[has_sublead][:, 1]
        # Reco_sublead_pho_phi = sorted_reco_photons.phi[has_sublead][:, 1]
        # Ensure both leading and subleading photons exist
        has_two_photons = ak.num(sorted_reco_photons.pt) > 1

        # Use the same mask for both
        Reco_lead_pho_pt     = sorted_reco_photons.pt[has_two_photons][:, 0]
        Reco_sublead_pho_pt  = sorted_reco_photons.pt[has_two_photons][:, 1]
        Reco_lead_pho_eta    = sorted_reco_photons.eta[has_two_photons][:, 0]
        Reco_sublead_pho_eta = sorted_reco_photons.eta[has_two_photons][:, 1]
        Reco_lead_pho_phi    = sorted_reco_photons.phi[has_two_photons][:, 0]
        Reco_sublead_pho_phi = sorted_reco_photons.phi[has_two_photons][:, 1]

        def dR(eta1, phi1, eta2, phi2):
            d_eta = eta1 - eta2
            d_phi = phi1 - phi2
            return np.sqrt(d_eta**2 + d_phi**2)

        dR_pho_1 = dR(pho_from_a_eta_1, pho_from_a_phi_1, Reco_pho_eta_uf, Reco_pho_phi_uf)
        dR_pho_2 = dR(pho_from_a_eta_2, pho_from_a_phi_2, Reco_pho_eta_uf, Reco_pho_phi_uf)

        min_idx_1 = ak.argmin(dR_pho_1, axis=1)
        min_idx_2 = ak.argmin(dR_pho_2, axis=1)

        photon_idx = ak.local_index(reco_photons)
        mask_idx_1 = photon_idx == min_idx_1[:, None]
        mask_idx_2 = photon_idx == min_idx_2[:, None]

        mask_dR_1 = dR_pho_1 < 0.1
        mask_dR_2 = dR_pho_2 < 0.1

        mask_sel_1 = mask_idx_1 & mask_dR_1
        mask_sel_2 = mask_idx_2 & mask_dR_2


        selected_photon_1 = ak.firsts(reco_photons[mask_sel_1])
        selected_photon_2 = ak.firsts(reco_photons[mask_sel_2])


        Gen_photon_pt_1 = ak.to_numpy(pho_from_a_pt_1[~ak.is_none(selected_photon_1.pt)])
        Gen_photon_pt_2 = ak.to_numpy(pho_from_a_pt_2[~ak.is_none(selected_photon_2.pt)])
        Genmatched_photon_1_clean = selected_photon_1[~ak.is_none(selected_photon_1)]
        Genmatched_photon_2_clean = selected_photon_2[~ak.is_none(selected_photon_2)]

        Genmatched_pho_1_pt  = ak.to_numpy(Genmatched_photon_1_clean.pt)
        Genmatched_pho_2_pt  = ak.to_numpy(Genmatched_photon_2_clean.pt)
        Genmatched_pho_1_eta = ak.to_numpy(Genmatched_photon_1_clean.eta)
        Genmatched_pho_2_eta = ak.to_numpy(Genmatched_photon_2_clean.eta)
        Genmatched_pho_1_phi = ak.to_numpy(Genmatched_photon_1_clean.phi)
        Genmatched_pho_2_phi = ak.to_numpy(Genmatched_photon_2_clean.phi)
        zeros = np.zeros_like(Genmatched_pho_2_phi, dtype=np.float32)
        Genmatched_pho_1_mass = np.zeros_like(Genmatched_pho_2_phi, dtype=np.float32)
        Genmatched_pho_2_mass = np.zeros_like(Genmatched_pho_2_phi, dtype=np.float32)


        bjets_all = Events.Jet[(Events.Jet.hadronFlavour == 5)]

        mask_bjets = ak.num(bjets_all) >= 2
        bjets = bjets_all[mask_bjets]

        gen_b = gen[(abs(gen.pdgId) == 5) & (gen.status == 23)]
        from_a_mask = abs(gen[gen_b.genPartIdxMother].pdgId) == 35
        gen_b_from_a = gen_b[from_a_mask]

        gen_b_from_a = gen_b_from_a[mask_bjets]

        # Sort gen_b_from_a by pt (descending)
        sorted_gen_b = gen_b_from_a[ak.argsort(gen_b_from_a.pt, axis=1, ascending=False)]

        # Extract leading and subleading gen b-quark kinematics
        gen_lead_b_pt  = sorted_gen_b.pt[:, 0]
        gen_lead_b_eta = sorted_gen_b.eta[:, 0]
        gen_lead_b_phi = sorted_gen_b.phi[:, 0]

        gen_sublead_b_pt  = sorted_gen_b.pt[:, 1]
        gen_sublead_b_eta = sorted_gen_b.eta[:, 1]
        gen_sublead_b_phi = sorted_gen_b.phi[:, 1]


        # Select leading two b-quarks
        b_pt_1 = gen_b_from_a.pt[:, 0]
        b_pt_2 = gen_b_from_a.pt[:, 1]
        b_eta_1 = gen_b_from_a.eta[:, 0]
        b_phi_1 = gen_b_from_a.phi[:, 0]
        b_eta_2 = gen_b_from_a.eta[:, 1]
        b_phi_2 = gen_b_from_a.phi[:, 1]

        jet_pt = bjets.pt
        jet_eta = bjets.eta
        jet_phi = bjets.phi

        # ΔR matching
        dR_b1 = dR(b_eta_1, b_phi_1, jet_eta, jet_phi)
        dR_b2 = dR(b_eta_2, b_phi_2, jet_eta, jet_phi)

        # Find closest reco jet
        min_idx_1 = ak.argmin(dR_b1, axis=1)
        min_idx_2 = ak.argmin(dR_b2, axis=1)

        # jet_idx = ak.local_index(bjets)

        # mask_1 = (jet_idx == min_idx_1[:, None]) & (dR_b1 < 0.4)
        # mask_2 = (jet_idx == min_idx_2[:, None]) & (dR_b2 < 0.4)

        jet_idx = ak.local_index(bjets, axis=1)

        mask_1 = (jet_idx == min_idx_1) & (dR_b1 < 0.4)
        mask_2 = (jet_idx == min_idx_2) & (dR_b2 < 0.4)


        # Select matched reco jets
        selected_bjet_1 = ak.firsts(bjets[mask_1])
        selected_bjet_2 = ak.firsts(bjets[mask_2])

        # Kinematics of selected reco b-jets
        bjet1_pt  = selected_bjet_1.pt
        bjet1_eta = selected_bjet_1.eta
        bjet1_phi = selected_bjet_1.phi

        bjet2_pt  = selected_bjet_2.pt
        bjet2_eta = selected_bjet_2.eta
        bjet2_phi = selected_bjet_2.phi

        sorted_reco_bjets = bjets[ak.argsort(bjets.pt, axis=1, ascending=False)]
        has_lead = ak.num(sorted_reco_bjets.pt) > 0
        Reco_lead_bjet_pt = sorted_reco_bjets.pt[has_lead][:, 0]
        Reco_lead_bjet_eta = sorted_reco_bjets.eta[has_lead][:, 0]
        Reco_sublead_bjet_phi = sorted_reco_bjets.phi[has_lead][:, 0]
        Reco_sublead_bjet_pt = sorted_reco_bjets.pt[has_lead][:, 1]
        Reco_sublead_bjet_eta = sorted_reco_bjets.eta[has_lead][:, 1]
        Reco_sublead_bjet_phi = sorted_reco_bjets.phi[has_lead][:, 1]

        gen_b1_pt = ak.to_numpy(gen_b_from_a.pt[:, 0])
        gen_b2_pt = ak.to_numpy(gen_b_from_a.pt[:, 1])
        reco_b1_pt = ak.to_numpy(bjet1_pt)
        reco_b2_pt = ak.to_numpy(bjet2_pt)
        Reco_bjet1_eta = ak.to_numpy(bjet1_eta)
        Reco_bjet1_phi = ak.to_numpy(bjet1_phi)
        Reco_bjet2_eta = ak.to_numpy(bjet2_eta)
        Reco_bjet2_phi = ak.to_numpy(bjet2_phi)

        # Assume these are Awkward Arrays of shape (N,)
        pt1, eta1, phi1, mass1 = selected_bjet_1.pt, selected_bjet_1.eta, selected_bjet_1.phi, selected_bjet_1.mass
        pt2, eta2, phi2, mass2 = selected_bjet_2.pt, selected_bjet_2.eta, selected_bjet_2.phi, selected_bjet_2.mass

        # Convert to NumPy arrays
        pt1 = ak.to_numpy(pt1)
        eta1 = ak.to_numpy(eta1)
        phi1 = ak.to_numpy(phi1)
        mass1 = ak.to_numpy(mass1)
        pt2 = ak.to_numpy(pt2)
        eta2 = ak.to_numpy(eta2)
        phi2 = ak.to_numpy(phi2)
        mass2 = ak.to_numpy(mass2)

        # Mask invalid entries (NaNs or None)
        valid_mask = ~np.isnan(pt1 + eta1 + phi1 + mass1 + pt2 + eta2 + phi2 + mass2)

        pt1, eta1, phi1, mass1 = pt1[valid_mask], eta1[valid_mask], phi1[valid_mask], mass1[valid_mask]
        pt2, eta2, phi2, mass2 = pt2[valid_mask], eta2[valid_mask], phi2[valid_mask], mass2[valid_mask]

        # Create TLorentzVectors and compute invariant mass
        invmasses_bb = np.empty(len(pt1))
        vec1 = TLorentzVector()
        vec2 = TLorentzVector()

        for i in range(len(pt1)):
            vec1.SetPtEtaPhiM(pt1[i], eta1[i], phi1[i], mass1[i])
            vec2.SetPtEtaPhiM(pt2[i], eta2[i], phi2[i], mass2[i])
            invmasses_bb[i] = (vec1 + vec2).M()

        both_matched_mask = ~ak.is_none(selected_photon_1.pt) & ~ak.is_none(selected_photon_2.pt)


        Genmatched_photon_1_clean_both = selected_photon_1[both_matched_mask]
        Genmatched_photon_2_clean_both = selected_photon_2[both_matched_mask]

        Genmatched_pho_1_pt_both  = ak.to_numpy(Genmatched_photon_1_clean_both.pt)
        Genmatched_pho_2_pt_both  = ak.to_numpy(Genmatched_photon_2_clean_both.pt)
        Genmatched_pho_1_eta_both = ak.to_numpy(Genmatched_photon_1_clean_both.eta)
        Genmatched_pho_2_eta_both = ak.to_numpy(Genmatched_photon_2_clean_both.eta)
        Genmatched_pho_1_phi_both = ak.to_numpy(Genmatched_photon_1_clean_both.phi)
        Genmatched_pho_2_phi_both = ak.to_numpy(Genmatched_photon_2_clean_both.phi)
        Genmatched_pho_1_mass_both = np.zeros_like(Genmatched_pho_2_phi_both, dtype=np.float32)
        Genmatched_pho_2_mass_both = np.zeros_like(Genmatched_pho_2_phi_both, dtype=np.float32)

        # Allocate output array
        invmasses_diphoton = np.empty(len(Genmatched_pho_1_pt_both), dtype=np.float32)

        # TLorentzVectors
        vec1 = TLorentzVector()
        vec2 = TLorentzVector()

        for i in range(len(Genmatched_pho_1_pt_both)):
            vec1.SetPtEtaPhiM(
                Genmatched_pho_1_pt_both[i],
                Genmatched_pho_1_eta_both[i],
                Genmatched_pho_1_phi_both[i],
                Genmatched_pho_1_mass_both[i],
            )
            vec2.SetPtEtaPhiM(
                Genmatched_pho_2_pt_both[i],
                Genmatched_pho_2_eta_both[i],
                Genmatched_pho_2_phi_both[i],
                Genmatched_pho_2_mass_both[i],
            )
            invmasses_diphoton[i] = (vec1 + vec2).M()

        gen_invmasses_diphoton = np.empty(len(pho_from_a_pt_1), dtype=np.float32)
        gen_invmasses_bb = np.empty(len(bquark_from_a_pt_1), dtype=np.float32)

        pho_from_a_mass_1 = np.zeros_like(pho_from_a_pt_1, dtype=np.float32)
        pho_from_a_mass_2 = np.zeros_like(pho_from_a_pt_1, dtype=np.float32)

        # bquark_from_a_mass_1 = ak.to_numpy(bquarks_from_a.mass[:, 0])
        # bquark_from_a_mass_2 = ak.to_numpy(bquarks_from_a.mass[:, 1])
        bquark_from_a_mass_1 = np.full_like(ak.to_numpy(bquarks_from_a.pt[:, 0]), 4.18)
        bquark_from_a_mass_2 = np.full_like(ak.to_numpy(bquarks_from_a.pt[:, 1]), 4.18)

        # Masks for valid (non-NaN) gen-level photons and b-quarks
        valid_mask_gg = ~np.isnan(pho_from_a_pt_1 + pho_from_a_eta_1 + pho_from_a_phi_1 + pho_from_a_mass_1 +
                                pho_from_a_pt_2 + pho_from_a_eta_2 + pho_from_a_phi_2 + pho_from_a_mass_2)

        valid_mask_bb = ~np.isnan(bquark_from_a_pt_1 + bquark_from_a_eta_1 + bquark_from_a_phi_1 + bquark_from_a_mass_1 +
                                bquark_from_a_pt_2 + bquark_from_a_eta_2 + bquark_from_a_phi_2 + bquark_from_a_mass_2)

        # Filter the arrays
        pho_1_pt = pho_from_a_pt_1[valid_mask_gg]
        pho_1_eta = pho_from_a_eta_1[valid_mask_gg]
        pho_1_phi = pho_from_a_phi_1[valid_mask_gg]
        pho_1_mass = pho_from_a_mass_1[valid_mask_gg]

        pho_2_pt = pho_from_a_pt_2[valid_mask_gg]
        pho_2_eta = pho_from_a_eta_2[valid_mask_gg]
        pho_2_phi = pho_from_a_phi_2[valid_mask_gg]
        pho_2_mass = pho_from_a_mass_2[valid_mask_gg]

        b1_pt = bquark_from_a_pt_1[valid_mask_bb]
        b1_eta = bquark_from_a_eta_1[valid_mask_bb]
        b1_phi = bquark_from_a_phi_1[valid_mask_bb]
        b1_mass = bquark_from_a_mass_1[valid_mask_bb]

        b2_pt = bquark_from_a_pt_2[valid_mask_bb]
        b2_eta = bquark_from_a_eta_2[valid_mask_bb]
        b2_phi = bquark_from_a_phi_2[valid_mask_bb]
        b2_mass = bquark_from_a_mass_2[valid_mask_bb]

        # Prepare output arrays
        gen_invmasses_diphoton = np.empty(len(pho_1_pt), dtype=np.float32)
        gen_invmasses_bb = np.empty(len(b1_pt), dtype=np.float32)

        p4_pho_1 = TLorentzVector()
        p4_pho_2 = TLorentzVector()
        p4_bq_1 = TLorentzVector()
        p4_bq_2 = TLorentzVector()


        for i, (PT1, ETA1, PHI1, M1, PT2, ETA2, PHI2, M2) in enumerate(zip(pho_1_pt, pho_1_eta, pho_1_phi, pho_1_mass,
                                                                        pho_2_pt, pho_2_eta, pho_2_phi, pho_2_mass)):
            p4_pho_1.SetPtEtaPhiM(PT1, ETA1, PHI1, M1)
            p4_pho_2.SetPtEtaPhiM(PT2, ETA2, PHI2, M2)
            gen_invmasses_diphoton[i] = (p4_pho_1 + p4_pho_2).M()


        for i, (PT1, ETA1, PHI1, M1, PT2, ETA2, PHI2, M2) in enumerate(zip(b1_pt, b1_eta, b1_phi, b1_mass,
                                                                        b2_pt, b2_eta, b2_phi, b2_mass)):
            p4_bq_1.SetPtEtaPhiM(PT1, ETA1, PHI1, M1)
            p4_bq_2.SetPtEtaPhiM(PT2, ETA2, PHI2, M2)
            gen_invmasses_bb[i] = (p4_bq_1 + p4_bq_2).M()

        output["pu_true"] += column_accumulator(ak.to_numpy(pu_true))
        output["nPhoton"] += column_accumulator(ak.to_numpy(nPhoton))
        output["higgs_pt"] += column_accumulator(ak.to_numpy(higgs_pt))
        output["higgs_eta"] += column_accumulator(ak.to_numpy(higgs_eta))
        output["higgs_phi"] += column_accumulator(ak.to_numpy(higgs_phi))

        output["A_pt_1"] += column_accumulator(ak.to_numpy(A_pt[:, 0]))
        output["A_pt_2"] += column_accumulator(ak.to_numpy(A_pt[:, 1]))
        output["A_eta_1"] += column_accumulator(ak.to_numpy(A_eta[:, 0]))
        output["A_eta_2"] += column_accumulator(ak.to_numpy(A_eta[:, 1]))
        output["A_phi_1"] += column_accumulator(ak.to_numpy(A_phi[:, 0]))
        output["A_phi_2"] += column_accumulator(ak.to_numpy(A_phi[:, 1]))

        output["leading_A_pt"] += column_accumulator(ak.to_numpy(leading_A_pt))
        output["subleading_A_pt"] += column_accumulator(ak.to_numpy(subleading_A_pt))
        output["leading_A_eta"] += column_accumulator(ak.to_numpy(leading_A_eta))
        output["subleading_A_eta"] += column_accumulator(ak.to_numpy(subleading_A_eta))
        output["leading_A_phi"] += column_accumulator(ak.to_numpy(leading_A_phi))
        output["subleading_A_phi"] += column_accumulator(ak.to_numpy(subleading_A_phi))

        output["pho_from_a_pt_1"] += column_accumulator(ak.to_numpy(pho_from_a_pt_1))
        output["pho_from_a_pt_2"] += column_accumulator(ak.to_numpy(pho_from_a_pt_2))
        output["pho_from_a_eta_1"] += column_accumulator(ak.to_numpy(pho_from_a_eta_1))
        output["pho_from_a_eta_2"] += column_accumulator(ak.to_numpy(pho_from_a_eta_2))
        output["pho_from_a_phi_1"] += column_accumulator(ak.to_numpy(pho_from_a_phi_1))
        output["pho_from_a_phi_2"] += column_accumulator(ak.to_numpy(pho_from_a_phi_2))

        output["lead_pt_pho_gen"] += column_accumulator(ak.to_numpy(lead_pt_pho_gen))
        output["sublead_pt_pho_gen"] += column_accumulator(ak.to_numpy(sublead_pt_pho_gen))
        output["lead_eta_pho_gen"] += column_accumulator(ak.to_numpy(lead_eta_pho_gen))
        output["sublead_eta_pho_gen"] += column_accumulator(ak.to_numpy(sublead_eta_pho_gen))
        output["lead_phi_pho_gen"] += column_accumulator(ak.to_numpy(lead_phi_pho_gen))
        output["sublead_phi_pho_gen"] += column_accumulator(ak.to_numpy(sublead_phi_pho_gen))

        output["bquark_from_a_pt_1"] += column_accumulator(ak.to_numpy(bquark_from_a_pt_1))
        output["bquark_from_a_pt_2"] += column_accumulator(ak.to_numpy(bquark_from_a_pt_2))
        output["bquark_from_a_eta_1"] += column_accumulator(ak.to_numpy(bquark_from_a_eta_1))
        output["bquark_from_a_eta_2"] += column_accumulator(ak.to_numpy(bquark_from_a_eta_2))
        output["bquark_from_a_phi_1"] += column_accumulator(ak.to_numpy(bquark_from_a_phi_1))
        output["bquark_from_a_phi_2"] += column_accumulator(ak.to_numpy(bquark_from_a_phi_2))

        output["lead_pt_bquark_gen"] += column_accumulator(ak.to_numpy(lead_pt_bquark_gen))
        output["sublead_pt_bquark_gen"] += column_accumulator(ak.to_numpy(sublead_pt_bquark_gen))
        output["lead_eta_bquark_gen"] += column_accumulator(ak.to_numpy(lead_eta_bquark_gen))
        output["sublead_eta_bquark_gen"] += column_accumulator(ak.to_numpy(sublead_eta_bquark_gen))
        output["lead_phi_bquark_gen"] += column_accumulator(ak.to_numpy(lead_phi_bquark_gen))
        output["sublead_phi_bquark_gen"] += column_accumulator(ak.to_numpy(sublead_phi_bquark_gen))

        output["Reco_pho_pt"] += column_accumulator(ak.to_numpy(Reco_pho_pt))
        output["Reco_pho_eta"] += column_accumulator(ak.to_numpy(Reco_pho_eta))
        output["Reco_pho_phi"] += column_accumulator(ak.to_numpy(Reco_pho_phi))

        output["Reco_photon_lead_pt_all"] += column_accumulator(ak.to_numpy(Reco_photon_lead_pt_all))
        output["Reco_photon_sublead_pt_all"] += column_accumulator(ak.to_numpy(Reco_photon_sublead_pt_all))
        output["Reco_photon_lead_eta_all"] += column_accumulator(ak.to_numpy(Reco_photon_lead_eta_all))
        output["Reco_photon_sublead_eta_all"] += column_accumulator(ak.to_numpy(Reco_photon_sublead_eta_all))
        output["Reco_photon_lead_phi_all"] += column_accumulator(ak.to_numpy(Reco_photon_lead_phi_all))
        output["Reco_photon_sublead_phi_all"] += column_accumulator(ak.to_numpy(Reco_photon_sublead_phi_all))

        output["Reco_lead_pho_pt"] += column_accumulator(ak.to_numpy(Reco_lead_pho_pt))
        output["Reco_sublead_pho_pt"] += column_accumulator(ak.to_numpy(Reco_sublead_pho_pt))
        output["Reco_lead_pho_eta"] += column_accumulator(ak.to_numpy(Reco_lead_pho_eta))
        output["Reco_sublead_pho_eta"] += column_accumulator(ak.to_numpy(Reco_sublead_pho_eta))
        output["Reco_lead_pho_phi"] += column_accumulator(ak.to_numpy(Reco_lead_pho_phi))
        output["Reco_sublead_pho_phi"] += column_accumulator(ak.to_numpy(Reco_sublead_pho_phi))

        output["gen_photon_from_a_1_pt"] += column_accumulator(ak.to_numpy(pho_from_a_pt_1))
        output["gen_photon_from_a_2_pt"] += column_accumulator(ak.to_numpy(pho_from_a_pt_2))

        output["Gen_photon_pt_1"] += column_accumulator(Gen_photon_pt_1)
        output["Gen_photon_pt_2"] += column_accumulator(Gen_photon_pt_2)

        output["Genmatched_pho_1_pt"] += column_accumulator(Genmatched_pho_1_pt)
        output["Genmatched_pho_2_pt"] += column_accumulator(Genmatched_pho_2_pt)
        output["Genmatched_pho_1_eta"] += column_accumulator(Genmatched_pho_1_eta)
        output["Genmatched_pho_2_eta"] += column_accumulator(Genmatched_pho_2_eta)
        output["Genmatched_pho_1_phi"] += column_accumulator(Genmatched_pho_1_phi)
        output["Genmatched_pho_2_phi"] += column_accumulator(Genmatched_pho_2_phi)

        # Gen-level b-quarks from A, pt-sorted
        output["gen_lead_b_pt"]  += column_accumulator(ak.to_numpy(gen_lead_b_pt))
        output["gen_lead_b_eta"] += column_accumulator(ak.to_numpy(gen_lead_b_eta))
        output["gen_lead_b_phi"] += column_accumulator(ak.to_numpy(gen_lead_b_phi))

        output["gen_sublead_b_pt"]  += column_accumulator(ak.to_numpy(gen_sublead_b_pt))
        output["gen_sublead_b_eta"] += column_accumulator(ak.to_numpy(gen_sublead_b_eta))
        output["gen_sublead_b_phi"] += column_accumulator(ak.to_numpy(gen_sublead_b_phi))

        # Gen b-quarks (unordered)
        output["gen_b1_pt"]  += column_accumulator(ak.to_numpy(b_pt_1))
        output["gen_b1_eta"] += column_accumulator(ak.to_numpy(b_eta_1))
        output["gen_b1_phi"] += column_accumulator(ak.to_numpy(b_phi_1))

        output["gen_b2_pt"]  += column_accumulator(ak.to_numpy(b_pt_2))
        output["gen_b2_eta"] += column_accumulator(ak.to_numpy(b_eta_2))
        output["gen_b2_phi"] += column_accumulator(ak.to_numpy(b_phi_2))

        # Reco b-jets (pt-sorted)
        output["reco_lead_bjet_pt"]  += column_accumulator(ak.to_numpy(Reco_lead_bjet_pt))
        output["reco_lead_bjet_eta"] += column_accumulator(ak.to_numpy(Reco_lead_bjet_eta))
        output["reco_lead_bjet_phi"] += column_accumulator(ak.to_numpy(Reco_sublead_bjet_phi))  # Note: typo? should be `Reco_lead_bjet_phi`?

        output["reco_sublead_bjet_pt"]  += column_accumulator(ak.to_numpy(Reco_sublead_bjet_pt))
        output["reco_sublead_bjet_eta"] += column_accumulator(ak.to_numpy(Reco_sublead_bjet_eta))
        output["reco_sublead_bjet_phi"] += column_accumulator(ak.to_numpy(Reco_sublead_bjet_phi))

        # Matched reco b-jets
        output["matched_bjet1_pt"]  += column_accumulator(ak.to_numpy(reco_b1_pt))
        output["matched_bjet1_eta"] += column_accumulator(ak.to_numpy(Reco_bjet1_eta))
        output["matched_bjet1_phi"] += column_accumulator(ak.to_numpy(Reco_bjet1_phi))

        output["matched_bjet2_pt"]  += column_accumulator(ak.to_numpy(reco_b2_pt))
        output["matched_bjet2_eta"] += column_accumulator(ak.to_numpy(Reco_bjet2_eta))
        output["matched_bjet2_phi"] += column_accumulator(ak.to_numpy(Reco_bjet2_phi))

        output["gen_invmasses_diphoton"] += column_accumulator(gen_invmasses_diphoton)
        output["gen_invmasses_bb"] += column_accumulator(gen_invmasses_bb)

        output["invmasses_diphoton"] += column_accumulator(invmasses_diphoton)
        output["invmasses_bb"] += column_accumulator(invmasses_bb)

        # Wrap output inside a dict keyed by dataset name!
        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator