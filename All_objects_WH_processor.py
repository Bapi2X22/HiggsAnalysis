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
            "nEvents": column_accumulator(np.array([], dtype=np.float32)),
            "pu_true": column_accumulator(np.array([], dtype=np.float32)),
            "nPhoton": column_accumulator(np.array([], dtype=np.float32)),
            "nPhoton_cut_10": column_accumulator(np.array([], dtype=np.float32)),
            "nPhoton_cut_18": column_accumulator(np.array([], dtype=np.float32)),
            "nPhoton_cut_30": column_accumulator(np.array([], dtype=np.float32)),
            "higgs_pt": column_accumulator(np.array([], dtype=np.float32)),
            "higgs_eta": column_accumulator(np.array([], dtype=np.float32)),
            "higgs_phi": column_accumulator(np.array([], dtype=np.float32)),
            "W_pt": column_accumulator(np.array([], dtype=np.float32)),
            "W_eta": column_accumulator(np.array([], dtype=np.float32)),
            "W_phi": column_accumulator(np.array([], dtype=np.float32)),

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

            "Reco_lead_pho_pt": column_accumulator(np.array([], dtype=np.float32)),
            "Reco_sublead_pho_pt": column_accumulator(np.array([], dtype=np.float32)),
            "Reco_lead_pho_eta": column_accumulator(np.array([], dtype=np.float32)),
            "Reco_sublead_pho_eta": column_accumulator(np.array([], dtype=np.float32)),
            "Reco_lead_pho_phi": column_accumulator(np.array([], dtype=np.float32)),
            "Reco_sublead_pho_phi": column_accumulator(np.array([], dtype=np.float32)),

            "gen_photon_from_a_1_pt": column_accumulator(np.array([], dtype=np.float32)),
            "gen_photon_from_a_2_pt": column_accumulator(np.array([], dtype=np.float32)),

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

            "gen_ele_pt": column_accumulator(np.array([], dtype=np.float32)),
            "gen_ele_eta": column_accumulator(np.array([], dtype=np.float32)),
            "gen_ele_phi": column_accumulator(np.array([], dtype=np.float32)),
            "matched_reco_ele_pt": column_accumulator(np.array([], dtype=np.float32)),
            "matched_reco_ele_eta": column_accumulator(np.array([], dtype=np.float32)),
            "matched_reco_ele_phi": column_accumulator(np.array([], dtype=np.float32)),
            "gen_mu_pt": column_accumulator(np.array([], dtype=np.float32)),
            "gen_mu_eta": column_accumulator(np.array([], dtype=np.float32)),
            "gen_mu_phi": column_accumulator(np.array([], dtype=np.float32)),
            "matched_reco_mu_pt": column_accumulator(np.array([], dtype=np.float32)),
            "matched_reco_mu_eta": column_accumulator(np.array([], dtype=np.float32)),
            "matched_reco_mu_phi": column_accumulator(np.array([], dtype=np.float32)),
            "gen_ele_tauW_pt": column_accumulator(np.array([], dtype=np.float32)),
            "gen_ele_tauW_eta": column_accumulator(np.array([], dtype=np.float32)),
            "gen_ele_tauW_phi": column_accumulator(np.array([], dtype=np.float32)),
            "matched_reco_ele_tau_pt": column_accumulator(np.array([], dtype=np.float32)),
            "matched_reco_ele_tau_eta": column_accumulator(np.array([], dtype=np.float32)),
            "matched_reco_ele_tau_phi": column_accumulator(np.array([], dtype=np.float32)),
            "gen_mu_tauW_pt": column_accumulator(np.array([], dtype=np.float32)),
            "gen_mu_tauW_eta": column_accumulator(np.array([], dtype=np.float32)),
            "gen_mu_tauW_phi": column_accumulator(np.array([], dtype=np.float32)),
            "matched_reco_mu_tau_pt": column_accumulator(np.array([], dtype=np.float32)),
            "matched_reco_mu_tau_eta": column_accumulator(np.array([], dtype=np.float32)),
            "matched_reco_mu_tau_phi": column_accumulator(np.array([], dtype=np.float32)),
        })

    def accumulator(self):
        # Return the accumulator prototype for coffea
        return self._accumulator

    def process(self, events):

        Events = events

        dataset = Events.metadata["dataset"]  # get dataset name dynamically

        output = self.accumulator()
        # nPhoton = Events["nPhoton"]

        nEvents = Events.event
        
        pu_true = Events.Pileup.nTrueInt
        # nPhoton = Events.Photon.nPhoton
        nPhoton = ak.num(events.Photon)

        # nPhoton = Events["nPhoton"]

        pt_cut_10 = Events.Photon.pt>10.0
        pt_cut_18 = Events.Photon.pt>18.0
        pt_cut_30 = Events.Photon.pt>30.0

        nPhoton_cut_10 = ak.num(Events.Photon[pt_cut_10].pt)
        nPhoton_cut_18 = ak.num(Events.Photon[pt_cut_18].pt)
        nPhoton_cut_30 = ak.num(Events.Photon[pt_cut_30].pt)

        gen = Events.GenPart
        higgs_mask = (gen.pdgId == 25) & (gen.status == 62)
        higgs = gen[higgs_mask]
        W_mask = (abs(gen.pdgId) == 24) & (gen.status == 62)
        W = gen[W_mask]

        higgs_pt = ak.flatten(higgs.pt)
        higgs_eta = ak.flatten(higgs.eta)
        higgs_phi = ak.flatten(higgs.phi)

        W_pt = ak.flatten(W.pt)
        W_eta = ak.flatten(W.eta)
        W_phi = ak.flatten(W.phi)

        is_A = (abs(gen.pdgId) == 35)

        A = gen[is_A]

        A_pt = A.pt
        A_eta = A.eta
        A_phi = A.phi

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

        # Pad to at least 2 photons per event (None if not available)
        photons_from_a_padded = ak.pad_none(photons_from_a, 2, axis=1, clip=True)

        # Extract pt, eta, phi
        pho_from_a_pt  = photons_from_a_padded.pt
        pho_from_a_eta = photons_from_a_padded.eta
        pho_from_a_phi = photons_from_a_padded.phi

        # Split leading and subleading, replacing None with NaN
        pho_from_a_pt_1  = ak.fill_none(pho_from_a_pt[:, 0], np.nan)
        pho_from_a_pt_2  = ak.fill_none(pho_from_a_pt[:, 1], np.nan)
        pho_from_a_eta_1 = ak.fill_none(pho_from_a_eta[:, 0], np.nan)
        pho_from_a_eta_2 = ak.fill_none(pho_from_a_eta[:, 1], np.nan)
        pho_from_a_phi_1 = ak.fill_none(pho_from_a_phi[:, 0], np.nan)
        pho_from_a_phi_2 = ak.fill_none(pho_from_a_phi[:, 1], np.nan)

        # Sort photons by pt
        sorted_photons = photons_from_a[ak.argsort(photons_from_a.pt, axis=1, ascending=False)]

        # Pad to at least 2 photons per event (None if missing)
        sorted_photons_padded = ak.pad_none(sorted_photons, 2, axis=1, clip=True)

        # Extract leading and subleading with NaN fill
        lead_pt_pho_gen  = ak.fill_none(sorted_photons_padded.pt[:, 0],  np.nan)
        sublead_pt_pho_gen  = ak.fill_none(sorted_photons_padded.pt[:, 1],  np.nan)

        lead_eta_pho_gen = ak.fill_none(sorted_photons_padded.eta[:, 0], np.nan)
        sublead_eta_pho_gen = ak.fill_none(sorted_photons_padded.eta[:, 1], np.nan)

        lead_phi_pho_gen = ak.fill_none(sorted_photons_padded.phi[:, 0], np.nan)
        sublead_phi_pho_gen = ak.fill_none(sorted_photons_padded.phi[:, 1], np.nan)

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

        # Don't drop events up front!
        reco_photons = Events.Photon  # keep all events

        Reco_pho_pt = ak.flatten(reco_photons.pt)
        Reco_pho_eta = ak.flatten(reco_photons.eta)
        Reco_pho_phi = ak.flatten(reco_photons.phi)

        # Sort photons by pt
        sorted_reco_photons = reco_photons[ak.argsort(reco_photons.pt, axis=1, ascending=False)]

        # Pad so every event has at least 2 photons (None if missing)
        sorted_reco_padded = ak.pad_none(sorted_reco_photons, 2, axis=1, clip=True)

        # Fill missing with NaN
        Reco_lead_pho_pt  = ak.fill_none(sorted_reco_padded.pt[:, 0],  np.nan)
        Reco_sublead_pho_pt  = ak.fill_none(sorted_reco_padded.pt[:, 1],  np.nan)

        Reco_lead_pho_eta = ak.fill_none(sorted_reco_padded.eta[:, 0], np.nan)
        Reco_sublead_pho_eta = ak.fill_none(sorted_reco_padded.eta[:, 1], np.nan)

        Reco_lead_pho_phi = ak.fill_none(sorted_reco_padded.phi[:, 0], np.nan)
        Reco_sublead_pho_phi = ak.fill_none(sorted_reco_padded.phi[:, 1], np.nan)

        def delta_r(eta1, phi1, eta2, phi2):
            dphi = (phi1 - phi2 + np.pi) % (2 * np.pi) - np.pi
            deta = eta1 - eta2
            return np.sqrt(deta**2 + dphi**2)

        # ---------------- ΔR between gen lead/sublead photons and all reco photons ----------------
        dr_pho_lead    = delta_r(lead_eta_pho_gen, lead_phi_pho_gen, reco_photons.eta, reco_photons.phi)
        dr_pho_sublead = delta_r(sublead_eta_pho_gen, sublead_phi_pho_gen, reco_photons.eta, reco_photons.phi)

        # --- select closest reco photon per gen photon ---
        min_idx_pho1 = ak.argmin(dr_pho_lead, axis=1, keepdims=False)
        min_idx_pho2 = ak.argmin(dr_pho_sublead, axis=1, keepdims=False)

        min_dr_pho1 = ak.min(dr_pho_lead, axis=1, initial=np.inf)
        min_dr_pho2 = ak.min(dr_pho_sublead, axis=1, initial=np.inf)

        # mask away if ΔR >= 0.1
        min_idx_pho1 = ak.mask(min_idx_pho1, min_dr_pho1 < 0.1)
        min_idx_pho2 = ak.mask(min_idx_pho2, min_dr_pho2 < 0.1)

        # --- conflict resolution (if needed, e.g., same reco photon matched to both gen photons) ---
        both_valid_pho = ~ak.is_none(min_idx_pho1) & ~ak.is_none(min_idx_pho2)
        indices_equal_pho = ak.where(both_valid_pho, min_idx_pho1 == min_idx_pho2, False)
        conflict_pho = both_valid_pho & indices_equal_pho
        valid_dr_pho = ~ak.is_none(dr_pho_lead) & ~ak.is_none(dr_pho_sublead)
        compare_mask_pho = conflict_pho & valid_dr_pho

        keep_pho1 = ak.where(compare_mask_pho, dr_pho_lead <= dr_pho_sublead, True)
        keep_pho2 = ak.where(compare_mask_pho, dr_pho_sublead < dr_pho_lead, True)

        min_idx_pho1 = ak.mask(min_idx_pho1, keep_pho1)
        min_idx_pho2 = ak.mask(min_idx_pho2, keep_pho2)

        photon_idx = ak.local_index(reco_photons)

        mask_idx_pho1 = photon_idx == min_idx_pho1
        mask_idx_pho2 = photon_idx == min_idx_pho2

        # --- now select reco photons ---
        selected_photon_1 = ak.firsts(reco_photons[mask_idx_pho1])
        selected_photon_2 = ak.firsts(reco_photons[mask_idx_pho2])

        # --- store reco photon kinematics (NaN if no match) ---
        Genmatched_pho_1_pt  = ak.to_numpy(ak.fill_none(selected_photon_1.pt,  np.nan))
        Genmatched_pho_2_pt  = ak.to_numpy(ak.fill_none(selected_photon_2.pt,  np.nan))
        Genmatched_pho_1_eta = ak.to_numpy(ak.fill_none(selected_photon_1.eta, np.nan))
        Genmatched_pho_2_eta = ak.to_numpy(ak.fill_none(selected_photon_2.eta, np.nan))
        Genmatched_pho_1_phi = ak.to_numpy(ak.fill_none(selected_photon_1.phi, np.nan))
        Genmatched_pho_2_phi = ak.to_numpy(ak.fill_none(selected_photon_2.phi, np.nan))

        # --- GEN b-quarks from A, with status==23; keep full event alignment ---

        # 1) Select b-quarks with status==23 (per-event filtering preserves event count)
        gen_b = gen[(abs(gen.pdgId) == 5) & (gen.status == 23)]

        # 2) Require valid mother and mother pdgId == 35 (A)
        mother_idx  = gen_b.genPartIdxMother
        valid_mom   = mother_idx >= 0
        idx_safe    = ak.where(valid_mom, mother_idx, 0)          # avoid negative indexing
        mom_pdgid   = gen[idx_safe].pdgId
        from_a_mask = valid_mom & (abs(mom_pdgid) == 35)

        # Keep only b’s from A (still preserves event count; events with no matches → empty lists)
        gen_b_from_a = gen_b[from_a_mask]

        # 3) Pad to at least 2 entries per event, then fill missing with NaN when extracting
        gen_b_padded = ak.pad_none(gen_b_from_a, 2, axis=1, clip=True)

        # --- As-stored (index 0 / 1) ---
        gen_b_pt_1  = ak.to_numpy(ak.fill_none(gen_b_padded.pt[:, 0],  np.nan))
        gen_b_eta_1 = ak.to_numpy(ak.fill_none(gen_b_padded.eta[:, 0], np.nan))
        gen_b_phi_1 = ak.to_numpy(ak.fill_none(gen_b_padded.phi[:, 0], np.nan))

        gen_b_pt_2  = ak.to_numpy(ak.fill_none(gen_b_padded.pt[:, 1],  np.nan))
        gen_b_eta_2 = ak.to_numpy(ak.fill_none(gen_b_padded.eta[:, 1], np.nan))
        gen_b_phi_2 = ak.to_numpy(ak.fill_none(gen_b_padded.phi[:, 1], np.nan))

        # --- Sorted by pT (leading / subleading) ---
        gen_b_sorted = gen_b_padded[ak.argsort(gen_b_padded.pt, axis=1, ascending=False)]

        gen_lead_b_pt  = ak.to_numpy(ak.fill_none(gen_b_sorted.pt[:, 0],  np.nan))
        gen_lead_b_eta = ak.to_numpy(ak.fill_none(gen_b_sorted.eta[:, 0], np.nan))
        gen_lead_b_phi = ak.to_numpy(ak.fill_none(gen_b_sorted.phi[:, 0], np.nan))

        gen_sublead_b_pt  = ak.to_numpy(ak.fill_none(gen_b_sorted.pt[:, 1],  np.nan))
        gen_sublead_b_eta = ak.to_numpy(ak.fill_none(gen_b_sorted.eta[:, 1], np.nan))
        gen_sublead_b_phi = ak.to_numpy(ak.fill_none(gen_b_sorted.phi[:, 1], np.nan))



        # --- Select reco b-jets ---
        bjets_all = Events.Jet[(Events.Jet.hadronFlavour == 5)]

        # Keep all events aligned, pad missing with NaN
        bjets_padded = ak.pad_none(bjets_all, 2, axis=1)  # at least 2 slots per event

        # ---------------- ΔR between gen b1/b2 and all reco b-jets ----------------
        # ---------------- ΔR between gen b1/b2 and all reco b-jets ----------------
        # dr_b1 = delta_r(gen_lead_b_eta, gen_lead_b_phi, bjets_all.eta, bjets_all.phi)
        # dr_b2 = delta_r(gen_sublead_b_eta, gen_sublead_b_phi, bjets_all.eta, bjets_all.phi)

        # # --- select closest reco jet per gen b ---
        # min_idx_b1 = ak.argmin(dr_b1, axis=1, keepdims=False)
        # min_idx_b2 = ak.argmin(dr_b2, axis=1, keepdims=False)

        # min_dr_b1 = ak.min(dr_b1, axis=1, initial=np.inf)
        # min_dr_b2 = ak.min(dr_b2, axis=1, initial=np.inf)

        # # mask away if ΔR >= 0.4
        # min_idx_b1 = ak.mask(min_idx_b1, min_dr_b1 < 0.4)
        # min_idx_b2 = ak.mask(min_idx_b2, min_dr_b2 < 0.4)

        # # Start with both_valid mask
        # both_valid = ~ak.is_none(min_idx_b1) & ~ak.is_none(min_idx_b2)

        # # Only compare min_idx_b1 == min_idx_b2 where both are valid
        # indices_equal = ak.where(both_valid, min_idx_b1 == min_idx_b2, False)

        # # True conflict: both valid AND indices equal
        # conflict = both_valid & indices_equal

        # # Only compare ΔR where conflict exists and both dr are valid
        # valid_dr = ~ak.is_none(dr_b1) & ~ak.is_none(dr_b2)
        # compare_mask = conflict & valid_dr

        # # Determine which index to keep
        # keep_b1 = ak.where(compare_mask, dr_b1 <= dr_b2, True)
        # keep_b2 = ak.where(compare_mask, dr_b2 < dr_b1, True)

        # # Apply masks
        # min_idx_b1 = ak.mask(min_idx_b1, keep_b1)
        # min_idx_b2 = ak.mask(min_idx_b2, keep_b2)
    
        # bjets_idx = ak.local_index(bjets_all)

        # mask_idx_b1 = bjets_idx == min_idx_b1
        # mask_idx_b2 = bjets_idx == min_idx_b2

        # # # --- now select reco jets ---
        # selected_bjet_1 = ak.firsts(bjets_all[mask_idx_b1])
        # selected_bjet_2 = ak.firsts(bjets_all[mask_idx_b2])

        # ΔR calculations (as you already have)
        # --- ΔR to all reco b-jets ---
        dr_b1 = delta_r(gen_lead_b_eta, gen_lead_b_phi, bjets_all.eta, bjets_all.phi)
        dr_b2 = delta_r(gen_sublead_b_eta, gen_sublead_b_phi, bjets_all.eta, bjets_all.phi)

        # --- closest reco jet per gen b ---
        min_idx_b1 = ak.argmin(dr_b1, axis=1)
        min_idx_b2 = ak.argmin(dr_b2, axis=1)

        min_dr_b1 = ak.min(dr_b1, axis=1, initial=np.inf)
        min_dr_b2 = ak.min(dr_b2, axis=1, initial=np.inf)

        # --- apply ΔR < 0.4 cut, else set to -1 (invalid) ---
        min_idx_b1 = ak.where(min_dr_b1 < 0.4, min_idx_b1, -1)
        min_idx_b2 = ak.where(min_dr_b2 < 0.4, min_idx_b2, -1)

        # --- conflict resolution ---
        conflict = (min_idx_b1 >= 0) & (min_idx_b2 >= 0) & (min_idx_b1 == min_idx_b2)

        # pick the closer one if conflict
        keep_b1 = ak.where(conflict, min_dr_b1 <= min_dr_b2, True)
        keep_b2 = ak.where(conflict, min_dr_b2 <  min_dr_b1, True)

        min_idx_b1 = ak.where(keep_b1, min_idx_b1, -1)
        min_idx_b2 = ak.where(keep_b2, min_idx_b2, -1)

        # # --- build masks without option types ---
        # bjets_idx = ak.local_index(bjets_all)

        # mask_idx_b1 = (min_idx_b1 >= 0) & (bjets_idx == min_idx_b1)
        # mask_idx_b2 = (min_idx_b2 >= 0) & (bjets_idx == min_idx_b2)

        bjets_idx = ak.local_index(bjets_all, axis=1)

        min_idx_b1 = ak.fill_none(min_idx_b1, -1)
        min_idx_b2 = ak.fill_none(min_idx_b2, -1)

        min_idx_b1_b, bjets_idx_b = ak.broadcast_arrays(min_idx_b1, bjets_idx)
        min_idx_b2_b, bjets_idx_b = ak.broadcast_arrays(min_idx_b2, bjets_idx)

        mask_idx_b1 = (min_idx_b1_b >= 0) & (bjets_idx_b == min_idx_b1_b)
        mask_idx_b2 = (min_idx_b2_b >= 0) & (bjets_idx_b == min_idx_b2_b)


        # --- select reco jets (None if no match) ---
        selected_bjet_1 = ak.firsts(bjets_all[mask_idx_b1])
        selected_bjet_2 = ak.firsts(bjets_all[mask_idx_b2])

        Genmatched_bjet1_pt  = ak.to_numpy(ak.fill_none(selected_bjet_1.pt,  np.nan))
        Genmatched_bjet1_eta = ak.to_numpy(ak.fill_none(selected_bjet_1.eta, np.nan))
        Genmatched_bjet1_phi = ak.to_numpy(ak.fill_none(selected_bjet_1.phi, np.nan))
        Genmatched_bjet1_mass  = ak.to_numpy(ak.fill_none(selected_bjet_1.mass,  np.nan))

        Genmatched_bjet2_pt  = ak.to_numpy(ak.fill_none(selected_bjet_2.pt,  np.nan))
        Genmatched_bjet2_eta = ak.to_numpy(ak.fill_none(selected_bjet_2.eta, np.nan))
        Genmatched_bjet2_phi = ak.to_numpy(ak.fill_none(selected_bjet_2.phi, np.nan))
        Genmatched_bjet2_mass = ak.to_numpy(ak.fill_none(selected_bjet_2.mass, np.nan))

        # ---------------- Sorted reco b-jets (leading/subleading) ----------------
        sorted_reco_bjets = bjets_padded[ak.argsort(bjets_padded.pt, axis=1, ascending=False)]
        sorted_reco_bjets_padded = ak.pad_none(sorted_reco_bjets, 2, axis=1, clip=True)

        Reco_lead_bjet_pt  = ak.to_numpy(ak.fill_none(sorted_reco_bjets_padded.pt[:, 0], np.nan))
        Reco_lead_bjet_eta = ak.to_numpy(ak.fill_none(sorted_reco_bjets_padded.eta[:, 0], np.nan))
        Reco_lead_bjet_phi = ak.to_numpy(ak.fill_none(sorted_reco_bjets_padded.phi[:, 0], np.nan))

        Reco_sublead_bjet_pt  = ak.to_numpy(ak.fill_none(sorted_reco_bjets_padded.pt[:, 1], np.nan))
        Reco_sublead_bjet_eta = ak.to_numpy(ak.fill_none(sorted_reco_bjets_padded.eta[:, 1], np.nan))
        Reco_sublead_bjet_phi = ak.to_numpy(ak.fill_none(sorted_reco_bjets_padded.phi[:, 1], np.nan))

        # Assume these are Awkward Arrays of shape (N,)
        pt1, eta1, phi1, mass1 = Genmatched_bjet1_pt  , Genmatched_bjet1_eta , Genmatched_bjet1_phi , Genmatched_bjet1_mass
        pt2, eta2, phi2, mass2 = Genmatched_bjet2_pt , Genmatched_bjet2_eta, Genmatched_bjet2_phi, Genmatched_bjet2_mass

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


        # Assume these are Awkward Arrays of shape (N,)
        pho1_pt, pho1_eta, pho1_phi, pho1_mass = Genmatched_pho_1_pt, Genmatched_pho_1_eta, Genmatched_pho_1_phi, np.zeros_like(Genmatched_pho_1_pt)
        pho2_pt, pho2_eta, pho2_phi, pho2_mass = Genmatched_pho_2_pt, Genmatched_pho_2_eta, Genmatched_pho_2_phi, np.zeros_like(Genmatched_pho_2_pt)

        # Convert to NumPy arrays
        pho1_pt   = ak.to_numpy(pho1_pt)
        pho1_eta  = ak.to_numpy(pho1_eta)
        pho1_phi  = ak.to_numpy(pho1_phi)
        pho1_mass = ak.to_numpy(pho1_mass)

        pho2_pt   = ak.to_numpy(pho2_pt)
        pho2_eta  = ak.to_numpy(pho2_eta)
        pho2_phi  = ak.to_numpy(pho2_phi)
        pho2_mass = ak.to_numpy(pho2_mass)

        # Mask invalid entries (NaNs or None)
        valid_mask_pho = ~np.isnan(pho1_pt + pho1_eta + pho1_phi + pho1_mass + pho2_pt + pho2_eta + pho2_phi + pho2_mass)

        pho1_pt, pho1_eta, pho1_phi, pho1_mass = pho1_pt[valid_mask_pho], pho1_eta[valid_mask_pho], pho1_phi[valid_mask_pho], pho1_mass[valid_mask_pho]
        pho2_pt, pho2_eta, pho2_phi, pho2_mass = pho2_pt[valid_mask_pho], pho2_eta[valid_mask_pho], pho2_phi[valid_mask_pho], pho2_mass[valid_mask_pho]

        # Create TLorentzVectors and compute invariant mass
        invmasses_photons = np.empty(len(pho1_pt))
        vec_pho1 = TLorentzVector()
        vec_pho2 = TLorentzVector()

        for i in range(len(pho1_pt)):
            vec_pho1.SetPtEtaPhiM(pho1_pt[i], pho1_eta[i], pho1_phi[i], pho1_mass[i])
            vec_pho2.SetPtEtaPhiM(pho2_pt[i], pho2_eta[i], pho2_phi[i], pho2_mass[i])
            invmasses_photons[i] = (vec_pho1 + vec_pho2).M()

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

        dr_cut = 0.1  # ΔR threshold
        pt_ratio_cut = 0.5  # optional (gen pt close to reco pt)

        is_ele1 = (abs(gen.pdgId) == 11) & (gen.status == 1)
        is_ele23 = (abs(gen.pdgId) == 11) & (gen.status == 23)
        is_W = abs(gen.pdgId) == 24

        ele1 = gen[is_ele1]
        ele23_from_W = is_ele23 & is_W[gen.genPartIdxMother]

        from_W = is_W[ele1.genPartIdxMother]
        from_ele23_from_W = ele23_from_W[ele1.genPartIdxMother]

        ele1_sel = ele1[from_W | from_ele23_from_W]   # final GEN electrons

        gen_ele_pt = ele1_sel.pt
        gen_ele_eta = ele1_sel.eta
        gen_ele_phi = ele1_sel.phi

        gen_ele_pt  = ak.pad_none(gen_ele_pt, 1)[:, 0]
        gen_ele_eta = ak.pad_none(gen_ele_eta, 1)[:, 0]
        gen_ele_phi = ak.pad_none(gen_ele_phi, 1)[:, 0]

        # --- Reco electrons ---
        reco_electrons = Events.Electron

        # --- Build pairwise ΔR ---
        def delta_r(eta1, phi1, eta2, phi2):
            dphi = (phi1 - phi2 + np.pi) % (2 * np.pi) - np.pi
            deta = eta1 - eta2
            return np.sqrt(deta**2 + dphi**2)

        # pairwise ΔR between each gen-ele1_sel and all reco electrons
        dr_ele = delta_r(
            ele1_sel.eta[:, :, None],   # (events, ngen, 1)
            ele1_sel.phi[:, :, None],
            reco_electrons.eta[:, None, :],  # (events, 1, nreco)
            reco_electrons.phi[:, None, :]
        )

        # --- Matching criterion ---

        # find best match per GEN electron
        best_reco_idx_ele = ak.argmin(dr_ele, axis=2)
        best_dr_ele = ak.min(dr_ele, axis=2)

        # mask with ΔR requirement
        matched_mask_ele = best_dr_ele < dr_cut

        # select matched reco electrons
        matched_reco_ele = reco_electrons[best_reco_idx_ele]
        matched_reco_ele = ak.mask(matched_reco_ele, matched_mask_ele)   # <-- cleaner

        matched_reco_ele_pt = matched_reco_ele.pt
        matched_reco_ele_eta = matched_reco_ele.eta
        matched_reco_ele_phi = matched_reco_ele.phi

        matched_reco_ele_pt = ak.pad_none(matched_reco_ele_pt, 1)[:, 0]
        matched_reco_ele_eta = ak.pad_none(matched_reco_ele_eta, 1)[:, 0]
        matched_reco_ele_phi = ak.pad_none(matched_reco_ele_phi, 1)[:, 0]

        is_mu1 = (abs(gen.pdgId) == 13) & (gen.status == 1)
        is_mu23 = (abs(gen.pdgId) == 13) & (gen.status == 23)
        is_W = abs(gen.pdgId) == 24

        mu1 = gen[is_mu1]
        mu23_from_W = is_mu23 & is_W[gen.genPartIdxMother]

        from_W = is_W[mu1.genPartIdxMother]
        from_mu23_from_W = mu23_from_W[mu1.genPartIdxMother]

        mu1_sel = mu1[from_W | from_mu23_from_W]   # final GEN muons

        gen_mu_pt = mu1_sel.pt
        gen_mu_eta = mu1_sel.eta
        gen_mu_phi = mu1_sel.phi

        gen_mu_pt  = ak.pad_none(gen_mu_pt, 1)[:, 0]
        gen_mu_eta = ak.pad_none(gen_mu_eta, 1)[:, 0]
        gen_mu_phi = ak.pad_none(gen_mu_phi, 1)[:, 0]

        # --- Reco muons ---
        reco_muons = Events.Muon

        # --- Pairwise ΔR between each gen-muon and all reco muons ---
        dr_mu = delta_r(
            mu1_sel.eta[:, :, None],   # (events, ngen_mu, 1)
            mu1_sel.phi[:, :, None],
            reco_muons.eta[:, None, :],  # (events, 1, nreco_mu)
            reco_muons.phi[:, None, :]
        )

        best_reco_mu_idx = ak.argmin(dr_mu, axis=2)
        best_dr_mu = ak.min(dr_mu, axis=2)

        # mask with ΔR requirement
        matched_mu_mask = best_dr_mu < dr_cut

        # select matched reco muons
        matched_reco_mu = reco_muons[best_reco_mu_idx]
        matched_reco_mu = ak.mask(matched_reco_mu, matched_mu_mask)

        matched_reco_mu_pt = matched_reco_mu.pt
        matched_reco_mu_eta = matched_reco_mu.eta
        matched_reco_mu_phi = matched_reco_mu.phi

        matched_reco_mu_pt = ak.pad_none(matched_reco_mu_pt, 1)[:, 0]
        matched_reco_mu_eta = ak.pad_none(matched_reco_mu_eta, 1)[:, 0]
        matched_reco_mu_phi = ak.pad_none(matched_reco_mu_phi, 1)[:, 0]

        # --- Identify relevant particles ---
        is_W      = abs(gen.pdgId) == 24
        is_tau    = abs(gen.pdgId) == 15
        is_ele1   = (abs(gen.pdgId) == 11) & (gen.status == 1)
        is_ele23  = (abs(gen.pdgId) == 11) & (gen.status == 23)

        mother = gen.genPartIdxMother
        has_mother = mother >= 0

        # --- Taus directly from W ---
        tau_from_W = is_tau & has_mother & is_W[mother]

        # --- Electrons directly from tau_from_W ---
        ele1_from_tauW = is_ele1 & has_mother & tau_from_W[mother]

        # --- Status=23 electrons directly from tau_from_W ---
        ele23_from_tauW = is_ele23 & has_mother & tau_from_W[mother]

        # --- Final state e± from those status=23 electrons ---
        ele1_from_ele23_tauW = (
            is_ele1 & has_mother & ele23_from_tauW[mother]
        )

        # --- Combine both cases ---
        ele_tauW = gen[ele1_from_tauW | ele1_from_ele23_tauW]

        gen_ele_tauW_pt = ele_tauW.pt
        gen_ele_tauW_eta = ele_tauW.eta
        gen_ele_tauW_phi = ele_tauW.phi

        gen_ele_tauW_pt  = ak.pad_none(gen_ele_tauW_pt, 1)[:, 0]
        gen_ele_tauW_eta = ak.pad_none(gen_ele_tauW_eta, 1)[:, 0]
        gen_ele_tauW_phi = ak.pad_none(gen_ele_tauW_phi, 1)[:, 0]

        # pairwise ΔR between each gen-ele1_sel and all reco electrons
        dr_ele_tau = delta_r(
            ele_tauW.eta[:, :, None],   # (events, ngen, 1)
            ele_tauW.phi[:, :, None],
            reco_electrons.eta[:, None, :],  # (events, 1, nreco)
            reco_electrons.phi[:, None, :]
        )

        # find best match per GEN electron
        best_reco_idx_ele_tau = ak.argmin(dr_ele_tau, axis=2)
        best_dr_ele_tau = ak.min(dr_ele_tau, axis=2)

        # mask with ΔR requirement
        matched_mask_ele_tau = best_dr_ele_tau < dr_cut

        # select matched reco electrons
        matched_reco_ele_tau = reco_electrons[best_reco_idx_ele_tau]
        matched_reco_ele_tau = ak.mask(matched_reco_ele_tau, matched_mask_ele_tau)   # <-- cleaner

        matched_reco_ele_tau_pt = matched_reco_ele_tau.pt
        matched_reco_ele_tau_eta = matched_reco_ele_tau.eta
        matched_reco_ele_tau_phi = matched_reco_ele_tau.phi

        matched_reco_ele_tau_pt = ak.pad_none(matched_reco_ele_tau_pt, 1)[:, 0]
        matched_reco_ele_tau_eta = ak.pad_none(matched_reco_ele_tau_eta, 1)[:, 0]
        matched_reco_ele_tau_phi = ak.pad_none(matched_reco_ele_tau_phi, 1)[:, 0]

        # --- Identify relevant particles ---
        is_mu1   = (abs(gen.pdgId) == 13) & (gen.status == 1)
        is_mu23  = (abs(gen.pdgId) == 13) & (gen.status == 23)

        # --- Muons directly from tau_from_W ---
        mu1_from_tauW = is_mu1 & has_mother & tau_from_W[mother]

        # --- Status=23 muons directly from tau_from_W ---
        mu23_from_tauW = is_mu23 & has_mother & tau_from_W[mother]

        # --- Final state μ± from those status=23 muons ---
        mu1_from_mu23_tauW = (
            is_mu1 & has_mother & mu23_from_tauW[mother]
        )

        # --- Combine both cases ---
        mu_tauW = gen[mu1_from_tauW | mu1_from_mu23_tauW]

        gen_mu_tauW_pt  = mu_tauW.pt
        gen_mu_tauW_eta = mu_tauW.eta
        gen_mu_tauW_phi = mu_tauW.phi

        gen_mu_tauW_pt  = ak.pad_none(gen_mu_tauW_pt, 1)[:, 0]
        gen_mu_tauW_eta = ak.pad_none(gen_mu_tauW_eta, 1)[:, 0]
        gen_mu_tauW_phi = ak.pad_none(gen_mu_tauW_phi, 1)[:, 0]

        # pairwise ΔR between each gen-muon and all reco muons
        dr_mu_tau = delta_r(
            mu_tauW.eta[:, :, None],   # (events, ngen, 1)
            mu_tauW.phi[:, :, None],
            reco_muons.eta[:, None, :],  # (events, 1, nreco)
            reco_muons.phi[:, None, :]
        )

        # find best match per GEN muon
        best_reco_idx_mu_tau = ak.argmin(dr_mu_tau, axis=2)
        best_dr_mu_tau = ak.min(dr_mu_tau, axis=2)

        # mask with ΔR requirement
        matched_mask_mu_tau = best_dr_mu_tau < dr_cut

        # select matched reco muons
        matched_reco_mu_tau = reco_muons[best_reco_idx_mu_tau]
        matched_reco_mu_tau = ak.mask(matched_reco_mu_tau, matched_mask_mu_tau)

        matched_reco_mu_tau_pt = matched_reco_mu_tau.pt
        matched_reco_mu_tau_eta = matched_reco_mu_tau.eta
        matched_reco_mu_tau_phi = matched_reco_mu_tau.phi   

        matched_reco_mu_tau_pt = ak.pad_none(matched_reco_mu_tau_pt, 1)[:, 0]
        matched_reco_mu_tau_eta = ak.pad_none(matched_reco_mu_tau_eta, 1)[:, 0]
        matched_reco_mu_tau_phi = ak.pad_none(matched_reco_mu_tau_phi, 1)[:, 0]

        output["nEvents"] += column_accumulator(ak.to_numpy(nEvents))
        output["pu_true"] += column_accumulator(ak.to_numpy(pu_true))
        output["nPhoton"] += column_accumulator(ak.to_numpy(nPhoton))
        output["nPhoton_cut_10"] += column_accumulator(ak.to_numpy(nPhoton_cut_10))
        output["nPhoton_cut_18"] += column_accumulator(ak.to_numpy(nPhoton_cut_18))
        output["nPhoton_cut_30"] += column_accumulator(ak.to_numpy(nPhoton_cut_30))
        output["higgs_pt"] += column_accumulator(ak.to_numpy(higgs_pt))
        output["higgs_eta"] += column_accumulator(ak.to_numpy(higgs_eta))
        output["higgs_phi"] += column_accumulator(ak.to_numpy(higgs_phi))
        output["W_pt"] += column_accumulator(ak.to_numpy(W_pt))
        output["W_eta"] += column_accumulator(ak.to_numpy(W_eta))
        output["W_phi"] += column_accumulator(ak.to_numpy(W_phi))

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

        output["Reco_lead_pho_pt"] += column_accumulator(ak.to_numpy(Reco_lead_pho_pt))
        output["Reco_sublead_pho_pt"] += column_accumulator(ak.to_numpy(Reco_sublead_pho_pt))
        output["Reco_lead_pho_eta"] += column_accumulator(ak.to_numpy(Reco_lead_pho_eta))
        output["Reco_sublead_pho_eta"] += column_accumulator(ak.to_numpy(Reco_sublead_pho_eta))
        output["Reco_lead_pho_phi"] += column_accumulator(ak.to_numpy(Reco_lead_pho_phi))
        output["Reco_sublead_pho_phi"] += column_accumulator(ak.to_numpy(Reco_sublead_pho_phi))

        output["gen_photon_from_a_1_pt"] += column_accumulator(ak.to_numpy(pho_from_a_pt_1))
        output["gen_photon_from_a_2_pt"] += column_accumulator(ak.to_numpy(pho_from_a_pt_2))

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
        output["gen_b1_pt"]  += column_accumulator(ak.to_numpy(gen_b_pt_1))
        output["gen_b1_eta"] += column_accumulator(ak.to_numpy(gen_b_eta_1))
        output["gen_b1_phi"] += column_accumulator(ak.to_numpy(gen_b_phi_1))

        output["gen_b2_pt"]  += column_accumulator(ak.to_numpy(gen_b_pt_2))
        output["gen_b2_eta"] += column_accumulator(ak.to_numpy(gen_b_eta_2))
        output["gen_b2_phi"] += column_accumulator(ak.to_numpy(gen_b_phi_2))

        # Reco b-jets (pt-sorted)
        output["reco_lead_bjet_pt"]  += column_accumulator(ak.to_numpy(Reco_lead_bjet_pt))
        output["reco_lead_bjet_eta"] += column_accumulator(ak.to_numpy(Reco_lead_bjet_eta))
        output["reco_lead_bjet_phi"] += column_accumulator(ak.to_numpy(Reco_lead_bjet_phi))  # Note: typo? should be `Reco_lead_bjet_phi`?

        output["reco_sublead_bjet_pt"]  += column_accumulator(ak.to_numpy(Reco_sublead_bjet_pt))
        output["reco_sublead_bjet_eta"] += column_accumulator(ak.to_numpy(Reco_sublead_bjet_eta))
        output["reco_sublead_bjet_phi"] += column_accumulator(ak.to_numpy(Reco_sublead_bjet_phi))

        # Matched reco b-jets
        output["matched_bjet1_pt"]  += column_accumulator(ak.to_numpy(Genmatched_bjet1_pt))
        output["matched_bjet1_eta"] += column_accumulator(ak.to_numpy(Genmatched_bjet1_eta))
        output["matched_bjet1_phi"] += column_accumulator(ak.to_numpy(Genmatched_bjet1_phi))

        output["matched_bjet2_pt"]  += column_accumulator(ak.to_numpy(Genmatched_bjet2_pt))
        output["matched_bjet2_eta"] += column_accumulator(ak.to_numpy(Genmatched_bjet2_eta))
        output["matched_bjet2_phi"] += column_accumulator(ak.to_numpy(Genmatched_bjet2_phi))

        output["gen_invmasses_diphoton"] += column_accumulator(gen_invmasses_diphoton)
        output["gen_invmasses_bb"] += column_accumulator(gen_invmasses_bb)

        output["invmasses_diphoton"] += column_accumulator(invmasses_photons)
        output["invmasses_bb"] += column_accumulator(invmasses_bb)

        output["gen_ele_pt"] += column_accumulator(ak.to_numpy(gen_ele_pt))
        output["gen_ele_eta"] += column_accumulator(ak.to_numpy(gen_ele_eta))
        output["gen_ele_phi"] += column_accumulator(ak.to_numpy(gen_ele_phi))
        output["matched_reco_ele_pt"] += column_accumulator(ak.to_numpy(matched_reco_ele_pt))
        output["matched_reco_ele_eta"] += column_accumulator(ak.to_numpy(matched_reco_ele_eta))
        output["matched_reco_ele_phi"] += column_accumulator(ak.to_numpy(matched_reco_ele_phi))
        output["gen_mu_pt"] += column_accumulator(ak.to_numpy(gen_mu_pt))
        output["gen_mu_eta"] += column_accumulator(ak.to_numpy(gen_mu_eta))
        output["gen_mu_phi"] += column_accumulator(ak.to_numpy(gen_mu_phi))
        output["matched_reco_mu_pt"] += column_accumulator(ak.to_numpy(matched_reco_mu_pt))
        output["matched_reco_mu_eta"] += column_accumulator(ak.to_numpy(matched_reco_mu_eta))
        output["matched_reco_mu_phi"] += column_accumulator(ak.to_numpy(matched_reco_mu_phi))
        output["gen_ele_tauW_pt"] += column_accumulator(ak.to_numpy(gen_ele_tauW_pt))
        output["gen_ele_tauW_eta"] += column_accumulator(ak.to_numpy(gen_ele_tauW_eta))
        output["gen_ele_tauW_phi"] += column_accumulator(ak.to_numpy(gen_ele_tauW_phi))
        output["matched_reco_ele_tau_pt"] += column_accumulator(ak.to_numpy(matched_reco_ele_tau_pt))
        output["matched_reco_ele_tau_eta"] += column_accumulator(ak.to_numpy(matched_reco_ele_tau_eta))
        output["matched_reco_ele_tau_phi"] += column_accumulator(ak.to_numpy(matched_reco_ele_tau_phi))
        output["gen_mu_tauW_pt"] += column_accumulator(ak.to_numpy(gen_mu_tauW_pt))
        output["gen_mu_tauW_eta"] += column_accumulator(ak.to_numpy(gen_mu_tauW_eta))
        output["gen_mu_tauW_phi"] += column_accumulator(ak.to_numpy(gen_mu_tauW_phi))
        output["matched_reco_mu_tau_pt"] += column_accumulator(ak.to_numpy(matched_reco_mu_tau_pt))
        output["matched_reco_mu_tau_eta"] += column_accumulator(ak.to_numpy(matched_reco_mu_tau_eta))
        output["matched_reco_mu_tau_phi"] += column_accumulator(ak.to_numpy(matched_reco_mu_tau_phi))

        # Wrap output inside a dict keyed by dataset name!
        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator