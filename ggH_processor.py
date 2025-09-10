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

            "Gen_photon_pt_1": column_accumulator(np.array([], dtype=np.float32)),
            "Gen_photon_pt_2": column_accumulator(np.array([], dtype=np.float32)),

            "Genmatched_pho_1_pt": column_accumulator(np.array([], dtype=np.float32)),
            "Genmatched_pho_2_pt": column_accumulator(np.array([], dtype=np.float32)),
            "Genmatched_pho_1_eta": column_accumulator(np.array([], dtype=np.float32)),
            "Genmatched_pho_2_eta": column_accumulator(np.array([], dtype=np.float32)),
            "Genmatched_pho_1_phi": column_accumulator(np.array([], dtype=np.float32)),
            "Genmatched_pho_2_phi": column_accumulator(np.array([], dtype=np.float32)),

            "gen_invmasses_diphoton": column_accumulator(np.array([], dtype=np.float32)),
            "invmasses_diphoton": column_accumulator(np.array([], dtype=np.float32)),
        })

    def accumulator(self):
        # Return the accumulator prototype for coffea
        return self._accumulator

    def process(self, events):

        # Events = events

        # dataset = Events.metadata["dataset"]  # get dataset name dynamically

        # output = self.accumulator()
        # # nPhoton = Events["nPhoton"]

        # nEvents = Events.event
        
        # pu_true = Events.Pileup.nTrueInt
        # # nPhoton = Events.Photon.nPhoton
        # nPhoton = ak.num(events.Photon)

        # pt_cut_10 = Events.Photon.pt>10.0
        # pt_cut_18 = Events.Photon.pt>18.0
        # pt_cut_30 = Events.Photon.pt>30.0

        # nPhoton_cut_10 = ak.num(Events.Photon[pt_cut_10].pt)
        # nPhoton_cut_18 = ak.num(Events.Photon[pt_cut_18].pt)
        # nPhoton_cut_30 = ak.num(Events.Photon[pt_cut_30].pt)

        # gen = Events.GenPart
        # higgs_mask = (gen.pdgId == 25) & (gen.status == 62)
        # higgs = gen[higgs_mask]

        # higgs_pt = ak.flatten(higgs.pt)
        # higgs_eta = ak.flatten(higgs.eta)
        # higgs_phi = ak.flatten(higgs.phi)

        # is_A = (abs(gen.pdgId) == 35)

        # A = gen[is_A]

        # A_pt = A.pt
        # A_eta = A.eta
        # A_phi = A.phi

        # A_pt_1 = A_pt[:,0]
        # A_pt_2 = A_pt[:,1]
        # A_eta_1 = A_eta[:,0]
        # A_eta_2 = A_eta[:,1]
        # A_phi_1 = A_phi[:,0]
        # A_phi_2 = A_phi[:,1]

        # sorted_As = A[ak.argsort(A_pt, axis=1, ascending=False)]

        # leading_A_pt = sorted_As.pt[:,0]
        # subleading_A_pt = sorted_As.pt[:,1]
        # leading_A_eta = sorted_As.eta[:,0]
        # subleading_A_eta = sorted_As.eta[:,1]
        # leading_A_phi = sorted_As.phi[:,0]
        # subleading_A_phi = sorted_As.phi[:,1]

        # photon = gen[gen.pdgId == 22]
        # mother_idx = photon.genPartIdxMother
        # mask = gen[mother_idx].pdgId == 35 
        # photons_from_a = photon[mask]
        # photons_from_a = photons_from_a[ak.num(A.pt)==2]

        # photons_from_a = photons_from_a[:, :2]

        # pho_from_a_pt = photons_from_a.pt
        # pho_from_a_eta = photons_from_a.eta
        # pho_from_a_phi = photons_from_a.phi

        # pho_from_a_pt_1 = pho_from_a_pt[:, 0]
        # pho_from_a_pt_2 = pho_from_a_pt[:, 1]
        # pho_from_a_eta_1 = pho_from_a_eta[:, 0]
        # pho_from_a_eta_2 = pho_from_a_eta[:, 1]
        # pho_from_a_phi_1 = pho_from_a_phi[:, 0]
        # pho_from_a_phi_2 = pho_from_a_phi[:, 1]

        # sorted_photons = photons_from_a[ak.argsort(photons_from_a.pt, axis=1, ascending=False)]
        # lead_pt_pho_gen = sorted_photons.pt[:, 0]
        # sublead_pt_pho_gen = sorted_photons.pt[:, 1]
        # lead_eta_pho_gen = sorted_photons.eta[:, 0]
        # sublead_eta_pho_gen = sorted_photons.eta[:, 1]
        # lead_phi_pho_gen = sorted_photons.phi[:, 0]
        # sublead_phi_pho_gen = sorted_photons.phi[:, 1]

        # reco_photons = Events.Photon[ak.num(A.pt)==2]
        # has_photon = ak.num(reco_photons) >= 1
        # Reco_photons_all = reco_photons[has_photon]
        # Reco_pho_pt = ak.flatten(reco_photons.pt)
        # Reco_pho_eta = ak.flatten(reco_photons.eta)
        # Reco_pho_phi = ak.flatten(reco_photons.phi)

        # Reco_pho_pt_uf = reco_photons.pt
        # Reco_pho_eta_uf = reco_photons.eta
        # Reco_pho_phi_uf = reco_photons.phi

        # sorted_reco_photons = reco_photons[ak.argsort(reco_photons.pt, axis=1, ascending=False)]
        # Sorted_reco_photons_all = Reco_photons_all[ak.argsort(Reco_photons_all.pt, axis=1, ascending=False)]

        # Reco_photon_lead_pt_all = Sorted_reco_photons_all.pt[:,0]
        # Reco_photon_sublead_pt_all = ak.pad_none(Sorted_reco_photons_all.pt, 2)[:, 1]
        # Reco_photon_lead_eta_all = Sorted_reco_photons_all.eta[:,0]
        # Reco_photon_sublead_eta_all = ak.pad_none(Sorted_reco_photons_all.eta, 2)[:, 1]
        # Reco_photon_lead_phi_all = Sorted_reco_photons_all.phi[:,0]
        # Reco_photon_sublead_phi_all = ak.pad_none(Sorted_reco_photons_all.phi, 2)[:, 1]

        # has_two_photons = ak.num(sorted_reco_photons.pt) > 1

        # # Use the same mask for both
        # Reco_lead_pho_pt     = sorted_reco_photons.pt[has_two_photons][:, 0]
        # Reco_sublead_pho_pt  = sorted_reco_photons.pt[has_two_photons][:, 1]
        # Reco_lead_pho_eta    = sorted_reco_photons.eta[has_two_photons][:, 0]
        # Reco_sublead_pho_eta = sorted_reco_photons.eta[has_two_photons][:, 1]
        # Reco_lead_pho_phi    = sorted_reco_photons.phi[has_two_photons][:, 0]
        # Reco_sublead_pho_phi = sorted_reco_photons.phi[has_two_photons][:, 1]

        # def dR(eta1, phi1, eta2, phi2):
        #     d_eta = eta1 - eta2
        #     d_phi = phi1 - phi2
        #     return np.sqrt(d_eta**2 + d_phi**2)

        # dR_pho_lead = dR(lead_eta_pho_gen, lead_phi_pho_gen, Reco_pho_eta_uf, Reco_pho_phi_uf)
        # dR_pho_sublead = dR(sublead_eta_pho_gen, sublead_phi_pho_gen, Reco_pho_eta_uf, Reco_pho_phi_uf)

        # min_idx_1 = ak.argmin(dR_pho_lead, axis=1)
        # min_idx_2 = ak.argmin(dR_pho_sublead, axis=1)

        # photon_idx = ak.local_index(reco_photons)
        # mask_idx_1 = photon_idx == min_idx_1[:, None]
        # mask_idx_2 = photon_idx == min_idx_2[:, None]

        # mask_dR_1 = dR_pho_lead < 0.1
        # mask_dR_2 = dR_pho_sublead < 0.1

        # mask_sel_1 = mask_idx_1 & mask_dR_1
        # mask_sel_2 = mask_idx_2 & mask_dR_2


        # selected_photon_1 = ak.firsts(reco_photons[mask_sel_1])
        # selected_photon_2 = ak.firsts(reco_photons[mask_sel_2])


        # Gen_photon_pt_1 = ak.to_numpy(pho_from_a_pt_1[~ak.is_none(selected_photon_1.pt)])  # These are not all the Gen photons
        # Gen_photon_pt_2 = ak.to_numpy(pho_from_a_pt_2[~ak.is_none(selected_photon_2.pt)])

        # Genmatched_pho_1_pt  = ak.to_numpy(ak.fill_none(selected_photon_1.pt,  np.nan))
        # Genmatched_pho_2_pt  = ak.to_numpy(ak.fill_none(selected_photon_2.pt,  np.nan))
        # Genmatched_pho_1_eta = ak.to_numpy(ak.fill_none(selected_photon_1.eta, np.nan))
        # Genmatched_pho_2_eta = ak.to_numpy(ak.fill_none(selected_photon_2.eta, np.nan))
        # Genmatched_pho_1_phi = ak.to_numpy(ak.fill_none(selected_photon_1.phi, np.nan))
        # Genmatched_pho_2_phi = ak.to_numpy(ak.fill_none(selected_photon_2.phi, np.nan))

        # both_matched_mask = ~ak.is_none(selected_photon_1.pt) & ~ak.is_none(selected_photon_2.pt)


        # Genmatched_photon_1_clean_both = selected_photon_1[both_matched_mask]
        # Genmatched_photon_2_clean_both = selected_photon_2[both_matched_mask]

        # Genmatched_pho_1_pt_both  = ak.to_numpy(Genmatched_photon_1_clean_both.pt)
        # Genmatched_pho_2_pt_both  = ak.to_numpy(Genmatched_photon_2_clean_both.pt)
        # Genmatched_pho_1_eta_both = ak.to_numpy(Genmatched_photon_1_clean_both.eta)
        # Genmatched_pho_2_eta_both = ak.to_numpy(Genmatched_photon_2_clean_both.eta)
        # Genmatched_pho_1_phi_both = ak.to_numpy(Genmatched_photon_1_clean_both.phi)
        # Genmatched_pho_2_phi_both = ak.to_numpy(Genmatched_photon_2_clean_both.phi)
        # Genmatched_pho_1_mass_both = np.zeros_like(Genmatched_pho_2_phi_both, dtype=np.float32)
        # Genmatched_pho_2_mass_both = np.zeros_like(Genmatched_pho_2_phi_both, dtype=np.float32)

        # # Allocate output array
        # invmasses_diphoton = np.empty(len(Genmatched_pho_1_pt_both), dtype=np.float32)

        # # TLorentzVectors
        # vec1 = TLorentzVector()
        # vec2 = TLorentzVector()

        # for i in range(len(Genmatched_pho_1_pt_both)):
        #     vec1.SetPtEtaPhiM(
        #         Genmatched_pho_1_pt_both[i],
        #         Genmatched_pho_1_eta_both[i],
        #         Genmatched_pho_1_phi_both[i],
        #         Genmatched_pho_1_mass_both[i],
        #     )
        #     vec2.SetPtEtaPhiM(
        #         Genmatched_pho_2_pt_both[i],
        #         Genmatched_pho_2_eta_both[i],
        #         Genmatched_pho_2_phi_both[i],
        #         Genmatched_pho_2_mass_both[i],
        #     )
        #     invmasses_diphoton[i] = (vec1 + vec2).M()

        # gen_invmasses_diphoton = np.empty(len(pho_from_a_pt_1), dtype=np.float32)

        # pho_from_a_mass_1 = np.zeros_like(pho_from_a_pt_1, dtype=np.float32)
        # pho_from_a_mass_2 = np.zeros_like(pho_from_a_pt_1, dtype=np.float32)

        # # Masks for valid (non-NaN) gen-level photons and b-quarks
        # valid_mask_gg = ~np.isnan(pho_from_a_pt_1 + pho_from_a_eta_1 + pho_from_a_phi_1 + pho_from_a_mass_1 +
        #                         pho_from_a_pt_2 + pho_from_a_eta_2 + pho_from_a_phi_2 + pho_from_a_mass_2)

        # # Filter the arrays
        # pho_1_pt = pho_from_a_pt_1[valid_mask_gg]
        # pho_1_eta = pho_from_a_eta_1[valid_mask_gg]
        # pho_1_phi = pho_from_a_phi_1[valid_mask_gg]
        # pho_1_mass = pho_from_a_mass_1[valid_mask_gg]

        # pho_2_pt = pho_from_a_pt_2[valid_mask_gg]
        # pho_2_eta = pho_from_a_eta_2[valid_mask_gg]
        # pho_2_phi = pho_from_a_phi_2[valid_mask_gg]
        # pho_2_mass = pho_from_a_mass_2[valid_mask_gg]

        # # Prepare output arrays
        # gen_invmasses_diphoton = np.empty(len(pho_1_pt), dtype=np.float32)

        # p4_pho_1 = TLorentzVector()
        # p4_pho_2 = TLorentzVector()


        # for i, (PT1, ETA1, PHI1, M1, PT2, ETA2, PHI2, M2) in enumerate(zip(pho_1_pt, pho_1_eta, pho_1_phi, pho_1_mass,
        #                                                                 pho_2_pt, pho_2_eta, pho_2_phi, pho_2_mass)):
        #     p4_pho_1.SetPtEtaPhiM(PT1, ETA1, PHI1, M1)
        #     p4_pho_2.SetPtEtaPhiM(PT2, ETA2, PHI2, M2)
        #     gen_invmasses_diphoton[i] = (p4_pho_1 + p4_pho_2).M()

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

        higgs_pt = ak.flatten(higgs.pt)
        higgs_eta = ak.flatten(higgs.eta)
        higgs_phi = ak.flatten(higgs.phi)

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
        photons_from_a = photons_from_a[:, :2]

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

        pho_from_a_mass_1 = np.zeros_like(pho_from_a_pt_1, dtype=np.float32)
        pho_from_a_mass_2 = np.zeros_like(pho_from_a_pt_1, dtype=np.float32)

        # Masks for valid (non-NaN) gen-level photons and b-quarks
        valid_mask_gg = ~np.isnan(pho_from_a_pt_1 + pho_from_a_eta_1 + pho_from_a_phi_1 + pho_from_a_mass_1 +
                                pho_from_a_pt_2 + pho_from_a_eta_2 + pho_from_a_phi_2 + pho_from_a_mass_2)

        # Filter the arrays
        pho_1_pt = pho_from_a_pt_1[valid_mask_gg]
        pho_1_eta = pho_from_a_eta_1[valid_mask_gg]
        pho_1_phi = pho_from_a_phi_1[valid_mask_gg]
        pho_1_mass = pho_from_a_mass_1[valid_mask_gg]

        pho_2_pt = pho_from_a_pt_2[valid_mask_gg]
        pho_2_eta = pho_from_a_eta_2[valid_mask_gg]
        pho_2_phi = pho_from_a_phi_2[valid_mask_gg]
        pho_2_mass = pho_from_a_mass_2[valid_mask_gg]

        # Prepare output arrays
        gen_invmasses_diphoton = np.empty(len(pho_1_pt), dtype=np.float32)

        p4_pho_1 = TLorentzVector()
        p4_pho_2 = TLorentzVector()


        for i, (PT1, ETA1, PHI1, M1, PT2, ETA2, PHI2, M2) in enumerate(zip(pho_1_pt, pho_1_eta, pho_1_phi, pho_1_mass,
                                                                        pho_2_pt, pho_2_eta, pho_2_phi, pho_2_mass)):
            p4_pho_1.SetPtEtaPhiM(PT1, ETA1, PHI1, M1)
            p4_pho_2.SetPtEtaPhiM(PT2, ETA2, PHI2, M2)
            gen_invmasses_diphoton[i] = (p4_pho_1 + p4_pho_2).M()

        output["nEvents"] += column_accumulator(ak.to_numpy(nEvents))
        output["pu_true"] += column_accumulator(ak.to_numpy(pu_true))
        output["nPhoton"] += column_accumulator(ak.to_numpy(nPhoton))
        output["nPhoton_cut_10"] += column_accumulator(ak.to_numpy(nPhoton_cut_10))
        output["nPhoton_cut_18"] += column_accumulator(ak.to_numpy(nPhoton_cut_18))
        output["nPhoton_cut_30"] += column_accumulator(ak.to_numpy(nPhoton_cut_30))
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

        output["gen_invmasses_diphoton"] += column_accumulator(gen_invmasses_diphoton)

        output["invmasses_diphoton"] += column_accumulator(invmasses_photons)

        # Wrap output inside a dict keyed by dataset name!
        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator