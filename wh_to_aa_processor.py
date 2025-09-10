from coffea import processor
from coffea.nanoevents import NanoAODSchema
import awkward as ak
import numpy as np


'''

class WHToAAProcessor(processor.ProcessorABC):
    def process(self, events):

        print("Started reading files.......")

        print("fields", events.fields)
        gen = events.GenPart

        # Higgs candidates
        higgs_mask = (gen.pdgId == 25) & (gen.status == 22)
        higgs = gen[higgs_mask]

        # Gen-level photons from a (pdgId=35)
        photons = gen[(gen.pdgId == 22) & (gen.status == 1)]
        mother_idx_pho = photons.genPartIdxMother
        from_a_mask_pho = gen[mother_idx_pho].pdgId == 35
        photons_from_a = photons[from_a_mask_pho]
        mask_two_photons = ak.num(photons_from_a.eta) >= 2
        photons_from_a = photons_from_a[mask_two_photons]
        sorted_photons = photons_from_a[ak.argsort(photons_from_a.pt, axis=1, ascending=False)]

        pho_from_a_pt_1 = sorted_photons.pt[:, 0]
        pho_from_a_pt_2 = sorted_photons.pt[:, 1]
        pho_from_a_eta_1 = sorted_photons.eta[:, 0]
        pho_from_a_eta_2 = sorted_photons.eta[:, 1]
        pho_from_a_phi_1 = sorted_photons.phi[:, 0]
        pho_from_a_phi_2 = sorted_photons.phi[:, 1]

        # Gen-level b-quarks from a
        bquarks = gen[abs(gen.pdgId) == 5]
        mother_idx_b = bquarks.genPartIdxMother
        from_a_mask_b = gen[mother_idx_b].pdgId == 35
        bquarks_from_a = bquarks[from_a_mask_b]
        sorted_bquarks = bquarks_from_a[ak.argsort(bquarks_from_a.pt, axis=1, ascending=False)]

        bquark_pt_1 = sorted_bquarks.pt[:, 0]
        bquark_pt_2 = sorted_bquarks.pt[:, 1]
        bquark_eta_1 = sorted_bquarks.eta[:, 0]
        bquark_eta_2 = sorted_bquarks.eta[:, 1]
        bquark_phi_1 = sorted_bquarks.phi[:, 0]
        bquark_phi_2 = sorted_bquarks.phi[:, 1]

        # Reco photons
        reco_photons = events.Photon[mask_two_photons]
        sorted_reco = reco_photons[ak.argsort(reco_photons.pt, axis=1, ascending=False)]

        has_lead = ak.num(sorted_reco.pt) > 0
        has_sublead = ak.num(sorted_reco.pt) > 1

        reco_lead_pt = sorted_reco.pt[has_lead][:, 0]
        reco_sublead_pt = sorted_reco.pt[has_sublead][:, 1]

        reco_eta = sorted_reco.eta
        reco_phi = sorted_reco.phi

        # ΔR matching
        def dR(eta1, phi1, eta2, phi2):
            d_eta = eta1 - eta2
            d_phi = phi1 - phi2
            return np.sqrt(d_eta**2 + d_phi**2)

        dR_pho_1 = dR(pho_from_a_eta_1, pho_from_a_phi_1, reco_eta, reco_phi)
        dR_pho_2 = dR(pho_from_a_eta_2, pho_from_a_phi_2, reco_eta, reco_phi)

        min_idx_1 = ak.argmin(dR_pho_1, axis=1)
        min_idx_2 = ak.argmin(dR_pho_2, axis=1)

        idxs = ak.local_index(reco_eta)
        mask_sel_1 = (idxs == min_idx_1[:, None]) & (dR_pho_1 < 0.1)
        mask_sel_2 = (idxs == min_idx_2[:, None]) & (dR_pho_2 < 0.1)

        selected_pho_1 = ak.firsts(reco_photons[mask_sel_1])
        selected_pho_2 = ak.firsts(reco_photons[mask_sel_2])

        # Output accumulators
        return {
            "gen_pho_pt_1": processor.column_accumulator(ak.to_numpy(pho_from_a_pt_1)),
            "gen_pho_pt_2": processor.column_accumulator(ak.to_numpy(pho_from_a_pt_2)),
            "gen_b_pt_1": processor.column_accumulator(ak.to_numpy(bquark_pt_1)),
            "gen_b_pt_2": processor.column_accumulator(ak.to_numpy(bquark_pt_2)),
            "matched_pho_pt_1": processor.column_accumulator(ak.to_numpy(selected_pho_1.pt[~ak.is_none(selected_pho_1.pt)])),
            "matched_pho_pt_2": processor.column_accumulator(ak.to_numpy(selected_pho_2.pt[~ak.is_none(selected_pho_2.pt)])),
        }

    def postprocess(self, accumulator):
        return accumulator

'''

class WHToAAProcessor(processor.ProcessorABC):
    def process(self, events):
        print("Reading a new chunk of events...")
        print("Fields available in this chunk:", events.fields)
        return {}

    def postprocess(self, accumulator):
        return accumulator
