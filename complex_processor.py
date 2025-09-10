# from coffea import processor
# import awkward as ak
# import numpy as np
# import pandas as pd
# from coffea.processor import column_accumulator

# class HiggsAnalysisProcessor(processor.ProcessorABC):
#     def process(self, events):
#         output = {}

#         gen = events.GenPart
#         photons = gen[(gen.pdgId == 22) & (gen.status == 1)]
#         bquarks = gen[abs(gen.pdgId) == 5]

#         mother_idx_pho = photons.genPartIdxMother
#         mother_idx_b = bquarks.genPartIdxMother

#         from_a_photons = photons[gen[mother_idx_pho].pdgId == 35]
#         from_a_bquarks = bquarks[gen[mother_idx_b].pdgId == 35]

#         # Filter events with >= 2 photons/b-quarks
#         photon_mask = ak.num(from_a_photons.eta) >= 2
#         bquark_mask = ak.num(from_a_bquarks.eta) >= 2
#         common_mask = photon_mask & bquark_mask

#         photons_from_a = from_a_photons[common_mask]
#         bquarks_from_a = from_a_bquarks[common_mask]
#         reco_photons = events.Photon[common_mask]

#         # Sort
#         sorted_pho = photons_from_a[ak.argsort(photons_from_a.pt, axis=1, ascending=False)]
#         sorted_b = bquarks_from_a[ak.argsort(bquarks_from_a.pt, axis=1, ascending=False)]
#         sorted_reco = reco_photons[ak.argsort(reco_photons.pt, axis=1, ascending=False)]

#         # Leading and subleading gen photons
#         lead_pho = sorted_pho[:, 0]
#         sublead_pho = sorted_pho[:, 1]

#         # Reco photon dR matching
#         def dR(eta1, phi1, eta2, phi2):
#             d_eta = eta1 - eta2
#             d_phi = np.arctan2(np.sin(phi1 - phi2), np.cos(phi1 - phi2))
#             return np.sqrt(d_eta**2 + d_phi**2)

#         dR_1 = dR(lead_pho.eta, lead_pho.phi, reco_photons.eta, reco_photons.phi)
#         dR_2 = dR(sublead_pho.eta, sublead_pho.phi, reco_photons.eta, reco_photons.phi)

#         min_idx_1 = ak.argmin(dR_1, axis=1)
#         min_idx_2 = ak.argmin(dR_2, axis=1)

#         idxs = ak.local_index(reco_photons)
#         mask1 = (idxs == min_idx_1[:, None]) & (dR_1 < 0.1)
#         mask2 = (idxs == min_idx_2[:, None]) & (dR_2 < 0.1)

#         sel_pho_1 = ak.firsts(reco_photons[mask1])
#         sel_pho_2 = ak.firsts(reco_photons[mask2])

#         not_none_1 = ~ak.is_none(sel_pho_1.pt)
#         not_none_2 = ~ak.is_none(sel_pho_2.pt)

#         # output["GenPhoton_pt_1"] = ak.to_numpy(lead_pho.pt[not_none_1])
#         # output["GenPhoton_pt_2"] = ak.to_numpy(sublead_pho.pt[not_none_2])
#         # output["RecoPhoton_pt_1"] = ak.to_numpy(sel_pho_1.pt[not_none_1])
#         # output["RecoPhoton_pt_2"] = ak.to_numpy(sel_pho_2.pt[not_none_2])
#         # output["RecoPhoton_eta_1"] = ak.to_numpy(sel_pho_1.eta[not_none_1])
#         # output["RecoPhoton_eta_2"] = ak.to_numpy(sel_pho_2.eta[not_none_2])
#         # output["RecoPhoton_phi_1"] = ak.to_numpy(sel_pho_1.phi[not_none_1])
#         # output["RecoPhoton_phi_2"] = ak.to_numpy(sel_pho_2.phi[not_none_2])

#         output["GenPhoton_pt_1"]     = column_accumulator(ak.to_numpy(lead_pho.pt[not_none_1]))
#         output["GenPhoton_pt_2"]     = column_accumulator(ak.to_numpy(sublead_pho.pt[not_none_2]))
#         output["RecoPhoton_pt_1"]    = column_accumulator(ak.to_numpy(sel_pho_1.pt[not_none_1]))
#         output["RecoPhoton_pt_2"]    = column_accumulator(ak.to_numpy(sel_pho_2.pt[not_none_2]))
#         output["RecoPhoton_eta_1"]   = column_accumulator(ak.to_numpy(sel_pho_1.eta[not_none_1]))
#         output["RecoPhoton_eta_2"]   = column_accumulator(ak.to_numpy(sel_pho_2.eta[not_none_2]))
#         output["RecoPhoton_phi_1"]   = column_accumulator(ak.to_numpy(sel_pho_1.phi[not_none_1]))
#         output["RecoPhoton_phi_2"]   = column_accumulator(ak.to_numpy(sel_pho_2.phi[not_none_2]))

#         return output

#     def postprocess(self, accumulator):
#         return accumulator


import awkward as ak
import numpy as np
from coffea.processor import ProcessorABC, dict_accumulator
# from coffea.processor.accumulator import column_accumulator
from coffea.processor.accumulator import column_accumulator, dict_accumulator

class HiggsAnalysisProcessor(ProcessorABC):
    def __init__(self):
        self._accumulator = self.accumulator()

    # def accumulator(self):
    #     return dict_accumulator({
    #         "GenPhoton_pt_1": column_accumulator([]),
    #         "GenPhoton_pt_2": column_accumulator([]),
    #         "RecoPhoton_pt_1": column_accumulator([]),
    #         "RecoPhoton_pt_2": column_accumulator([]),
    #         "RecoPhoton_eta_1": column_accumulator([]),
    #         "RecoPhoton_eta_2": column_accumulator([]),
    #         "RecoPhoton_phi_1": column_accumulator([]),
    #         "RecoPhoton_phi_2": column_accumulator([]),
    #     })
    def accumulator(self):
        return {
            "GenPhoton_pt_1": column_accumulator(np.array([], dtype=np.float32)),
            "GenPhoton_pt_2": column_accumulator(np.array([], dtype=np.float32)),
            "RecoPhoton_pt_1": column_accumulator(np.array([], dtype=np.float32)),
            "RecoPhoton_pt_2": column_accumulator(np.array([], dtype=np.float32)),
            "RecoPhoton_eta_1": column_accumulator(np.array([], dtype=np.float32)),
            "RecoPhoton_eta_2": column_accumulator(np.array([], dtype=np.float32)),
            "RecoPhoton_phi_1": column_accumulator(np.array([], dtype=np.float32)),
            "RecoPhoton_phi_2": column_accumulator(np.array([], dtype=np.float32)),
            # Add more fields here if needed
        }

    def process(self, events):
        output = self.accumulator()

        # Gen photon pt (sorted by pt)
        gen_photons = ak.pad_none(events.GenPart[(events.GenPart.pdgId == 22) & (events.GenPart.hasFlags(["isHardProcess"]) == True)], 2, axis=1)
        # gen_photons = ak.sort(gen_photons, axis=1, ascending=False, stable=True)
        gen_photons = ak.sort(gen_photons, axis=-1, ascending=False, stable=True)


        gen_pho_pt_1 = ak.to_numpy(ak.fill_none(gen_photons[:, 0].pt, 0.0))
        gen_pho_pt_2 = ak.to_numpy(ak.fill_none(gen_photons[:, 1].pt, 0.0))

        # Reco photons (sorted by pt)
        reco_photons = ak.pad_none(events.Photon, 2, axis=1)
        reco_photons = ak.sort(reco_photons, axis=1, ascending=False, stable=True)

        reco_pt_1  = ak.to_numpy(ak.fill_none(reco_photons[:, 0].pt, 0.0))
        reco_pt_2  = ak.to_numpy(ak.fill_none(reco_photons[:, 1].pt, 0.0))
        reco_eta_1 = ak.to_numpy(ak.fill_none(reco_photons[:, 0].eta, 0.0))
        reco_eta_2 = ak.to_numpy(ak.fill_none(reco_photons[:, 1].eta, 0.0))
        reco_phi_1 = ak.to_numpy(ak.fill_none(reco_photons[:, 0].phi, 0.0))
        reco_phi_2 = ak.to_numpy(ak.fill_none(reco_photons[:, 1].phi, 0.0))

        # Fill accumulator
        output["GenPhoton_pt_1"] += column_accumulator(gen_pho_pt_1)
        output["GenPhoton_pt_2"] += column_accumulator(gen_pho_pt_2)
        output["RecoPhoton_pt_1"] += column_accumulator(reco_pt_1)
        output["RecoPhoton_pt_2"] += column_accumulator(reco_pt_2)
        output["RecoPhoton_eta_1"] += column_accumulator(reco_eta_1)
        output["RecoPhoton_eta_2"] += column_accumulator(reco_eta_2)
        output["RecoPhoton_phi_1"] += column_accumulator(reco_phi_1)
        output["RecoPhoton_phi_2"] += column_accumulator(reco_phi_2)

        return output

    def postprocess(self, accumulator):
        return accumulator
