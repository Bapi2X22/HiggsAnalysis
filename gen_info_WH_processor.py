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
            "higgs_pt": column_accumulator(np.array([], dtype=np.float32)),
            "higgs_eta": column_accumulator(np.array([], dtype=np.float32)),
            "higgs_phi": column_accumulator(np.array([], dtype=np.float32)),

            "A_pt_1": column_accumulator(np.array([], dtype=np.float32)),
            "A_pt_2": column_accumulator(np.array([], dtype=np.float32)),
            "A_eta_1": column_accumulator(np.array([], dtype=np.float32)),
            "A_eta_2": column_accumulator(np.array([], dtype=np.float32)),
            "A_phi_1": column_accumulator(np.array([], dtype=np.float32)),
            "A_phi_2": column_accumulator(np.array([], dtype=np.float32)),

            # Prompt photons from a (unsorted)
            "pho_from_a_pt_1": column_accumulator(np.array([], dtype=np.float32)),
            "pho_from_a_pt_2": column_accumulator(np.array([], dtype=np.float32)),
            "pho_from_a_eta_1": column_accumulator(np.array([], dtype=np.float32)),
            "pho_from_a_eta_2": column_accumulator(np.array([], dtype=np.float32)),
            "pho_from_a_phi_1": column_accumulator(np.array([], dtype=np.float32)),
            "pho_from_a_phi_2": column_accumulator(np.array([], dtype=np.float32)),

            # Prompt b-quarks from a (unsorted)
            "bquark_from_a_pt_1": column_accumulator(np.array([], dtype=np.float32)),
            "bquark_from_a_pt_2": column_accumulator(np.array([], dtype=np.float32)),
            "bquark_from_a_eta_1": column_accumulator(np.array([], dtype=np.float32)),
            "bquark_from_a_eta_2": column_accumulator(np.array([], dtype=np.float32)),
            "bquark_from_a_phi_1": column_accumulator(np.array([], dtype=np.float32)),
            "bquark_from_a_phi_2": column_accumulator(np.array([], dtype=np.float32)),
        })

    def accumulator(self):
        # Return the accumulator prototype for coffea
        return self._accumulator

    def process(self, events):

        Events = events

        dataset = Events.metadata["dataset"]  # get dataset name dynamically

        output = self.accumulator()

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

        A1_pt = A_pt[:, 0]
        A2_pt = A_pt[:, 1]
        A1_eta = A_eta[:, 0]
        A2_eta = A_eta[:, 1]
        A1_phi = A_phi[:, 0]
        A2_phi = A_phi[:, 1]

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

        output["higgs_pt"] += column_accumulator(ak.to_numpy(higgs_pt))
        output["higgs_eta"] += column_accumulator(ak.to_numpy(higgs_eta))
        output["higgs_phi"] += column_accumulator(ak.to_numpy(higgs_phi))

        output["A_pt_1"] += column_accumulator(ak.to_numpy(A1_pt))
        output["A_pt_2"] += column_accumulator(ak.to_numpy(A2_pt))
        output["A_eta_1"] += column_accumulator(ak.to_numpy(A1_eta))
        output["A_eta_2"] += column_accumulator(ak.to_numpy(A2_eta))
        output["A_phi_1"] += column_accumulator(ak.to_numpy(A1_phi))
        output["A_phi_2"] += column_accumulator(ak.to_numpy(A2_phi))

        output["pho_from_a_pt_1"] += column_accumulator(ak.to_numpy(pho_from_a_pt_1))
        output["pho_from_a_pt_2"] += column_accumulator(ak.to_numpy(pho_from_a_pt_2))
        output["pho_from_a_eta_1"] += column_accumulator(ak.to_numpy(pho_from_a_eta_1))
        output["pho_from_a_eta_2"] += column_accumulator(ak.to_numpy(pho_from_a_eta_2))
        output["pho_from_a_phi_1"] += column_accumulator(ak.to_numpy(pho_from_a_phi_1))
        output["pho_from_a_phi_2"] += column_accumulator(ak.to_numpy(pho_from_a_phi_2))

        output["bquark_from_a_pt_1"] += column_accumulator(ak.to_numpy(bquark_from_a_pt_1))
        output["bquark_from_a_pt_2"] += column_accumulator(ak.to_numpy(bquark_from_a_pt_2))
        output["bquark_from_a_eta_1"] += column_accumulator(ak.to_numpy(bquark_from_a_eta_1))
        output["bquark_from_a_eta_2"] += column_accumulator(ak.to_numpy(bquark_from_a_eta_2))
        output["bquark_from_a_phi_1"] += column_accumulator(ak.to_numpy(bquark_from_a_phi_1))
        output["bquark_from_a_phi_2"] += column_accumulator(ak.to_numpy(bquark_from_a_phi_2))

        # Wrap output inside a dict keyed by dataset name!
        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator