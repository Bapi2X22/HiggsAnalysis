from coffea import processor
import awkward as ak
import numpy as np
from coffea.nanoevents import NanoAODSchema
'''

class MyAnalysisProcessor(processor.ProcessorABC):
    def __init__(self):
        # Here you would define your accumulator (output objects)
        # self._accumulator = processor.dict_accumulator({
        #     "lead_pho_pt": processor.column_accumulator([]),
        #     "sublead_pho_pt": processor.column_accumulator([]),
        #     "lead_bquark_pt": processor.column_accumulator([]),
        #     "sublead_bquark_pt": processor.column_accumulator([]),
        #     # Add more arrays as needed for all outputs (eta, phi, etc.)
        # })
        self._accumulator = processor.dict_accumulator({
            "lead_pho_pt": processor.column_accumulator(np.array([])),
            "sublead_pho_pt": processor.column_accumulator(np.array([])),
            "lead_bquark_pt": processor.column_accumulator(np.array([])),
            "sublead_bquark_pt": processor.column_accumulator(np.array([])),
            # ...and so on
        })


    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.accumulator.identity()

        # ========== Begin your selection and calculation code: ==========

        gen = events.GenPart
        photons = gen[(gen.pdgId == 22) & (gen.status == 1)]
        mother_idx = photons.genPartIdxMother
        from_a_mask = gen[mother_idx].pdgId == 35
        photons_from_a = photons[from_a_mask]
        mask_two_photons = ak.num(photons_from_a.eta) >= 2
        photons_from_a = photons_from_a[mask_two_photons]

        sorted_photons = photons_from_a[ak.argsort(photons_from_a.pt, axis=1, ascending=False)]
        lead_pt_pho_gen  = sorted_photons.pt[:, 0]
        sublead_pt_pho_gen = sorted_photons.pt[:, 1]

        bquarks = gen[abs(gen.pdgId) == 5]
        mother_idx_b = bquarks.genPartIdxMother
        from_a_mask_b = gen[mother_idx_b].pdgId == 35
        bquarks_from_a = bquarks[from_a_mask_b]
        sorted_bquarks = bquarks_from_a[ak.argsort(bquarks_from_a.pt, axis=1, ascending=False)]
        lead_pt_bquark_gen  = sorted_bquarks.pt[:, 0]
        sublead_pt_bquark_gen = sorted_bquarks.pt[:, 1]

        # Fill accumulator (extend with .tolist() or ak.to_numpy as needed)
        # output["lead_pho_pt"] += processor.column_accumulator(ak.to_numpy(lead_pt_pho_gen).tolist())
        # output["sublead_pho_pt"] += processor.column_accumulator(ak.to_numpy(sublead_pt_pho_gen).tolist())
        # output["lead_bquark_pt"] += processor.column_accumulator(ak.to_numpy(lead_pt_bquark_gen).tolist())
        # output["sublead_bquark_pt"] += processor.column_accumulator(ak.to_numpy(sublead_pt_bquark_gen).tolist())
        output["lead_pho_pt"].add(ak.to_numpy(lead_pt_pho_gen))
        output["sublead_pho_pt"].add(ak.to_numpy(sublead_pt_pho_gen))
        output["lead_bquark_pt"].add(ak.to_numpy(lead_pt_bquark_gen))
        output["sublead_bquark_pt"].add(ak.to_numpy(sublead_pt_bquark_gen))

        # Add your reco and matching analysis similarly

        return output

    def postprocess(self, accumulator):
        return accumulator
'''
'''

class TestProcessor(processor.ProcessorABC):
    def __init__(self):
        self._accumulator = processor.dict_accumulator({})

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        print("Fields in events:", events.fields)
        print("Events type:", type(events))
        print("GenPart count:", ak.num(events.GenPart))
        return self.accumulator.identity()

    def postprocess(self, accumulator):
        return accumulator
'''

from coffea import processor
import awkward as ak

class TestProcessor:
    def __init__(self):
        pass  # no accumulator needed unless you're returning something

    def process(self, events):
        print("Fields in events:", events.fields)
        print("Events type:", type(events))
        print("GenPart count:", ak.num(events.GenPart))
        return {"nEvents": len(events)}

    def postprocess(self, accumulator):
        return accumulator
