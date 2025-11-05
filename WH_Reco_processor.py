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
        })

    def accumulator(self):
        # Return the accumulator prototype for coffea
        return self._accumulator

    def process(self, events):

        eta_rho_corr = 1.5
        low_eta_rho_corr = 0.16544
        high_eta_rho_corr = 0.13212
        max_pho_iso_EB_low_r9 = 4.0
        max_pho_iso_EE_low_r9 = 4.0
        min_full5x5_r9_EB_high_r9 = 0.85
        min_full5x5_r9_EE_high_r9 = 0.9
        min_full5x5_r9_EB_low_r9 = 0.5
        min_full5x5_r9_EE_low_r9 = 0.8
        max_trkSumPtHollowConeDR03_EB_low_r9 = (6.0) 
        max_trkSumPtHollowConeDR03_EE_low_r9 = 6.0
        max_sieie_EB_low_r9 = 0.015
        max_sieie_EE_low_r9 = 0.035
        min_pt_photon = 25.0
        min_mvaid = -0.9
        max_hovere = 0.08
        min_full5x5_r9 = 0.8
        max_chad_iso = 20.0
        max_chad_rel_iso = 0.3

        Events = events

        dataset = Events.metadata["dataset"]  # get dataset name dynamically

        output = self.accumulator()
        # nPhoton = Events["nPhoton"]

        nEvents = Events.event

        # Don't drop events up front!
        photons = Events.Photon  # keep all events

        rho = events.Rho.fixedGridRhoAll  * ak.ones_like(photons.pt)
        photon_abs_eta = np.abs(photons.eta)

        # Run 2, use standard photon preselection
        pass_phoIso_rho_corr_EB = (
            (photon_abs_eta < eta_rho_corr)
            & (
                # photons.pfPhoIso03 - rho * self.low_eta_rho_corr
                photons.pfPhoIso03 - rho * low_eta_rho_corr
                < max_pho_iso_EB_low_r9
            )
        ) | (
            # should almost never happen because of the requirement of (photons.isScEtaEB) earlier, thus might be slightly redundant
            (photon_abs_eta > eta_rho_corr)
            & (
                # photons.pfPhoIso03 - rho * self.high_eta_rho_corr
                photons.pfPhoIso03 - rho * high_eta_rho_corr
                < max_pho_iso_EB_low_r9
            )
        )

        pass_phoIso_rho_corr_EE = (
            (photon_abs_eta < eta_rho_corr)
            & (
                photons.pfPhoIso03 - rho * low_eta_rho_corr
                < max_pho_iso_EE_low_r9
            )
        ) | (
            (photon_abs_eta > eta_rho_corr)
            & (
                photons.pfPhoIso03 - rho * high_eta_rho_corr
                < max_pho_iso_EE_low_r9
            )
        )

        isEB_high_r9 = (photons.isScEtaEB) & (photons.r9 > min_full5x5_r9_EB_high_r9)
        isEE_high_r9 = (photons.isScEtaEE) & (photons.r9 > min_full5x5_r9_EE_high_r9)
        iso = photons.trkSumPtHollowConeDR03 if hasattr(photons, "trkSumPtHollowConeDR03") else photons.pfChargedIsoPFPV  # photons.pfChargedIsoPFPV for v11, photons.trkSumPtHollowConeDR03 v12 and above
        rel_iso = photons.pfRelIso03_chg if hasattr(photons, "pfRelIso03_chg") else photons.pfRelIso03_chg_quadratic  # photons.pfChargedIsoPFPV for v1?, photons.pfChargedIsoPFPV_quadratic v12 and above
        isEB_low_r9 = (
            (photons.isScEtaEB)
            & (photons.r9 > min_full5x5_r9_EB_low_r9)
            & (photons.r9 < min_full5x5_r9_EB_high_r9)
            & (
                iso
                < max_trkSumPtHollowConeDR03_EB_low_r9
            )
            & (photons.sieie < max_sieie_EB_low_r9)
            & (pass_phoIso_rho_corr_EB)
        )
        isEE_low_r9 = (
            (photons.isScEtaEE)
            & (photons.r9 > min_full5x5_r9_EE_low_r9)
            & (photons.r9 < min_full5x5_r9_EE_high_r9)
            & (
                iso
                < max_trkSumPtHollowConeDR03_EE_low_r9
            )
            & (photons.sieie < max_sieie_EE_low_r9)
            & (pass_phoIso_rho_corr_EE)
        )

        e_veto_cut = (photons.electronVeto == 1)
        photons = photons[
            e_veto_cut
            & (photons.pt > min_pt_photon)
            & (photons.isScEtaEB | photons.isScEtaEE)
            & (photons.mvaID > min_mvaid)
            & (photons.hoe < max_hovere)
            & (
                (photons.r9 > min_full5x5_r9)
                | (
                    rel_iso * photons.pt < max_chad_iso
                )
                | (rel_iso < max_chad_rel_iso)
            )
            & (isEB_high_r9 | isEB_low_r9 | isEE_high_r9 | isEE_low_r9)
        ]

        Reco_pho_pt = ak.flatten(photons.pt)
        Reco_pho_eta = ak.flatten(photons.eta)
        Reco_pho_phi = ak.flatten(photons.phi)

        # Sort photons by pt
        sorted_reco_photons = photons[ak.argsort(photons.pt, axis=1, ascending=False)]

        # Pad so every event has at least 2 photons (None if missing)
        sorted_reco_padded = ak.pad_none(sorted_reco_photons, 2, axis=1, clip=True)

        # Fill missing with NaN
        Reco_lead_pho_pt  = ak.fill_none(sorted_reco_padded.pt[:, 0],  np.nan)
        Reco_sublead_pho_pt  = ak.fill_none(sorted_reco_padded.pt[:, 1],  np.nan)

        Reco_lead_pho_eta = ak.fill_none(sorted_reco_padded.eta[:, 0], np.nan)
        Reco_sublead_pho_eta = ak.fill_none(sorted_reco_padded.eta[:, 1], np.nan)

        Reco_lead_pho_phi = ak.fill_none(sorted_reco_padded.phi[:, 0], np.nan)
        Reco_sublead_pho_phi = ak.fill_none(sorted_reco_padded.phi[:, 1], np.nan)

        output["nEvents"] += column_accumulator(ak.to_numpy(nEvents))

        output["Reco_pho_pt"] += column_accumulator(ak.to_numpy(Reco_pho_pt))
        output["Reco_pho_eta"] += column_accumulator(ak.to_numpy(Reco_pho_eta))
        output["Reco_pho_phi"] += column_accumulator(ak.to_numpy(Reco_pho_phi))

        output["Reco_lead_pho_pt"] += column_accumulator(ak.to_numpy(Reco_lead_pho_pt))
        output["Reco_sublead_pho_pt"] += column_accumulator(ak.to_numpy(Reco_sublead_pho_pt))
        output["Reco_lead_pho_eta"] += column_accumulator(ak.to_numpy(Reco_lead_pho_eta))
        output["Reco_sublead_pho_eta"] += column_accumulator(ak.to_numpy(Reco_sublead_pho_eta))
        output["Reco_lead_pho_phi"] += column_accumulator(ak.to_numpy(Reco_lead_pho_phi))
        output["Reco_sublead_pho_phi"] += column_accumulator(ak.to_numpy(Reco_sublead_pho_phi))

        # Wrap output inside a dict keyed by dataset name!
        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator