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

        dr_cut = 0.1  # ΔR threshold
        pt_ratio_cut = 0.5  # optional (gen pt close to reco pt)
        
        # --- Select GEN electrons ---
        gen = Events.GenPart

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