# import uproot
# import awkward as ak
# from particle import Particle
# import argparse
# from tqdm import tqdm

# def pdg_to_name(pdgid: int) -> str:
#     """Convert PDG ID to particle name if possible."""
#     try:
#         p = Particle.from_pdgid(pdgid)
#         return p.name
#     except Exception:
#         return str(pdgid)  # fallback if not found

# def build_children_map(mothers):
#     """Build a dict: parent_idx -> list of child indices for a single event."""
#     children = {i: [] for i in range(len(mothers))}
#     for i, mom in enumerate(mothers):
#         if mom >= 0:
#             children[mom].append(i)
#     return children

# def get_decay_tree(pdgIds, status, pts, etas, phis, masses, children_map, idx, level=0):
#     """Return tree-like decay structure with names, pdgId, status and 4-vector info."""
#     pdgId = int(pdgIds[idx])
#     st = int(status[idx])
#     name = pdg_to_name(pdgId)
#     indent = "  " * level

#     # kinematics
#     pt = float(pts[idx])
#     eta = float(etas[idx])
#     phi = float(phis[idx])
#     mass = float(masses[idx])

#     line = (
#         f"{indent}{name} [{pdgId}] (status={st}) "
#         f"pt={pt:.2f}, eta={eta:.2f}, phi={phi:.2f}, mass={mass:.2f}\n"
#     )

#     if idx in children_map:
#         for c in children_map[idx]:
#             line += get_decay_tree(
#                 pdgIds, status, pts, etas, phis, masses, children_map, c, level + 1
#             )

#     return line

# def dump_photon_chains(events, output_file="decay_chains.txt"):
#     pdgIds = events["GenPart_pdgId"]
#     status = events["GenPart_status"]
#     mothers = events["GenPart_genPartIdxMother"]
#     pts = events["GenPart_pt"]
#     etas = events["GenPart_eta"]
#     phis = events["GenPart_phi"]
#     masses = events["GenPart_mass"]

#     with open(output_file, "w") as f:
#         for evt_idx in tqdm(range(len(pdgIds)), desc="Processing events"):
#             pdg_evt = pdgIds[evt_idx]
#             st_evt = status[evt_idx]
#             mom_evt = mothers[evt_idx]
#             pt_evt = pts[evt_idx]
#             eta_evt = etas[evt_idx]
#             phi_evt = phis[evt_idx]
#             mass_evt = masses[evt_idx]

#             children_map = build_children_map(mom_evt)

#             # select photons in this event
#             photons = ak.where(pdg_evt == 22)[0]
#             for pho in photons:
#                 mom_idx = mom_evt[pho]
#                # if mom_idx >= 0 and pdg_evt[mom_idx] == 35 and st_evt[pho] == 23:
#                 if mom_idx >= 0 and pdg_evt[mom_idx] == 35:
#                     # climb up to top ancestor
#                     root = mom_idx
#                     while root >= 0 and mom_evt[root] >= 0:
#                         root = mom_evt[root]

#                     tree = get_decay_tree(
#                         pdg_evt, st_evt, pt_evt, eta_evt, phi_evt, mass_evt,
#                         children_map, root
#                     )
#                     f.write(f"\nEvent {evt_idx}\n")
#                     f.write(tree)

# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description="Dump photon decay chains from ROOT file")
#     parser.add_argument("input_file", help="Path to input ROOT file")
#     parser.add_argument("-o", "--output", default="decay_chains.txt", help="Output text file name (default: decay_chains.txt)")
#     args = parser.parse_args()

#     events = uproot.open(args.input_file)["Events"].arrays(
#         [
#             "GenPart_pdgId", "GenPart_status", "GenPart_genPartIdxMother",
#             "GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass"
#         ],
#         library="ak"
#     )

#     dump_photon_chains(events, output_file=args.output)

import uproot
import awkward as ak
from particle import Particle
import argparse
from tqdm import tqdm

def pdg_to_name(pdgid: int) -> str:
    """Convert PDG ID to particle name if possible."""
    try:
        p = Particle.from_pdgid(pdgid)
        return p.name
    except Exception:
        return str(pdgid)  # fallback if not found

def build_children_map(mothers):
    """Build a dict: parent_idx -> list of child indices for a single event."""
    children = {i: [] for i in range(len(mothers))}
    for i, mom in enumerate(mothers):
        if mom >= 0:
            children[mom].append(i)
    return children

def get_decay_tree(pdgIds, status, pts, etas, phis, masses, children_map, idx, level=0):
    """Return tree-like decay structure with names, pdgId, status and 4-vector info."""
    pdgId = int(pdgIds[idx])
    st = int(status[idx])
    name = pdg_to_name(pdgId)
    indent = "  " * level

    # kinematics
    pt = float(pts[idx])
    eta = float(etas[idx])
    phi = float(phis[idx])
    mass = float(masses[idx])

    line = (
        f"{indent}{name} [{pdgId}] (status={st}) "
        f"pt={pt:.2f}, eta={eta:.2f}, phi={phi:.2f}, mass={mass:.2f}\n"
    )

    if idx in children_map:
        for c in children_map[idx]:
            line += get_decay_tree(
                pdgIds, status, pts, etas, phis, masses, children_map, c, level + 1
            )

    return line

def dump_photon_chains(events, output_file="decay_chains.txt", require_status23=False, maxevents=None):
    pdgIds = events["GenPart_pdgId"]
    status = events["GenPart_status"]
    mothers = events["GenPart_genPartIdxMother"]
    pts = events["GenPart_pt"]
    etas = events["GenPart_eta"]
    phis = events["GenPart_phi"]
    masses = events["GenPart_mass"]

    n_events = len(pdgIds) if maxevents is None else min(maxevents, len(pdgIds))

    with open(output_file, "w") as f:
        for evt_idx in tqdm(range(n_events), desc="Processing events"):
            pdg_evt = pdgIds[evt_idx]
            st_evt = status[evt_idx]
            mom_evt = mothers[evt_idx]
            pt_evt = pts[evt_idx]
            eta_evt = etas[evt_idx]
            phi_evt = phis[evt_idx]
            mass_evt = masses[evt_idx]

            children_map = build_children_map(mom_evt)

            # select photons in this event
            photons = ak.where(pdg_evt == 22)[0]
            for pho in photons:
                mom_idx = mom_evt[pho]

                # if mom_idx >= 0 and pdg_evt[mom_idx] == 35:
                if mom_idx >= 0:
                    if require_status23 and st_evt[pho] != 23:
                        continue  # skip photons that are not status 23

                    # climb up to top ancestor
                    root = mom_idx
                    while root >= 0 and mom_evt[root] >= 0:
                        root = mom_evt[root]

                    tree = get_decay_tree(
                        pdg_evt, st_evt, pt_evt, eta_evt, phi_evt, mass_evt,
                        children_map, root
                    )
                    f.write(f"\nEvent {evt_idx}\n")
                    f.write(tree)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Dump photon decay chains from ROOT file")
    parser.add_argument("input_file", help="Path to input ROOT file")
    parser.add_argument("-o", "--output", default="decay_chains.txt", help="Output text file name (default: decay_chains.txt)")
    parser.add_argument("--require-status23", action="store_true", help="Require photons to have status == 23")
    parser.add_argument("--maxevents", type=int, default=None, help="Maximum number of events to process")

    args = parser.parse_args()

    events = uproot.open(args.input_file)["Events"].arrays(
        [
            "GenPart_pdgId", "GenPart_status", "GenPart_genPartIdxMother",
            "GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass"
        ],
        library="ak"
    )

    dump_photon_chains(events, output_file=args.output,
                       require_status23=args.require_status23,
                       maxevents=args.maxevents)


