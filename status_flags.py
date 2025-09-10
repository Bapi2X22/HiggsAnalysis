import uproot
import awkward as ak
from ROOT import TLorentzVector
import numpy as np


fname = "86D3D6A3-EE02-5849-8FC0-C89B9BE7D930.root"
with uproot.open(fname, timeout=120) as Hfile:
    print(Hfile.keys())
    Tree = Hfile["Events"]  # can access TTrees by name
    Events = Tree.arrays(library="ak", how="zip")

from tqdm import tqdm   # progress bar

gen = Events.GenPart
photons = gen[gen.pdgId == 22]

# Helper: check if photon comes from A
def find_status_if_from_A(event_gen, idx, target_pdgid=35):
    """Return statusFlags if photon comes from A, else None"""
    mother = event_gen[idx].genPartIdxMother
    while mother >= 0:
        if event_gen[mother].pdgId == target_pdgid:
            return int(event_gen[idx].statusFlags)
        mother = event_gen[mother].genPartIdxMother
    return None

# Loop per event with tqdm
status_from_A = []
for genev in tqdm(gen, desc="Processing events"):
    for i, part in enumerate(genev):
        if part.pdgId == 22:  # photon
            sflag = find_status_if_from_A(genev, i, 35)
            if sflag is not None:
                status_from_A.append(sflag)

# Convert to numpy and get unique
status_from_A = np.array(status_from_A)
unique_flags = np.unique(status_from_A)

print("Unique statusFlags for photons from A:", unique_flags)
