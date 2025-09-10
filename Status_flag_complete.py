from tqdm import tqdm   # progress bar
import uproot
import awkward as ak
import numpy as np

fname = "86D3D6A3-EE02-5849-8FC0-C89B9BE7D930.root"
with uproot.open(fname, timeout=120) as Hfile:
    print(Hfile.keys())
    Tree = Hfile["Events"]  # can access TTrees by name
    Events = Tree.arrays(library="ak", how="zip")

gen = Events.GenPart
photons = gen[gen.pdgId == 22]



def find_status_if_from_A_from_Higgs(event_gen, idx, A_pdgid=35, H_pdgid=25):
    """
    Return statusFlags if photon comes from an A whose immediate mother is a Higgs.
    Else return None.
    """
    mother = event_gen[idx].genPartIdxMother
    while mother >= 0:
        if event_gen[mother].pdgId == A_pdgid:
            # Check only the direct parent of A
            mother_of_A = event_gen[mother].genPartIdxMother
            if mother_of_A >= 0 and event_gen[mother_of_A].pdgId == H_pdgid:
                return int(event_gen[idx].statusFlags)
            return None
        mother = event_gen[mother].genPartIdxMother
    return None
# Define bit meanings
bit_meanings = {
    0: "isPrompt",
    1: "isDecayedLeptonHadron",
    2: "isTauDecayProduct",
    3: "isPromptTauDecayProduct",
    4: "isDirectTauDecayProduct",
    5: "isDirectPromptTauDecayProduct",
    6: "isDirectHadronDecayProduct",
    7: "isHardProcess",
    8: "fromHardProcess",
    9: "isHardProcessTauDecayProduct",
    10: "isDirectHardProcessTauDecayProduct",
    11: "fromHardProcessBeforeFSR",
    12: "isFirstCopy",
    13: "isLastCopy",
    14: "isLastCopyBeforeFSR",
}

# Example input list of decimal numbers
decimal_numbers = [1, 5, 42, 511, 8192, 16383]

# Write results to file
with open("gen_status_flags.txt", "w") as f:
    for num in decimal_numbers:
        binary_str = format(num, "015b")  # 15-bit binary
        f.write(f"Decimal: {num}\n")
        f.write(f"Binary : {binary_str}\n")
        f.write("Properties:\n")
        
        for bit, meaning in bit_meanings.items():
            if num & (1 << bit):
                f.write(f"  - {meaning}\n")
        f.write("\n")


status_from_AH = []
for genev in tqdm(gen, desc="Processing events"):
    for i, part in enumerate(genev):
        if part.pdgId == 22:  # photon
            sflag = find_status_if_from_A_from_Higgs(genev, i)
            if sflag is not None:
                status_from_AH.append(sflag)

# Convert to numpy and get unique
status_from_AH = np.array(status_from_AH)
unique_flags = np.unique(status_from_AH)

print("Unique statusFlags for photons from A->H:", unique_flags)

