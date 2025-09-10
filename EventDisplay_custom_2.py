import uproot
import awkward as ak
import numpy as np
import sys

def dump_events(input_file, output_file, max_events=10000):
    try:
        events = uproot.open(input_file)["Events"]
    except Exception as e:
        print(f"Error opening file: {e}")
        return

    # Load GenPart branches
    branches = events.arrays(
        [
            "GenPart_pdgId", "GenPart_status", "GenPart_statusFlags",
            "GenPart_genPartIdxMother", "GenPart_eta", "GenPart_mass",
            "GenPart_phi", "GenPart_pt"
        ],
        entry_stop=max_events
    )

    with open(output_file, "w") as f:
        print(f"# File used : {input_file}\n", file=f)

        # Loop over events
        for i in range(len(branches["GenPart_pdgId"])):
            pdgid   = branches["GenPart_pdgId"][i]
            status  = branches["GenPart_status"][i]
            mother  = branches["GenPart_genPartIdxMother"][i]

            # Select photons
            photon_mask = pdgid == 22
            photons_idx = ak.where(photon_mask)[0]

            if len(photons_idx) == 0:
                continue

            # Check mothers
            mother_idx = mother[photons_idx]
            # mother_pdg = pdgid[mother_idx]  # lookup mother pdgId
            statusFlags = branches["GenPart_statusFlags"][i]

            # Apply condition: from A (pdgId=35) and status==23
            # good_photons = (mother_pdg == 35) & (status[photons_idx] == 23)
            good_photons = (statusFlags[photons_idx] == 12308)

            if not ak.any(good_photons):
                continue  # skip events with no matching photon

            # --- If we reach here, event has at least one photon from A with status 23 ---
            print(f"\nEvent Number : {i}", file=f)
            print("Idx  PdgId  Status  StatusFlag  genPartIdxMother     eta    mass     phi      pt", file=f)
            print("---- ------ ------- ----------- ------------------ -------- ------- -------- -------", file=f)

            n_particles = len(pdgid)
            for j in range(n_particles):
                print(f"{j:>4} {pdgid[j]:>6} {status[j]:>7} "
                      f"{branches['GenPart_statusFlags'][i][j]:>11} "
                      f"{mother[j]:>18} "
                      f"{branches['GenPart_eta'][i][j]:>8.4f} "
                      f"{branches['GenPart_mass'][i][j]:>7.3f} "
                      f"{branches['GenPart_phi'][i][j]:>8.4f} "
                      f"{branches['GenPart_pt'][i][j]:>7.5f}",
                      file=f)

    print(f"\nEvent listing saved to: {output_file}")


# Main logic to parse command-line arguments
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python eventlisting_custom.py <input_root_file> <output_txt_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    dump_events(input_file, output_file)