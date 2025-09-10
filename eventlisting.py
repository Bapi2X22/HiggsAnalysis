# import uproot
# import awkward as ak
# import numpy as np
# import sys

# def dump_events(input_file, output_file, max_events=20):
#     try:
#         events = uproot.open(input_file)["Events"]
#     except Exception as e:
#         print(f"Error opening file: {e}")
#         return

#     # Load GenPart branches
#     branches = events.arrays(
#         [
#             "GenPart_pdgId", "GenPart_status", "GenPart_statusFlags",
#             "GenPart_genPartIdxMother", "GenPart_eta", "GenPart_mass",
#             "GenPart_phi", "GenPart_pt"
#         ],
#         entry_stop=max_events
#     )

#     # Write to output file
#     with open(output_file, "w") as f:
#         for i in range(len(branches["GenPart_pdgId"])):
#             print(f"\nEvent Number : {i}", file=f)
#             print("Idx  PdgId  Status  StatusFlag  genPartIdxMother     eta    mass     phi      pt", file=f)
#             print("---- ------ ------- ----------- ------------------ -------- ------- -------- -------", file=f)

#             n_particles = len(branches["GenPart_pdgId"][i])
#             for j in range(n_particles):
#                 pdgid = branches["GenPart_pdgId"][i][j]
#                 status = branches["GenPart_status"][i][j]
#                 status_flag = branches["GenPart_statusFlags"][i][j]
#                 mother = branches["GenPart_genPartIdxMother"][i][j]
#                 eta = branches["GenPart_eta"][i][j]
#                 mass = branches["GenPart_mass"][i][j]
#                 phi = branches["GenPart_phi"][i][j]
#                 pt = branches["GenPart_pt"][i][j]

#                 print(f"{j:>4} {pdgid:>6} {status:>7} {status_flag:>11} {mother:>18} {eta:>8.4f} {mass:>7.3f} {phi:>8.4f} {pt:>7.3f}", file=f)

#     print(f"\nEvent listing saved to: {output_file}")


# # Main logic to parse command-line arguments
# if __name__ == "__main__":
#     if len(sys.argv) != 3:
#         print("Usage: python eventlisting_custom.py <input_root_file> <output_txt_file>")
#         sys.exit(1)

#     input_file = sys.argv[1]
#     output_file = sys.argv[2]

#     dump_events(input_file, output_file)

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

    # Write to output file
    with open(output_file, "w") as f:
        print(f"# File used : {input_file}\n", file=f)
        for i in range(len(branches["GenPart_pdgId"])):
            print(f"\nEvent Number : {i}", file=f)
            print("Idx  PdgId  Status  StatusFlag  genPartIdxMother     eta    mass     phi      pt", file=f)
            print("---- ------ ------- ----------- ------------------ -------- ------- -------- -------", file=f)

            n_particles = len(branches["GenPart_pdgId"][i])
            for j in range(n_particles):
                pdgid = branches["GenPart_pdgId"][i][j]
                status = branches["GenPart_status"][i][j]
                status_flag = branches["GenPart_statusFlags"][i][j]
                mother = branches["GenPart_genPartIdxMother"][i][j]
                eta = branches["GenPart_eta"][i][j]
                mass = branches["GenPart_mass"][i][j]
                phi = branches["GenPart_phi"][i][j]
                pt = branches["GenPart_pt"][i][j]

                print(f"{j:>4} {pdgid:>6} {status:>7} {status_flag:>11} {mother:>18} {eta:>8.4f} {mass:>7.3f} {phi:>8.4f} {pt:>7.5f}", file=f)

    print(f"\nEvent listing saved to: {output_file}")


# Main logic to parse command-line arguments
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python eventlisting_custom.py <input_root_file> <output_txt_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    dump_events(input_file, output_file)

