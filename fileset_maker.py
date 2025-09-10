# import subprocess

# dataset_file = "datasets_ZH.txt"
# output_file = "fileset_ZH.txt"

# # Read dataset list from a text file
# with open(dataset_file) as f:
#     dataset_list = [line.strip() for line in f if line.strip()]

# fileset = {}

# for dataset in dataset_list:
#     # Step 1: Extract mass from parent dataset
#     mass = dataset.split("_M-")[1].split("_")[0]

#     # Step 2: Find child NanoAOD dataset
#     child_query = f"child dataset={dataset}"
#     child_result = subprocess.run(
#         ["dasgoclient", "-query", child_query],
#         stdout=subprocess.PIPE, text=True
#     )
#     child_datasets = [line.strip() for line in child_result.stdout.splitlines() if line.strip()]
#     if not child_datasets:
#         print(f"No NanoAOD child found for {dataset}")
#         continue

#     nano_dataset = [d for d in child_datasets if "NANOAODSIM" in d][0]

#     # Step 3: Extract NanoAOD campaign string
#     campaign = nano_dataset.split("/")[2].split("-")[0]

#     # Step 4: Query files from NanoAOD dataset
#     file_query = f"file dataset={nano_dataset}"
#     file_result = subprocess.run(
#         ["dasgoclient", "-query", file_query],
#         stdout=subprocess.PIPE, text=True
#     )
#     files = [
#         f"root://cms-xrd-global.cern.ch/{line.strip()}"
#         for line in file_result.stdout.splitlines() if line.strip()
#     ]

#     # Step 5: Add to fileset dictionary
#     key = f"M{mass}_{campaign}"
#     fileset[key] = files

# # Step 6: Write output to a file
# with open(output_file, "w") as out:
#     out.write("fileset = {\n")
#     for k, v in fileset.items():
#         out.write(f'    "{k}": [\n')
#         for file in v:
#             out.write(f'        "{file}",\n')
#         out.write("    ],\n")
#     out.write("}\n")

# print(f"fileset saved to {output_file}")

import subprocess

dataset_file = "datasets_ggH.txt"
output_file = "fileset_ggH.txt"

# Read dataset list from a text file
with open(dataset_file) as f:
    dataset_list = [line.strip() for line in f if line.strip()]

# Prepare separate dictionaries
fileset_global = {}
fileset_fnal   = {}
fileset_cnaf   = {}

for dataset in dataset_list:
    # Step 1: Extract mass from parent dataset
    mass = dataset.split("_M-")[1].split("_")[0]

    # Step 2: Find child NanoAOD dataset
    child_query = f"child dataset={dataset}"
    child_result = subprocess.run(
        ["dasgoclient", "-query", child_query],
        stdout=subprocess.PIPE, text=True
    )
    child_datasets = [line.strip() for line in child_result.stdout.splitlines() if line.strip()]
    if not child_datasets:
        print(f"No NanoAOD child found for {dataset}")
        continue

    nano_dataset = [d for d in child_datasets if "NANOAODSIM" in d][0]

    # Step 3: Extract NanoAOD campaign string
    campaign = nano_dataset.split("/")[2].split("-")[0]

    # Step 4: Query files from NanoAOD dataset
    file_query = f"file dataset={nano_dataset}"
    file_result = subprocess.run(
        ["dasgoclient", "-query", file_query],
        stdout=subprocess.PIPE, text=True
    )
    files = [line.strip() for line in file_result.stdout.splitlines() if line.strip()]

    # Step 5: Add to fileset dictionaries with different redirectors
    key = f"M{mass}_{campaign}"

    fileset_global[key] = [f"root://cms-xrd-global.cern.ch/{f}" for f in files]
    fileset_fnal[key]   = [f"root://cmsxrootd.fnal.gov/{f}" for f in files]
    fileset_cnaf[key]   = [f"root://xrootd-cms.infn.it/{f}" for f in files]

# Step 6: Write output to a Python file
with open(output_file, "w") as out:
    out.write("fileset_global = {\n")
    for k, v in fileset_global.items():
        out.write(f'    "{k}": [\n')
        for file in v:
            out.write(f'        "{file}",\n')
        out.write("    ],\n")
    out.write("}\n\n")

    out.write("fileset_fnal = {\n")
    for k, v in fileset_fnal.items():
        out.write(f'    "{k}": [\n')
        for file in v:
            out.write(f'        "{file}",\n')
        out.write("    ],\n")
    out.write("}\n\n")

    out.write("fileset_cnaf = {\n")
    for k, v in fileset_cnaf.items():
        out.write(f'    "{k}": [\n')
        for file in v:
            out.write(f'        "{file}",\n')
        out.write("    ],\n")
    out.write("}\n")

print(f"filesets saved to {output_file}")
