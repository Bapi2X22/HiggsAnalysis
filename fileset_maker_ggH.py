import subprocess

dataset_file = "datasets_ggH.txt"
output_file = "fileset_ggH.txt"

# Read dataset list from a text file (NanoAOD datasets directly)
with open(dataset_file) as f:
    dataset_list = [line.strip() for line in f if line.strip()]

# Prepare separate dictionaries
fileset_global = {}
fileset_fnal   = {}
fileset_cnaf   = {}

for dataset in dataset_list:
    # Step 1: Extract mass from dataset
    try:
        mass = dataset.split("_MA-")[1].split("_")[0]
    except IndexError:
        print(f"Could not extract mass from dataset: {dataset}")
        continue

    # Step 2: Extract campaign string
    try:
        campaign = dataset.split("/")[2].split("-")[0]
    except IndexError:
        print(f"Could not extract campaign from dataset: {dataset}")
        continue

    # Step 3: Query files from NanoAOD dataset
    file_query = f"file dataset={dataset}"
    file_result = subprocess.run(
        ["dasgoclient", "-query", file_query],
        stdout=subprocess.PIPE, text=True
    )
    files = [line.strip() for line in file_result.stdout.splitlines() if line.strip()]
    if not files:
        print(f"No files found for {dataset}")
        continue

    # Step 4: Add to fileset dictionaries with different redirectors
    key = f"M{mass}_{campaign}"

    fileset_global[key] = [f"root://cms-xrd-global.cern.ch/{f}" for f in files]
    fileset_fnal[key]   = [f"root://cmsxrootd.fnal.gov/{f}" for f in files]
    fileset_cnaf[key]   = [f"root://xrootd-cms.infn.it/{f}" for f in files]

# Step 5: Write output to a Python file
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

print(f"Filesets saved to {output_file}")
