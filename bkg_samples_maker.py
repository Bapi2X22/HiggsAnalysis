import json
import subprocess
import argparse

# def strip_dataset_name(ds_name):
#     """Strip dataset name like TTGJets_20UL16APV"""
#     process = ds_name.split("_Tune")[0].lstrip("/")
#     if "UL16NanoAODAPV" in ds_name:
#         year = "20UL16APV"
#     elif "UL16NanoAOD" in ds_name:
#         year = "20UL16"
#     elif "UL17NanoAOD" in ds_name:
#         year = "20UL17"
#     elif "UL18NanoAOD" in ds_name:
#         year = "20UL18"
#     else:
#         year = "Unknown"
#     return f"{process}_{year}"

def strip_dataset_name(ds_name):
    """Strip dataset name like TTGJets_20UL16APV or TTGJets_23BPix"""
    process = ds_name.split("_Tune")[0].lstrip("/")

    # Run 2 (UL) campaigns
    if "UL16NanoAODAPV" in ds_name:
        year = "20UL16APV"
    elif "UL16NanoAOD" in ds_name:
        year = "20UL16"
    elif "UL17NanoAOD" in ds_name:
        year = "20UL17"
    elif "UL18NanoAOD" in ds_name:
        year = "20UL18"

    # Run 3 campaigns
    elif "Run3Winter22NanoAOD" in ds_name:
        year = "22WinterRun3"
    elif "Run3Summer22NanoAOD" in ds_name:
        year = "22SummerRun3"
    elif "Run3Summer22EENanoAOD" in ds_name:
        year = "22EESummerRun3"
    elif "Run3Summer23NanoAOD" in ds_name:
        year = "23SummerRun3"
    elif "Run3Summer23BPixNanoAOD" in ds_name:
        year = "23BPixSummerRun3"
    elif "RunIII2024Summer24NanoAOD" in ds_name:
        year = "24SummerRun3"
    else:
        year = "Unknown"

    return f"{process}_{year}"


def get_files_from_das(dataset):
    """Return list of root files from dasgoclient for a dataset"""
    cmd = ["dasgoclient", "-query", f"file dataset={dataset}"]
    try:
        output = subprocess.check_output(cmd, universal_newlines=True)
        files = output.strip().split("\n")
        # Add xrootd prefix for HiggsDNA style
        files = [f"root://xrootd-cms.infn.it//{f}" for f in files if f]
        return files
    except subprocess.CalledProcessError as e:
        print(f"Error querying DAS for {dataset}: {e}")
        return []

def main():
    parser = argparse.ArgumentParser(description="Create HiggsDNA-style JSON for NanoAOD datasets.")
    parser.add_argument("-i", "--input", required=True, help="Input text file with dataset names (one per line).")
    parser.add_argument("-o", "--output", default="nanoaod_higgsdna.json", help="Output JSON file name.")
    args = parser.parse_args()

    # Read datasets from input file
    with open(args.input, "r") as f:
        datasets = [line.strip() for line in f if line.strip()]

    # Build the JSON dictionary
    dataset_dict = {}
    for ds in datasets:
        key = strip_dataset_name(ds)
        print(f"Processing dataset: {ds} -> {key}")
        files = get_files_from_das(ds)
        dataset_dict[key] = files

    # Save to JSON
    with open(args.output, "w") as f:
        json.dump(dataset_dict, f, indent=4)

    print(f"JSON file '{args.output}' created successfully!")

if __name__ == "__main__":
    main()

