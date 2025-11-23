import pandas as pd

# Mass points
mass_points = [12, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]

# Data construction
rows = []
for M in mass_points:
    row = {
        "dataset": f"/ZHTo2A_ATo2B_ATo2G_M-{M}_TuneCP5_13p6TeV_madgraph-pythia8",
        "fragment": f"Fragments/ZH/Zh_M125_ToAA_{M}_Tobbgg_TuneCP5_PSweights_13p6TeV-madgraph_pythia8-fragment.py",
        "events": 1000000,
        "generator": "madgraph",
        "gridpack": f"GridPacks/gridpack_Zh_M125_ToAA_Tobbgg/Zh_M125_ToAA_{M}_Tobbgg_el9_amd64_gcc11_CMSSW_13_2_9_tarball.tar.xz",
        "card location": f"Gen_Cards/Cards/ZH/Inputcards_M{M}",
        "time per event": "",
        "size per event": "",
        "dataset name in previous campaign (if any)": f"/ZHToAA_AToBB_AToGG_M-{M}_TuneCP5_13TeV_madgraph_pythia8"

    }
    rows.append(row)

# Create DataFrame
df = pd.DataFrame(rows)

# Save to CSV
csv_path = "./ZH_gen_sample_summary.csv"
df.to_csv(csv_path, index=False)
