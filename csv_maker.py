import pandas as pd

# Mass points
mass_points = [12, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]

# Data construction
rows = []
for M in mass_points:
    row = {
        "dataset": f"WHTo2A_ATo2B_ATo2G_M-{M}_TuneCP5_13p6TeV_madgraph-pythia8",
        "fragment": f"https://gitlab.cern.ch/bbapi/Bapi2X22-Run3-HtoAAtobbgg-Samples/-/tree/main/Fragments/WH/Wh_M125_ToAA_{M}_Tobbgg_TuneCP5_PSweights_13p6TeV-madgraph_pythia8-fragment.py",
        "events": 1000000,
        "generator": "madgraph",
        "gridpack": f"https://gitlab.cern.ch/bbapi/Bapi2X22-Run3-HtoAAtobbgg-Samples/-/tree/main/GridPacks/gridpack_Wh_M125_ToAA_Tobbgg/Wh_M125_ToAA_{M}_Tobbgg_el9_amd64_gcc11_CMSSW_13_2_9_tarball.tar.xz",
        "card location": f"https://gitlab.cern.ch/bbapi/Bapi2X22-Run3-HtoAAtobbgg-Samples/-/tree/main/Gen_Cards/Cards/WH/Inputcards_M{M}",
        "time per event": "",
        "size per event": "",
        "dataset name in previous campaign (if any)": f"/WHToAA_AToBB_AToGG_M-{M}_TuneCP5_13TeV_madgraph_pythia8"

    }
    rows.append(row)

# Create DataFrame
df = pd.DataFrame(rows)

# Save to CSV
csv_path = "./WH_gen_sample_summary.csv"
df.to_csv(csv_path, index=False)


# import pandas as pd

# # Mass points
# mass_points = [12, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]

# # Data construction
# rows = []
# for M in mass_points:
#     row = {
#         "dataset": f"HAHMHTo2A_ATo2B_ATo2G_M-{M}_TuneCP5_13p6TeV_madgraphMLM-pythia8",
#         "fragment": f"https://gitlab.cern.ch/bbapi/Bapi2X22-Run3-HtoAAtobbgg-Samples/-/tree/main/hToaaTo2gamma2b_ma{M}GeV_TuneCP5_PSweights_13p6TeV-madgraphMLM-pythia8-fragment.py",
#         "events": 1000000,
#         "generator": "madgraph",
#         "gridpack": f"https://gitlab.cern.ch/bbapi/Bapi2X22-Run3-HtoAAtobbgg-Samples/-/tree/main/GridPacks/gridpack_hToaaTo2gamma2b/hToaaTo2gamma2b_ma{M}GeV_MLM_4f_max1j_el9_amd64_gcc11_CMSSW_13_2_9_tarball.tar.xz",
#         "card location": f"https://gitlab.cern.ch/bbapi/Bapi2X22-Run3-HtoAAtobbgg-Samples/-/tree/main/Gen_Cards/Cards/ggH/Inputcards_M{M}",
#         "time per event": "",
#         "size per event": "",
#         "dataset name in previous campaign (if any)": f"--"

#     }
#     rows.append(row)

# # Create DataFrame
# df = pd.DataFrame(rows)

# # Save to CSV
# csv_path = "./ggH_gen_sample_summary.csv"
# df.to_csv(csv_path, index=False)


