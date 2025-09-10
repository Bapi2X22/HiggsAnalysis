
variables = [
    "h_higgs_pt",
    "h_A_pt_1", "h_A_pt_2",
    "h_leading_A_pt", "h_subleading_A_pt",
    "h_pho_from_a_pt_1", "h_pho_from_a_pt_2",
    "h_lead_pt_pho_gen", "h_sublead_pt_pho_gen",
    "h_Reco_pho_pt",
    "h_Reco_lead_pho_pt", "h_Reco_sublead_pho_pt",
    "h_gen_photon_from_a_1_pt", "h_gen_photon_from_a_2_pt",
    "h_Gen_photon_pt_1", "h_Gen_photon_pt_2",
    "h_Genmatched_pho_1_pt", "h_Genmatched_pho_2_pt",
]


datasets = {
    "M15_22EE": "M15_Run3Summer22EENanoAODv13",
    "M15_22":   "M15_Run3Summer22NanoAODv13",

    "M20_22EE": "M20_Run3Summer22EENanoAODv13",
    "M20_22":   "M20_Run3Summer22NanoAODv13",

    "M25_22EE": "M25_Run3Summer22EENanoAODv13",
    "M25_22":   "M25_Run3Summer22NanoAODv13",

    "M30_22EE": "M30_Run3Summer22EENanoAODv13",
    "M30_22":   "M30_Run3Summer22NanoAODv13",

    "M35_22EE": "M35_Run3Summer22EENanoAODv13",
    "M35_22":   "M35_Run3Summer22NanoAODv13",

    "M40_22EE": "M40_Run3Summer22EENanoAODv13",
    "M40_22":   "M40_Run3Summer22NanoAODv13",

    "M45_22EE": "M45_Run3Summer22EENanoAODv13",
    "M45_22":   "M45_Run3Summer22NanoAODv13",

    "M50_22EE": "M50_Run3Summer22EENanoAODv13",
    "M50_22":   "M50_Run3Summer22NanoAODv13",

    "M55_22EE": "M55_Run3Summer22EENanoAODv13",
    "M55_22":   "M55_Run3Summer22NanoAODv13",

    "M60_22EE": "M60_Run3Summer22EENanoAODv13",
    "M60_22":   "M60_Run3Summer22NanoAODv13",
}

datasets_by_year = {
    "2022": ["M15_22","M15_22","M20_22","M25_22","M30_22","M35_22","M40_22","M45_22","M50_22","M55_22","M60_22"],
    "2022EE": ["M15_22EE","M15_22EE","M20_22EE","M25_22EE","M30_22EE","M35_22EE","M40_22EE","M45_22EE","M50_22EE","M55_22EE","M60_22EE"],
}


# def get_pair(var):
#     if "_1" in var:
#         return var.replace("_1", "_2")
#     elif "_2" in var:
#         return var.replace("_2", "_1")
#     elif "lead" in var and "sublead" not in var:
#         return var.replace("lead", "sublead")
#     elif "sublead" in var:
#         return var.replace("sublead", "lead")
#     return None

# done_pairs = set()

# for var in variables:
#     # Only dirs mode for lead/sublead/_1/_2
#     if any(x in var for x in ["lead", "_1", "_2"]):
#         for year, keys in datasets_by_year.items():
#             title = f"{var} distribution {year}"
#             args = [
#                 "Overlay_hist.py", "hist_output_ggH.root", "dirs",
#                 *keys, var,
#                 "--rebin", "3", "--ratio", "--xrange", "0", "200",
#                 "--title", title, "--save_name", title
#             ]
#             print(args,",")

#     # Pair plots (only if not already done)
#     pair_var = get_pair(var)
#     if pair_var and pair_var in variables and (tuple(sorted([var, pair_var])) not in done_pairs):
#         for key in datasets:
#             title = f"{var} vs {pair_var} {key.replace('_', ' ')}"
#             args = [
#                 "Overlay_hist.py", "hist_output_ggH.root", "branches",
#                 key, var, pair_var,
#                 "--rebin", "3", "--ratio", "--xrange", "0", "200",
#                 "--title", title, "--save_name", title
#             ]
#             print(args,",")
#         done_pairs.add(tuple(sorted([var, pair_var])))

def get_pair(var):
    if "_1" in var:
        return var.replace("_1", "_2")
    elif "_2" in var:
        return var.replace("_2", "_1")
    elif "lead" in var and "sublead" not in var:
        return var.replace("lead", "sublead")
    elif "sublead" in var:
        return var.replace("sublead", "lead")
    return None


done_pairs = set()

for var in variables:
    # Always make per-year plots
    for year, keys in datasets_by_year.items():
        title = f"{var} distribution {year}"
        args = [
            "Overlay_hist.py", "hist_output_ggH.root", "dirs",
            *keys, var,
            "--rebin", "3", "--ratio", "--xrange", "0", "200",
            "--title", title, "--save_name", title
        ]
        print(args, ",")

    # Make pair plots only if var has a valid partner
    pair_var = get_pair(var)
    if pair_var and pair_var in variables and (tuple(sorted([var, pair_var])) not in done_pairs):
        for key in datasets:
            title = f"{var} vs {pair_var} {key.replace('_', ' ')}"
            args = [
                "Overlay_hist.py", "hist_output_ggH.root", "branches",
                key, var, pair_var,
                "--rebin", "3", "--ratio", "--xrange", "0", "200",
                "--title", title, "--save_name", title
            ]
            print(args, ",")
        done_pairs.add(tuple(sorted([var, pair_var])))

