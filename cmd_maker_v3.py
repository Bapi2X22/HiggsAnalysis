
variables = [
    "h_higgs_pt",
    "h_A_pt_1", "h_A_pt_2",
    "h_leading_A_pt", "h_subleading_A_pt",
    "h_pho_from_a_pt_1", "h_pho_from_a_pt_2",
    "h_lead_pt_pho_gen", "h_sublead_pt_pho_gen",
    "h_bquark_from_a_pt_1", "h_bquark_from_a_pt_2",
    "h_lead_pt_bquark_gen", "h_sublead_pt_bquark_gen",
    "h_Reco_pho_pt",
    "h_Reco_lead_pho_pt", "h_Reco_sublead_pho_pt",
    "h_gen_photon_from_a_1_pt", "h_gen_photon_from_a_2_pt",
    "h_Gen_photon_pt_1", "h_Gen_photon_pt_2",
    "h_Genmatched_pho_1_pt", "h_Genmatched_pho_2_pt",
    "h_gen_lead_b_pt", "h_gen_sublead_b_pt",
    "h_gen_b1_pt", "h_gen_b2_pt",
    "h_reco_lead_bjet_pt", "h_reco_sublead_bjet_pt",
    "h_matched_bjet1_pt", "h_matched_bjet2_pt"
]


datasets = {
    "M55_18": "M55_RunIISummer20UL18NanoAODv9", "M60_18": "M60_RunIISummer20UL18NanoAODv9", 
    "M50_18": "M50_RunIISummer20UL18NanoAODv9", "M45_18": "M45_RunIISummer20UL18NanoAODv9", 
    "M40_18": "M40_RunIISummer20UL18NanoAODv9", "M35_18": "M35_RunIISummer20UL18NanoAODv9", 
    "M30_18": "M30_RunIISummer20UL18NanoAODv9", "M25_18": "M25_RunIISummer20UL18NanoAODv9", 
    "M20_18": "M20_RunIISummer20UL18NanoAODv9", "M60_17": "M60_RunIISummer20UL17NanoAODv9",
    "M55_17": "M55_RunIISummer20UL17NanoAODv9", "M50_17": "M50_RunIISummer20UL17NanoAODv9",
    "M45_17": "M45_RunIISummer20UL17NanoAODv9", "M40_17": "M40_RunIISummer20UL17NanoAODv9",
    "M35_17": "M35_RunIISummer20UL17NanoAODv9", "M30_17": "M30_RunIISummer20UL17NanoAODv9",
    "M25_17": "M25_RunIISummer20UL17NanoAODv9", "M20_17": "M20_RunIISummer20UL17NanoAODv9",
    "M60_16": "M60_RunIISummer20UL16NanoAODv9", "M55_16": "M55_RunIISummer20UL16NanoAODv9",
    "M45_16": "M45_RunIISummer20UL16NanoAODv9", "M50_16": "M50_RunIISummer20UL16NanoAODv9",
    "M40_16": "M40_RunIISummer20UL16NanoAODv9", "M35_16": "M35_RunIISummer20UL16NanoAODv9",
    "M30_16": "M30_RunIISummer20UL16NanoAODv9", "M25_16": "M25_RunIISummer20UL16NanoAODv9",
    "M20_16": "M20_RunIISummer20UL16NanoAODv9", "M60_16APV": "M60_RunIISummer20UL16NanoAODAPVv9",
    "M55_16APV": "M55_RunIISummer20UL16NanoAODAPVv9", "M50_16APV": "M50_RunIISummer20UL16NanoAODAPVv9",
    "M45_16APV": "M45_RunIISummer20UL16NanoAODAPVv9", "M40_16APV": "M40_RunIISummer20UL16NanoAODAPVv9",
    "M35_16APV": "M35_RunIISummer20UL16NanoAODAPVv9", "M20_16APV": "M20_RunIISummer20UL16NanoAODAPVv9",
    "M30_16APV": "M30_RunIISummer20UL16NanoAODAPVv9", "M25_16APV": "M25_RunIISummer20UL16NanoAODAPVv9"
}

datasets_by_year = {
    "2018": ["M20_18","M25_18","M30_18","M35_18","M40_18","M45_18","M50_18","M55_18","M60_18"],
    "2017": ["M20_17","M25_17","M30_17","M35_17","M40_17","M45_17","M50_17","M55_17","M60_17"],
    "2016": ["M20_16","M25_16","M30_16","M35_16","M40_16","M45_16","M50_16","M55_16","M60_16"],
    "2016APV": ["M20_16APV","M25_16APV","M30_16APV","M35_16APV","M40_16APV","M45_16APV","M50_16APV","M55_16APV","M60_16APV"]
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
#     # Always make per-year plots
#     for year, keys in datasets_by_year.items():
#         title = f"{var} distribution {year}"
#         args = [
#             "Overlay_hist.py", "hist_output_WH.root", "dirs",
#             *keys, var,
#             "--rebin", "3", "--ratio", "--xrange", "0", "200",
#             "--title", title, "--save_name", title
#         ]
#         print(args, ",")

#     # Make pair plots only if var has a valid partner
#     pair_var = get_pair(var)
#     if pair_var and pair_var in variables and (tuple(sorted([var, pair_var])) not in done_pairs):
#         for key in datasets:
#             title = f"{var} vs {pair_var} {key.replace('_', ' ')}"
#             args = [
#                 "Overlay_hist.py", "hist_output_WH.root", "branches",
#                 key, var, pair_var,
#                 "--rebin", "3", "--ratio", "--xrange", "0", "200",
#                 "--title", title, "--save_name", title
#             ]
#             print(args, ",")
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

# --- Configurable option ---
skip_mass_points = True   # set False to include all, True to skip every other

for var in variables:
    # Always make per-year plots
    for year, keys in datasets_by_year.items():
        if skip_mass_points:
            keys_to_use = keys[::2]  # take every second point (M20, M30, M40, ...)
        else:
            keys_to_use = keys

        title = f"{var} distribution {year}"
        args = [
            "Overlay_hist.py", "hist_output_ZH.root", "dirs",
            *keys_to_use, var,
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
                "Overlay_hist.py", "hist_output_ZH.root", "branches",
                key, var, pair_var,
                "--rebin", "3", "--ratio", "--xrange", "0", "200",
                "--title", title, "--save_name", title
            ]
            print(args, ",")
        done_pairs.add(tuple(sorted([var, pair_var])))
