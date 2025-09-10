datasets = {
    "M55_18": "M55_RunIISummer20UL18NanoAODv9",
    "M60_18": "M60_RunIISummer20UL18NanoAODv9",
    "M50_18": "M50_RunIISummer20UL18NanoAODv9",
    "M45_18": "M45_RunIISummer20UL18NanoAODv9",
    "M40_18": "M40_RunIISummer20UL18NanoAODv9",
    "M35_18": "M35_RunIISummer20UL18NanoAODv9",
    "M30_18": "M30_RunIISummer20UL18NanoAODv9",
    "M25_18": "M25_RunIISummer20UL18NanoAODv9",
    "M20_18": "M20_RunIISummer20UL18NanoAODv9",
    "M60_17": "M60_RunIISummer20UL17NanoAODv9",
    "M55_17": "M55_RunIISummer20UL17NanoAODv9",
    "M50_17": "M50_RunIISummer20UL17NanoAODv9",
    "M45_17": "M45_RunIISummer20UL17NanoAODv9",
    "M40_17": "M40_RunIISummer20UL17NanoAODv9",
    "M35_17": "M35_RunIISummer20UL17NanoAODv9",
    "M30_17": "M30_RunIISummer20UL17NanoAODv9",
    "M25_17": "M25_RunIISummer20UL17NanoAODv9",
    "M20_17": "M20_RunIISummer20UL17NanoAODv9",
    "M60_16": "M60_RunIISummer20UL16NanoAODv9",
    "M55_16": "M55_RunIISummer20UL16NanoAODv9",
    "M45_16": "M45_RunIISummer20UL16NanoAODv9",
    "M50_16": "M50_RunIISummer20UL16NanoAODv9",
    "M40_16": "M40_RunIISummer20UL16NanoAODv9",
    "M35_16": "M35_RunIISummer20UL16NanoAODv9",
    "M30_16": "M30_RunIISummer20UL16NanoAODv9",
    "M25_16": "M25_RunIISummer20UL16NanoAODv9",
    "M20_16": "M20_RunIISummer20UL16NanoAODv9",
    "M60_16APV": "M60_RunIISummer20UL16NanoAODAPVv9",
    "M55_16APV": "M55_RunIISummer20UL16NanoAODAPVv9",
    "M50_16APV": "M50_RunIISummer20UL16NanoAODAPVv9",
    "M45_16APV": "M45_RunIISummer20UL16NanoAODAPVv9",
    "M40_16APV": "M40_RunIISummer20UL16NanoAODAPVv9",
    "M35_16APV": "M35_RunIISummer20UL16NanoAODAPVv9",
    "M20_16APV": "M20_RunIISummer20UL16NanoAODAPVv9",
    "M30_16APV": "M30_RunIISummer20UL16NanoAODAPVv9",
    "M25_16APV": "M25_RunIISummer20UL16NanoAODAPVv9"
}

datasets_by_year = {
    "2018": ["M20_18","M25_18","M30_18","M35_18","M40_18","M45_18","M50_18","M55_18","M60_18"],
    "2017": ["M20_17","M25_17","M30_17","M35_17","M40_17","M45_17","M50_17","M55_17","M60_17"],
    "2016": ["M20_16","M25_16","M30_16","M35_16","M40_16","M45_16","M50_16","M55_16","M60_16"],
    "2016APV": ["M20_16APV","M25_16APV","M30_16APV","M35_16APV","M40_16APV","M45_16APV","M50_16APV","M55_16APV","M60_16APV"]
}

base_cmd = [
    "Overlay_hist.py",
    "hist_output_WH.root",
    "branches",
    None,  # placeholder for dataset key
    "h_lead_pt_pho_gen",
    "h_sublead_pt_pho_gen",
    "--ratio",
    "--rebin", "3",
    "--xrange", "0", "200",
    "--title", None,  # placeholder for title
    "--save_name", None  # placeholder for save name
]

for key in datasets:
    # Extract readable title from key (e.g., "M20_16APV" -> "M20 2016APV")
    title_part = key.replace("_", " ")
    cmd = base_cmd.copy()
    cmd[3] = key
    cmd[13] = f"Gen lead vs sublead Photon pT distribution {title_part}"
    cmd[15] = f"Gen lead vs sublead Photon pT distribution {title_part}"
    print(cmd, ",")


base_cmd = [
    "Overlay_hist.py",
    "hist_output_WH.root",
    "dirs",
    None,  # placeholder for dataset keys
    "h_leading_A_pt",
    "--rebin", "3",
    "--xrange", "0", "200",
    "--title", None,  # placeholder for title
    "--save_name", None  # placeholder for save name
]

for year, keys in datasets_by_year.items():
    cmd = base_cmd.copy()
    cmd[3] = None  # not used here
    full_cmd = ["Overlay_hist.py", "hist_output_WH.root", "dirs"] + keys + [
        "h_leading_A_pt",
        "--rebin", "3",
        "--ratio",
        "--xrange", "0", "200",
        "--title", f"Leading A pT distribution {year}",
        "--save_name", f"Leading A pT distribution {year}"
    ]
    print(full_cmd,",")


