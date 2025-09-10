import os

variables = [
    "h_leading_A_pt", "h_subleading_A_pt",
    "h_lead_pt_pho_gen", "h_sublead_pt_pho_gen",
    "h_Reco_pho_pt",
    "h_Reco_lead_pho_pt", "h_Reco_sublead_pho_pt",
    "h_Genmatched_pho_1_pt", "h_Genmatched_pho_2_pt",
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

# lumi per year
lumi_map = {
    "16": 16900,      # 2016 postVFP UL
    "16APV": 19500,   # 2016 preVFP UL
    "17": 41500,      # 2017 UL
    "18": 59700,      # 2018 UL
    "22": 34700,      # Run3 2022
}

commands = []

for var in variables:
    for tag, dirname in datasets.items():
        # detect year
        if "16APV" in tag:
            year = "16APV"
        elif "_16" in tag:
            year = "16"
        elif "_17" in tag:
            year = "17"
        elif "_18" in tag:
            year = "18"
        else:
            year = "22"

        lumi = lumi_map[year]

        for do_scale in [False, True]:
            # dirname is like "M55_RunIISummer20UL18NanoAODv9"
            mass = dirname.split("_")[0]   # --> "M55"

            # Extract campaign, e.g. "20UL18"
            parts = dirname.split("_")
            campaign = parts[1].replace("RunIISummer", "").replace("NanoAODv9", "")

            cmd = f'''compare_hists(
                histname="{var}",
                file1="hist_output_WH.root",
                dir1="{dirname}",
                xsec1=1.373,
                lumi1={lumi},

                file2="hist_output_ggH.root",
                dir2="{mass}_Run3Summer22NanoAODv13",
                xsec2=48.58,
                lumi2={lumi_map["22"]},

                file3="hist_output_ZH.root",
                dir3="{dirname}",
                xsec3=0.8839,
                lumi3={lumi},

                leg1="WH ({mass}, {campaign})",
                leg2="ggH ({mass}, 2022 Run3)",
                leg3="ZH ({mass}, {campaign})",

                xlabel="{var}",
                rebin=2,
                outprefix="{var}_WH_vs_ggH_vs_ZH_{tag}_{'scaled' if do_scale else 'raw'}",
                outroot="Histo_comparison.root",
                do_scale={str(do_scale)}
            )'''

            commands.append(cmd)



# save all generated compare_hists calls
with open("generate_commands.py", "w") as f:
    f.write("\n\n".join(commands))

print(f"Generated {len(commands)} commands into generate_commands.py")


