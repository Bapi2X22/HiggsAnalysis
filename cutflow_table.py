# import pandas as pd
# import re

# def parse_dataset_name(name):
#     # Extract mass
#     mass_match = re.match(r"^(M\d+)", name)
#     mass = mass_match.group(1) if mass_match else "Unknown"

#     # Extract year
#     year_match = re.search(r"UL(\d+)(NanoAODAPVv9|NanoAODv9)", name)
#     if year_match:
#         year = "20" + year_match.group(1)
#         if "APV" in year_match.group(2):
#             year += "APV"
#     else:
#         year = "Unknown"
#     return mass, year


# # ---- Get the stage column order from the *first dict* ----
# first_key = next(iter(count_stages))
# stage_order = list(count_stages[first_key].keys())

# # ---- Flatten into rows ----
# rows = []
# for dset, stages in count_stages.items():
#     mass, year = parse_dataset_name(dset)
#     row = {"Year": year, "Mass": mass}
#     # preserve order of stage keys
#     for col in stage_order:
#         row[col] = stages.get(col, None)
#     # add efficiency (last / first)
#     first_val = row[stage_order[0]]
#     last_val = row[stage_order[-1]]
#     row["Efficiency"] = last_val / first_val if first_val else None
#     rows.append(row)

# # ---- Build DataFrame with correct order ----
# df = pd.DataFrame(rows)
# df = df[["Year", "Mass"] + stage_order + ["Efficiency"]]

# # ---- Sorting ----
# df["Mass_num"] = df["Mass"].str.extract(r"M(\d+)").astype(int)
# df["Year_sort"] = df["Year"].replace({"2016APV": 2015.5}).replace("Unknown", 9999)
# df = df.sort_values(by=["Year_sort", "Mass_num"]).reset_index(drop=True)
# df = df.drop(columns=["Mass_num", "Year_sort"])

# # ---- Create grouped table for Excel ----
# excel_rows = []
# for year, block in df.groupby("Year"):
#     excel_rows.append([f"Year: {year}"] + [""] * (len(stage_order) + 1))
#     excel_rows.extend(block.drop(columns=["Year"]).values.tolist())
#     excel_rows.append([""] * (len(stage_order) + 2))

# df_excel = pd.DataFrame(
#     excel_rows,
#     columns=["Year/Mass"] + stage_order + ["Efficiency"]
# )

# df_excel.to_excel("photon_cutflow_table.xlsx", index=False)

# print("Saved photon_cutflow_table.xlsx with preserved column order and efficiency column.")

import pandas as pd
import re

def parse_dataset_name(name):
    """Extract mass and year info from dataset name."""
    mass_match = re.match(r"^(M\d+)", name)
    mass = mass_match.group(1) if mass_match else "Unknown"

    year_match = re.search(r"UL(\d+)(NanoAODAPVv9|NanoAODv9)", name)
    if year_match:
        year = "20" + year_match.group(1)
        if "APV" in year_match.group(2):
            year += "APV"
    else:
        year = "Unknown"
    return mass, year


def make_cutflow_table(count_stages: dict, output_file: str):
    """
    Build a photon cutflow table from count_stages and save to Excel.
    
    Parameters
    ----------
    count_stages : dict
        Nested dict of {dataset_name: {stage_name: count}}
    output_file : str
        Path to the Excel file to save
    """

    # ---- Get the stage column order from the *first dict* ----
    first_key = next(iter(count_stages))
    stage_order = list(count_stages[first_key].keys())

    # ---- Flatten into rows ----
    rows = []
    for dset, stages in count_stages.items():
        mass, year = parse_dataset_name(dset)
        row = {"Year": year, "Mass": mass}
        # preserve order of stage keys
        for col in stage_order:
            row[col] = stages.get(col, None)
        # add efficiency (last / first)
        first_val = row[stage_order[0]]
        last_val = row[stage_order[-1]]
        # row["Efficiency"] = last_val / first_val if first_val else None
        row["Efficiency"] = round(last_val / first_val, 3) if first_val else None
        rows.append(row)

    # ---- Build DataFrame with correct order ----
    df = pd.DataFrame(rows)
    df = df[["Year", "Mass"] + stage_order + ["Efficiency"]]

    # ---- Sorting ----
    df["Mass_num"] = df["Mass"].str.extract(r"M(\d+)").astype(int)
    df["Year_sort"] = df["Year"].replace({"2016APV": 2015.5}).replace("Unknown", 9999)
    df = df.sort_values(by=["Year_sort", "Mass_num"]).reset_index(drop=True)
    df = df.drop(columns=["Mass_num", "Year_sort"])

    # ---- Create grouped table for Excel ----
    excel_rows = []
    for year, block in df.groupby("Year"):
        excel_rows.append([f"Year: {year}"] + [""] * (len(stage_order) + 1))
        excel_rows.extend(block.drop(columns=["Year"]).values.tolist())
        excel_rows.append([""] * (len(stage_order) + 2))

    df_excel = pd.DataFrame(
        excel_rows,
        columns=["Year/Mass"] + stage_order + ["Efficiency"]
    )

    # ---- Save ----
    # df_excel.to_excel(output_file, index=False)
    with pd.ExcelWriter(output_file, engine="xlsxwriter") as writer:
        df_excel.to_excel(writer, index=False, sheet_name="Cutflow")

        # Apply 3 decimal format to Efficiency column
        workbook  = writer.book
        worksheet = writer.sheets["Cutflow"]
        format3dec = workbook.add_format({"num_format": "0.000"})

        # Efficiency column index (last col)
        eff_col_idx = len(df_excel.columns) - 1
        worksheet.set_column(eff_col_idx, eff_col_idx, None, format3dec)
    print(f"Saved {output_file} with preserved column order and efficiency column.")
