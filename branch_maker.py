import pandas as pd

# Read the Parquet file
Events = pd.read_parquet("Mass_20_2018_WH.parquet")

# Get list of all branch names
branches = Events.columns

# Print branch_arrays structure
print("branch_arrays = {")
for b in branches:
    print(f'    "{b}": {b},')
print("}\n")

# Print bin_settings structure with empty tuples
print("bin_settings = {")
for b in branches:
    print(f'    "{b}": (),')
print("}")
