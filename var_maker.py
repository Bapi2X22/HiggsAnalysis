import awkward as ak

# Load events from the Parquet file
Events_ak = ak.from_parquet("Mass_20_2018_WH.parquet")

# Iterate over all top-level keys and print the variable assignment line
for key in Events_ak.fields:
    print(f'{key} = ak.to_numpy(Events_ak["{key}"])')

