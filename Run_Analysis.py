from coffea.processor import Runner, FuturesExecutor
from coffea.nanoevents import NanoAODSchema
from coffea.processor.accumulator import column_accumulator
from WH_processor import HiggsAnalysisProcessor  # Your custom processor
import pandas as pd
import awkward as ak

# # Define the fileset
# fileset = {
#     "WHToAA": [
#         "root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/WHToAA_AToBB_AToGG_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2530000/FB6A4557-3B3E-5E4A-900D-45A77C107EA2.root"
#     ]
# }

fileset = {
    "WHToAA": [
        "root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/WHToAA_AToBB_AToGG_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2530000/FB6A4557-3B3E-5E4A-900D-45A77C107EA2.root",
        "root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/WHToAA_AToBB_AToGG_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2530000/52967538-6671-C748-9CEC-C21D276D640B.root",
    ]
}

# NanoAODSchema.config["include"] = None
# Set up the runner
runner = Runner(
    executor=FuturesExecutor(compression=None),
    schema=NanoAODSchema,
    savemetrics=True,
)

# Run the processor
result, metrics = runner(
    fileset=fileset,
    treename="Events",
    processor_instance=HiggsAnalysisProcessor(),
)

# Print result
print(type(result))
print(result)

# # Convert accumulators to raw arrays
# if isinstance(result["pho_from_a_pt_1"], column_accumulator):
#     result = {k: v.value for k, v in result.items()}

# # Optional: Save to pandas for inspection
# df = pd.DataFrame(result)
# print(df.head())

# Convert all column_accumulator values to NumPy or Awkward arrays
ak_arrays = {
    key: ak.Array(value.value) if isinstance(value, column_accumulator) else ak.Array(value)
    for key, value in result.items()
}

# Save to CSV
# df.to_csv("output.csv", index=False)
# df.to_parquet("full_result.parquet")
# ak.to_parquet(ak_arrays, "Full_output.parquet")
ak.to_parquet(ak_arrays, "Mass_20_2018_WH.parquet")



