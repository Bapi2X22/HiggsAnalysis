# from coffea.processor import Runner, FuturesExecutor
# from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
# from coffea.processor import run_uproot_job
# from complex_processor import HiggsAnalysisProcessor

# fileset = {
#     "WHToAA": [
#         "root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/WHToAA_AToBB_AToGG_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2530000/FB6A4557-3B3E-5E4A-900D-45A77C107EA2.root"
#     ]
# }

# runner = Runner(
#     executor=FuturesExecutor(compression=None),
#     schema=NanoAODSchema,
#     savemetrics=True,
# )

# result, metrics = runner(
#     fileset=fileset,
#     treename="Events",
#     processor_instance=HiggsAnalysisProcessor(),
# )

# # Save output
# import pandas as pd
# df = pd.DataFrame(result)
# df.to_parquet("higgs_analysis_output.parquet")

# from coffea.processor import Runner, FuturesExecutor
# from coffea.nanoevents import NanoAODSchema
# from complex_processor import HiggsAnalysisProcessor  # your custom processor
# from coffea.processor.accumulator import column_accumulator

# import pandas as pd

# # Define fileset
# fileset = {
#     "WHToAA": [
#         "root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/WHToAA_AToBB_AToGG_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2530000/FB6A4557-3B3E-5E4A-900D-45A77C107EA2.root"
#     ]
# }

# # Set up runner
# runner = Runner(
#     executor=FuturesExecutor(compression=None),
#     schema=NanoAODSchema,
#     savemetrics=True,
# )

# # Run processor
# result, metrics = runner(
#     fileset=fileset,
#     treename="Events",
#     processor_instance=HiggsAnalysisProcessor(),
# )

# # Check what 'result' contains
# print(type(result))
# print(result)

# # If using column_accumulators:
# if isinstance(result["GenPhoton_pt_1"], column_accumulator):
#     result = {k: v.value for k, v in result.items()}

# # # Handle result (assuming result is an accumulator with a 'output' key or a dict)
# # if isinstance(result, dict):
# #     df = pd.DataFrame(result)
# # elif hasattr(result, "items"):
# #     # You may need to adapt this to your processor's output format
# #     df = pd.DataFrame({k: v.value if hasattr(v, "value") else v for k, v in result.items()})
# # else:
# #     raise TypeError("Unexpected result format from processor")

# # # Save to parquet
# # df.to_parquet("higgs_analysis_output.parquet")

# # Optional: Save to pandas for inspection
# df = pd.DataFrame(result)
# print(df.head())
# df.to_parquet("higgs_analysis_output.parquet")

from coffea.processor import Runner, FuturesExecutor
from coffea.nanoevents import NanoAODSchema
from coffea.processor.accumulator import column_accumulator
from complex_processor import HiggsAnalysisProcessor  # Your custom processor
import pandas as pd

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

# Convert accumulators to raw arrays
if isinstance(result["GenPhoton_pt_1"], column_accumulator):
    result = {k: v.value for k, v in result.items()}

# Optional: Save to pandas for inspection
df = pd.DataFrame(result)
print(df.head())

# Save to CSV
# df.to_csv("output.csv", index=False)
df.to_parquet("output_full.parquet")



