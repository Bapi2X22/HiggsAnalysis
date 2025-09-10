import awkward as ak
from coffea.processor import Runner, FuturesExecutor
from coffea.nanoevents import NanoAODSchema
from WH_processor import HiggsAnalysisProcessor

fileset = {
    "M20": [
        "root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/WHToAA_AToBB_AToGG_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2530000/FB6A4557-3B3E-5E4A-900D-45A77C107EA2.root",
        "root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/WHToAA_AToBB_AToGG_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2530000/52967538-6671-C748-9CEC-C21D276D640B.root",
    ],
    "M25": [
        "root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/WHToAA_AToBB_AToGG_M-25_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2520000/4531AE5A-5C3B-F446-A0E3-B9DAA19B87C5.root",
        "root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/WHToAA_AToBB_AToGG_M-25_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2520000/7E11D93E-1C81-2849-B3E5-BC7B90404C3D.root",
    ]
}

runner = Runner(
    executor=FuturesExecutor(compression=None),
    schema=NanoAODSchema,
    savemetrics=True,
)

all_results = {}

for ds_name, files in fileset.items():
    print(f"Processing dataset: {ds_name}")
    result, metrics = runner(
        fileset={ds_name: files},
        treename="Events",
        processor_instance=HiggsAnalysisProcessor(),
    )
    # Convert column_accumulator to awkward arrays
    processed = {key: ak.Array(value.value) for key, value in result.items()}
    all_results[ds_name] = processed

# Convert dictionary of dicts to one nested awkward array
nested_ak = ak.Array(all_results)

# Save everything to a single Parquet file
ak.to_parquet(nested_ak, "all_datasets_output.parquet")
print("Saved all datasets in one Parquet file: all_datasets_output.parquet")
