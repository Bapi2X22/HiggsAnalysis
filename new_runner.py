# run_processor.py

from coffea.processor import Runner, FuturesExecutor
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from vanilla_processor import PrintFieldsProcessor

fileset = {
    "WHToAA": [
        "root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/WHToAA_AToBB_AToGG_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2530000/FB6A4557-3B3E-5E4A-900D-45A77C107EA2.root",
        "root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/WHToAA_AToBB_AToGG_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2530000/52967538-6671-C748-9CEC-C21D276D640B.root",
    ]
}

runner = Runner(
    executor=FuturesExecutor(compression=None),
    schema=NanoAODSchema,
    savemetrics=True,
)

out, metrics = runner(
    fileset,
    treename="Events",
    processor_instance=PrintFieldsProcessor(),
)
