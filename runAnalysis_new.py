# from coffea.nanoevents import NanoEventsFactory, BaseSchema
# from coffea.processor import Runner, FuturesExecutor
# # from MyProcessor import MyAnalysisProcessor
# from MyProcessor import TestProcessor
# from coffea.nanoevents import NanoAODSchema

# # fileset = {"SampleName": ["file1.root", "file2.root", "..."]}
# fileset = {
#     "WHToAA": [
#         "root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/WHToAA_AToBB_AToGG_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2530000/FB6A4557-3B3E-5E4A-900D-45A77C107EA2.root",
#         "root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/WHToAA_AToBB_AToGG_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2530000/52967538-6671-C748-9CEC-C21D276D640B.root"
#         # Add more files here
#     ]
# }

# runner = Runner(
#     executor=FuturesExecutor(workers=4),  # Adjust for cores or use DaskExecutor for a cluster
#     schema=NanoAODSchema,
#     chunksize=100_000,
#     savemetrics=True,
# )
# output, metrics = runner(
#     fileset,
#     treename="Events",  # or whatever your tree is
#     processor_instance=TestProcessor(),
# )


'''
from coffea import processor
import awkward as ak
import numpy as np
from coffea.nanoevents import NanoAODSchema
# from coffea.processor import run_uproot_job
from coffea.processor.executor import run_uproot_job



class TestProcessor(processor.ProcessorBase):
    def __init__(self):
        self._accumulator = processor.dict_accumulator({})

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        print("Fields in events:", events.fields)
        print("Events type:", type(events))
        print("GenPart count:", ak.num(events.GenPart))
        return self.accumulator.identity()

    def postprocess(self, accumulator):
        return accumulator


fileset = {
    "WHToAA": [
        "root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/WHToAA_AToBB_AToGG_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2530000/FB6A4557-3B3E-5E4A-900D-45A77C107EA2.root",
        "root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/WHToAA_AToBB_AToGG_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2530000/52967538-6671-C748-9CEC-C21D276D640B.root"
    ]
}

result = run_uproot_job(
    fileset=fileset,
    treename="Events",
    processor_instance=TestProcessor(),
    executor="futures",
    executor_args={"workers": 4},
    schema=NanoAODSchema,
    chunksize=100000,
)

print("Result:", result)
'''

from coffea import processor
from coffea.processor import Runner, FuturesExecutor
from coffea.nanoevents import NanoAODSchema
from MyProcessor import TestProcessor

fileset = {
    "WHToAA": [
        "root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/WHToAA_AToBB_AToGG_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2530000/FB6A4557-3B3E-5E4A-900D-45A77C107EA2.root",
        "root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/WHToAA_AToBB_AToGG_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2530000/52967538-6671-C748-9CEC-C21D276D640B.root",
    ]
}

runner = Runner(
    executor=FuturesExecutor(workers=2),
    schema=NanoAODSchema,
    chunksize=100_000,
    savemetrics=True,
)

output, metrics = runner(
    fileset,
    treename="Events",
    processor_instance=TestProcessor(),
)

print("Output:", output)
print("Metrics:", metrics)
