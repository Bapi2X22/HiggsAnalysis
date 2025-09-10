# from coffea import processor
# from coffea.nanoevents import NanoAODSchema
# from dask.distributed import Client, LocalCluster
# import pickle
# from wh_to_aa_processor import WHToAAProcessor

# fileset = {
#     "WHToAA": [
#         "root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/WHToAA_AToBB_AToGG_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2530000/FB6A4557-3B3E-5E4A-900D-45A77C107EA2.root"
#         "root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/WHToAA_AToBB_AToGG_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2530000/52967538-6671-C748-9CEC-C21D276D640B.root"
#         # Add more files here
#     ]
# }

# # Setup Dask
# cluster = LocalCluster()
# client = Client(cluster)

# # Run the processor
# output = processor.run_uproot_job(
#     fileset,
#     "Events",
#     processor_instance=WHToAAProcessor(),
#     executor=processor.dask_executor,
#     executor_args={"client": client, "schema": NanoAODSchema},
# )

# # Save merged result
# with open("merged_output.pkl", "wb") as fout:
#     pickle.dump(output, fout)

from coffea import processor
from coffea.nanoevents import NanoAODSchema
from dask.distributed import Client, LocalCluster
import pickle
from wh_to_aa_processor import WHToAAProcessor
import multiprocessing

def main():
    fileset = {
        "WHToAA": [
            "root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/WHToAA_AToBB_AToGG_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2530000/FB6A4557-3B3E-5E4A-900D-45A77C107EA2.root",
            "root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/WHToAA_AToBB_AToGG_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2530000/52967538-6671-C748-9CEC-C21D276D640B.root"
            # Add more files here
        ]
    }

    # Setup Dask
    cluster = LocalCluster()
    client = Client(cluster)

    # # Run the processor
    # output = processor.run_uproot_job(
    #     fileset,
    #     "Events",
    #     processor_instance=WHToAAProcessor(),
    #     executor=processor.dask_executor,
    #     executor_args={"client": client, "schema": NanoAODSchema},
    # )
    print("Debugging...")

    print("Number of files:", len(fileset["WHToAA"]))

    output = processor.run_uproot_job(
    fileset,
    "Events",
    processor_instance=WHToAAProcessor(),
    executor=processor.futures_executor,
    executor_args={"workers": 4, "schema": NanoAODSchema},
)

    # Save merged result
    with open("merged_output.pkl", "wb") as fout:
        pickle.dump(output, fout)

if __name__ == "__main__":
    multiprocessing.set_start_method("fork", force=True)  # Fixes multiprocessing error
    main()

