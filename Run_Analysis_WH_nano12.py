from coffea.processor import Runner, FuturesExecutor
from coffea.nanoevents import NanoAODSchema
# from Custom_Processor import HiggsAnalysisProcessor
from All_objects_WH_processor import HiggsAnalysisProcessor
import awkward as ak
import pyarrow.parquet as pq
import pyarrow as pa
import os
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from coffea import processor
import glob


with open("samples_local.json") as f:
    fileset = json.load(f)

parquet_dir = "parquet_files_WH_all_nano12"

def process_and_save(dataset_name, files):
    try:
        runner = Runner(
            executor=FuturesExecutor(compression=None, workers=20),
            schema=NanoAODSchema,
            # schemaclass=NanoAODSchema,
            # schemaargs={"version": "v12"},
            savemetrics=True,
        )

        output, _ = runner(
            fileset={dataset_name: files},
            treename="Events",
            processor_instance=HiggsAnalysisProcessor(),
        )

        events = {
            key: ak.Array(val.value) if hasattr(val, "value") else val
            for key, val in output[dataset_name].items()
        }

        if not events or any(len(v) == 0 for v in events.values()):
            print(f"Skipping {dataset_name}: No events")
            return None

        num_events = len(next(iter(events.values())))
        events["dataset"] = ak.Array([dataset_name] * num_events)

        os.makedirs(parquet_dir, exist_ok=True)
        out_file = os.path.join(parquet_dir, f"{dataset_name}.parquet")
        table = ak.to_arrow_table(events)
        pq.write_table(table, out_file, compression=None)

        return out_file

    except Exception as e:
        print(f"Error processing {dataset_name}: {e}")
        return None


def is_parquet_valid(file_path, delete_if_invalid=False):
    """Check if parquet file exists, is non-empty, and can be opened."""
    try:
        if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
            pq.read_table(file_path)
            return True
    except Exception:
        pass

    if delete_if_invalid and os.path.exists(file_path):
        print(f"Removing broken file: {file_path}")
        try:
            os.remove(file_path)
        except Exception as e:
            print(f"Could not delete {file_path}: {e}")

    return False

fileset_sources = [fileset]

def run_all_datasets(fileset_sources, max_workers=10, max_retries=3):
    failed_datasets = list(fileset_sources[0].keys())  # assume same keys across filesets
    successful = []
    skipped = []

    for attempt in range(1, max_retries + 1):
        if not failed_datasets:
            break

        print(f"\nAttempt {attempt} with {len(failed_datasets)} datasets...")

        # Choose which redirector to use on this attempt
        fileset = fileset_sources[(attempt - 1) % len(fileset_sources)]
        print(f"   → Using redirector set: {['global','fnal','cnaf'][(attempt - 1) % len(fileset_sources)]}")

        remaining = []
        for ds in failed_datasets:
            out_path = os.path.join(parquet_dir, f"{ds}.parquet")
            if is_parquet_valid(out_path, delete_if_invalid=True):
                print(f"Skipping {ds} (already processed)")
                skipped.append(ds)
            else:
                remaining.append(ds)

        if not remaining:
            break

        new_failures = []
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(process_and_save, name, fileset[name]): name
                for name in remaining if name in fileset
            }

            for future in as_completed(futures):
                dataset = futures[future]
                try:
                    out_file = future.result()
                    if out_file and is_parquet_valid(out_file, delete_if_invalid=True):
                        print(f"Finished {dataset}")
                        successful.append(dataset)
                    else:
                        print(f"No valid file for {dataset}")
                        new_failures.append(dataset)
                except Exception as e:
                    print(f"Error processing {dataset}: {e}")
                    new_failures.append(dataset)

        failed_datasets = new_failures

    if failed_datasets:
        print(f"\nStill failed after retries: {failed_datasets}")

    # Final Summary
    print("\nFinal Summary:")
    print(f"Successful: {len(successful)} → {successful}")
    print(f"Skipped (already valid): {len(skipped)} → {skipped}")
    print(f"Failed after all retries: {len(failed_datasets)} → {failed_datasets}")

    return successful, skipped, failed_datasets


def load_with_label(file):
    try:
        label = os.path.basename(file).split(".")[0]
        arr = ak.from_parquet(file)
        return ak.with_field(arr, label, "dataset")
    except Exception as e:
        print(f"Could not load {file}: {e}")
        return None


# Main execution
# successful, skipped, failed = run_all_datasets(fileset, max_workers=10, max_retries=2)

successful, skipped, failed = run_all_datasets(fileset_sources, max_workers=10, max_retries=3)

# Merge only valid parquet files
files = [f for f in glob.glob(os.path.join(parquet_dir, "*.parquet")) if is_parquet_valid(f)]

if files:
    print(f"\nMerging {len(files)} parquet files...")
    labeled_arrays = [arr for arr in (load_with_label(f) for f in files) if arr is not None]

    if labeled_arrays:
        Events = ak.concatenate(labeled_arrays, axis=0)
        table = ak.to_arrow_table(Events)
        pq.write_table(table, "WH_all_nano12.parquet", compression=None)
        print("Done! All datasets saved and merged.")
    else:
        print("No valid parquet files to merge.")
else:
    print("No parquet files found for merging.")
