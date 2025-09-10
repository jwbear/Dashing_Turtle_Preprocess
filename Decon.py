import sys
import os
import pandas as pd
from sklearn.preprocessing import RobustScaler
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
import uuid
from datetime import datetime
import numpy as np

cols_to_scale = ['event_level_mean', 'event_length', 'event_stdv']


def process_read_chunk(df_subset, read_list):
    """Process a chunk of reads and return concatenated DataFrame."""
    result = []
    for read in read_list:
        dft = df_subset[df_subset['read_index'] == read].copy()
        dft[cols_to_scale] = RobustScaler().fit_transform(dft[cols_to_scale])
        dft = dft.sort_values(by=['contig', 'read_index', 'position'])
        result.append(dft)
    return pd.concat(result, ignore_index=True)


def chunk_reads(reads, chunk_size):
    """Split a list of reads into chunks."""
    for i in range(0, len(reads), chunk_size):
        yield reads[i:i + chunk_size]


def preprocess_df(df, mod, max_workers=5, chunk_size=100):
    """Parallel preprocessing of reads."""
    reads = df['read_index'].unique()
    results = []

    # Split reads into chunks
    read_chunks = list(chunk_reads(reads, chunk_size))

    # Use only the subset of DataFrame relevant for each chunk
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(process_read_chunk, df[df['read_index'].isin(chunk)], chunk): chunk
            for chunk in read_chunks
        }
        for future in as_completed(futures):
            results.append(future.result())

    # Combine all processed reads
    dfn = pd.concat(results, ignore_index=True)

    # Save per contig
    for seq in dfn['contig'].unique():
        run_id = datetime.now().strftime("%Y%m%d_%H%M%S") + "_" + uuid.uuid4().hex[:8]
        dft = dfn[dfn['contig'] == seq]
        fn = f"{mod}_{seq}_{run_id}_decon_signal_preprocess.csv"
        dft.to_csv(fn, index=False)

    return


def run_decon(input_path, modification, cpus):
    cpu = min(multiprocessing.cpu_count(), cpus)
    multiprocessing.freeze_support()  # Required on macOS/Windows
    df = pd.read_csv(input_path, sep='\t')
    preprocess_df(df=df, mod=modification, max_workers=cpu)
    print("Processing complete.")
    return


if __name__ == "__main__":
    input_path = "/Users/timshel/transcripts/ACIM_422.txt"
    mod = "ACIM"
    CPUS = min(multiprocessing.cpu_count(), 10)
    run_decon(input_path, modification=mod, cpus=CPUS)
