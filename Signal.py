#!/usr/bin/env python3
import os
import subprocess
from pathlib import Path
import shutil
import pandas as pd

# -------------------------------
# Count sequences from signal file
# -------------------------------
def count_sequences(signal_file: str):
    signal_path = Path(signal_file)
    if not signal_path.exists():
        raise FileNotFoundError(f"Signal file not found: {signal_file}")

    print(f"[INFO] Counting sequences in signal file: {signal_file}")

    df = pd.read_csv(signal_file, sep="\t", usecols=['contig', 'read_index'])
    counts = df.groupby('contig')['read_index'].nunique()
    counts_dict = counts.to_dict()

    for contig, count in counts_dict.items():
        print(f"{contig}\t{count}")

    return counts_dict


# -------------------------------
# Dependency check
# -------------------------------
def check_dependency(tool: str):
    if shutil.which(tool) is None:
        raise EnvironmentError(f"Required tool '{tool}' is not installed or not in PATH.")


# -------------------------------
# Run shell command helper
# -------------------------------
def run_shell_command(cmd, cwd=None):
    print(f"\n[RUNNING] {cmd}")
    result = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr)
        raise RuntimeError(f"Command failed: {cmd}")


# -------------------------------
# Verify required files
# -------------------------------
def verify_required_files(alignment_dir: Path):
    fastq_file = alignment_dir / "single.fastq"
    fasta_file = alignment_dir / "library.fasta"

    missing = []
    if not fastq_file.exists():
        missing.append(str(fastq_file))
    if not fasta_file.exists():
        missing.append(str(fasta_file))
    if missing:
        raise FileNotFoundError(f"Missing required files:\n" + "\n".join(missing))

    return fastq_file, fasta_file


# -------------------------------
# Prepare compressed fast5 input
# -------------------------------
def prepare_compressed_fast5(fast5_path: Path):
    fast5_name = fast5_path.stem if fast5_path.is_file() else fast5_path.name

    # Compressed dir inside fast5_path if dir, or fast5_path.parent if file
    if fast5_path.is_file():
        compressed_dir = fast5_path.parent / f"{fast5_name}_gzip"
    else:
        compressed_dir = fast5_path / f"{fast5_name}_gzip"

    if compressed_dir.exists():
        print(f"[INFO] Removing old compressed directory: {compressed_dir}")
        shutil.rmtree(compressed_dir)

    print(f"[INFO] Creating compressed directory: {compressed_dir}")

    cmd = f"compress_fast5 -i {fast5_path} -s {compressed_dir} -c gzip"
    run_shell_command(cmd)

    return compressed_dir


# -------------------------------
# Ensure unique filename
# -------------------------------
def get_unique_filename(base_path: Path, base_name: str, extension: str = ".txt") -> Path:
    candidate = base_path / f"{base_name}{extension}"
    counter = 1
    while candidate.exists():
        candidate = base_path / f"{base_name}_{counter}{extension}"
        counter += 1
    return candidate


# -------------------------------
# Process single fast5 directory or file
# -------------------------------
def process_single_fast5(fast5_path: Path, alignment_dir: Path):
    if not fast5_path.exists():
        raise FileNotFoundError(f"{fast5_path} does not exist")

    # Step 1: compress fast5
    compressed_dir = prepare_compressed_fast5(fast5_path)

    # Step 2: Verify required files exist in alignment directory
    fastq_file, fasta_file = verify_required_files(alignment_dir)

    # Step 3: Nanopolish indexing
    run_shell_command(f"nanopolish index -d {compressed_dir} {fastq_file}", cwd=alignment_dir)

    # Step 4: Mapping and sorting
    run_shell_command(
        f"minimap2 -ax map-ont -t 8 {fasta_file} {fastq_file} | "
        f"samtools sort -o reads-ref.sorted.bam -T reads.tmp", cwd=alignment_dir
    )
    run_shell_command("samtools index reads-ref.sorted.bam", cwd=alignment_dir)

    # Step 5: Nanopolish eventalign with unique filename
    output_txt = get_unique_filename(alignment_dir, "aligned_signal")
    run_shell_command(
        f"nanopolish eventalign --reads {fastq_file} --bam reads-ref.sorted.bam "
        f"--genome {fasta_file} --scale-events > {output_txt}", cwd=alignment_dir
    )

    print(f"[INFO] Output txt file: {output_txt}")
    return output_txt


# -------------------------------
# Run nanopolish pipeline (single input only)
# -------------------------------
def run_nanopolish_pipeline(fast5_path: str, alignment_dir: str):
    fast5_path = Path(fast5_path)
    alignment_dir = Path(alignment_dir)

    for tool in ["compress_fast5", "nanopolish", "minimap2", "samtools", "python3"]:
        check_dependency(tool)

    if fast5_path.is_file() and fast5_path.suffix == ".fast5":
        print("[INFO] Detected single fast5 file")
        output_txt = process_single_fast5(fast5_path, alignment_dir)

    elif fast5_path.is_dir():
        print("[INFO] Detected single fast5 directory")
        output_txt = process_single_fast5(fast5_path, alignment_dir)

    else:
        raise FileNotFoundError(f"{fast5_path} does not exist or is not a valid .fast5 file/directory")

    # Count sequences in final output
    count_sequences(output_txt)
    # Print the path to the generated output txt
    print(f"[INFO] Generated output file: {output_txt}")

    return output_txt


# -------------------------------
# Example usage
# -------------------------------
if __name__ == "__main__":
    fast5_dir = "fast5_422"
    alignment_dir = "Alignment"

    run_nanopolish_pipeline(fast5_dir, alignment_dir)
