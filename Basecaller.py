#!/usr/bin/env python3
import os
import subprocess
import argparse
import shutil
from pathlib import Path


def find_basecaller():
    """Attempt to find guppy_basecaller in PATH."""
    guppy_path = shutil.which("guppy_basecaller")
    if guppy_path is None:
        raise FileNotFoundError(
            "Guppy not found. Ensure guppy_basecaller is installed and in your PATH.\n"
            "Dorado is not yet supported.\n"
            "Download Guppy: https://nanoporetech.com/software/other/guppy"
        )
    print(f"Guppy Path Found: {guppy_path}")
    return Path(guppy_path)


def find_fast5_inputs(input_path):
    """Detect .fast5 input: single file or directories containing .fast5 files."""
    input_path = Path(input_path)

    if input_path.is_file() and input_path.suffix == ".fast5":
        return [str(input_path)]

    elif input_path.is_dir():
        fast5_dirs = []
        for root, _, files in os.walk(input_path):
            if any(f.endswith(".fast5") for f in files):
                fast5_dirs.append(root)
        if not fast5_dirs:
            raise FileNotFoundError(f"No .fast5 files found in {input_path}")
        return fast5_dirs

    else:
        raise FileNotFoundError(f"{input_path} is not a valid .fast5 file or directory")


def validate_device(device):
    """Basic validation of device string."""
    if device is None:
        return None
    device = str(device).lower()
    if device != "cpu" and not device.startswith("cuda:"):
        raise ValueError(f"Invalid device string '{device}'. Use 'cpu' or 'cuda:N'.")
    return device


def unique_output_dir(base_path):
    """Ensure output directory is unique if *_compressed already exists."""
    base = Path(base_path)
    if not base.exists():
        return base
    i = 1
    while True:
        candidate = Path(f"{base}_{i}")
        if not candidate.exists():
            return candidate
        i += 1


def run_guppy(guppy_path, fast5_input, config, device=None, output_override=None):
    """Run guppy_basecaller and return output directory."""
    fast5_path = Path(fast5_input)

    if output_override:
        if len(str(fast5_path).split(os.sep)) > 1:
            output_dir = Path(output_override) / fast5_path.name
        else:
            output_dir = Path(output_override)
    else:
        if fast5_path.is_file():
            output_dir = unique_output_dir(fast5_path.parent / f"{fast5_path.stem}_compressed")
        else:
            output_dir = unique_output_dir(Path(f"{fast5_path}_compressed"))

    output_dir.mkdir(parents=True, exist_ok=True)
    input_arg = str(fast5_path)

    print(f"Processing {fast5_path} -> {output_dir}")

    cmd = [
        str(guppy_path),
        "--compress_fastq",
        "-i", input_arg,
        "-s", str(output_dir),
        "-c", str(config),
        "--disable_qscore_filtering"
    ]
    if device is not None:
        cmd.extend(["--device", str(device)])

    # Run Guppy with live stdout (no log)
    subprocess.run(cmd, check=True)
    print(f"✅ Completed basecalling for {fast5_path}")
    return output_dir


def run_basecaller(input_path, config="rna_r9.4.1_70bps_hac.cfg", device=None, output=None):
    """Wrapper for Guppy basecalling."""
    input_path = Path(input_path)
    if not input_path.exists():
        raise FileNotFoundError(f"{input_path} does not exist")

    guppy_path = find_basecaller()
    device = validate_device(device)

    fast5_inputs = find_fast5_inputs(input_path)
    print(f"Found {len(fast5_inputs)} fast5 inputs.")

    output_dirs = []
    for f in fast5_inputs:
        out = run_guppy(guppy_path, f, config, device=device, output_override=output)
        output_dirs.append(out)

    return output_dirs


def main():
    parser = argparse.ArgumentParser(description="Run Guppy basecaller on fast5 files.")
    parser.add_argument("input", help="Path to .fast5 file or directory")
    parser.add_argument(
        "-c", "--config",
        default="rna_r9.4.1_70bps_hac.cfg",
        help="Guppy config file (default: rna_r9.4.1_70bps_hac.cfg)."
    )
    parser.add_argument(
        "-d", "--device",
        help="Device for Guppy (e.g. 'cuda:0' or 'cpu')."
    )
    parser.add_argument(
        "-o", "--output",
        help="Optional output directory. Defaults to *_compressed next to input."
    )
    args = parser.parse_args()

    outputs = run_basecaller(args.input, config=args.config, device=args.device, output=args.output)

    print("\nAll outputs written to:")
    for o in outputs:
        print(f" - {o}")

if __name__ == "__main__":
    #main()
    run_basecaller("/Users/timshel/transcripts/fast5_422")
