```python
; List of Steps Inputs and Outputs required for preprocessing
; Advised to run separately


;========================= REQUIREMENTS =========================
; REQUIRES APPLICATIONS (command line installation):  ont-guppy basecaller
; - Can be found here: https://nanoporetech.com/software/other/guppy (login required)
; REQUIRED FILES in Input Directory:
; - fast5 files
; returns: compressed fast5 files for Alignment processing
;================================================================
input_path = ""
basecall_dirs = Basecaller.run_guppybasecall(
    input_path, config="rna_r9.4.1_70bps_hac.cfg"
)
print(basecall_dirs)


;========================= REQUIREMENTS =========================
; REQUIRES APPLICATIONS (command line installation):  minimap2, samtools, gzip
; REQUIRED FILES in Input Directory:
; - Compressed fastq files: fastq.gz
; - library.fasta / library.fa of reference sequences
; returns: Alignment Directory of sam, bam, baq, and single.fastq files required for Nanopolish processing
;================================================================
basecall_dir = ""
alignment_dirs = Alignment.run_alignment(basecall_dir)
print(alignment_dirs)


;========================= REQUIREMENTS =========================
; REQUIRED ENVIRONMENT: LINUX preferred
; "compress_fast5", "nanopolish", "minimap2", "samtools", "python3"
; nanopolish is older and can be difficult to install
; REQUIRES APPLICATIONS (command line installation):  nanopolish
; REQUIRED FILES in Input Directory:
; - directory/single fast5 files
; - single.fastq created by Alignment
; - alignment files from Alignment
; - library.fasta / library.fa of reference sequences
; INPUT: Path to Directory of Fast5 files, Path to Directory of Alignment files (output by Alignment)
; returns: tab delimited text file of signal events
;================================================================
fast5_path = ""
alignment_dir = ""
output_file = Signal.run_signal(fast5_path, alignment_dir)
print(output_file)


;========================= REQUIREMENTS =========================
; REQUIRES APPLICATIONS (command line installation):  python, sklearn, multiprocesssing, concurrent
; REQUIRED FILES:
; - nanopolish event align tab delimited text file
; - modification used for naming output
; - number of cpus based on available cores for preprocessing
; returns: robust scaling of event_level_mean and event_length columns for ML processing
;================================================================
input_path = "/Users/timshel/transcripts/ACIM/ACIM_422.txt"
mod = "ACIM"
CPUS = 10
Decon.run_decon(input_path, modification=mod, cpus=CPUS)

flowchart TD
    A[FAST5 Files] --> B[Guppy Basecaller<br/>(Basecaller.run_guppybasecall)]
    B -->|Compressed FASTQ| C[Alignment<br/>(Alignment.run_alignment)]
    C -->|single.fastq + BAM| D[Nanopolish Signal Extraction<br/>(Signal.run_signal)]
    D -->|Event Align TXT| E[Decon Preprocessing<br/>(Decon.run_decon)]
    E --> F[ML-Ready Scaled Features]

    classDef step fill:#f9f9f9,stroke:#333,stroke-width:1px,color:#111;
    class A,B,C,D,E,F step;
