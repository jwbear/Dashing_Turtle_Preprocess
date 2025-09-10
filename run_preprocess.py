import Basecaller
import Alignment
import Signal
import Decon


#### List of Steps Inputs and Outputs required for preprocessing
#### Advised to run separately


########################## REQUIREMENTS ##################
# REQUIRES APPLICATIONS (command line installation):  ont-guppy basecaller
# - Can be found her https://nanoporetech.com/software/other/guppy (login required)
# REQUIRED FILES in Input Directory:
# - fast5 files
# returns: compressed fast5 files for Alignment processing
###########################################################
input_path=""
output_file=""
basecall_dirs = Basecaller.run_guppybasecall(input_path, config="rna_r9.4.1_70bps_hac.cfg")
print(basecall_dirs)
########################## REQUIREMENTS ##################
# REQUIRES APPLICATIONS (command line installation):  minimap2, samtools, gzip
# REQUIRED FILES in Input Directory:
# - Compressed fastq files: fastq.gz
# - library.fasta/ library.fa of reference sequences
# returns: Alignment Directory of sam, bam, baq, and single.fastq files required for Nanoplish processing
###########################################################
bascall_dir =""
alignment_dirs = Alignment.run_alignment(basecall_dir)
print(alignment_dirs)
########################## REQUIREMENTS ##################
# REQUIRED ENVIRONMENT: LINUX preferred
# "compress_fast5", "nanopolish", "minimap2", "samtools", "python3"
# nanopolish is older and can be difficult to install
# REQUIRES APPLICATIONS (command line installation):  nanopolish
# REQUIRED FILES in Input Directory:
# - directory/single fast5 files
# - single.fastq created by Alignment
# - alignment files from Alignment
# - library.fasta/ library.fa of reference sequences
# INPUT: Path to Directory of Fast 5 files, Path to Directory of Alignment files (output by Alignment)
# returns: tab delimited text file of signal events
###########################################################
fast5_path = ""
alignment_dir = ""
output_file = Signal.run_signal(fast5_path, alignment_dir)
print(output_file)
########################## REQUIREMENTS ##################
# REQUIRES APPLICATIONS (command line installation):  python, sklearn, multiprocesssing, concurrent
# REQUIRED FILES:
# - nanopolish event align tab delimited text file
# - modification used for naming output
# - number of cpus based on available cores for preprocessing
# returns: robust scaling of event_level_mean and event_length columns for ML processing
###########################################################
input_path = "/Users/timshel/transcripts/ACIM/ACIM_422.txt"
mod = "ACIM"
CPUS = 10
Decon.run_decon(input_path, modification=mod, cpus=CPUS)
