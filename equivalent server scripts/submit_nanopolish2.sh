#!/bin/bash
#SBATCH --account=def-jeromew
#SBATCH --mem=300G                # memory (per node)
#SBATCH --time=0-06:00            # time (DD-HH:MM)

module purge
module load cuda
module load python/3.10
module load StdEnv/2020
module load gcc/9.3.0
module load minimap2
module load samtools
module load nanopolish/0.13.2
virtualenv --no-download $SLURM_TMPDIR/nanopore_env
source $SLURM_TMPDIR/nanopore_env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index -r requirements.txt
pip install --no-index ont-fast5-api

#input start directory $1 
#input output file name <file>.txt $2
if [ $# -lt 3 ];then
  echo "$0: Missing arguments"
  exit 1
elif [ $# -gt 3 ];then
  echo "$0: Too many arguments: $@"
  exit 1
elif [ $# -eq 3 ];then
  echo "Processing $1"
  echo "==========================="
  echo "Number of arguments.: $#"
  echo "List of arguments...: $@"
  echo "Path to Directory ..............: $1"
  echo "Output File ..............: $2"
  echo "Single (1) or Multiple Directories (2) ..............: $3"
  echo "==========================="
fi

myhome=$"$1"

# input is fast5 directory not compressed directory
# single directory fast5
# TODO: add check compressed directory with Alignment and library must exist
if [[ "$3" -eq 1 ]]; then
	f=$"$1"
	if [ -d "$f" ]; then
		echo $f
		fn=$(basename $f)
		echo $fn
		rm -rf "${f}_compressed/Alignment/${f}_gzip"
		echo "${f}_compressed/Alignment/${fn}_gzip"
		compress_fast5 -i "$f" -s "${f}_compressed/Alignment/${fn}_gzip" -c gzip
		cd "${f}_compressed/Alignment/"
		nanopolish index -d "${fn}_gzip" single.fastq
		minimap2 -ax map-ont -t 8 library.fasta single.fastq | samtools sort -o reads-ref.sorted.bam -T reads.tmp
		samtools index reads-ref.sorted.bam
		nanopolish eventalign --reads single.fastq --bam reads-ref.sorted.bam --genome library.fasta --scale-events > "${2}.txt"
		python3 /home/jwbear/projects/def-jeromew/jwbear/count_nanopolish_reads.py "${2}.txt"
		cd $myhome
	fi
elif [[ "$3" -eq 2 ]]; then
	# directory of multiple fast5 directories
	cd $myhome
	for f in *; do
		echo "$f"
		if [ -d "$f" ]; then
			compressed=$(find $f -name "*compressed" | wc -c)
			if [ $compressed -eq 0 ]; then
				echo "$f"
				rm -rf "${f}_gzip"
				rm -rf "${f}_compressed/Alignment/${f}_gzip"
				compress_fast5 -i "$f" -s "${f}_gzip" -c gzip
				mv "${f}_gzip" "${f}_compressed/Alignment/"
				cd "${f}_compressed/Alignment/"
				nanopolish index -d "${f}_gzip" single.fastq
				minimap2 -ax map-ont -t 8 library.fasta single.fastq | samtools sort -o reads-ref.sorted.bam -T reads.tmp
				samtools index reads-ref.sorted.bam
				nanopolish eventalign --reads single.fastq --bam reads-ref.sorted.bam --genome library.fasta --scale-events > "${2}.txt"
				python3 /home/jwbear/projects/def-jeromew/jwbear/count_nanopolish_reads.py "${2}.txt"
				cd $myhome
			fi
		fi
	done
else
	echo "Invalid Option: $3"
	exit 1
fi
