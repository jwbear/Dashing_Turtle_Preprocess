#!/bin/bash
#SBATCH --account=def-jeromew
#SBATCH --gres=gpu:1              # Number of GPUs (per node)
#SBATCH --mem=100G                # memory (per node)
#SBATCH --time=0-02:00            # time (DD-HH:MM)

module purge
module load cuda
module load python/3.10
module load minimap2
module load samtools
virtualenv --no-download $SLURM_TMPDIR/nanopore_env
source $SLURM_TMPDIR/nanopore_env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index -r requirements.txt


# input is compressed directory, or main directory containing compressed files
#input start directory $1 
if [ $# -lt 2 ];then
  echo "$0: Missing arguments"
  exit 1
elif [ $# -gt 2 ];then
  echo "$0: Too many arguments: $@"
  exit 1
elif [ $# -eq 2 ];then
  echo "Processing $1"
  echo "==========================="
  echo "Number of arguments.: $#"
  echo "List of arguments...: $@"
  echo "Path to Directory ..............: $1"
  echo "Single (1) or Multiple Directories (2) ..............: $2"
  echo "==========================="
fi

myhome=$"$1"
# single directory 
if [[ "$2" -eq 1 ]]; then
        f=$"$1"
        if [ -d "$f" ]; then
        		cd $f
			echo $f
			rm -r *.fastq
			mkdir -p Alignment
			gzip -dk *.gz
			cat *.fastq > single.fastq
			mv single.fastq Alignment/single.fastq
			if [ ! -f Alignment/library.fasta ]; then
					mv library.fasta Alignment/library.fasta
			fi
			cd Alignment
			if [ ! -f library.fasta ]; then
					echo "library.fasta file not found."
					exit 1
			fi
			minimap2 -ax map-ont library.fasta single.fastq > aln.sam 
			samtools view -bT library.fasta -F 16 aln.sam > aln.bam
			samtools sort aln.bam > aln.s.bam
			samtools index aln.s.bam 
			samtools calmd -bAr aln.s.bam library.fasta > aln.baq.bam
			samtools index aln.baq.bam
			cd $myhome
		fi
elif [[ "$2" -eq 2 ]]; then
# multiple directory
	cd $myhome
	for f in *; do
		if [ -d "$f" ]; then
			compressed=$(find $f -name "*compressed" | wc -c)
			if [ ! $compressed -eq 0 ]; then
					cd $f
					echo $f
					rm -r *.fastq
					mkdir -p Alignment
					gzip -dk *.gz
					cat *.fastq > single.fastq
					mv single.fastq Alignment/single.fastq
					if [ ! -f Alignment/library.fasta ]; then
							mv library.fasta Alignment/library.fasta
					fi
					cd Alignment
					if [ ! -f library.fasta ]; then
							echo "library.fasta file not found."
							exit 1
					fi
					minimap2 -ax map-ont library.fasta single.fastq > aln.sam 
					samtools view -bT library.fasta -F 16 aln.sam > aln.bam
					samtools sort aln.bam > aln.s.bam
					samtools index aln.s.bam 
					samtools calmd -bAr aln.s.bam library.fasta > aln.baq.bam
					samtools index aln.baq.bam
					cd $myhome
			fi
		fi
	done
else
        echo "Invalid Option: $2"
        exit 1
fi


