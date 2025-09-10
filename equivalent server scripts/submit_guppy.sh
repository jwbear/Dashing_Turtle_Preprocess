#!/bin/bash
#SBATCH --account=def-jeromew
#SBATCH --gres=gpu:1              # Number of GPUs (per node)
#SBATCH --mem=100G                # memory (per node)
#SBATCH --time=0-01:00            # time (DD-HH:MM)

module purge
module load cuda
module load python/3.10
virtualenv --no-download $SLURM_TMPDIR/nanopore_env
source $SLURM_TMPDIR/nanopore_env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index -r requirements.txt

#input is directory of fast5 files, compressed fast5 directory is created and required by submit_aln
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
# single directory fast5
if [[ "$2" -eq 1 ]]; then
        f=$"$1"
        if [ -d "$f" ]; then
                echo $f
                /home/jwbear/projects/def-jeromew/jwbear/ont-guppy/bin/guppy_basecaller --compress_fastq -i "$f" -s "$f"_compressed -c rna_r9.4.1_70bps_hac.cfg --device cuda:0 --disable_qscore_filtering
                cd $myhome
        fi
elif [[ "$2" -eq 2 ]]; then
        # directory of multiple fast5 directories
        cd $myhome
        for f in *; do
                echo "$f"
                if [ -d "$f" ]; then
                        compressed=$(find $f -name "*compressed" | wc -c)
                        if [ $compressed -eq 0 ]; then
				/home/jwbear/projects/def-jeromew/jwbear/ont-guppy/bin/guppy_basecaller --compress_fastq -i "$f" -s "$f"_compressed -c rna_r9.4.1_70bps_hac.cfg --device cuda:0 --disable_qscore_filtering
                                cd $myhome
                        fi
                fi
        done
else
        echo "Invalid Option: $2"
        exit 1
fi        
        
