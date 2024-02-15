#!/bin/bash
#SBATCH -J deinter
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -t 24:00:00
#SBATCH -o deinterleave.o%j
#SBATCH -e deinterleave.e%j
#SBATCH --mail-type=all
#SBATCH --mail-user=rohwer@utexas.edu
module load launcher
export LAUNCHER_JOB_FILE=deinter.sh
export LAUNCHER_SCHED=interleaved
${LAUNCHER_DIR}/paramrun