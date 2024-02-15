#!/bin/bash
#SBATCH -J targzing
#SBATCH -p normal
#SBATCH -N 4
#SBATCH -n 472
#SBATCH -t 1:00:00
#SBATCH -o targzing.o%j
#SBATCH -e targzing.e%j
#SBATCH --mail-type=all
#SBATCH --mail-user=rohwer@utexas.edu
module load launcher
export LAUNCHER_JOB_FILE=tar_gz_all_sample_folders.launcher
export LAUNCHER_SCHED=interleaved
${LAUNCHER_DIR}/paramrun