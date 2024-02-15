#!/bin/bash
#SBATCH -J tarring
#SBATCH -p normal
#SBATCH -N 4
#SBATCH -n 472
#SBATCH -t 1:00:00
#SBATCH -o tarring.o%j
#SBATCH -e tarring.e%j
#SBATCH --mail-type=all
#SBATCH --mail-user=rohwer@utexas.edu
module load launcher
export LAUNCHER_JOB_FILE=tar_all_sample_folders.launcher
export LAUNCHER_SCHED=interleaved
${LAUNCHER_DIR}/paramrun