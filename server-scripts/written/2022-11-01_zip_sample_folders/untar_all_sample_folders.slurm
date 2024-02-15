#!/bin/bash
#SBATCH -J untarring
#SBATCH -p normal
#SBATCH -N 4
#SBATCH -n 472
#SBATCH -t 1:00:00
#SBATCH -o untarring.o%j
#SBATCH -e untarring.e%j
#SBATCH --mail-type=all
#SBATCH --mail-user=rohwer@utexas.edu
module load launcher
export LAUNCHER_JOB_FILE=untar_all_sample_folders.launcher
export LAUNCHER_SCHED=interleaved
${LAUNCHER_DIR}/paramrun