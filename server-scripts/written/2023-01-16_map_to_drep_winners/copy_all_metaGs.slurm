#!/bin/bash
#SBATCH -J copying
#SBATCH -p normal
#SBATCH -N 4
#SBATCH -n 472
#SBATCH -t 6:00:00
#SBATCH -o copying.o%j
#SBATCH -e copying.e%j
#SBATCH --mail-type=all
#SBATCH --mail-user=rohwer@utexas.edu
module load launcher
export LAUNCHER_JOB_FILE=copy_all_metaGs.launcher
export LAUNCHER_SCHED=interleaved
${LAUNCHER_DIR}/paramrun