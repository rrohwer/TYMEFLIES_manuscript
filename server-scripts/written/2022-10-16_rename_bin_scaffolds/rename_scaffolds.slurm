#!/bin/bash
#SBATCH -J renaming
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:0:00
#SBATCH -o renaming.o%j
#SBATCH -e renaming.e%j
#SBATCH --mail-type=all
#SBATCH --mail-user=rohwer@utexas.edu
module load launcher
export LAUNCHER_JOB_FILE=rename_scaffolds.launcher
export LAUNCHER_SCHED=interleaved
${LAUNCHER_DIR}/paramrun
