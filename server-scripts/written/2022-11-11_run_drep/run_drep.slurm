#!/bin/bash
#SBATCH -J drep
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:0:00
#SBATCH -o drep.o%j
#SBATCH -e drep.e%j
#SBATCH --mail-type=all
#SBATCH --mail-user=rohwer@utexas.edu
module load launcher
export LAUNCHER_JOB_FILE=run_drep.launcher
export LAUNCHER_SCHED=interleaved
${LAUNCHER_DIR}/paramrun
