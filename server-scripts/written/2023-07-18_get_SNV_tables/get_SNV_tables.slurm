#!/bin/bash
#SBATCH -J SNVs
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -t 00:02:00
#SBATCH -o snvs.o%j
#SBATCH -e snvs.e%j
#SBATCH --mail-type=all
#SBATCH --mail-user=rohwer@utexas.edu
module load launcher
export LAUNCHER_JOB_FILE=get_SNV_tables.launcher
export LAUNCHER_SCHED=interleaved
${LAUNCHER_DIR}/paramrun
