#!/bin/bash
#SBATCH -J testplots
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 10:0:00
#SBATCH -o instrain.o%j
#SBATCH -e instrain.e%j
#SBATCH --mail-type=all
#SBATCH --mail-user=rohwer@utexas.edu
module load launcher
export LAUNCHER_JOB_FILE=test_instrain_profile_plots.launcher
export LAUNCHER_SCHED=interleaved
${LAUNCHER_DIR}/paramrun
