#!/bin/bash
#SBATCH -J singlem
#SBATCH -p normal
#SBATCH -N 5
#SBATCH -n 20
#SBATCH -t 36:00:00
#SBATCH -o singlem.o%j
#SBATCH -e singlem.e%j
#SBATCH --mail-type=all
#SBATCH --mail-user=rohwer@utexas.edu
module load launcher
export LAUNCHER_JOB_FILE=singlem.sh
export LAUNCHER_SCHED=interleaved
${LAUNCHER_DIR}/paramrun
