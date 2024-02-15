#!/bin/bash

# Goal is to run each command on a different node

#SBATCH -J testsubmit
#SBATCH -p development
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -t 00:30:00
#SBATCH -o myjob.o%j
#SBATCH -e myjob.e%j
#SBATCH --mail-type=all
#SBATCH --mail-user=rrr3485@tacc.utexas.edu

module load biocontainers
module load bbmap

echo Hello Kitty 1
echo Hello Kitty 2
echo Hello Kitty 3
echo Hello Kitty 4