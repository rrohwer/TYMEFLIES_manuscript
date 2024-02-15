#!/bin/bash
#SBATCH -N 5
#SBATCH -n 5
#SBATCH --output out_%j.txt

# replace these with relevant module loads
module purge
module load Python/intel/3.6.9/intel-19.0.5

# mpiexec may need to be replaced with ibrun(?)
mpiexec -n 5 python pyscript.py

