#!/bin/bash

# be in the bin.42 directory with ref already generated
# source as: sbatch ./slurm_hello.sh
# type sqs to check on the job's status

#SBATCH --qos=debug
#SBATCH --time=00:30:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=64
#SBATCH --constraint=haswell

export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_NUM_THREADS=16

srun ./map_to_bin_42.sh 


# Notes:
# --qos   what priority the job has
# --time  maximum walltime (time I wait, not counting paralelle processes)
# --nodes there are 9,000 total nodes 
# --cpus-per-task each node has 32 cpus-per-task (aka cores)
# --ntasks-per-node 
# --constraint which server to use (haswell is faster is not as paralellized)
# export OMP things: don't understand but those are the cori defaults