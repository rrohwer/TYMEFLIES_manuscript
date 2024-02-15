#!/bin/bash

# have say_hello.sh in same directory
# source as: sbatch ./slurm_hello.sh
# type sqs to check on the job's status

#SBATCH --qos=debug
#SBATCH --time=00:05:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=32
#SBATCH --constraint=haswell

export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_NUM_THREADS=16

srun ./say_hello.sh 1 &
srun ./say_hello.sh 2 &
srun ./say_hello.sh 3 &
srun ./say_hello.sh 4 &
wait


# Notes:
# --qos   what priority the job has
# --time  maximum walltime (time I wait, not counting paralelle processes)
# --nodes there are 9,000 total nodes 
# --cpus-per-task each node has 32 cpus-per-task (aka cores)
# --ntasks-per-node 
# --constraint which server to use (haswell is faster is not as paralellized)
# export OMP things: don't understand but those are the cori defaults