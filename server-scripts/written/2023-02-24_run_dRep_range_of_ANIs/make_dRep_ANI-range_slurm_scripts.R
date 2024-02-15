# RRR
# Realizing that the dRep ANI choice is really important to downstream interpretation/analysis
# Last time I just did the default- 95%
# This time I will do a range, and look at how the number of representative genomes changes
# Ideally, there will be an ANI threshold where number rep genomes changes more sharply, or levels off
# For the plot, need to do ANI's that are also coarser than 95% to see the curve
# so do 90-99% ANI, that will be 10 runs
# each run took 10 hrs on one node before, using all the threads
# I can't find in my notes if it actually used all the threads, but I'm going to assume it did and not do additional testing

# the command from before:
# source /work/08922/rrohwer/ls6/miniforge3/etc/profile.d/conda.sh ; 
# conda activate drep ; 
# dRep dereplicate drep_output 
#   -p 128 
#   -g bin_locations_for_dRep.txt 
#   -comp 50 
#   -con 10 
#   --genomeInfo bin_info_for_dRep.csv 
#   --S_algorithm fastANI 
#   --P_ani 0.9 
#   --S_ani 0.95 
#   --multiround_primary_clustering 
#   --skip_plots

# run from within "rundrep" folder that contains:
#   bin info file
#   bin locations file
#   metabat2_bins_renamed/

# takes 14 hrs to run

ani.range <- seq(from = .90, to = .99, by = .01)

drep.slurm <- c("#!/bin/bash",
                "#SBATCH -J dRep",
                "#SBATCH -p normal",
                "#SBATCH -N 10",
                "#SBATCH -n 10",
                "#SBATCH -t 20:00:00",
                "#SBATCH -o drep.o%j",
                "#SBATCH -e drep.e%j",
                "#SBATCH --mail-type=all",
                "#SBATCH --mail-user=rohwer@utexas.edu",
                # "module load biocontainers",
                "module load launcher",
                "export LAUNCHER_JOB_FILE=run_drep_ANIrange.launcher",
                "export LAUNCHER_SCHED=interleaved",
                "${LAUNCHER_DIR}/paramrun")

drep.launcher <- paste0("source /work/08922/rrohwer/ls6/miniforge3/etc/profile.d/conda.sh ; conda activate drep ; ",
                        "dRep dereplicate drep_output_",ani.range," -p 128 -g bin_locations_for_dRep.txt -comp 50 -con 10 --genomeInfo bin_info_for_dRep.csv ",
                        "--S_algorithm fastANI --P_ani 0.9 --S_ani ",ani.range," --multiround_primary_clustering --skip_plots")

write.table(x = drep.slurm, file = "server-scripts/generated/2023-02-24_run_dRep_range_of_ANIs/run_drep_ANIrange.slurm",
            row.names = F, col.names = F, quote = F)
write.table(x = drep.launcher, file = "server-scripts/generated/2023-02-24_run_dRep_range_of_ANIs/run_drep_ANIrange.launcher",
            row.names = F, col.names = F, quote = F)



