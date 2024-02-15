# RRR
# Need to run prodigal before doing inStrain
# better to run on individual genomes, not concatenated
# it does some machine learning on a per-genome level that way
# For draft genomes use normal mode

# from inStrain example command: (-e 1 means no genes in gaps)
# prodigal -i genome.fna -o genome.genes -a genome.faa -d genome.genes.fna -m -p single
# it's so fast, seconds

all.bins <- readRDS("data/2022-11-15_bin_stats/all_bins.rds")
drep.bins <- all.bins[all.bins$winner, ]
bin.paths <- paste0("../metabat2_bins_renamed/",drep.bins$tymeflies.name,"/",drep.bins$bin.filename)
genes.files <- paste0(drep.bins$bin.full.name,".genes")
aa.files <- paste0(drep.bins$bin.full.name,".genes.faa")
gene.fna <- paste0(drep.bins$bin.full.name,".genes.fna")

prodigal.launcher <- paste0("source /work/08922/rrohwer/ls6/miniforge3/etc/profile.d/conda.sh ; conda activate instrain ; prodigal -i ",bin.paths," -o ",genes.files," -a ",aa.files," -d ",gene.fna," -m -p single")

prodigal.slurm <- c("#!/bin/bash",
                    "#SBATCH -J prodigal",
                    "#SBATCH -p development", # fast enough for development queue
                    paste("#SBATCH -N 1"),
                    paste("#SBATCH -n 128"),
                    paste("#SBATCH -t 00:15:00"),
                    "#SBATCH -o prodigal.o%j",
                    "#SBATCH -e prodigal.e%j",
                    "#SBATCH --mail-type=all",
                    "#SBATCH --mail-user=rohwer@utexas.edu",
                    "module load launcher",
                    "export LAUNCHER_JOB_FILE=run_prodigal.launcher",
                    "export LAUNCHER_SCHED=interleaved",
                    "${LAUNCHER_DIR}/paramrun")

write.table(x = prodigal.launcher, file = "server-scripts/generated/2023-02-11_run_prodigal/run_prodigal.launcher", 
            row.names = F, col.names = F, quote = F)
write.table(x = prodigal.slurm, file = "server-scripts/generated/2023-02-11_run_prodigal/run_prodigal.slurm", 
            row.names = F, col.names = F, quote = F)

# Slurm Job_id=694767 Name=prodigal Ended, Run time 00:04:56, COMPLETED, ExitCode 0

# after creating all the files, concatenate them for input to instrain:
# cat *.genes.fna > drep_winners_concat_genes.fna   **** this one used by inStrain ****
# cat *.faa > drep_winners_concat.faa