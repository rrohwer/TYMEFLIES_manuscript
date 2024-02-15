# RRR
# as submit scripts, can copy in paralelle and much faster

sample.names <- readRDS("data/2022-11-15_bin_stats/all_bins.rds")
sample.names <- unique(sample.names$tymeflies.name)

sample.names.tar <- paste0(sample.names,".tar")
sample.names.tar.gz <- paste0(sample.names,".tar.gz")

sample.names.fastq <- paste0(sample.names, "_filtered-reads.fastq.gz")
sample.names.assemblies <- paste0(sample.names, "_assembled_contigs.fasta")
sample.names.instrain.folders <- paste0(sample.names, ".IS")

# move fastqs ----

fastqs.slurm <- c("#!/bin/bash",
                  "#SBATCH -J copyfastqs",
                  "#SBATCH -p development",
                  "#SBATCH -N 4",
                  "#SBATCH -n 472",
                  "#SBATCH -t 1:00:00",
                  "#SBATCH -o copyfastqs.o%j",
                  "#SBATCH -e copyfastqs.e%j",
                  "#SBATCH --mail-type=all",
                  "#SBATCH --mail-user=rohwer@utexas.edu",
                  "module load launcher",
                  "export LAUNCHER_JOB_FILE=copy_fastqs.launcher",
                  "export LAUNCHER_SCHED=interleaved",
                  "${LAUNCHER_DIR}/paramrun")

fastqs.launcher <- paste0("cp /corral/utexas/MCB22035/tymeflies/filtered_reads/", sample.names.fastq, " ./")

write.table(x = fastqs.slurm, file = "server-scripts/generated/2023-02-24_copying_files_in_paralelle/copy_fastqs.slurm",
            row.names = F, col.names = F, quote = F)

write.table(x = fastqs.launcher, file = "server-scripts/generated/2023-02-24_copying_files_in_paralelle/copy_fastqs.launcher",
            row.names = F, col.names = F, quote = F)

# move assemblies ----

assemblies.slurm <- c("#!/bin/bash",
                  "#SBATCH -J copyassemblies",
                  "#SBATCH -p development",
                  "#SBATCH -N 4",
                  "#SBATCH -n 472",
                  "#SBATCH -t 1:00:00",
                  "#SBATCH -o copyassemblies.o%j",
                  "#SBATCH -e copyassemblies.e%j",
                  "#SBATCH --mail-type=all",
                  "#SBATCH --mail-user=rohwer@utexas.edu",
                  "module load launcher",
                  "export LAUNCHER_JOB_FILE=copy_assemblies.launcher",
                  "export LAUNCHER_SCHED=interleaved",
                  "${LAUNCHER_DIR}/paramrun")

assemblies.launcher <- paste0("cp /corral/utexas/MCB22035/tymeflies/assembled_contigs/", sample.names.assemblies, " ./")

write.table(x = assemblies.slurm, file = "server-scripts/generated/2023-02-24_copying_files_in_paralelle/copy_assemblies.slurm",
            row.names = F, col.names = F, quote = F)

write.table(x = assemblies.launcher, file = "server-scripts/generated/2023-02-24_copying_files_in_paralelle/copy_assemblies.launcher",
            row.names = F, col.names = F, quote = F)

# move zipped bins ----

bins.slurm <- c("#!/bin/bash",
                "#SBATCH -J copybins",
                "#SBATCH -p development",
                "#SBATCH -N 4",
                "#SBATCH -n 472",
                "#SBATCH -t 1:00:00",
                "#SBATCH -o copybins.o%j",
                "#SBATCH -e copybins.e%j",
                "#SBATCH --mail-type=all",
                "#SBATCH --mail-user=rohwer@utexas.edu",
                "module load launcher",
                "export LAUNCHER_JOB_FILE=copy_bins.launcher",
                "export LAUNCHER_SCHED=interleaved",
                "${LAUNCHER_DIR}/paramrun")

bins.launcher <- paste0("cp /corral/utexas/MCB22035/tymeflies/metabat2_bins_renamed/", sample.names.tar.gz, " ./")

write.table(x = bins.slurm, file = "server-scripts/generated/2023-02-24_copying_files_in_paralelle/copy_bins.slurm",
            row.names = F, col.names = F, quote = F)

write.table(x = bins.launcher, file = "server-scripts/generated/2023-02-24_copying_files_in_paralelle/copy_bins.launcher",
            row.names = F, col.names = F, quote = F)

# move zipped ... ----

# copy to make more moving scripts

# untar.gz zipped sample-folders ----

untar.gz.slurm <- c("#!/bin/bash",
                    "#SBATCH -J untargz",
                    "#SBATCH -p development",
                    "#SBATCH -N 4",
                    "#SBATCH -n 472",
                    "#SBATCH -t 1:00:00",
                    "#SBATCH -o untargz.o%j",
                    "#SBATCH -e untargz.e%j",
                    "#SBATCH --mail-type=all",
                    "#SBATCH --mail-user=rohwer@utexas.edu",
                    "module load launcher",
                    "export LAUNCHER_JOB_FILE=untargz_sample_folders.launcher",
                    "export LAUNCHER_SCHED=interleaved",
                    "${LAUNCHER_DIR}/paramrun")

untar.gz.launcher <- paste0("tar xvfz ", sample.names.tar.gz)

write.table(x = untar.gz.slurm, file = "server-scripts/generated/2023-02-24_copying_files_in_paralelle/untargz_sample_folders.slurm",
            row.names = F, col.names = F, quote = F)

write.table(x = untar.gz.launcher, file = "server-scripts/generated/2023-02-24_copying_files_in_paralelle/untargz_sample_folders.launcher",
            row.names = F, col.names = F, quote = F)

# untar zipped sample-folders ----

untar.slurm <- c("#!/bin/bash",
                    "#SBATCH -J untar",
                    "#SBATCH -p development",
                    "#SBATCH -N 4",
                    "#SBATCH -n 472",
                    "#SBATCH -t 1:00:00",
                    "#SBATCH -o untar.o%j",
                    "#SBATCH -e untar.e%j",
                    "#SBATCH --mail-type=all",
                    "#SBATCH --mail-user=rohwer@utexas.edu",
                    "module load launcher",
                    "export LAUNCHER_JOB_FILE=untar_sample_folders.launcher",
                    "export LAUNCHER_SCHED=interleaved",
                    "${LAUNCHER_DIR}/paramrun")

untar.launcher <- paste0("tar xvf ", sample.names.tar)

write.table(x = untar.slurm, file = "server-scripts/generated/2023-02-24_copying_files_in_paralelle/untar_sample_folders.slurm",
            row.names = F, col.names = F, quote = F)

write.table(x = untar.launcher, file = "server-scripts/generated/2023-02-24_copying_files_in_paralelle/untar_sample_folders.launcher",
            row.names = F, col.names = F, quote = F)


# touch the instrain output ----

touch.IS.slurm <- c("#!/bin/bash",
                 "#SBATCH -J touchIS",
                 "#SBATCH -p development",
                 "#SBATCH -N 4",
                 "#SBATCH -n 472",
                 "#SBATCH -t 2:00:00",
                 "#SBATCH -o touch.o%j",
                 "#SBATCH -e touch.e%j",
                 "#SBATCH --mail-type=all",
                 "#SBATCH --mail-user=rohwer@utexas.edu",
                 "module load launcher",
                 "export LAUNCHER_JOB_FILE=touch_instrain_sample_folders.launcher",
                 "export LAUNCHER_SCHED=interleaved",
                 "${LAUNCHER_DIR}/paramrun")

touch.IS.launcher <- paste0("find ",sample.names.instrain.folders," -exec touch {} \\;")

write.table(x = touch.IS.slurm, file = "server-scripts/generated/2023-02-24_copying_files_in_paralelle/touch_instrain_sample_folders.slurm",
            row.names = F, col.names = F, quote = F)

write.table(x = touch.IS.launcher, file = "server-scripts/generated/2023-02-24_copying_files_in_paralelle/touch_instrain_sample_folders.launcher",
            row.names = F, col.names = F, quote = F)


