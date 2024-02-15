# RRR

# inStrain profile
#   sample.bam
#   drepbins.fasta
#   -o sample.IS
#   -p 128                        *** fewer threads saves RAM ***
#   -g drep.genes.fna
#   -s instrain_stb.txt           *** the scaffold name \t genome name file ***
#   --database_mode
#   --skip_mm                     *** I think default with --database_mode anyway ***
#   --min_read_ani 93             *** recommends 3 less than used for drep, and used drep 96 ***
#   --skip_plot_generation        *** it seems like it broke on the plotting step ***

# So profiling runs once per each bam, and takes an hour or two and uses all the node's memory with 1
# running out of memory is an issue, can use less memory if use fewer processors.
# -p 128 reached memory limit on largest bam test file (14GB)
# -p 112 ran out of memory
# -p 96 worked, took 03:30:39
# -p 64 worked, took 03:37:59
# -p 18 worked, took 05:52:02 (but with plots it timed out at 10 hrs)

# ---- install inStrain ----

# conda create --name instrain_py310 python=3.10  # in Feb 25, 2023 github response, Olm says updated to python 3.10
# pip install instrain                            # Downloading inStrain-1.7.1.tar.gz
# conda install -c bioconda prodigal              # bioconda/linux-64::prodigal-2.6.3-hec16e2b_4
# conda install -c bioconda samtools              # bioconda/linux-64::samtools-1.6-h3f2fef4_8

# ---- functions ----

get.timestamp.from.minutes <- function(m){
  if (m > 48 * 60){cat("\n\nOh No! This is gonna take too long for 1 job\n\n")}
  h <- 0
  if (m >= 60){
    h <- floor(m / 60)
    m <- m - (h * 60)
  }
  t <- paste0(h,":",m,":00")
  return(t)
}

get.slum.specs <- function(N, simultaneous.tasks.per.node, min.per.task, tot.tasks){
  
  tasks.per.node <- floor( tot.tasks / N )
  if (tasks.per.node == 0){
    tasks.per.node <- 1
  }
  
  extra.tasks <- tot.tasks - (tasks.per.node * N)
  if (extra.tasks < 0){
    extra.tasks <- 0
  }
  extra.tasks.per.node <- ceiling(extra.tasks / N)
  
  n <- simultaneous.tasks.per.node * N
  
  threads <- floor(128 / simultaneous.tasks.per.node)
  
  m <- ceiling(min.per.task * tasks.per.node) + ceiling(min.per.task * extra.tasks.per.node) # 20 min / task * n tasks / node
  t <- get.timestamp.from.minutes(m = m)
  
  cat("N = ",N," n = ",n, " total = ",tot.tasks,"\ntime = ",t," (",min.per.task," min each)\n",simultaneous.tasks.per.node," at a time on each node" , sep = "")
  return(list("n"=n,"N"=N,"t"=t,"threads" = threads))
}

# ---- make scripts ----

all.bins <- readRDS(file = "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")
all.samples <- unique(all.bins$tymeflies.name)
bam.files <- paste0("../maptodrep96/",all.samples,".bam")
drep.fasta <- "../maptodrep96/drep_96_winners_concat.fna"
output.names <- paste0(all.samples,".IS")
genes.file <- "../runprodigal96/drep_96_winners_concat_genes.fna"
stb.file <- "instrain_stb.txt"


slurm <- get.slum.specs(N = 60, simultaneous.tasks.per.node = 1, min.per.task = 285, tot.tasks = length(all.samples)) # memory limit- 1/node, time ~ 3:30hrs, give it 4:45 hrs each.

instrain.slurm <- c("#!/bin/bash",
                  "#SBATCH -J instrain",
                  "#SBATCH -p normal",
                  paste("#SBATCH -N", slurm$N),
                  paste("#SBATCH -n", slurm$n),
                  paste("#SBATCH -t", slurm$t),
                  "#SBATCH -o instrain.o%j",
                  "#SBATCH -e instrain.e%j",
                  "#SBATCH --mail-type=all",
                  "#SBATCH --mail-user=rohwer@utexas.edu",
                  # "module load biocontainers",
                  "module load launcher",
                  "export LAUNCHER_JOB_FILE=run_instrain_profile.launcher",
                  "export LAUNCHER_SCHED=interleaved",
                  "${LAUNCHER_DIR}/paramrun")

instrain.launcher <- paste0("source /work/08922/rrohwer/ls6/miniforge3/etc/profile.d/conda.sh ; ",
                            "conda activate instrain_py310 ; ",
                            "inStrain profile ",
                            bam.files," ",
                            drep.fasta,
                            " -o ",output.names,
                            " -g ",genes.file,
                            " -s ", stb.file,
                            " --database_mode --skip_mm --min_read_ani 93 -p 96 --skip_plot_generation")

write.table(x = instrain.slurm, file = "server-scripts/generated/2023-03-13_run_inStrain_on_drep96_winners/run_instrain_profile.slurm",
            row.names = F, col.names = F, quote = F)
write.table(x = instrain.launcher, file = "server-scripts/generated/2023-03-13_run_inStrain_on_drep96_winners/run_instrain_profile.launcher",
            row.names = F, col.names = F, quote = F)

# Slurm Job_id=751529 Name=instrain Failed, Run time 1-14:00:16, TIMEOUT, ExitCode 0

# ---- find samples finished running ----

# cat instrain.o751529 | grep Job > Job_completed_list.txt
# cat instrain.o751529 | grep job > job_commands_list.txt

completed <- read.table("data/2023-03-13_inStrain_on_drep96/errors/Job_completed_list.txt", sep = " ")
commands <- read.table("data/2023-03-13_inStrain_on_drep96/errors/job_commands_list.txt", sep = " ", skip = 1)

head(completed)
completed <- completed[ ,c(3,6,4)]
colnames(completed) <- c("Job","seconds","status")
unique(completed$status) # these 395 all finished

head(commands)
commands <- commands[ ,c(3,6,21)]
colnames(commands) <- c("Task","Job","Sample")
commands$Sample <- sub("\\.IS", "",commands$Sample)
nrow(commands) # these 410 were launched
410 - 395 # so 15 were running when it timed out and didn't finish ***why not 60 ??***
471 - 395 # and 76 are left to run

# find those 76 to re-run:

completed <- merge(x = completed, y = commands, by = "Job", all.x = TRUE, all.y = FALSE)
head(completed)

unfinished <- setdiff(x = all.samples, y = completed$Sample) # these are all to re-run
partial.run <- setdiff(x = commands$Sample, y = completed$Sample) # these have folders to delete

# ---- delete partial-run output folders ----

delete.script <- c("#!/bin/bash",
                   paste0("rm -Rf ",partial.run,".IS"))

write.table(x = delete.script, file = "server-scripts/generated/2023-03-13_run_inStrain_on_drep96_winners/delete_partial-run_folders.sh",
            row.names = F, col.names = F, quote = F)

# ---- make new submit scripts ----

all.samples <- unfinished
bam.files <- paste0("../maptodrep96/",all.samples,".bam")
drep.fasta <- "../maptodrep96/drep_96_winners_concat.fna"
output.names <- paste0(all.samples,".IS")
genes.file <- "../runprodigal96/drep_96_winners_concat_genes.fna"
stb.file <- "instrain_stb.txt"


slurm <- get.slum.specs(N = 60, simultaneous.tasks.per.node = 1, min.per.task = 900, tot.tasks = length(all.samples)) # memory limit- 1/node, longest time was ~ 9:30hrs, give it 15 hrs each.

instrain.slurm <- c("#!/bin/bash",
                    "#SBATCH -J unfinishedinstrain",
                    "#SBATCH -p normal",
                    paste("#SBATCH -N", slurm$N),
                    paste("#SBATCH -n", slurm$n),
                    paste("#SBATCH -t", slurm$t),
                    "#SBATCH -o instrain_repeat_unfinished.o%j",
                    "#SBATCH -e instrain_repeat_unfinished.e%j",
                    "#SBATCH --mail-type=all",
                    "#SBATCH --mail-user=rohwer@utexas.edu",
                    # "module load biocontainers",
                    "module load launcher",
                    "export LAUNCHER_JOB_FILE=run_instrain_profile_repeat_unfinished.launcher",
                    "export LAUNCHER_SCHED=interleaved",
                    "${LAUNCHER_DIR}/paramrun")

instrain.launcher <- paste0("source /work/08922/rrohwer/ls6/miniforge3/etc/profile.d/conda.sh ; ",
                            "conda activate instrain_py310 ; ",
                            "inStrain profile ",
                            bam.files," ",
                            drep.fasta,
                            " -o ",output.names,
                            " -g ",genes.file,
                            " -s ", stb.file,
                            " --database_mode --skip_mm --min_read_ani 93 -p 96 --skip_plot_generation")

write.table(x = instrain.slurm, file = "server-scripts/generated/2023-03-13_run_inStrain_on_drep96_winners/run_instrain_profile_repeat_unfinished.slurm",
            row.names = F, col.names = F, quote = F)
write.table(x = instrain.launcher, file = "server-scripts/generated/2023-03-13_run_inStrain_on_drep96_winners/run_instrain_profile_repeat_unfinished.launcher",
            row.names = F, col.names = F, quote = F)

# Slurm Job_id=759701 Name=unfinishedinstrain Failed, Run time 1-06:00:12, TIMEOUT, ExitCode 0


# ---- find samples that finished running from the first re-run ----

# cat instrain_repeat_unfinished.o759701 | grep Job > Job_complete_list_unfinished.txt
# cat instrain_repeat_unfinished.o759701 | grep job > job_commands_list_unfinished.txt

completed <- read.table("data/2023-03-13_inStrain_on_drep96/errors_unfinished/Job_complete_list_unfinished.txt", sep = " ")
commands <- read.table("data/2023-03-13_inStrain_on_drep96/errors_unfinished/job_commands_list_unfinished.txt", sep = " ", skip = 1)

head(completed)
completed <- completed[ ,c(3,6,4)]
colnames(completed) <- c("Job","seconds","status")
unique(completed$status) # these 57 all finished

head(commands)
commands <- commands[ ,c(3,6,21)]
colnames(commands) <- c("Task","Job","Sample")
commands$Sample <- sub("\\.IS", "",commands$Sample)
nrow(commands) # these 67 were launched
67 - 57 # so 10 were running when it timed out and didn't finish 
76 - 57 # and 19 are left to run

# find those 19 to re-run:

completed <- merge(x = completed, y = commands, by = "Job", all.x = TRUE, all.y = FALSE)
head(completed)

unfinished <- setdiff(x = all.samples, y = completed$Sample) # these are all to re-run
partial.run <- setdiff(x = commands$Sample, y = completed$Sample) # these have folders to delete

# I noticed that there are a lot of repeats between what didn't run this time and what didn't run last time
# Maybe they are running out of RAM
# decrease the processors this time too from 96 to 64

# ---- delete partial-run output folders ----

delete.script <- c("#!/bin/bash",
                   paste0("rm -Rf ",partial.run,".IS"))

write.table(x = delete.script, file = "server-scripts/generated/2023-03-13_run_inStrain_on_drep96_winners/delete_partial-run_folders_unfinished.sh",
            row.names = F, col.names = F, quote = F)

# ---- make new submit scripts ----

all.samples <- unfinished
bam.files <- paste0("../maptodrep96/",all.samples,".bam")
drep.fasta <- "../maptodrep96/drep_96_winners_concat.fna"
output.names <- paste0(all.samples,".IS")
genes.file <- "../runprodigal96/drep_96_winners_concat_genes.fna"
stb.file <- "instrain_stb.txt"


slurm <- get.slum.specs(N = 19, simultaneous.tasks.per.node = 1, min.per.task = 1500, tot.tasks = length(all.samples)) # memory limit- 1/node, longest time was ~ 9:30hrs, give it 25 hrs each I just want to be done with this mess.

instrain.slurm <- c("#!/bin/bash",
                    "#SBATCH -J unfinishedinstrain2",
                    "#SBATCH -p normal",
                    paste("#SBATCH -N", slurm$N),
                    paste("#SBATCH -n", slurm$n),
                    paste("#SBATCH -t", slurm$t),
                    "#SBATCH -o instrain_repeat_unfinished_2.o%j",
                    "#SBATCH -e instrain_repeat_unfinished_2.e%j",
                    "#SBATCH --mail-type=all",
                    "#SBATCH --mail-user=rohwer@utexas.edu",
                    # "module load biocontainers",
                    "module load launcher",
                    "export LAUNCHER_JOB_FILE=run_instrain_profile_repeat_unfinished_2.launcher",
                    "export LAUNCHER_SCHED=interleaved",
                    "${LAUNCHER_DIR}/paramrun")

instrain.launcher <- paste0("source /work/08922/rrohwer/ls6/miniforge3/etc/profile.d/conda.sh ; ",
                            "conda activate instrain_py310 ; ",
                            "inStrain profile ",
                            bam.files," ",
                            drep.fasta,
                            " -o ",output.names,
                            " -g ",genes.file,
                            " -s ", stb.file,
                            " --database_mode --skip_mm --min_read_ani 93 -p 64 --skip_plot_generation")

write.table(x = instrain.slurm, file = "server-scripts/generated/2023-03-13_run_inStrain_on_drep96_winners/run_instrain_profile_repeat_unfinished_2.slurm",
            row.names = F, col.names = F, quote = F)
write.table(x = instrain.launcher, file = "server-scripts/generated/2023-03-13_run_inStrain_on_drep96_winners/run_instrain_profile_repeat_unfinished_2.launcher",
            row.names = F, col.names = F, quote = F)

# Slurm Job_id=764406 Name=unfinishedinstrain2 Ended, Run time 07:41:12, COMPLETED, ExitCode 0
