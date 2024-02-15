# RRR
# run checkM2 on all the genomes downloaded from GTDB
# before I fed it folders, and parallelized the bin folders
# but this time I don't have folders, so one launcher line per bin file
# It should be super fast

bin.names <- read.table("data/2023-04-16_GTDB_accessions/3a_gtdb_downloaded_bin_names.txt")
bin.names <- bin.names$V1
folder.names <- sub(pattern = "\\.fna\\.gz$", replacement = "", x = bin.names)

# ---- funcitons ----

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

get.slum.specs.single.threaded <- function(N, simultaneous.tasks.per.node, min.per.task, tot.tasks){
  
  tasks.per.node <- floor( tot.tasks / N )
  
  if (N == 0){
    N <- 1
  }
  
  n <- simultaneous.tasks.per.node * N
  
  extra.tasks <- tot.tasks - (tasks.per.node * N)
  if (extra.tasks < 0){
    extra.tasks <- 0
  }
  extra.tasks.per.node <- ceiling(extra.tasks / N)
  threads <- floor(128 / simultaneous.tasks.per.node)
  launcher.rounds <- tot.tasks / n
  m <- ceiling(min.per.task) * ceiling(launcher.rounds) 
  t <- get.timestamp.from.minutes(m = m)
  cat("N = ",N," n = ",n, " total = ",tot.tasks,"\ntime = ",t," (",min.per.task," min each)\n",simultaneous.tasks.per.node," at a time on each node" , sep = "")
  return(list("n"=n,"N"=N,"t"=t,"threads" = threads))
}


# ---- generate scripts- didn't work, just stalled out somehow ----

# slurm <- get.slum.specs.single.threaded(N = 6, simultaneous.tasks.per.node = 128, min.per.task = 600, tot.tasks = length(bin.names)) # damn it timed out at 30 minutes- how???
# 
# checkm2.slurm <- c("#!/bin/bash",
#                    "#SBATCH -J checkm2",
#                    "#SBATCH -p normal",
#                    paste("#SBATCH -N", slurm$N),
#                    paste("#SBATCH -n", slurm$n),
#                    paste0("#SBATCH -t ",slurm$t),
#                    "#SBATCH -o checkm.o%j",
#                    "#SBATCH -e checkm.e%j",
#                    "#SBATCH --mail-type=all",
#                    "#SBATCH --mail-user=rohwer@utexas.edu",
#                    "module load launcher",
#                    # paste0("export OMP_NUM_THREADS=",slurm$threads),
#                    "export LAUNCHER_JOB_FILE=checkm2.launcher",
#                    "export LAUNCHER_SCHED=interleaved",
#                    "${LAUNCHER_DIR}/paramrun")
# 
# checkm2.launcher <- paste0("mkdir output/",folder.names," ; source /work/08922/rrohwer/ls6/miniforge3/etc/profile.d/conda.sh ; conda activate checkm2 ; checkm2 predict --threads 1 --input ", bin.names, " --output-directory output/", folder.names)

# ---- try single node, multi-threaded instead ----

checkm2.slurm <- c("#!/bin/bash",
                  "#SBATCH -J checkm2",
                  "#SBATCH -p normal",
                  "#SBATCH -N 1",
                  "#SBATCH -n 1",
                  "#SBATCH -t 48:00:00",
                  "#SBATCH -o checkm.o%j",
                  "#SBATCH -e checkm.e%j",
                  "#SBATCH --mail-type=all",
                  "#SBATCH --mail-user=rohwer@utexas.edu",
                  "module load launcher",
                  "export LAUNCHER_JOB_FILE=checkm2.launcher",
                  "export LAUNCHER_SCHED=interleaved",
                  "${LAUNCHER_DIR}/paramrun")

checkm2.launcher <- "source /work/08922/rrohwer/ls6/miniforge3/etc/profile.d/conda.sh ; conda activate checkm2 ; checkm2 predict --threads 127 -x gz --input GTDB_refs --output-directory checkM2_on_GTDB_Actinomycetia"

# did 127 threads so could check on it with htop. Using all threads. but in waves. Max I saw was 20 GB
# Slurm Job_id=811481 Name=checkm2 Ended, Run time 00:09:38, COMPLETED, ExitCode 0
# wtf that was so fast. why didn't it finish before???

# ---- export ----

write.table(x = checkm2.slurm, file = "server-scripts/generated/2023-04-16_get_Actinomycetia_backbone/checkm2.slurm", quote = F, row.names = F, col.names = F, append = F, sep = "\n")
write.table(x = checkm2.launcher, file = "server-scripts/generated/2023-04-16_get_Actinomycetia_backbone/checkm2.launcher", quote = F, row.names = F, col.names = F, append = F, sep = "\n")
