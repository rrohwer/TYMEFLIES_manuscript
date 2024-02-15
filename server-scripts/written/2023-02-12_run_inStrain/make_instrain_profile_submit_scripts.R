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
#   --min_read_ani 92             *** default for 95% ani, but an important choice for results ***
#   --skip_plot_generation        *** it seems like it broke on the plotting step ***

# So profiling runs once per each bam, and takes an hour or two and uses all the node's memory with 1

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

all.bins <- readRDS(file = "data/2022-11-15_bin_stats/all_bins.rds")
all.samples <- unique(all.bins$tymeflies.name)
bam.files <- paste0("../maptodrep/",all.samples,".bam")
drep.fasta <- "../maptodrep/drep_winners_concat.fna"
output.names <- paste0(all.samples,".IS")
genes.file <- "../runprodigal/drep_winners_concat_genes.fna"
stb.file <- "instrain_stb.txt"


slurm <- get.slum.specs(N = 60, simultaneous.tasks.per.node = 1, min.per.task = 180, tot.tasks = length(all.samples)) # memory limit- 1/node, time ~ 2hrs, give it 3 hrs each.

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
                            "conda activate instrain ; ",
                            "inStrain profile ",
                            bam.files," ",
                            drep.fasta,
                            " -o ",output.names,
                            " -g ",genes.file,
                            " -s ", stb.file,
                            " --database_mode --skip_mm --min_read_ani 92 -p 128 --skip_plot_generation")

write.table(x = instrain.slurm, file = "server-scripts/generated/2023-02-12_run_inStrain/run_instrain_profile.slurm",
            row.names = F, col.names = F, quote = F)
write.table(x = instrain.launcher, file = "server-scripts/generated/2023-02-12_run_inStrain/run_instrain_profile.launcher",
            row.names = F, col.names = F, quote = F)
