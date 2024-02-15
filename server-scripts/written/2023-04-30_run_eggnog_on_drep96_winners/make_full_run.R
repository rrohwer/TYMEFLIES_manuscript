# RRR

# notes are in the test runs script

all.bins <- readRDS("data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")
all.bins <- all.bins[all.bins$winner, "bin.full.name"]

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
  
  threads <- floor(127 / simultaneous.tasks.per.node) # made 127 so I can ssh in to htop
  
  m <- ceiling(min.per.task * tasks.per.node) + ceiling(min.per.task * extra.tasks.per.node) # 20 min / task * n tasks / node
  t <- get.timestamp.from.minutes(m = m)
  
  cat("N = ",N," n = ",n, " total = ",tot.tasks,"\ntime = ",t," (",min.per.task," min each)\n",simultaneous.tasks.per.node," at a time on each node" , sep = "")
  return(list("n"=n,"N"=N,"t"=t,"threads" = threads))
}

# ---- full run ----

slurm <- get.slum.specs(N = 20, simultaneous.tasks.per.node = 2, min.per.task = 10, tot.tasks = length(all.bins)) 

slurm.file <- c("#!/bin/bash",
                "#SBATCH -J eggnog",
                "#SBATCH -p normal",
                paste("#SBATCH -N", slurm$N),
                paste("#SBATCH -n", slurm$n),
                paste("#SBATCH -t", slurm$t),
                "#SBATCH -o eggnog.o%j",
                "#SBATCH -e eggnog.e%j",
                "#SBATCH --mail-type=all",
                "#SBATCH --mail-user=rohwer@utexas.edu",
                "module load launcher",
                "export LAUNCHER_JOB_FILE=run_eggnog.launcher",
                "export LAUNCHER_SCHED=interleaved",
                "${LAUNCHER_DIR}/paramrun")

launcher.file <- paste0("source /work/08922/rrohwer/ls6/miniforge3/etc/profile.d/conda.sh ; ",
                        "conda activate eggnog ; ",
                        "emapper.py -i ../runprodigal96_take2/",all.bins,".genes.faa -o ", all.bins, " --cpu ",slurm$threads," --data_dir ../eggnog_databases --dbmem")

write.table(x = launcher.file, file = "server-scripts/generated/2023-04-30_run_eggnog_on_drep96_winners/run_eggnog.launcher", 
            quote = F, row.names = F, col.names = F)

write.table(x = slurm.file, file = "server-scripts/generated/2023-04-30_run_eggnog_on_drep96_winners/run_eggnog.slurm", 
            quote = F, row.names = F, col.names = F)

# Slurm Job_id=833229 Name=eggnog Ended, Run time 20:28:54, COMPLETED, ExitCode 0

