# RRR

# eggnog is in a conda environment called eggnog
# the database is in a scratch folder called eggnog_database
# the protein files are in runprodigal96_take2

# $ emapper.py --version --data_dir /scratch/08922/rrohwer/eggnog_databases/
# emapper-2.1.10 / Expected eggNOG DB version: 5.0.2 / Installed eggNOG DB version: 5.0.2 / Diamond version found: diamond version 2.1.6 / MMseqs2 version found: 13.45111

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

# ---- test one at a time ----

test.bins <- all.bins[1]

slurm <- get.slum.specs(N = 1, simultaneous.tasks.per.node = 1, min.per.task = 400, tot.tasks = length(test.bins)) # memory limit- 3/node, time ~ 300 min / line

slurm.file <- c("#!/bin/bash",
                "#SBATCH -J testeggnog1",
                "#SBATCH -p normal",
                paste("#SBATCH -N", slurm$N),
                paste("#SBATCH -n", slurm$n),
                paste("#SBATCH -t", slurm$t),
                "#SBATCH -o eggnog.o%j",
                "#SBATCH -e eggnog.e%j",
                "#SBATCH --mail-type=all",
                "#SBATCH --mail-user=rohwer@utexas.edu",
                "module load launcher",
                "export LAUNCHER_JOB_FILE=testeggnog1.launcher",
                "export LAUNCHER_SCHED=interleaved",
                "${LAUNCHER_DIR}/paramrun")

launcher.file <- paste0("source /work/08922/rrohwer/ls6/miniforge3/etc/profile.d/conda.sh ; ",
                        "conda activate eggnog ; ",
                        "emapper.py -i ../runprodigal96_take2/",test.bins,".genes.faa -o ", test.bins, " --cpu ",slurm$threads," --data_dir ../eggnog_databases --dbmem")

write.table(x = launcher.file, file = "server-scripts/generated/2023-04-30_run_eggnog_on_drep96_winners/testeggnog1.launcher", 
            quote = F, row.names = F, col.names = F)

write.table(x = slurm.file, file = "server-scripts/generated/2023-04-30_run_eggnog_on_drep96_winners/testeggnog1.slurm", 
            quote = F, row.names = F, col.names = F)


# ---- test 2 at a time ----

test.bins <- all.bins[1:2]

slurm <- get.slum.specs(N = 1, simultaneous.tasks.per.node = 2, min.per.task = 400, tot.tasks = length(test.bins)) # memory limit- 3/node, time ~ 300 min / line

slurm.file <- c("#!/bin/bash",
                "#SBATCH -J testeggnog2",
                "#SBATCH -p normal",
                paste("#SBATCH -N", slurm$N),
                paste("#SBATCH -n", slurm$n),
                paste("#SBATCH -t", slurm$t),
                "#SBATCH -o eggnog.o%j",
                "#SBATCH -e eggnog.e%j",
                "#SBATCH --mail-type=all",
                "#SBATCH --mail-user=rohwer@utexas.edu",
                "module load launcher",
                "export LAUNCHER_JOB_FILE=testeggnog2.launcher",
                "export LAUNCHER_SCHED=interleaved",
                "${LAUNCHER_DIR}/paramrun")

launcher.file <- paste0("source /work/08922/rrohwer/ls6/miniforge3/etc/profile.d/conda.sh ; ",
                        "conda activate eggnog ; ",
                        "emapper.py -i ../runprodigal96_take2/",test.bins,".genes.faa -o ", test.bins, " --cpu ",slurm$threads," --data_dir ../eggnog_databases --dbmem")

write.table(x = launcher.file, file = "server-scripts/generated/2023-04-30_run_eggnog_on_drep96_winners/testeggnog2.launcher", 
            quote = F, row.names = F, col.names = F)

write.table(x = slurm.file, file = "server-scripts/generated/2023-04-30_run_eggnog_on_drep96_winners/testeggnog2.slurm", 
            quote = F, row.names = F, col.names = F)

# ---- test 4 at a time ----

test.bins <- all.bins[1:4]

slurm <- get.slum.specs(N = 1, simultaneous.tasks.per.node = 4, min.per.task = 400, tot.tasks = length(test.bins)) # memory limit- 3/node, time ~ 300 min / line

slurm.file <- c("#!/bin/bash",
                "#SBATCH -J testeggnog4",
                "#SBATCH -p normal",
                paste("#SBATCH -N", slurm$N),
                paste("#SBATCH -n", slurm$n),
                paste("#SBATCH -t", slurm$t),
                "#SBATCH -o eggnog.o%j",
                "#SBATCH -e eggnog.e%j",
                "#SBATCH --mail-type=all",
                "#SBATCH --mail-user=rohwer@utexas.edu",
                "module load launcher",
                "export LAUNCHER_JOB_FILE=testeggnog4.launcher",
                "export LAUNCHER_SCHED=interleaved",
                "${LAUNCHER_DIR}/paramrun")

launcher.file <- paste0("source /work/08922/rrohwer/ls6/miniforge3/etc/profile.d/conda.sh ; ",
                        "conda activate eggnog ; ",
                        "emapper.py -i ../runprodigal96_take2/",test.bins,".genes.faa -o ", test.bins, " --cpu ",slurm$threads," --data_dir ../eggnog_databases --dbmem")

write.table(x = launcher.file, file = "server-scripts/generated/2023-04-30_run_eggnog_on_drep96_winners/testeggnog4.launcher", 
            quote = F, row.names = F, col.names = F)

write.table(x = slurm.file, file = "server-scripts/generated/2023-04-30_run_eggnog_on_drep96_winners/testeggnog4.slurm", 
            quote = F, row.names = F, col.names = F)

# ---- results ----

# Slurm Job_id=832649 Name=testeggnog1 Ended, Run time 00:07:57, COMPLETED, ExitCode 0
# Slurm Job_id=832650 Name=testeggnog2 Ended, Run time 00:08:41, COMPLETED, ExitCode 0
# Slurm Job_id=832651 Name=testeggnog4 Ended, Run time 00:19:03, COMPLETED, ExitCode 0

((57 / 60) + 7 ) / 1 # 7.95 per job with 1/node
((41 / 60) + 8 ) / 2 # 4.34 per job with 2/node
((3 / 60) + 19 ) / 4 # 4.76 per job with 4/node

# so run it with 2/node, expecting it to take 4:20 minutes each

((2855 / 2 ) * 4.5) / 60 # 107 hours total
107 / 20 # 5 hours on 20 nodes!





