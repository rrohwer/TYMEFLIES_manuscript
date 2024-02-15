# RRR

# ---- set-up ----

x <- readRDS(file = "data/2022-09-14_list_of_all_bins/bin_key.rds")
head(x)

sample.names <- unique(x$sample)

generated.script.folder <- "server-scripts/generated/2022-09-15_run_checkm2/"

get.timestamp.from.minutes <- function(m){
  if (m >= 45 * 60){cat("\n\nOh No! This is gonna take too long for 1 job\n\n")}
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
  m <- ceiling(min.per.task * tasks.per.node) + ceiling(min.per.task * extra.tasks.per.node) # 20 min / task * n tasks / node
  t <- get.timestamp.from.minutes(m = m)
  cat("N = ",N," n = ",n, " total = ",tot.tasks,"\ntime = ",t," (",min.per.task," min each)\n",simultaneous.tasks.per.node," at a time on each node" , sep = "")
  return(list("n"=n,"N"=N,"t"=t,"threads" = threads))
}


# ---- make test version ----

# x <- x[x$sample %in% sample.names[1:5], ]
# sample.names <- sample.names[1:5]

# ---- get full scripts ----

slurm <- get.slum.specs(N = 20, simultaneous.tasks.per.node = 2, min.per.task = 10, tot.tasks = length(sample.names))

checkm2.slurm <- c("#!/bin/bash",
                  "#SBATCH -J checkm2",
                  "#SBATCH -p normal",
                  paste("#SBATCH -N", slurm$N),
                  paste("#SBATCH -n", slurm$n),
                  paste0("#SBATCH -t ",slurm$t),
                  "#SBATCH -o checkm.o%j",
                  "#SBATCH -e checkm.e%j",
                  "#SBATCH --mail-type=all",
                  "#SBATCH --mail-user=rohwer@utexas.edu",
                  "module load launcher",
                  # paste0("export OMP_NUM_THREADS=",slurm$threads),
                  "export LAUNCHER_JOB_FILE=checkm2.launcher",
                  "export LAUNCHER_SCHED=interleaved",
                  "${LAUNCHER_DIR}/paramrun")

checkm2.launcher <- NULL
for (s in sample.names){
  input.folder <- paste0("metabat2_bins/", s)
  output.folder <- paste0("checkm2_output/", s)
  checkm2.launcher <- c(checkm2.launcher,
                       paste0("mkdir ", output.folder, " ; ",
                              "module unload python3 ; ",
                              "source /work/08922/rrohwer/ls6/miniforge3/etc/profile.d/conda.sh ; ",
                              "conda activate checkm2 ; ",
                              "checkm2 predict --threads ", slurm$threads, " --input ", input.folder, " --output-directory ", output.folder))
}


write.table(x = checkm2.slurm, file = file.path(generated.script.folder, "checkm2.slurm"), quote = F, row.names = F, col.names = F, append = F, sep = "\n")
write.table(x = checkm2.launcher, file = file.path(generated.script.folder, "checkm2.launcher"), quote = F, row.names = F, col.names = F, append = F, sep = "\n")

