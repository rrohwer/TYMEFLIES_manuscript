# RRR

# ---- set-up ----

x <- readRDS(file = "data/2022-09-14_list_of_all_bins/bin_key.rds")
head(x)

sample.names <- unique(x$sample)

generated.script.folder <- "server-scripts/generated/2022-09-22_run_gtdb-tk/"

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

slurm <- get.slum.specs(N = 20, simultaneous.tasks.per.node = 1, min.per.task = 120, tot.tasks = length(sample.names))

gtdbtk.slurm <- c("#!/bin/bash",
                  "#SBATCH -J gtdbtk",
                  "#SBATCH -p normal",
                  paste("#SBATCH -N", slurm$N),
                  paste("#SBATCH -n", slurm$n),
                  paste0("#SBATCH -t ",slurm$t),
                  "#SBATCH -o gtdbtk.o%j",
                  "#SBATCH -e gtdbtk.e%j",
                  "#SBATCH --mail-type=all",
                  "#SBATCH --mail-user=rohwer@utexas.edu",
                  "module load launcher",
                  "export LAUNCHER_JOB_FILE=gtdbtk.launcher",
                  "export LAUNCHER_SCHED=interleaved",
                  "${LAUNCHER_DIR}/paramrun")

gtdbtk.launcher <- NULL
for (s in sample.names){
  input.folder <- paste0("metabat2_bins/", s)
  output.folder <- paste0("gtdbtk_output/", s)
  gtdbtk.launcher <- c(gtdbtk.launcher,
                       paste0("mkdir ", output.folder, " ; ",
                              "source /work/08922/rrohwer/ls6/miniforge3/etc/profile.d/conda.sh ; ",
                              "conda activate gtdbtk-2.1.1 ; ",
                              "gtdbtk classify_wf --genome_dir ", input.folder, " --out_dir ", output.folder, " --cpus 15 --pplacer_cpus 3 "))
}


write.table(x = gtdbtk.slurm, file = file.path(generated.script.folder, "gtdbtk.slurm"), quote = F, row.names = F, col.names = F, append = F, sep = "\n")
write.table(x = gtdbtk.launcher, file = file.path(generated.script.folder, "gtdbtk.launcher"), quote = F, row.names = F, col.names = F, append = F, sep = "\n")

