# RRR
# make the final coverM table
# exclude 1000 bp from contig ends
# exclude top and bottom 5% coverage sections
# return all the metrics in long-format
# calculate all the possible abundance metrics to have on hand

# coverm --version  # coverm 0.6.1

# don't forget to add the genome reference file to the coverM run folder on scratch

# ---- set up ----

sample.names <- readRDS(file = "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")
sample.names <- unique(sample.names$tymeflies.name)

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

all.bams <- paste0(sample.names,".bam")
all.bams <- paste0("../maptodrep96/", all.bams)
all.bams <- paste0(all.bams, collapse = " ")

exclude.ends <- 1000
min.ANI <- c(0,93) # inStrain will use a 93% cutoff with drep at 96

slurm <- get.slum.specs(N = 1, simultaneous.tasks.per.node = 2, min.per.task = 600, tot.tasks = length(min.ANI)) # memory limit- 3/node, time ~ 300 min / line

coverm.slurm <- c("#!/bin/bash",
                  "#SBATCH -J coverm",
                  "#SBATCH -p normal",
                  paste("#SBATCH -N", slurm$N),
                  paste("#SBATCH -n", slurm$n),
                  paste("#SBATCH -t", slurm$t),
                  "#SBATCH -o coverm.o%j",
                  "#SBATCH -e coverm.e%j",
                  "#SBATCH --mail-type=all",
                  "#SBATCH --mail-user=rohwer@utexas.edu",
                  # "module load biocontainers",
                  "module load launcher",
                  "export LAUNCHER_JOB_FILE=get_coverm_stats.launcher",
                  "export LAUNCHER_SCHED=interleaved",
                  "${LAUNCHER_DIR}/paramrun")


coverm.launcher <- paste0("source /work/08922/rrohwer/ls6/miniforge3/etc/profile.d/conda.sh ; ",
                          "conda activate coverm ; ",
                          "coverm genome --bam-files ",all.bams,
                          " --genome-definition coverM_genome_reference.txt ",
                          "--min-read-percent-identity-pair ",min.ANI," ",
                          "--proper-pairs-only ",
                          "--methods relative_abundance mean rpkm tpm ",
                          "--contig-end-exclusion ",exclude.ends," ",
                          "--output-file coverM_drep96_minANI",min.ANI,".txt ",
                          "--output-format sparse ",
                          "--threads ",slurm$threads," ",
                          "&> coverM_drep96_terminal_output_minANI",min.ANI,".txt") # save terminal output



# ---- export ----

write.table(x = coverm.launcher, file = "server-scripts/generated/2023-03-15_run_coverM_on_drep96_winners/get_coverm_stats.launcher", 
            quote = F, row.names = F, col.names = F)

write.table(x = coverm.slurm, file = "server-scripts/generated/2023-03-15_run_coverM_on_drep96_winners/get_coverm_stats.slurm", 
            quote = F, row.names = F, col.names = F)

# Slurm Job_id=755210 Name=coverm Ended, Run time 14:42:47, COMPLETED, ExitCode 0