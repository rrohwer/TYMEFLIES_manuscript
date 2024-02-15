# RRR
# installed java, but need to add it to the path at the start of each compute node
# export PATH="/work/08922/rrohwer/ls6/jdk-17.0.6+10/bin:$PATH"
# export PATH="/work/08922/rrohwer/ls6/my_interproscan/interproscan-5.61-93.0:$PATH"
# command options:
# interproscan.sh 
# -i /scratch/08922/rrohwer/runprodigal96_take2/BIN-NAME.genes.faa
# -f tsv gff3
# -goterms also lookup GO terms (needs lookup)
# -cpu
# -b BIN-NAME (basename of output files)
# -dp (disable pre-calc. try without doing this first?)
# -pa also lookup pathways (needs lookup)

bins <- readRDS("data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")

bins <- bins[bins$winner, ]
head(bins)

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

# ---- test runs ----

bins <- bins[1:4, ]

slurm <- get.slum.specs(N = 1, simultaneous.tasks.per.node = 4, min.per.task = 400, tot.tasks = length(bins$bin.full.name)) # memory limit- 3/node, time ~ 300 min / line

slurm.file <- c("#!/bin/bash",
                "#SBATCH -J test4Nlookup",
                "#SBATCH -p normal",
                paste("#SBATCH -N", slurm$N),
                paste("#SBATCH -n", slurm$n),
                paste("#SBATCH -t", slurm$t),
                "#SBATCH -o interproscan.o%j",
                "#SBATCH -e interproscan.e%j",
                "#SBATCH --mail-type=all",
                "#SBATCH --mail-user=rohwer@utexas.edu",
                "module load launcher",
                "export LAUNCHER_JOB_FILE=test4Nlookup.launcher",
                "export LAUNCHER_SCHED=interleaved",
                "${LAUNCHER_DIR}/paramrun")

launcher.file <- paste0("export PATH=\"/work/08922/rrohwer/ls6/jdk-17.0.6+10/bin:$PATH\" ; export PATH=\"/work/08922/rrohwer/ls6/my_interproscan/interproscan-5.61-93.0:$PATH\" ; ",
                        "interproscan.sh -i ","/scratch/08922/rrohwer/runprodigal96_take2/",bins$bin.full.name,".genes.faa",
                        " -f tsv gff3 -goterms -cpu ",127," -b ",bins$bin.full.name," -pa")

write.table(x = slurm.file, file = "server-scripts/generated/2023-04-19_run_interproscan_on_drep96_winners/test4Nlookup.slurm", quote = F, row.names = F, col.names = F)
write.table(x = launcher.file, file = "server-scripts/generated/2023-04-19_run_interproscan_on_drep96_winners/test4Nlookup.launcher", quote = F, row.names = F, col.names = F)

slurm.file <- c("#!/bin/bash",
                "#SBATCH -J test4Ndp",
                "#SBATCH -p normal",
                paste("#SBATCH -N", slurm$N),
                paste("#SBATCH -n", slurm$n),
                paste("#SBATCH -t", slurm$t),
                "#SBATCH -o interproscan.o%j",
                "#SBATCH -e interproscan.e%j",
                "#SBATCH --mail-type=all",
                "#SBATCH --mail-user=rohwer@utexas.edu",
                "module load launcher",
                "export LAUNCHER_JOB_FILE=test4Ndp.launcher",
                "export LAUNCHER_SCHED=interleaved",
                "${LAUNCHER_DIR}/paramrun")

launcher.file <- paste0("export PATH=\"/work/08922/rrohwer/ls6/jdk-17.0.6+10/bin:$PATH\" ; export PATH=\"/work/08922/rrohwer/ls6/my_interproscan/interproscan-5.61-93.0:$PATH\" ; ",
                        "interproscan.sh -i ","/scratch/08922/rrohwer/runprodigal96_take2/",bins$bin.full.name,".genes.faa",
                        " -f tsv gff3 -cpu ",127," -b ",bins$bin.full.name," -dp")

write.table(x = slurm.file, file = "server-scripts/generated/2023-04-19_run_interproscan_on_drep96_winners/test4Ndp.slurm", quote = F, row.names = F, col.names = F)
write.table(x = launcher.file, file = "server-scripts/generated/2023-04-19_run_interproscan_on_drep96_winners/test4Ndp.launcher", quote = F, row.names = F, col.names = F)

bins <- bins[1, ]

slurm <- get.slum.specs(N = 1, simultaneous.tasks.per.node = 1, min.per.task = 300, tot.tasks = length(bins$bin.full.name)) # memory limit- 3/node, time ~ 300 min / line

slurm.file <- c("#!/bin/bash",
                "#SBATCH -J test1Nlookup",
                "#SBATCH -p normal",
                paste("#SBATCH -N", slurm$N),
                paste("#SBATCH -n", slurm$n),
                paste("#SBATCH -t", slurm$t),
                "#SBATCH -o interproscan.o%j",
                "#SBATCH -e interproscan.e%j",
                "#SBATCH --mail-type=all",
                "#SBATCH --mail-user=rohwer@utexas.edu",
                "module load launcher",
                "export LAUNCHER_JOB_FILE=test1Nlookup.launcher",
                "export LAUNCHER_SCHED=interleaved",
                "${LAUNCHER_DIR}/paramrun")

launcher.file <- paste0("export PATH=\"/work/08922/rrohwer/ls6/jdk-17.0.6+10/bin:$PATH\" ; export PATH=\"/work/08922/rrohwer/ls6/my_interproscan/interproscan-5.61-93.0:$PATH\" ; ",
                        "interproscan.sh -i ","/scratch/08922/rrohwer/runprodigal96_take2/",bins$bin.full.name,".genes.faa",
                        " -f tsv gff3 -goterms -cpu ",127," -b ",bins$bin.full.name," -pa")

write.table(x = slurm.file, file = "server-scripts/generated/2023-04-19_run_interproscan_on_drep96_winners/test1Nlookup.slurm", quote = F, row.names = F, col.names = F)
write.table(x = launcher.file, file = "server-scripts/generated/2023-04-19_run_interproscan_on_drep96_winners/test1Nlookup.launcher", quote = F, row.names = F, col.names = F)

slurm.file <- c("#!/bin/bash",
                "#SBATCH -J test1Ndp",
                "#SBATCH -p normal",
                paste("#SBATCH -N", slurm$N),
                paste("#SBATCH -n", slurm$n),
                paste("#SBATCH -t", slurm$t),
                "#SBATCH -o interproscan.o%j",
                "#SBATCH -e interproscan.e%j",
                "#SBATCH --mail-type=all",
                "#SBATCH --mail-user=rohwer@utexas.edu",
                "module load launcher",
                "export LAUNCHER_JOB_FILE=test1Ndp.launcher",
                "export LAUNCHER_SCHED=interleaved",
                "${LAUNCHER_DIR}/paramrun")

launcher.file <- paste0("export PATH=\"/work/08922/rrohwer/ls6/jdk-17.0.6+10/bin:$PATH\" ; export PATH=\"/work/08922/rrohwer/ls6/my_interproscan/interproscan-5.61-93.0:$PATH\" ; ",
                        "interproscan.sh -i ","/scratch/08922/rrohwer/runprodigal96_take2/",bins$bin.full.name,".genes.faa",
                        " -f tsv gff3 -cpu ",127," -b ",bins$bin.full.name," -dp")

write.table(x = slurm.file, file = "server-scripts/generated/2023-04-19_run_interproscan_on_drep96_winners/test1Ndp.slurm", quote = F, row.names = F, col.names = F)
write.table(x = launcher.file, file = "server-scripts/generated/2023-04-19_run_interproscan_on_drep96_winners/test1Ndp.launcher", quote = F, row.names = F, col.names = F)

slurm.file <- c("#!/bin/bash",
                "#SBATCH -J test1NdpGO",
                "#SBATCH -p normal",
                paste("#SBATCH -N", slurm$N),
                paste("#SBATCH -n", slurm$n),
                paste("#SBATCH -t", slurm$t),
                "#SBATCH -o interproscan.o%j",
                "#SBATCH -e interproscan.e%j",
                "#SBATCH --mail-type=all",
                "#SBATCH --mail-user=rohwer@utexas.edu",
                "module load launcher",
                "export LAUNCHER_JOB_FILE=test1NdpGO.launcher",
                "export LAUNCHER_SCHED=interleaved",
                "${LAUNCHER_DIR}/paramrun")

launcher.file <- paste0("export PATH=\"/work/08922/rrohwer/ls6/jdk-17.0.6+10/bin:$PATH\" ; export PATH=\"/work/08922/rrohwer/ls6/my_interproscan/interproscan-5.61-93.0:$PATH\" ; ",
                        "interproscan.sh -i ","/scratch/08922/rrohwer/runprodigal96_take2/",bins$bin.full.name,".genes.faa",
                        " -f tsv gff3 -cpu ",127," -b ",bins$bin.full.name," -dp -iprlookup -goterms")

write.table(x = slurm.file, file = "server-scripts/generated/2023-04-19_run_interproscan_on_drep96_winners/test1NdpGO.slurm", quote = F, row.names = F, col.names = F)
write.table(x = launcher.file, file = "server-scripts/generated/2023-04-19_run_interproscan_on_drep96_winners/test1NdpGO.launcher", quote = F, row.names = F, col.names = F)



# ---- timing results ----

# Slurm Job_id=816310 Name=test4Nlookup Ended, Run time 05:04:23, COMPLETED, ExitCode 0
# Slurm Job_id=816315 Name=test4Ndp Ended, Run time 04:51:13, COMPLETED, ExitCode 0
# Slurm Job_id=816320 Name=test1Nlookup Ended, Run time 01:27:58, COMPLETED, ExitCode 0
# Slurm Job_id=816324 Name=test1Ndp Ended, Run time 01:22:11, COMPLETED, ExitCode 0
# Slurm Job_id=826817 Name=test1NdpGO Ended, Run time 01:27:51, COMPLETED, ExitCode 0






