# RRR
# installed java, but need to add it to the path at the start of each compute node
# same with interproscan- it's not a conda environment, just a shell script to source
# export PATH="/work/08922/rrohwer/ls6/jdk-17.0.6+10/bin:$PATH"
# export PATH="/work/08922/rrohwer/ls6/interproscan-5.63-95.0:$PATH"

# command options:
# interproscan.sh 
#   -i input.faa 	(with no asterisks)
#   -f tsv gff3 	(output a tsv format and a gff3 format)
#   -cpu 128 		  (use all the cores - did 127 for testing so could htop in)
#   -b bin39	 	  (basename of output files)
#   -dp 			    (don't use the web-lookup for pre-calculated gene calls)
# 	-iprlookup 		(do use the ipr-lookup for GO terms and pathways)
# 	-goterms		  (report the GO terms in the output)
# 	-pa				    (also look up pathways)

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

# ---- first run on my example actinos ----

cherries <- bins[bins$bin.full.name == "ME2011-09-21_3300043464_group3_bin69" |
               bins$bin.full.name == "ME2016-07-20_3300033996_group7_bin32" |
               bins$bin.full.name == "ME2011-09-04_3300044729_group3_bin142", ]

slurm <- get.slum.specs(N = 1, simultaneous.tasks.per.node = 3, min.per.task = 180, tot.tasks = length(cherries$bin.full.name)) # memory limit- 3/node, time ~ 1:30 / job, double to 180 min

slurm.file <- c("#!/bin/bash",
                "#SBATCH -J interproactinos",
                "#SBATCH -p normal",
                paste("#SBATCH -N", slurm$N),
                paste("#SBATCH -n", slurm$n),
                paste("#SBATCH -t", slurm$t),
                "#SBATCH -o interproscan.o%j",
                "#SBATCH -e interproscan.e%j",
                "#SBATCH --mail-type=all",
                "#SBATCH --mail-user=rohwer@utexas.edu",
                "module load launcher",
                "export LAUNCHER_JOB_FILE=interpro_actinos.launcher",
                "export LAUNCHER_SCHED=interleaved",
                "${LAUNCHER_DIR}/paramrun")

launcher.file <- paste0("export PATH=\"/work/08922/rrohwer/ls6/jdk-17.0.6+10/bin:$PATH\" ; export PATH=\"/work/08922/rrohwer/ls6/interproscan-5.63-95.0:$PATH\" ; ",
                        "interproscan.sh -i /scratch/08922/rrohwer/runprodigal96_take2/",cherries$bin.full.name,".genes.faa",
                        " -f tsv gff3 -cpu ",128," -b ",cherries$bin.full.name," -dp -iprlookup -goterms -pa")

write.table(x = slurm.file, file = "server-scripts/generated/2023-04-19_run_interproscan_on_drep96_winners/interpro_actinos.slurm", quote = F, row.names = F, col.names = F)
write.table(x = launcher.file, file = "server-scripts/generated/2023-04-19_run_interproscan_on_drep96_winners/interpro_actinos.launcher", quote = F, row.names = F, col.names = F)

# ---- full run of all bins ----

# need to make folders for each run
# run 3 per node, 20 nodes at a time = 60 at a time = 47 jobs of 9 hrs.
# run 3 per node, 30 nodes at a time = 90 at a time = 32 jobs of 9 hrs.
# run 3 at once, 20 nodes at a time, but 3 in a row = 180 per job = 16 jobs of 27 hrs ***do this***
# run 3 at once, 30 nodes at a time, but 3 in a row = 270 per job = 11 jobs of 27 hrs

# ---- get bin groups for each run ----

bins$interpro.group <- rep(x = 1:16, each = 180)[1:2855]

# write these into a file structure that I can zip and transfer to tacc
dir.create("server-scripts/generated/2023-04-19_run_interproscan_on_drep96_winners/run_interproscan", showWarnings = F)

for (group in unique(bins$interpro.group)){
  
  slurm <- get.slum.specs(N = 20, simultaneous.tasks.per.node = 3, min.per.task = 180, tot.tasks = sum(bins$interpro.group == group)) # memory limit- 3/node, time ~ 1:30 / job, double to 180 min
  
  slurm.file <- c("#!/bin/bash",
                  "#SBATCH -J interpro",
                  "#SBATCH -p normal",
                  paste("#SBATCH -N", slurm$N),
                  paste("#SBATCH -n", slurm$n),
                  paste("#SBATCH -t", slurm$t),
                  "#SBATCH -o interproscan.o%j",
                  "#SBATCH -e interproscan.e%j",
                  "#SBATCH --mail-type=all",
                  "#SBATCH --mail-user=rohwer@utexas.edu",
                  "module load launcher",
                  "export LAUNCHER_JOB_FILE=runinterproscan.launcher",
                  "export LAUNCHER_SCHED=interleaved",
                  "${LAUNCHER_DIR}/paramrun")
  
  launcher.file <- paste0("export PATH=\"/work/08922/rrohwer/ls6/jdk-17.0.6+10/bin:$PATH\" ; export PATH=\"/work/08922/rrohwer/ls6/interproscan-5.63-95.0:$PATH\" ; ",
                          "interproscan.sh -i /scratch/08922/rrohwer/runprodigal96_take2/",bins$bin.full.name[bins$interpro.group == group],".genes.faa",
                          " -f tsv gff3 -cpu ",128," -b ",bins$bin.full.name[bins$interpro.group == group]," -dp -iprlookup -goterms -pa")
  
  dir.create(paste0("server-scripts/generated/2023-04-19_run_interproscan_on_drep96_winners/run_interproscan/group",group), showWarnings = F)
  write.table(x = slurm.file, file = paste0("server-scripts/generated/2023-04-19_run_interproscan_on_drep96_winners/run_interproscan/group",group,"/runinterproscan.slurm"), quote = F, row.names = F, col.names = F)
  write.table(x = launcher.file, file = paste0("server-scripts/generated/2023-04-19_run_interproscan_on_drep96_winners/run_interproscan/group",group,"/runinterproscan.launcher"), quote = F, row.names = F, col.names = F)
  
}

# ---- end ----
