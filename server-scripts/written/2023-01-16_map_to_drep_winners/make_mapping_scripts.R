# RRR
# Make a script that maps all the metaGs to the concatenated drep bins
# this will run a ~1 hr mapping job 471 times = 20 days, so I will run in parallel
# say on 20 nodes, it will take one day, on 40 nodes 1/2 day
# since the mapping parallelizes fully, can run 1 at a time on a node. 
# this will also save memory since I'm not using usemodulo anymore.
# 1/node at a time means using launcher's scheduler differently than before
# It will run like
# N=20
# n=471
# and when the first 20 finish, start the next in line
# No need for bonus lines anymore, but do still need N to divide n evenly

# Example command:
# bbmap.sh in=ME2000-03-30D8pf_3300042899_filtered-reads.fastq.gz out=ME2000-03-30D8pf_3300042899.sam 32bit=t
#   notes: no need to save covstats since I am saving the BAMS this time
#   notes: needed 32bit=t because got an error without it in my test run, doesn't slow it down
#   notes: not using usemodulo, that only saves memory doesn't speed it up
#   notes: didn't specify threads before, so no need to start now


# ---- set-up ----

sample.names <- readRDS(file = "data/2022-11-15_bin_stats/all_bins.rds")
sample.names <- unique(sample.names$tymeflies.name)


file.map.slurm <- "server-scripts/generated/2023-01-16_map_to_drep_winners/map_to_drep.slurm"
file.map.launcher <- "server-scripts/generated/2023-01-16_map_to_drep_winners/map_to_drep.launcher"

file.sort.slurm <- "server-scripts/generated/2023-01-16_map_to_drep_winners/sort_drep_sams.slurm"
file.sort.launcher <- "server-scripts/generated/2023-01-16_map_to_drep_winners/sort_drep_sams.launcher"

file.delete.sams <- "server-scripts/generated/2023-01-16_map_to_drep_winners/delete_sams.sh"

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

# ---- bbmap ref -----

# Run in idev, single time and goes fast:

# idev -t 2:00:00
# module load biocontainers
# module load bbmap
# bbmap.sh ref=drep_winners_concat.fna
# exit

# Total time:     	224.565 seconds.

# ---- bbmap map ----

slurm <- get.slum.specs(N = 40, simultaneous.tasks.per.node = 1, min.per.task = 180, tot.tasks = 471)

map.slurm <- c("#!/bin/bash",
               "#SBATCH -J mapping",
               "#SBATCH -p normal",
               paste("#SBATCH -N", slurm$N),
               paste("#SBATCH -n", slurm$n),
               paste("#SBATCH -t", slurm$t),
               "#SBATCH -o mapping.o%j",
               "#SBATCH -e mapping.e%j",
               "#SBATCH --mail-type=all",
               "#SBATCH --mail-user=rohwer@utexas.edu",
               "module load biocontainers",
               "module load launcher",
               "export LAUNCHER_JOB_FILE=map_to_drep.launcher",
               "export LAUNCHER_SCHED=interleaved",
               "${LAUNCHER_DIR}/paramrun")

map.launcher <- NULL
for(s in sample.names){
  map.launcher <- c(map.launcher,
                    paste0("module load bbmap; bbmap.sh in=",s,"_filtered-reads.fastq.gz out=",s,".sam 32bit=t"))
}

write.table(x = map.slurm, file = file.map.slurm, quote = F, row.names = F, col.names = F, append = F, sep = "\n")
write.table(x = map.launcher, file = file.map.launcher, quote = F, row.names = F, col.names = F, append = F, sep = "\n")

# with 40 nodes, took 10:32:33

# ---- samtools sort ----

slurm <- get.slum.specs(N = 40, simultaneous.tasks.per.node = 3, min.per.task = 60, tot.tasks = 471)

sort.slurm <- c("#!/bin/bash",
               "#SBATCH -J mapping",
               "#SBATCH -p normal",
               paste("#SBATCH -N", slurm$N),
               paste("#SBATCH -n", slurm$n),
               paste("#SBATCH -t", slurm$t),
               "#SBATCH -o mapping.o%j",
               "#SBATCH -e mapping.e%j",
               "#SBATCH --mail-type=all",
               "#SBATCH --mail-user=rohwer@utexas.edu",
               "module load biocontainers",
               "module load launcher",
               "export LAUNCHER_JOB_FILE=sort_drep_sams.launcher",
               "export LAUNCHER_SCHED=interleaved",
               "${LAUNCHER_DIR}/paramrun")

sort.launcher <- NULL
for(s in sample.names){
  sort.launcher <- c(sort.launcher,
                    paste0("module load samtools; samtools sort -o ",s,".bam -O bam -@ ",slurm$threads - 1," ",s,".sam")) 
}

write.table(x = sort.slurm, file = file.sort.slurm, quote = F, row.names = F, col.names = F, append = F, sep = "\n")
write.table(x = sort.launcher, file = file.sort.launcher, quote = F, row.names = F, col.names = F, append = F, sep = "\n")

# with 40 nodes, took 00:52:22

# ---- delete sams ----

delete.sams <- "#!/bin/bash"
for (s in sample.names){
  delete.sams <- c(delete.sams,
                   paste0("rm ",s,".sam"))
}

write.table(x = delete.sams, file = file.delete.sams, quote = F, row.names = F, col.names = F, append = F, sep = "\n")

