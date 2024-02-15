# RRR

x <- readRDS(file = "data/2022-06-27_binning_groups/mapping_groups_list.rds")
output.folder <- "server-scripts/generated/2022-06-29_slurm_scripts/"

# ---- functions ----

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

get.slum.specs.multi.threaded <- function(tasks.per.node, min.per.task, sample.names){
  # tasks.per.node is the highest tasks per node allowed
  tot.tasks <- length(sample.names)
  N <- ceiling(tot.tasks / tasks.per.node)
  if (N == 0){
    N <- 1
  }
  adj.tasks.per.node <- ceiling(tot.tasks / N)
  n <- adj.tasks.per.node * N
  extra.tasks <- n - tot.tasks
  threads <- floor(128 / adj.tasks.per.node)
  m <- ceiling(min.per.task * adj.tasks.per.node)
  t <- get.timestamp.from.minutes(m = m)
  cat("N = ",N," n = ",n, " total = ",tot.tasks,"\ntime = ",t," (",min.per.task," min each)","\nextra filler lines = ",extra.tasks , sep = "", "\nmax tasks/node = ", tasks.per.node, ", adj. tasks per node = ", adj.tasks.per.node)
  return(list("n"=n,"N"=N,"t"=t,"threads" = threads,"filler.lines" = extra.tasks))
}

get.slum.specs.single.threaded <- function(tasks.per.node, min.per.task, sample.names){
  tot.tasks <- length(sample.names)
  N <- ceiling(tot.tasks / tasks.per.node)
  if (N == 0){
    N <- 1
  }
  adj.tasks.per.node <- ceiling(tot.tasks / N)
  n <- adj.tasks.per.node * N
  extra.tasks <- n - tot.tasks
  threads <- floor(128 / adj.tasks.per.node)
  m <- min.per.task 
  t <- get.timestamp.from.minutes(m = m)
  cat("N = ",N," n = ",n, " total = ",tot.tasks,"\ntime = ",t," (",min.per.task," min each)","\nextra filler lines = ",extra.tasks , sep = "", "\nmax tasks/node = ", tasks.per.node, ", adj. tasks per node = ", adj.tasks.per.node)
  return(list("n"=n,"N"=N,"t"=t,"threads" = threads,"filler.lines" = extra.tasks))
}

# -----

group <- names(x)[9]
for (group in names(x)){
  if( !dir.exists(paste0(output.folder,group))){
    dir.create(paste0(output.folder,group))
  }
  
  sample.names <- x[[group]]
  
  # ---- 0. setup directories ----
  min.per.task <- 5 * length(sample.names)
  slurm <- get.slum.specs.single.threaded(tasks.per.node = 32, min.per.task = min.per.task, sample.names = sample.names) # uses ~5 GB and 1 thread
  
  step.0.slurm <- c("#!/bin/bash",
                    "#SBATCH -J setup",
                    "#SBATCH -p normal",
                    paste("#SBATCH -N", slurm$N),
                    paste("#SBATCH -n", slurm$n),
                    paste0("#SBATCH -t ",slurm$t),
                    "#SBATCH -o setup.o%j",
                    "#SBATCH -e setup.e%j",
                    "#SBATCH --mail-type=all",
                    "#SBATCH --mail-user=rohwer@utexas.edu",
                    "module load biocontainers",
                    "module load launcher",
                    paste0("export OMP_NUM_THREADS=",slurm$threads),
                    "export LAUNCHER_JOB_FILE=step.0.launcher",
                    "export LAUNCHER_SCHED=interleaved",
                    "${LAUNCHER_DIR}/paramrun")
  
  step.0.launcher <- NULL
  for (assembly in sample.names){
    launcher.line <- NULL
    fasta <- paste0(assembly, "_assembled_contigs.fasta")
    launcher.line <- c(launcher.line, 
                paste("mkdir", assembly),
                paste("cd", assembly),
                paste0("cp /work/08922/rrohwer/stampede2/tymeflies/assembled_contigs/", fasta, " ./"))
    for (s in sample.names){
      fastq <- paste0(s, "_filtered-reads.fastq.gz")
      launcher.line <- c(launcher.line, 
                  paste0("cp /work/08922/rrohwer/stampede2/tymeflies/filtered_reads/", fastq, " ./"))
    }
    launcher.line <- c(launcher.line, 
                "cp ../step.2.slurm ./",
                "cp ../step.2.launcher ./",
                "cp ../step.3.slurm ./",
                "cp ../step.3.launcher ./",
                "cd ../")
    launcher.line <- paste(launcher.line, collapse = " ; ")
    step.0.launcher <- c(step.0.launcher, launcher.line)
  }
  
  write.table(x = step.0.slurm, file = file.path(output.folder, group, "step.0.slurm"), quote = F, row.names = F, col.names = F, append = F, sep = "\n")
  write.table(x = step.0.launcher, file = file.path(output.folder, group, "step.0.launcher"), quote = F, row.names = F, col.names = F, append = F, sep = "\n")
  
  # ---- 1. bbmap refs ----
  slurm <- get.slum.specs.single.threaded(tasks.per.node = 10, min.per.task = 10, sample.names = sample.names) # uses ~5 threads and ~15GB memory each
  
  step.1.slurm <- c("#!/bin/bash",
                    "#SBATCH -J bbmap_refs",
                    "#SBATCH -p normal",
                    paste("#SBATCH -N", slurm$N),
                    paste("#SBATCH -n", slurm$n),
                    paste0("#SBATCH -t ",slurm$t),
                    "#SBATCH -o bbmap_refs.o%j",
                    "#SBATCH -e bbmap_refs.e%j",
                    "#SBATCH --mail-type=all",
                    "#SBATCH --mail-user=rohwer@utexas.edu",
                    "module load biocontainers",
                    "module load launcher",
                    paste0("export OMP_NUM_THREADS=",slurm$threads),
                    "export LAUNCHER_JOB_FILE=step.1.launcher",
                    "export LAUNCHER_SCHED=interleaved",
                    "${LAUNCHER_DIR}/paramrun")
  
  step.1.launcher <- NULL
  for (s in sample.names){
    assembly <- paste0(s, "_assembled_contigs.fasta")
    step.1.launcher <- c(step.1.launcher,
                         paste0("module load bbmap; cd ", s, "; bbmap.sh ref=", assembly, " t=", slurm$threads, " usemodulo"))
  }
  if (slurm$filler.lines > 0){
    for (l in 1:slurm$filler.lines){
      step.1.launcher <- c(step.1.launcher, "echo bonus line")
    }
  }
    
  write.table(x = step.1.slurm, file = file.path(output.folder, group, "step.1.slurm"), quote = F, row.names = F, col.names = F, append = F, sep = "\n")
  write.table(x = step.1.launcher, file = file.path(output.folder, group, "step.1.launcher"), quote = F, row.names = F, col.names = F, append = F, sep = "\n")
  
  # ---- 2. mapping ----
  # same script in each folder
  slurm <- get.slum.specs.multi.threaded(tasks.per.node = 4, min.per.task = 60, sample.names = sample.names)
  
  step.2.slurm <- c("#!/bin/bash",
                    "#SBATCH -J mapping",
                    "#SBATCH -p normal",
                    paste("#SBATCH -N", slurm$N),
                    paste("#SBATCH -n", slurm$n),
                    paste0("#SBATCH -t ",slurm$t),
                    "#SBATCH -o mapping.o%j",
                    "#SBATCH -e mapping.e%j",
                    "#SBATCH --mail-type=all",
                    "#SBATCH --mail-user=rohwer@utexas.edu",
                    "module load biocontainers",
                    "module load launcher",
                    paste0("export OMP_NUM_THREADS=",slurm$threads),
                    "export LAUNCHER_JOB_FILE=step.2.launcher",
                    "export LAUNCHER_SCHED=interleaved",
                    "${LAUNCHER_DIR}/paramrun")
  
  step.2.launcher <- NULL
  for (s in sample.names){
    fastq <- paste0(s, "_filtered-reads.fastq.gz")
    sam <- paste0(s, ".sam")
    covstats <- paste0(s, ".cov")
    step.2.launcher <- c(step.2.launcher,
                         paste0("module load bbmap; bbmap.sh in=", fastq, " out=", sam," covstats=", covstats," usemodulo"))
  }
  if (slurm$filler.lines > 0){
    for (l in 1:slurm$filler.lines){
      step.2.launcher <- c(step.2.launcher, "echo bonus line")
    }
  }
  
  write.table(x = step.2.slurm, file = file.path(output.folder, group, "step.2.slurm"), quote = F, row.names = F, col.names = F, append = F, sep = "\n")
  write.table(x = step.2.launcher, file = file.path(output.folder, group, "step.2.launcher"), quote = F, row.names = F, col.names = F, append = F, sep = "\n")
  
  # ---- 3. sorting -----
  # same script in each assembly folder
  slurm <- get.slum.specs.single.threaded(tasks.per.node = 3, min.per.task = 60, sample.names = sample.names) # 78GB per task, view single-threaded, sort multi-threaded, so is single then multi
  
  step.3.slurm <- c("#!/bin/bash",
                    "#SBATCH -J sorting",
                    "#SBATCH -p normal",
                    paste("#SBATCH -N", slurm$N),
                    paste("#SBATCH -n", slurm$n),
                    paste0("#SBATCH -t ",slurm$t),
                    "#SBATCH -o sorting.o%j",
                    "#SBATCH -e sorting.e%j",
                    "#SBATCH --mail-type=all",
                    "#SBATCH --mail-user=rohwer@utexas.edu",
                    "module load biocontainers",
                    "module load launcher",
                    paste0("export OMP_NUM_THREADS=",slurm$threads),
                    "export LAUNCHER_JOB_FILE=step.3.launcher",
                    "export LAUNCHER_SCHED=interleaved",
                    "${LAUNCHER_DIR}/paramrun")
  
  step.3.launcher <- NULL
  for (s in sample.names){
    sam <- paste0(s, ".sam")
    bam <- paste0(s, ".bam")
    fastq <- paste0(s, "_filtered-reads.fastq.gz")
    step.3.launcher <- c(step.3.launcher,
                         paste0("module load samtools; samtools sort -o ", bam, " -O bam -@ ", slurm$threads, " ", sam, "; rm ", sam, "; rm ", fastq))
  }
  if (slurm$filler.lines > 0){
    for (l in 1:slurm$filler.lines){
      step.3.launcher <- c(step.3.launcher, "echo bonus line")
    }
  }
  
  write.table(x = step.3.slurm, file = file.path(output.folder, group, "step.3.slurm"), quote = F, row.names = F, col.names = F, append = F, sep = "\n")
  write.table(x = step.3.launcher, file = file.path(output.folder, group, "step.3.launcher"), quote = F, row.names = F, col.names = F, append = F, sep = "\n")
  
  # ---- 4. depth ----
  # one script for the whole group
  # single-threaded task
  min.per.task = 40 * length(sample.names)
  slurm <- get.slum.specs.single.threaded(tasks.per.node = 80, min.per.task = min.per.task, sample.names = sample.names) # single threaded, 5 GB/3 samples, so max/node is ~ 80
  
  step.4.slurm <- c("#!/bin/bash",
                    "#SBATCH -J depths",
                    "#SBATCH -p normal",
                    paste("#SBATCH -N", slurm$N),
                    paste("#SBATCH -n", slurm$n),
                    paste0("#SBATCH -t ",slurm$t),
                    "#SBATCH -o depths.o%j",
                    "#SBATCH -e depths.e%j",
                    "#SBATCH --mail-type=all",
                    "#SBATCH --mail-user=rohwer@utexas.edu",
                    "module load biocontainers",
                    "module load launcher",
                    paste0("export OMP_NUM_THREADS=",slurm$threads),
                    "export LAUNCHER_JOB_FILE=step.4.launcher",
                    "export LAUNCHER_SCHED=interleaved",
                    "${LAUNCHER_DIR}/paramrun")
  
  step.4.launcher <- NULL
  for (assembly in sample.names){
    all.bams <- paste0(sample.names, ".bam")
    all.bams <- paste(all.bams, collapse = " ")
    step.4.launcher <- c(step.4.launcher,
                         paste0("module load metabat2; cd ", assembly,"; jgi_summarize_bam_contig_depths --outputDepth ", assembly, ".depth ", all.bams))
  }
  if (slurm$filler.lines > 0){
    for (l in 1:slurm$filler.lines){
      step.4.launcher <- c(step.4.launcher, "echo bonus line")
    }
  }
  write.table(x = step.4.slurm, file = file.path(output.folder, group, "step.4.slurm"), quote = F, row.names = F, col.names = F, append = F, sep = "\n")
  write.table(x = step.4.launcher, file = file.path(output.folder, group, "step.4.launcher"), quote = F, row.names = F, col.names = F, append = F, sep = "\n")
  
  # ---- 5. metabat2 ----
  # one script for whole group
  # single-threaded task
  min.per.task = 20 * length(sample.names)
  slurm <- get.slum.specs.single.threaded(tasks.per.node = 31, min.per.task = min.per.task, sample.names = sample.names)
  
  step.5.slurm <- c("#!/bin/bash",
                    "#SBATCH -J metabat",
                    "#SBATCH -p normal",
                    paste("#SBATCH -N", slurm$N),
                    paste("#SBATCH -n", slurm$n),
                    paste0("#SBATCH -t ",slurm$t),
                    "#SBATCH -o metabat.o%j",
                    "#SBATCH -e metabat.e%j",
                    "#SBATCH --mail-type=all",
                    "#SBATCH --mail-user=rohwer@utexas.edu",
                    "module load biocontainers",
                    "module load launcher",
                    paste0("export OMP_NUM_THREADS=",slurm$threads),
                    "export LAUNCHER_JOB_FILE=step.5.launcher",
                    "export LAUNCHER_SCHED=interleaved",
                    "${LAUNCHER_DIR}/paramrun")
  
  step.5.launcher <- NULL
  for (assembly in sample.names){
    fasta <- paste0(assembly,"_assembled_contigs.fasta")
    depth <- paste0(assembly, ".depth")
    bin.prefix <- paste0(assembly,"/bin")
    step.5.launcher <- c(step.5.launcher,
                         paste0("module load metabat2; cd ", assembly,"; metabat2 -i ",fasta, " -a ", depth, " -o ", bin.prefix, " -t ", slurm$threads ))
  }
  if (slurm$filler.lines > 0){
    for (l in 1:slurm$filler.lines){
      step.5.launcher <- c(step.5.launcher, "echo bonus line")
    }
  }
  
  write.table(x = step.5.slurm, file = file.path(output.folder, group, "step.5.slurm"), quote = F, row.names = F, col.names = F, append = F, sep = "\n")
  write.table(x = step.5.launcher, file = file.path(output.folder, group, "step.5.launcher"), quote = F, row.names = F, col.names = F, append = F, sep = "\n")
  
  # ---- 6. delete bams and other ----
  # one script for whole group
  # can run from login node
  
  step.6 <- "#!/usr/bin/bash"
  for (assembly in sample.names){
    step.6 <- c(step.6, 
                paste("cd", assembly))
    fasta <- paste0(assembly, "_assembled_contigs.fasta")
    step.6 <- c(step.6, 
                paste0("rm ", fasta),
                "rm -Rf ref")
    for (s in sample.names){
      bam <- paste0(s, ".bam")
      step.6 <- c(step.6, 
                  paste0("rm ", bam))
    }
    step.6 <- c(step.6, "cd ../")
  }
  
  write.table(x = step.6, file = file.path(output.folder, group, "step.6.sh"), quote = F, row.names = F, col.names = F, append = F, sep = "\n")
  
  # ---- 7. save rest of data, transfer Corral ----
  # zip folders into convenient single file
  min.per.task <- 5 * length(sample.names)
  slurm <- get.slum.specs.single.threaded(tasks.per.node = 32, min.per.task = min.per.task, sample.names = sample.names) # uses ~5 GB and 1 thread
  
  step.7.slurm <- c("#!/bin/bash",
                    "#SBATCH -J cleanup",
                    "#SBATCH -p normal",
                    paste("#SBATCH -N", slurm$N),
                    paste("#SBATCH -n", slurm$n),
                    paste0("#SBATCH -t ",slurm$t),
                    "#SBATCH -o cleanup.o%j",
                    "#SBATCH -e cleanup.e%j",
                    "#SBATCH --mail-type=all",
                    "#SBATCH --mail-user=rohwer@utexas.edu",
                    "module load biocontainers",
                    "module load launcher",
                    paste0("export OMP_NUM_THREADS=",slurm$threads),
                    "export LAUNCHER_JOB_FILE=step.7.launcher",
                    "export LAUNCHER_SCHED=interleaved",
                    "${LAUNCHER_DIR}/paramrun")
  
  step.7.launcher <-  NULL
  for (assembly in sample.names){
    launcher.line <- NULL
    bin.folder <- assembly
    tarball.bins <- paste0(bin.folder,".tar.gz")
    tarball.covs <- paste0(assembly,".tar.gz")
    launcher.line <- c(launcher.line, 
                       paste0("cd ", assembly),
                       paste0("tar cvfz ",tarball.bins," ",bin.folder),
                       paste0("cp ",tarball.bins," /corral/utexas/MCB22035/tymeflies/metabat2_bins"),
                       paste0("mv ",tarball.bins," /scratch/08922/rrohwer/robinning/metabat2_bins"),
                       paste0("rm -Rf ", bin.folder),
                       "cd ../",
                       paste0("tar cvfz ",tarball.covs," ",assembly),
                       paste0("cp ",tarball.covs," /corral/utexas/MCB22035/tymeflies/mapping_coverages"),
                       paste0("mv ",tarball.covs," /scratch/08922/rrohwer/robinning/mapping_coverages"),
                       paste0("rm -Rf ", assembly))
    launcher.line <- paste(launcher.line, collapse = " ; ")
    step.7.launcher <- c(step.7.launcher, launcher.line)
  }
  
  write.table(x = step.7.slurm, file = file.path(output.folder, group, "step.7.slurm"), quote = F, row.names = F, col.names = F, append = F, sep = "\n")
  write.table(x = step.7.launcher, file = file.path(output.folder, group, "step.7.launcher"), quote = F, row.names = F, col.names = F, append = F, sep = "\n")
  
}


# ---- make tracking excel sheet ----

all.sample.names <- NULL
all.mapping.groups <- NULL
for (group in names(x)){
  all.sample.names <- c(all.sample.names, x[[group]])
  all.mapping.groups <- c(all.mapping.groups, rep(group,length(x[[group]])))
}

tracking.sheet <- data.frame("Sample.Name" = all.sample.names, "Mapping.Group" = all.mapping.groups, 
                             "step.0" = "", "0.run.t" = "",
                             "step.1" = "", "1.N" = "", "1.n" = "", "1.Done" = "", "1.wait.t" = "", "1.run.t" = "", 
                             "step.2" = "", "2.N" = "", "2.n" = "", "2.Done" = "", "2.wait.t" = "", "2.run.t" = "", 
                             "step.3" = "", "3.N" = "", "3.n" = "", "3.Done" = "", "3.wait.t" = "", "3.run.t" = "", 
                             "step.4" = "", "4.N" = "", "4.n" = "", "4.Done" = "", "4.wait.t" = "", "4.run.t" = "", 
                             "step.5" = "", "5.N" = "", "5.n" = "", "5.Done" = "", "5.wait.t" = "", "5.run.t" = "", 
                             "step.6" = "", "6.run.t" = "",  
                             "step.7" = "", "7.run.t" = "" )
write.table(x = tracking.sheet, file = "data/2022-07-07_tracking_mapping/tracking_template.csv", append = F, quote = F, sep = ",",row.names = F, col.names = T)





