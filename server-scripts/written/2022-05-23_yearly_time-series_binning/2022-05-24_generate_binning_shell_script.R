# RRR

created.shell.script <- "server-scripts/bash/2022-05-23_yearly_time-series_binning/bin-on-helheim-"

groups <- readRDS(file = "data/2022-05-24_metabat_by_year/samplenames_by_year_subsampled.rds")

for (g in 1:length(groups)){
  my.commands <- c("#!/bin/bash","# RRR","\n","set -ex")
  my.g <- groups[[g]]
  generated.file <- paste0(created.shell.script,names(groups)[g],".sh")
  
  group.message <- paste0("\n\n# starting on group ", g, "- ",names(groups)[g])
  
  starting.fasta <- paste0("/home/yggshare/current_projects/TYMEFLIES/tymeflies/assembled_contigs/",my.g, "_assembled_contigs.fasta")
  copy.fasta <- paste0("cp ",starting.fasta," ./ &")
  starting.fastq <- paste0("/home/yggshare/current_projects/TYMEFLIES/tymeflies/filtered_reads/",my.g, "_filtered-reads.fastq.gz")
  copy.fastq <- paste0("cp ",starting.fastq," ./ &")
  
  my.commands <- c(my.commands, group.message, copy.fasta, copy.fastq, "wait")
  
  for (a in 1:length(my.g)){
    assembly <- paste0("../", my.g[a], "_assembled_contigs.fasta")
    
    assembly.message <- paste0("\n# starting on assembly ", a, "- ", my.g[a])
    make.assembly.folder <- paste0("mkdir ", my.g[a])
    cd.assembly.folder <- paste0("cd ", my.g[a])
    activate.bbmap <- "source /home/rrohwer/miniconda3/bin/activate bbmap_robin"
    
    bbmap.ref <- paste0("bbmap.sh ref=", assembly, " t=118 usemodulo")
    
    my.commands <- c(my.commands, assembly.message, make.assembly.folder, cd.assembly.folder, activate.bbmap, "wait", bbmap.ref)
    
    created.bams <- NULL
    for (r in 1:length(my.g)){
      reads <- paste0("../", my.g[r], "_filtered-reads.fastq.gz")
      
      bbmap.reads <- paste0("bbmap.sh in=", reads, " out=", my.g[r],".sam covstats=",my.g[r],".cov t=118 usemodulo bamscript=bs.sh; sh bs.sh")
      
      delete.sam <- paste0("rm ",my.g[r],".sam &")
      
      my.commands <- c(my.commands, bbmap.reads, delete.sam)
      
      created.bams <- c(created.bams, paste0(my.g[r],"_sorted.bam"))
    }
    
    activate.metabat <- "source /home/rrohwer/miniconda3/bin/activate metabat2_robin"
    
    depth <- paste0(my.g[a],".depth ")
    all.bams <- paste(created.bams, collapse = " ")
    metabat.depth <- paste0("jgi_summarize_bam_contig_depths --outputDepth ", depth,all.bams)
    
    bin.folder <- paste0(my.g[a],"/",my.g[a],".bin")
    metabat.bin <- paste0("metabat2 -i ", assembly, " -a ", depth, " -o ", bin.folder, " -t 118")
    
    leave.folder <- "cd ../"
    
    my.commands <- c(my.commands, activate.metabat, "wait", metabat.depth, metabat.bin, "wait", leave.folder)
    
  }
  
  finished.message <- "\n# finished group"
  working.fasta <- paste0(my.g, "_assembled_contigs.fasta")
  delete.fasta <- paste0("rm ",working.fasta," &")
  working.fastq <- paste0(my.g, "_filtered-reads.fastq.gz")
  delete.fastq <- paste0("rm ",working.fastq," &")
  
  my.commands <- c(my.commands, finished.message, delete.fasta, delete.fastq)
  
  write.table(x = my.commands, file = generated.file, quote = F, sep = "\n", append = F, row.names = F, col.names = F)  
}

  
