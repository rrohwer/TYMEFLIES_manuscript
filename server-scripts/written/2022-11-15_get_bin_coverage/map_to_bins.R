# RRR
# 85,684 bins * 471 samples * 2 files (cov & bam) = 80,714,328 file
# 80 million is too many files to create, so can't run step by step, have to run whole workflow in each step
# for example, when I mapped the reads to the assemblies, I made all the references, then I did all the mapping, then I did the sorting, etc
# but this time I need to combine steps that do mapping, sorting, and then combining files so that file amts don't build up

bins <- readRDS("data/2022-11-15_bin_stats/all_bins.rds")

generated.scripts.folder <- "server-scripts/generated/2022-11-15_get_bin_coverage/"

# ----

index <- order(bins$bin.full.name)
bins <- bins[index, ]

samples <- unique(bins$tymeflies.name)

sample.paths <- paste0("filtered_reads/",samples,"_filtered-reads.fastq.gz")

bin.paths <- paste0("metabat2_bins_renamed/",bins$tymeflies.name,"/",bins$bin.filename)

results.paths <- paste0("mapping_wrkdir/",bins$tymeflies.name,"/",bins$bin.full.name,"/")




# ---- set up directories for each bin ----

setup.folders.launcher <- paste0("mkdir -p ", results.paths)

write.table(x = setup.folders.launcher, file = file.path(generated.scripts.folder,"setup_folders.launcher"), row.names = F, col.names = F, quote = F)

# ---- make mapping references ----

threads.per.process <- 128

make.ref.launcher <- paste0("module load bbmap; bbmap.sh ref=", bin.paths, " path=",results.paths, " t=", threads.per.process)

write.table(x = make.ref.launcher, file = file.path(generated.scripts.folder,"make_refs.launcher"), row.names = F, col.names = F, quote = F)

# ---- map, sort, and combine files ----

map.portion <- paste0("module load bbmap; bbmap.sh in=",sample.paths, " out=", paste0(results.paths, samples,".sam")," covstats=", paste0(results.paths, samples,".cov"), " path=", results.paths, " t=", threads.per.process)

sort.portion <- paste0("module load samtools; samtools sort -o ", paste0(results.paths, samples,".bam"), " -O bam -@ ", threads.per.process, " ", paste0(results.paths, samples,".sam"), "; rm ", paste0(results.paths, samples,".sam"))

combining.bams <- paste0("module load samtools; samtools merge")

combining.covs <- paste0("module load Rstats ; Rscript myscript")

map.sort.combine.launcher <- paste0(map.portion, "; ",sort.portion)

write.table(x = map.sort.combine.launcher, file = file.path(generated.scripts.folder, "map_sort_combine.launcher"), quote = F, row.names = F, col.names = F)

# ---- 



