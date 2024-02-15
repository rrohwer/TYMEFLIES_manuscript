# RRR
# my test run of kofanscan on one genome was super fast, so I tried to just run it in one big go:
# exec_annotation drep_96_winners_concat.faa -o drep_96_winners_concat_kofamscan.tsv -p ../kofam_databases/profiles -k ../kofam_databases/ko_list --cpu 100 --format detail-tsv
# This started out running hmmer on 100 threads, but now it's been running ruby on a single thread for the whole day.
# So... I don't think it is fully parallelized and that unparallelized step is taking forever.
# I'm going to kill that job and instead run one where each genome gets 1 thread but they run in parallel on 100 threads using parallel
# one genome on one thread took about 9 minutes, so 10 min/2855 genomes/100 threads ~ 5 hrs

# oh but the tmp folder is going to get overwritten as they go in parallel. can create different ones with a flag

key <- readRDS(file = "data/2023-02-24_dRep_with_range_ANIs/output_processed/drep_results_all_genomes_0.96.rds")
key[1:5, ]
key <- key[key$winner, ]

bin.faas <- paste0(key$bin.full.name,".genes.faa")
bin.output <- paste0(key$bin.full.name,".kofamscan.tsv")
bin.tmp.dir <- paste0("temp/temp_",key$bin.full.name)

my.commands <- paste0("exec_annotation ",
                      "../../yggshare/current_projects/TYMEFLIES/tymeflies/runprodigal96_take2/",bin.faas," ",
                      "-o ",bin.output," ",
                      "-p ../kofam_databases/profiles -k ../kofam_databases/ko_list --cpu 1 --format detail-tsv ",
                      "--tmp-dir=",bin.tmp.dir)

write.table(x = my.commands, file = "server-scripts/generated/2023-09-13_run_kofamscan/kofamscan_by_genome.jobfile", quote = F, sep = "\n", row.names = F, col.names = F)

# and run this with
# make a working directory in home on Vanaheim
# add the jobfile
# add a directory called "temp"
# conda activate kofamscan
# cat kofamscan_by_genome.jobfile | parallel -j 115

# and this is kofamscan version 1.3.0 and databases last updated 2023-09-06