#!/bin/bash

module load python
conda activate robinning

bbmap.sh in=../reads/ME2000-06-06pf_3300042549_filtered-reads.fastq.gz out=ME2000-06-06pf_3300042549.bam covstats=ME2000-06-06pf_3300042549.cov &
bbmap.sh in=../reads/ME2000-06-20pf_3300042868_filtered-reads.fastq.gz out=ME2000-06-20pf_3300042868.bam covstats=ME2000-06-20pf_3300042868.cov &
bbmap.sh in=../reads/ME2001-07-02D8pf_3300042397_filtered-reads.fastq.gz out=ME2001-07-02D8pf_3300042397.bam covstats=ME2001-07-02D8pf_3300042397.cov &
bbmap.sh in=../reads/ME2008-09-15D10_3300042539_filtered-reads.fastq.gz out=ME2008-09-15D10_3300042539.bam covstats=ME2008-09-15D10_3300042539.cov &
bbmap.sh in=../reads/ME2011-06-28_3300042360_filtered-reads.fastq.gz out=ME2011-06-28_3300042360.bam covstats=ME2011-06-28_3300042360.cov &
bbmap.sh in=../reads/ME2013-05-21_3300042460_filtered-reads.fastq.gz out=ME2013-05-21_3300042460.bam covstats=ME2013-05-21_3300042460.cov &
bbmap.sh in=../reads/ME2013-08-23_3300034106_filtered-reads.fastq.gz out=ME2013-08-23_3300034106.bam covstats=ME2013-08-23_3300034106.cov &
bbmap.sh in=../reads/ME2015-06-10_3300042356_filtered-reads.fastq.gz out=ME2015-06-10_3300042356.bam covstats=ME2015-06-10_3300042356.cov &
bbmap.sh in=../reads/ME2015-06-27_3300042558_filtered-reads.fastq.gz out=ME2015-06-27_3300042558.bam covstats=ME2015-06-27_3300042558.cov &
bbmap.sh in=../reads/ME2017-05-25_3300042863_filtered-reads.fastq.gz out=ME2017-05-25_3300042863.bam covstats=ME2017-05-25_3300042863.cov &
bbmap.sh in=../reads/ME2017-06-16_3300042330_filtered-reads.fastq.gz out=ME2017-06-16_3300042330.bam covstats=ME2017-06-16_3300042330.cov &
bbmap.sh in=../reads/ME2018-06-07_3300034013_filtered-reads.fastq.gz out=ME2018-06-07_3300034013.bam covstats=ME2018-06-07_3300034013.cov &

wait