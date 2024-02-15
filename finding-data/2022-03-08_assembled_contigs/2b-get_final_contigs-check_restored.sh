#!/bin/bash

# RRR

# ---- set-up ----

module load jamo

SPIDVECT=`awk '{ print $1 }' SPIDs_with_different_filenames.txt`
SPIDARRAY=($SPIDVECT)

# echo "FileStatus" > jamo_status.txt


# ---- restore with JAMO, and save filenames ----

for SPID in ${SPIDARRAY[@]}
do
	echo $SPID
	FILESTATUS=`jamo info all spid $SPID | grep -i final.contigs.fasta | awk '{ print $3 }'`	
	echo $FILESTATUS >> jamo_status.txt
done

# ---- end ----

exit