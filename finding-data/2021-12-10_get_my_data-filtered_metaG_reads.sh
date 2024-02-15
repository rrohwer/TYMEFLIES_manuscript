#!/bin/bash

# RRR
# Run each section alone with others commented out:
# terminal input syntax includes sample key file as an argument,
# sample key should be col 1 = SPIDs, col 2 = TYMEFLIES names:
# 		./2021_12-10_get_my_data_filtered_metaG_reads.txt samplekey.txt

# ---- set-up ----

FOLDERPATH=/global/cfs/cdirs/pkmeco/TYMEFLIES/filtered_reads/

SPIDVECT=`awk '{ print $1 }' SPID-NAME_test.txt`
NAMESVECT=`awk '{ print $2 }' SPID-NAME_test.txt`

SPIDARRAY=($SPIDVECT)
NAMESARRAY=($NAMESARRAY)

module load jamo

echo $SPIDLIST
echo $MYNAMES


# ---- restore with JAMO ----

for SPID in $SPIDARRAY\.fastq\.gz
do
	echo $SPID
	FILEID=`jamo info all spid $SPID | grep -i filter.*\.fastq.gz | awk '{ print $4 }'`
	echo $FILEID
# 	jamo fetch all id $FILEID
done




# ---- check that all are restored ----
# 
# for SPID in $SPIDLIST
# do
# 	FILESTATUS=`jamo info all spid $SPID | grep -i filter.*\.fastq.gz | awk '{ print $3 }'`
# 	if [ $FILESTATUS != "restored"]; then
# 		echo $SPID is still $FILESTATUS \n
# 	fi
# done
# 
# 
# 
# 
# ---- copy files (run from dtn, use tmux) ----
# 
# touch original_filenames.txt
# echo "Name_In_Jamo,SPID,Transfer_Date,Name_In_filtered_reads" > $FOLDERPATH/original_filenames.txt
# 
# echo There are ${#SPIDLIST[@]} samples total
# 
# INDEX=0
# while [$INDEX -le ${#SPIDLIST[@]}
# do
# 	FILENAME=`jamo info all spid ${SPIDLIST[${INDEX}] | grep -i filter.*\.fastq.gz | awk '{ print $2 }'`
# 	cp -i $FILENAME $FOLDERPATH/${MYNAMES[$INDEX]}.fastq.gz # -i makes interactive if filename already exists
# 	curr_time=$(date "+%c")
# 	echo $FILENAME,${SPIDLIST[$INDEX]},${MYNAMES[$INDEX]},$curr_time >> $FOLDERPATH/original_filenames.txt
# done
# 
# 
# exit