#!/bin/bash

# RRR
# Run each section alone with others commented out:
# terminal input syntax includes sample key file as an argument,
# sample key should be col 1 = SPIDs, col 2 = TYMEFLIES names:
# 		./2021_12-10_get_my_data_filtered_metaG_reads.txt samplekey.txt

# ---- set-up ----

module load jamo

SPIDVECT=`awk '{ print $1 }' SPID.txt`
SPIDARRAY=($SPIDVECT)

echo "SPID FILENAME" > uniqueSPIDS.txt
echo "SPID FILENAME" > dupSPIDS.txt


# ---- restore with JAMO, and save filenames ----

for SPID in ${SPIDARRAY[@]}
do
	echo $SPID
	FILEID=`jamo info all spid $SPID | grep -i assembly.contigs | awk '{ print $4 }'`	
	FILENAME=`jamo info all spid $SPID | grep -i assembly.contigs | awk '{ print $2 }'`	
	
	jamo fetch all id $FILEID

	FILENAME=($FILENAME)
	if [ ${#FILENAME[@]} -gt 1 ]
	then
		for EACHFILE in ${FILENAME[@]}
		do
			echo $SPID $EACHFILE >> dupSPIDS.txt
		done
	else
		echo $SPID $FILENAME >> uniqueSPIDS.txt
	fi

done

# ---- end ----

exit