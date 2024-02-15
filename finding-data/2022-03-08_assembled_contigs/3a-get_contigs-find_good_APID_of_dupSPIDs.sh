#!/bin/bash

# RRR
# This is from Alicia
# Gets the ITS APID that is good for each SPID

# source as ./this.sh > output.txt


for i in `awk '{print $1}' dupSPIDS.txt` 
do 
	jat report select metadata.sequencing_project_id,publish,metadata.analysis_project_id,metadata.analysis_task_id where outputs.label=contigs and publish=true and metadata.sequencing_project_id=$i
done