#!/bin/bash

# RRR
# This is from Alicia
# Gets the contigs filename for the APID that was good

# source as ./this.sh > output.txt



for i in `cat ats`
do
	jamo report select metadata.sequencing_project_id,metadata.analysis_project_id,file_path,file_name where metadata.analysis_task_id=$i and file_type=contigs
done