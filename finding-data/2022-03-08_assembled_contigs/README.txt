README
RRR 3/8/22
These are the steps for moving and renaming the assembly files

1- fetch the assembly.contigs.fasta file for each SPID
	use 
		SPIDs.txt
	call
		./1-get_assembled_contigs-find_and_fetch.sh
	create
		dupSPIDS.txt
		uniqueSPIDS.txt
		
1b- A few SPIDs have a different filename: final.contigs.fasta
	use 
		SPIDs_with_different_filenames.txt
	call
		./1-get_assembled_contigs-find_and_fetch.sh
	add to
		dupSPIDS.txt
		uniqueSPIDS.txt
		
2- Check that all the fetched files have been restored
	use 
		SPID.txt
	call
		./2-get_assembled_contigs-check_restored.sh
	create
		jamo_status.txt

2b- check status on the oddly named ones too
	use 
		SPIDs_with_different_filenames.txt
	call
		./2b-get_final_contigs-check_restored.sh
	add to
		jamo_status.txt

3a- find the good APID for each duplicated SPID
	* This is copied from Alicia's email, and I got a permission error so did not use this *
3b- find the filenames of the assembly file for the APID
	* Instead, just used the output table Alicia had pasted into the email *

3c- choose the correct filenames
	manually create and then use
		dupSPIDS_chosen.xlsx
	call
		3c-choose_correct_filenames.R
	create
		dupSPIDS_chosen.txt

4- generate a shell script that will copy over the files, and rename them
	use either
		uniqueSPIDS.txt
		BIGTYME.xlsx
	OR
		dupSPIDS_chosen.txt
		BIGTYME.xlsx
	run
		4-generate_bash_file_to_move_and_rename_files.R
	create either
		4b-get_contigs-copy_and_rename-dupSPIDs.sh
	OR
		4a-get_contigs-copy_and_rename-uniqueSPIDs.sh 
	call both 4a and 4b scripts on the dtn node
	
5- add the filenames used into the BIGTYME.xlsx reference sheet
	use
		BIGTYME.xlsx
		filename_key_unique.rds
		filename_key_dupSPIDs.rds
	run
		5-add_filenames_to_BIGTYME.R
	create
		2022-03-08_bigtyme_plus_assembly_filenames.tsv
	and paste in the new column manually (checking that SPIDs match)

6- add in the SPID I accidentally left out!
	* this is one that had the different filename, but I accidentally missed it *
	* just do all steps manually for this one file *
		SPID: 1229871
	commands
		module load jamo
		jamo info all spid 1229871 | grep -i final.contigs.fasta
			[1229871] /global/dna/dm_archive/rqc/analyses/AUTO-228479/final.contigs.fasta PURGED 5c916a9d46d1e61fc468a1b4 
		jamo fetch all id 5c916a9d46d1e61fc468a1b4
		jamo info all spid 1229871 | grep -i final.contigs.fasta
			[1229871] /global/dna/dm_archive/rqc/analyses/AUTO-228479/final.contigs.fasta RESTORED 5c916a9d46d1e61fc468a1b4
		manually create filename using bigtyme sheet: ME2014-08-24_3300033981_assembled_contigs.fasta
		cp -i /global/dna/dm_archive/rqc/analyses/AUTO-228479/final.contigs.fasta /global/cfs/cdirs/pkmeco/tymeflies/assembled_contigs/ME2014-08-24_3300033981_assembled_contigs.fasta
		manually add filename to Bigtyme sheet: /global/dna/dm_archive/rqc/analyses/AUTO-228479/final.contigs.fasta
