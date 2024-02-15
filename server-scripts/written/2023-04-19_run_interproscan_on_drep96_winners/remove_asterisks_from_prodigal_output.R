# RRR

# you need to remove the asterisks that prodigal puts at the end of sequences
# they are stop codons, which are usually not included in the amino acid sequence

# also need to remove any _ , but I checked and there are none:
# $ cat runprodigal96/drep_96_winners_concat.faa | grep ">" | wc -l
# 7607541
# $ cat runprodigal96/drep_96_winners_concat.faa | grep "_" | wc -l
# 7607541
# since _ is in all the sequence names, and number sequences eauals number underscors, no underscores in sequence part.
# refs:
# https://github.com/milaboratory/mixcr/issues/204
# https://github.com/Gaius-Augustus/BRAKER/issues/56

# oh but there are no asterisks in any sequence name
# $ cat runprodigal96/drep_96_winners_concat.faa | grep ">" | grep "\*" | wc -l
# 0
# so I can just remove all asterisks, no need for a script!!
# sed -i 's/\*//g'
# ref on what -i means (replace within the file, avoids making a new file and then renaming it)
# https://stackoverflow.com/questions/18527365/what-does-sed-i-option-do

# also, the first time I ran prodigal I only saved the concatenated output, which is what fed into inStrain
# but it would be better to split up the file to run interproscan in parallel
# interproscan can do multithreading, across multiple nodes the same command
# but this seems very complicated to figure out. Instead I'll just run multiple instances with launcher.
# I re-ran prodigal (it's very fast), and found that it will give identical output each time 
# so I can use the new by-genome prodigal output for this, and have separate runs per sample.
# ref on prodigal being the same each time: https://groups.google.com/g/prodigal-discuss/c/RQ-BGFbC2xM/m/pJVvPyQxCQAJ

# so remove asterisks from every genome file in the take_2 folder
# $ ls *\.faa | wc -l
# 2855
# so one .faa file per bin
# RUN:
# find *\.faa -exec sed -i 's/\*//g' {} \;