#!/bin/bash -l
# #$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N grabtaxaID
#$ -l mem_per_core=4G
#$ -l h_rt=12:00:00
#$ -P incrna
#$ -M james.e.fifer@gmail.com
##run blast.sh before running the following commands
##download taxonkit from https://github.com/shenwei356/taxonkit

##We want to take the NCBI acc numbers from the blast output find the taxa ID and then use taxon kit 
##
##first take NCBI accession numbers
#cut -f2 nt.blast.out --output-delimiter='|' |rev| cut -f2 -d '|' |rev > NCBI_acc.input

## 
#module load perl
#results=();cat NCBI_acc.input | while read p ;do results=($(efetch -db nuccore -id "$p" -format docsum |\
# xtract -pattern DocumentSummary -element TaxId)); echo "$p" "${results[@]}" >> NCBI_acc.output; done

#cut -f2 NCBI_acc.output -d " ">taxa.ID

##get the taxonomy ID files 
#wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.Z
##unzip file to get .dmp files
##directory below is location of .dmp files
taxonkit lineage --data-dir "/projectnb/davieslab/jfifer/Fouha/Analysis/"  ./taxa.ID > taxa.ID.output
cut -f1 nt.blast.out>taxa.seq.ID; paste taxa.seq.ID taxa.ID.output > taxa.ID.output.withseq
grep "Eukaryota" taxa.ID.output.withseq | cut -f1 | sort | uniq >euk_contam_asvs
