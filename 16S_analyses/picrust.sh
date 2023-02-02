#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N newpipeline # job name, anything you want
#$ -l h_rt=12:00:00
#$ -M james.e.fifer@gmail.com #your email
#$ -m be

#follow this pipeline https://github.com/picrust/picrust2/wiki/PICRUSt2-Tutorial-(v2.3.0-beta)

#below can be run without submitting job (~5 mins)
module load miniconda
conda activate picrust2
#place_seqs.py -s ./all.rare.trim.fa -o out.tre -p 1 \
#              --intermediate intermediate/place_seqs

#add -n to get nsti 
#hsp.py -i EC -t out.tre -o EC_predicted.tsv.gz -p 1 -n 
#hsp.py -i KO -t out.tre -o KO_predicted.tsv.gz -p 1 -n
#hsp.py -i 16S -t out.tre -o marker_predicted_and_nsti.tsv.gz -p 1 -n

#edit the .tsv so it looks like this
#OTU ID  A10_S10 A11_S11 A12_S12 A2_S2   A3_S3   A4_S4   A5_S5   A6_S6   A7_S7   A8_S8   A9_S9   B10_S22 B11_S23 B12_S24 B1_S13  B2_S14  B3_S15  B4_S16
#ASV_2   18953    9350    6892    8573    1723   10514   15760   13914   11931   10381   12470   10306   19810   13322    3863    4441    8854   15720
#ASV_4      0       7      39       4       3       0       0       1       0       0       0       0       0       0      17       5       0       0
#ASV_7   6279       0       0       0       0       0       0       0       0       0       0       0       0     360       0       0       0       0
#then convert to biom
#you'll need to choose a utf
#LC_ALL=aa_DJ.utf8 biom convert -i ASVs_counts.all.rare.trim.tsv -o converted_table.biom --to-hdf5
#make sure the biom looks good
#biom head -i converted_table.biom
# Constructed from biom file
#OTU ID A10_S10 A11_S11 A12_S12 A2_S2   A3_S3
#ASV_2   18953.0 9350.0  6892.0  8573.0  1723.0
#ASV_4   0.0     7.0     39.0    4.0     3.0
#ASV_7   6279.0  0.0     0.0     0.0     0.0
#ASV_8   0.0     1.0     15.0    3.0     0.0
#ASV_9   3.0     2.0     11.0    13.0    3.0
metagenome_pipeline.py -i converted_table.biom -m marker_predicted_and_nsti.tsv.gz -f EC_predicted.tsv.gz -o EC_metagenome_out --strat_out
metagenome_pipeline.py -i converted_table.biom -m marker_predicted_and_nsti.tsv.gz -f KO_predicted.tsv.gz -o KO_metagenome_out --strat_out
#can only do tis for EC
#this last one has to be run as a job
pathway_pipeline.py -i EC_metagenome_out/pred_metagenome_contrib.tsv.gz \
                    -o EC_pathways_out -p 1

#swith to aldex.R

#I didn't do this but you should otherwise a bunch will just be called "pathways"

add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
                    -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
                    -o pathways_out/path_abun_unstrat_descrip.tsv.gz

#run aldex

module load R
Rscript ./aldex2.R

#run through post-processing aldex script 
