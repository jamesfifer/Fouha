#!/bin/bash -l
# #$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N Blast
#$ -pe omp 28
#$ -l mem_per_core=4G
#$ -l h_rt=12:00:00 
#$ -P incrna 
#$ -M james.e.fifer@gmail.com
#$ -t 1-841:100 #I'd run this 1-10:1 first as a test and then run for the full number of .fa files

##First run the following command on fasta file to split fasta file into many .fa files
#awk '/^>/ {if(x>0) close(outname); x++; outname=sprintf("_%d.fa",x); print > outname;next;} {if(x>0) print >> outname;}' ./foo.fasta

module load blast+/2.7.1 

nthreads=4 # how many threads to use for each indivdual execution of blastx

cntrtop=6 # $(( $NSLOTS / $nthreads )) # max number of concurrent executions (integer div of ompCores/num_threads)
cntr=0


for (( i=$SGE_TASK_ID; i<=$SGE_TASK_ID+540; i+=1 )); do
    # if modulo is not equal to 0
    if [[ $(( cntr % cntrtop )) -ne 0 ]]; then
            cntr=$(( cntr + 1 ))
            blastn -query /projectnb/davieslab/jfifer/Fouha/Analysis_16S/Blast/_"$i".fa -db /projectnb/incrna/ncbi/projectnb/incrna/ncbi/NT/nt -outfmt 6 -num_threads $nthreads -evalue 1e-5 -max_target_seqs 10 -out /projectnb/davieslab/jfifer/Fouha/Analysis_16S/Blast/"$i".out &
            # cheap hack to have multiple background jobs,
    else
            cntr=$(( cntr + 1 ))
            nice blastn -query /projectnb/davieslab/jfifer/Fouha/Analysis_16S/Blast/_"$i".fa -db /projectnb/incrna/ncbi/projectnb/incrna/ncbi/NT/nt -outfmt 6 -num_threads $nthreads -evalue 1e-5 -max_target_seqs 10 -out /projectnb/davieslab/jfifer/Fouha/Analysis_16S/Blast/"$i".out
    fi
done
wait


