#!/bin/bash                               
#$ -l rmem=8G                          
#$ -M n.bulgakova@sheffield.ac.uk
#$ -m bea

cd /mnt/fastdata/md1nbu/PRJNA369618/fastq

module load apps/R/3.3.1               

R CMD BATCH count_reads_cluster.r count_reads_cluster.r.o$JOB_ID
