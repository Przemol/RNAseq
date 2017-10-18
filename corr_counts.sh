#!/bin/bash                               
#$ -l rmem=48G                      
#$ -M n.bulgakova@sheffield.ac.uk
#$ -m bea

module load compilers/gcc/4.8.2

module load apps/R/3.3.1               

R CMD BATCH correlate_counts_qsub.R correlate_counts_qsub.R.o$JOB_ID
