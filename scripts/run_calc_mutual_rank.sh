#!/bin/sh -l
#PBS -M user@purdue.edu
#PBS -m e
#PBS -q standby
#PBS -l naccesspolicy=singlejob
#PBS -l nodes=1:ppn=1
#PBS -l walltime=48:00:00 
#PBS -N pcc2mr

cd $PBS_O_WORKDIR

perl scripts/calc_mutual_rank.pl -i rlog_pcc -o rlog_mutual_ranks.txt
