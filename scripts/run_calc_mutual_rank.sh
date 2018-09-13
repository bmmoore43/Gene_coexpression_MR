#!/bin/sh -l
#PBS -M jwisecav@gmail.com
#PBS -m e
#PBS -q jwisecav 
#PBS -l naccesspolicy=singlejob
#PBS -l nodes=1:ppn=1
#PBS -l walltime=48:00:00 
#PBS -N pcc2mr
#PBS -l epilogue=scripts/epilogue.sh

cd $PBS_O_WORKDIR

perl scripts/calc_mutual_rank.pl -i rlog_pcc -o rlog_mutual_ranks.txt
