#!/bin/sh -l
#PBS -M user@purdue.edu
#PBS -m e
#PBS -q standby 
#PBS -l naccesspolicy=singlejob
#PBS -l nodes=1:ppn=20 
#PBS -l walltime=96:00:00 
#PBS -N pcc

cd $PBS_O_WORKDIR

perl scripts/calc_pearsons_correlation.pl -i rlog_transformed.matrix -o rlog_pcc
