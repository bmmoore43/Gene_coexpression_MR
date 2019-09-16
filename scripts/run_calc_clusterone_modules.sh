#!/bin/sh -l
#PBS -M user@purdue.edu
#PBS -m e
#PBS -q standby
#PBS -l naccesspolicy=singleuser
#PBS -l nodes=1:ppn=1
#PBS -l walltime=48:00:00 
#PBS -N clusterone

cd $PBS_O_WORKDIR

perl scripts/calc_clusterone_modules.pl -i rlog_mutual_ranks.txt -c $CLUSTERONE -d 5 -p 0.1 -q 0.1
perl scripts/calc_clusterone_modules.pl -i rlog_mutual_ranks.txt -c $CLUSTERONE -d 10 -p 0.1 -q 0.1
perl scripts/calc_clusterone_modules.pl -i rlog_mutual_ranks.txt -c $CLUSTERONE -d 25 -p 0.1 -q 0.1
perl scripts/calc_clusterone_modules.pl -i rlog_mutual_ranks.txt -c $CLUSTERONE -d 50 -p 0.1 -q 0.1
perl scripts/calc_clusterone_modules.pl -i rlog_mutual_ranks.txt -c $CLUSTERONE -d 100 -p 0.1 -q 0.1

