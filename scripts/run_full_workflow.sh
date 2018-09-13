#!/bin/sh -l
#PBS -M jwisecav@gmail.com
#PBS -m e
#PBS -q jwisecav 
#PBS -l naccesspolicy=singleuser
#PBS -l nodes=1:ppn=20
#PBS -l walltime=48:00:00 
#PBS -N calcexpmodules
#PBS -l epilogue=scripts/epilogue.sh

cd $PBS_O_WORKDIR

# Initialize variables
INFILE='rlog_transformed.matrix'	#matrix of vst or rlog transformed counts
PCCDIR='rlog_pcc'					#name of directory to save pearsons correlations
THREADS=20
MRFILE='rlog_mutual_ranks.txt'		#name of output file to save mutual ranks
PCC_THRESHOLD=0.3					#retained edges must have a PCC >= this value
WEIGHT_THRESHOLD=0.01				#retained edges must have a weight (transformed MR) >= this value
#CLUSTERONE=''						#path to clusterone jar file
PVAL=1								#retained modules must have p value <= this value
QUAL=0								#retained modules must have weight >= this value 


# calculate pearsons correlation of co-expression
perl scripts/calc_pearsons_correlation.pl -i $INFILE -o $PCCDIR -t $THREADS

# calculate mutual ranks and edge weights
perl scripts/calc_mutual_rank.pl -i  $PCCDIR -o $MRFILE -c $PCC_THRESHOLD -w $WEIGHT_THRESHOLD

# call modules
perl scripts/calc_clusterone_modules.pl -i $MRFILE -c $CLUSTERONE -d 5 -p $PVAL -q $QUAL
perl scripts/calc_clusterone_modules.pl -i $MRFILE -c $CLUSTERONE -d 10 -p $PVAL -q $QUAL
perl scripts/calc_clusterone_modules.pl -i $MRFILE -c $CLUSTERONE -d 25 -p $PVAL -q $QUAL
perl scripts/calc_clusterone_modules.pl -i $MRFILE -c $CLUSTERONE -d 50 -p $PVAL -q $QUAL
perl scripts/calc_clusterone_modules.pl -i $MRFILE -c $CLUSTERONE -d 100 -p $PVAL -q $QUAL

