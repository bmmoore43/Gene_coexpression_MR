#!/usr/bin/env python
import os
import os.path
import subprocess
import getopt
import sys
import math
from numpy import mean
from numpy import std
from scipy.stats import pearsonr
import operator

def runCMD(cmd):
    val = subprocess.Popen(cmd, shell=True).wait()
    if val == 0:
        pass
    else:
        print ('command failed')
        print (cmd)
        sys.exit(1)

def usage():
    print('\n  Usage: '+sys.argv[0]+' -i <input.matrix> -o <output directory> [-t <threads>]')
    print("    -i: <FILENAME> INFILE: path to input file: matrix of gene counts")
    print("    -o: <DIRECTORY> OUTDIR: path to output directory. Default creates output directory 'pcc' in working directory.")
    print("    -t: <INTEGER> THREADS: multithreaded using perl fork (default = 1)\n\n")
  
# Set default variables
output_dirname = 'pcc'
threads = 1

# Read in command line arguments
try:
    options, remainder = getopt.getopt(sys.argv[1:], 'hi:o:t:', ['help', 'infile=', 'outdir=', 'threads=',])
except getopt.GetoptError as err:
    # print help information and exit:
    print(err)  # will print something like "option -a not recognized"
    usage()
    sys.exit(2)

for opt, arg in options:
    if opt in ('-h', '--help'):
        usage()
        sys.exit()
    elif opt in ('-o', '--output'):
        output_dirname = arg
    elif opt in ('-i', '--infile'):
        infile = arg
    elif opt in ('-t', '--threads'):
        threads = int(arg)
command = " ".join(sys.argv)
print("\nCOMMAND:", command, '\n')

if os.path.exists(output_dirname) == 'false':
    cmd = 'mkdir ' + output_dirname
    runCMD(cmd)
    print('Creating output directory:', output_dirname, '...\n')
else:
    print('Writing to existing output directory:', output_dirname, '...\n')


# Read in matrix file
fi = open(infile)
print('Reading in gene expression values from', infile, '...\n')

geneDict = {}
skip = 0
for line in fi:
    if skip == 0:
        skip += 1
        continue
    
    col = line.rstrip().split('\t')
    gene = col.pop(0)
    col = [float(i) for i in col]
    #print(gene, len(col))
    
    geneDict[gene] = col
    
fi.close()

# Split gene expression vectors into n=threads subgroups 
subgeneDict = {}
num_genes = len(geneDict)
subcount = math.ceil(num_genes/threads)

gct = 0
tct = 1
for gene in geneDict:
    gct += 1 
    
    if gct == 1:
        subgeneDict[tct] = []
    
    if gct > subcount:
        gct = 0
        tct += 1
        subgeneDict[tct] = []
        
    subgeneDict[tct].append(gene)

# Define function for calculating PCCs
def sub_fork(subgeneList):
    
    for gene1 in subgeneList:
        pccs = {}
        for gene2 in geneDict:
            corr,_ = pearsonr(geneDict[gene1], geneDict[gene2])
            #print(gene1,gene2,corr)
            pccs[gene2] = corr
        sorted_d = sorted(pccs.items(), key=operator.itemgetter(1), reverse = True)
        
        outfile = output_dirname + '/' + gene1
        fo = open(outfile, 'w')
        for gene2corr in sorted_d:
            gene2 = gene2corr[0]
            corr = gene2corr[1]
            line = gene2 + '\t' + str(corr) + '\n'
            fo.write(line)
            
        fo.close()


forks = threads
for i in range(forks):
    try:
        pid = os.fork()
    except OSError:
        sys.stderr.write("Could not create a child process\n")
        continue
    
    if pid == 0:
        tct = i + 1
        subgeneList = list(subgeneDict[tct])
        print("  thread", tct, ": Calculating PCCs for", len(subgeneList), 'genes...') 
        sub_fork(subgeneList)
        exit()

for i in range(forks):
    tct = i + 1
    finished = os.waitpid(0, 0)
    print("  thread", tct, ": Finished calculating PCCs!") 

