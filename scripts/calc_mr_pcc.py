#!/usr/bin/env python
import os
import os.path
import subprocess
import getopt
import sys
import math
from scipy.stats import pearsonr
import operator
import time
import glob
import numpy as np

def runCMD(cmd):
    val = subprocess.Popen(cmd, shell=True).wait()
    if val == 0:
        pass
    else:
        print ('command failed')
        print (cmd)
        sys.exit(1)


def runMR(id1,id2):
    mranks = ranka[id1][id2] * ranka[id2][id1]
    mr = np.round(np.sqrt(mranks), 3)

    mra[id1][id2] = mr


def usage():
    print('\n  Usage: '+sys.argv[0]+' -i <input.matrix> -o <output directory> [-t <threads> -k <skipfile>]')
    print("    -i|--infile  <FILENAME>  INFILE: path to input file: matrix of gene counts")
    print("    -o|--outdir  <DIRECTORY> OUTDIR: path to output directory")
    print("    -t|--threads <INTEGER>   THREADS: multithreaded using perl fork (default = 1)")
    print("    -k|--skip    <FILENAME>  SKIPFILE: path to file containing list of gene to skip\n\n")
  
# Set default variables
threads = 1

# Read in command line arguments
try:
    options, remainder = getopt.getopt(sys.argv[1:], 'hi:o:t:k:', ['help', 'infile=', 'outdir=', 'threads=', 'skip='])
except getopt.GetoptError as err:
    # print help information and exit:
    print(err)  # will print something like "option -a not recognized"
    usage()
    sys.exit(2)

for opt, arg in options:
    if opt in ('-h', '--help'):
        usage()
        sys.exit()
    elif opt in ('-o', '--outdir'):
        outdir = arg
    elif opt in ('-i', '--infile'):
        infile = arg
    elif opt in ('-t', '--threads'):
        threads = int(arg)
    elif opt in ('-k', '--skip'):
        skipfile = arg
command = " ".join(sys.argv)

try:
    infile
    outdir
except NameError:
    usage()
    sys.exit(2)
else:
    print("\nCOMMAND:", command, '\n')

if os.path.exists(outdir) == False:
    cmd = 'mkdir ' + outdir
    runCMD(cmd)
    print('Creating output directory:', outdir, '...\n')
else:
    print('Writing to existing output directory:', outdir, '...\n')


# Read in matrix file
expDict = {}
header = 0

fi = open(infile)
print('Reading in gene expression values from', infile, '...\n')
for line in fi:
    if header == 0:
        header += 1
        continue
    
    exp = line.rstrip().split('\t')
    gene = exp.pop(0)
    exp = [float(i) for i in exp]
    #print(gene, len(col))

    expDict[gene] = exp
    
fi.close()

try:
    skipfile
except NameError:
    print('Skipping 0 genes...\n')
else:
    fi = open(skipfile)
    skipcount = 0
    for gene in fi:
        gene = gene.rstrip()
        if gene in expDict:
            del expDict[gene]
            skipcount += 1
    print('Skipping', skipcount, 'genes...\n')

    fi.close()

print('Creating gene lookup tables...\n')
idlookup = {}
genelookup = {}
gid = 0
for gene in expDict:
    #print(gene,gid)
    idlookup[gene] = gid
    genelookup[gid] = gene
    gid = gid + 1 

for gene in idlookup:
    gid = idlookup[gene]
    expDict[gid] = expDict.pop(gene)

#print(expDict.keys())

    
    
# Split gene expression vectors into n=threads subgroups 
subDict = {}
num_genes = len(expDict)
subcount = math.ceil(num_genes/threads)

gct = 0
tct = 1
for gid in expDict:
    gct += 1 
    
    if gct == 1:
        subDict[tct] = []
    
    if gct > subcount:
        gct = 1
        tct += 1
        subDict[tct] = []
        
    subDict[tct].append(gid)

sum = 0
for tct in subDict:
    sublen = len(list(subDict[tct]))
    sum = sum + sublen
if num_genes != sum:
    print('WARNING')
    print('total genes:', num_genes)
    print('total split genes:', sum)
    sys.exit(2)
        
numgenes = len(expDict)
mra = np.zeros(shape=(numgenes,numgenes))
ranka = np.zeros(shape=(numgenes,numgenes))

    
# Define function for calculating PCC and MR
def fork_pcc(mra,geneList):
    for gid1 in geneList:
        for gid2 in range(0, numgenes):
            corr,_ = pearsonr(expDict[gid1], expDict[gid2])
            #print(gene1,gene2,corr)
            mra[gid1][gid2] = corr
            mra[gid2][gid1] = corr

def fork_rank(geneList):
    for gid1 in geneList:
        count = 0
        d = dict(enumerate(mra[gid1,]))
        sorted_d = sorted(d.items(), key=operator.itemgetter(1), reverse = True)
    
        for gc in sorted_d:
            gid2 = gc[0]
            corr = gc[1]    
            ranka[gid1][gid2] = count
            count += 1

def fork_mr(geneList):
    for gid1 in geneList:
        for gid2 in range(0, numgenes):
            if gid2 > gid1:
                runMR(gid1,gid2)

def fork_out(geneList):
    for gid1 in geneList:
        mrDict = {}
        pccDict = {}

        for gid2 in range(0, numgenes):
            if gid1 == gid2:
                continue

            mr = ''
            pcc = ''
            if gid2 > gid1:
                mr = mra[gid1][gid2]
                pcc = mra[gid2][gid1]        
            elif gid1 > gid2:
                pcc = mra[gid1][gid2]
                mr = mra[gid2][gid1]
        
            mrDict[gid2] = mr
            pccDict[gid2] = pcc
    
        gene = genelookup[gid1]
        outfile = outdir + '/' + gene
        fo = open(outfile, 'w')

        sorted_d = sorted(mrDict.items(), key=operator.itemgetter(1), reverse = False)
        for gm in sorted_d:
            [gid2, mr] = gm
            gene2 = genelookup[int(gid2)]

            line = gene2 + '\t' + str(np.round(mr,3)) + '\t' + str(pccDict[gid2]) + '\n'
            fo.write(line)
    
        fo.close()

from multiprocessing import Process, Array


# Create, start, and finish the child process
geneList = list(subDict[1])
p = Process(target=fork_pcc, args=(mra,geneList))
p.start()
p.join()


# geneList = list(subDict[1])
# print(geneList)
# 
# forks = threads
# for i in range(forks):
#     try:
#         pid = os.fork()
#     except OSError:
#         sys.stderr.write("Could not create a child process\n")
#         continue
#     if pid == 0:
#         tct = i + 1
#         geneList = list(subDict[tct])
#         print(geneList)
#         print("  thread", tct, ": Calculating PCCs for", len(geneList), 'genes...') 
#         fork_pcc(geneList)
#         exit()
# for i in range(forks):
#     os.waitpid(0, 0)
#     
# print("Finished calculating PCCs!\n\n") 

#fork_pcc(geneList)
fork_rank(geneList)
fork_mr(geneList)
fork_out(geneList)



# 
# 
# for i in range(forks):
#     try:
#         pid = os.fork()
#     except OSError:
#         sys.stderr.write("Could not create a child process\n")
#         continue
#     if pid == 0:
#         tct = i + 1
#         geneList = list(subDict[tct])
#         print("  thread", tct, ": Calculating ranks for", len(geneList), 'genes...') 
#         fork_rank(geneList)
#         exit()
# for i in range(forks):
#     finished = os.waitpid(0, 0)
# print("Finished calculating ranks!\n\n") 
# 
# 
# for i in range(forks):
#     try:
#         pid = os.fork()
#     except OSError:
#         sys.stderr.write("Could not create a child process\n")
#         continue
#     if pid == 0:
#         tct = i + 1
#         geneList = list(subDict[tct])
#         print("  thread", tct, ": Calculating MRs for", len(geneList), 'genes...') 
#         fork_mr(geneList)
#         exit()
# for i in range(forks):
#     finished = os.waitpid(0, 0)
# print("Finished calculating MRs!\n\n") 
# 
# 
# 
# for i in range(forks):
#     try:
#         pid = os.fork()
#     except OSError:
#         sys.stderr.write("Could not create a child process\n")
#         continue
#     if pid == 0:
#         tct = i + 1
#         geneList = list(subDict[tct])
#         print("  thread", tct, ": Writing outfiles for", len(geneList), 'genes...') 
#         fork_out(geneList)
#         exit()
# for i in range(forks):
#     finished = os.waitpid(0, 0)
# print("Finished writing outfiles!\n\n") 
# 



