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
import multiprocessing


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

t = time.time()
tmpdir = '/tmp/' + 'COEXP' +str(t)

if os.path.exists(tmpdir) == False:
    cmd = 'mkdir ' + tmpdir
    runCMD(cmd)
    print('Creating output directory:', tmpdir, '...\n')
else:
    print('Writing to existing output directory:', tmpdir, '...\n')
    
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
print("Calculating PCCs...") 
def fork_pcc(geneList):
    for gid1 in geneList:
        pccs = {}
        for gid2 in range(0, numgenes):
            corr,_ = pearsonr(expDict[gid1], expDict[gid2])
            #print(gene1,gene2,corr)
            pccs[gid2] = corr
        sorted_d = sorted(pccs.items(), key=operator.itemgetter(1), reverse = True)

        outfile = tmpdir + '/' + str(gid1)
        fo = open(outfile, 'w')
        for g2c in sorted_d:
            gid2 = g2c[0]
            corr = g2c[1]
            line = str(gid2) + '\t' + str(corr) + '\n'
            fo.write(line)
            
        fo.close()



for i in range(threads):
    try:
        pid = os.fork()
    except OSError:
        sys.stderr.write("Could not create a child process\n")
        continue
    
    if pid == 0:
        tct = i + 1
        geneList = list(subDict[tct])
        print("  thread", tct, ": Calculating PCCs for", len(geneList), 'genes...') 
        fork_pcc(geneList)
        exit()

for i in range(threads):
    finished = os.waitpid(0, 0)
print("Finished calculating PCCs!\n\n") 


# read in pccs and store pccs and ranks in np arrays
files = glob.glob(tmpdir + '/*')

print('Calculating ranks...\n')
for file in files:
    #print(infile)
    gid1 = file.split('/')[-1]
    #print("\tReading", gid1)
    rank = 0
    
    fi = open(file)
    for line in fi:
        [gid2, pcc] = line.rstrip().split('\t')
        gid1 = int(gid1)
        gid2 = int(gid2)
        pcc = np.round(float(pcc),3)
        
        mra[gid1][gid2] = pcc
        mra[gid2][gid1] = pcc
                    
        ranka[gid1][gid2] = rank
        rank += 1
    fi.close()

print('Creating mutual ranks...\n')
#calculate MR
for gid1 in range(0, numgenes):
    for gid2 in range(0, numgenes):
        if gid2 > gid1:
            runMR(gid1,gid2)

# write outfiles
print("Writing final correlations to file...") 
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


for i in range(threads):
    try:
        pid = os.fork()
    except OSError:
        sys.stderr.write("Could not create a child process\n")
        continue
    
    if pid == 0:
        tct = i + 1
        geneList = list(subDict[tct])
        print("  thread", tct, ": Writing correlations for", len(geneList), 'genes...') 
        fork_out(geneList)
        exit()

for i in range(threads):
    finished = os.waitpid(0, 0)
print("Finished writing final correlations!\n\n") 


