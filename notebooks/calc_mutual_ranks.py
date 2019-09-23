#!/usr/bin/env python
import os
import os.path
import sys
import subprocess
import getopt
import glob
import operator
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
    mranks = bigDict[id1][id2] * bigDict[id2][id1]
    mr = np.round(np.sqrt(mranks), 3)

    if id1 not in mrDict:
        mrDict[id1] = {}
    mrDict[id1][id2] = mr

    if id2 not in mrDict:
        mrDict[id2] = {}
    mrDict[id2][id1] = mr

def usage():
    print('\n  Usage: '+sys.argv[0]+' -i <pcc directory> -o <mr directory> [-k <skipfile> -c <minimum PCC>]')
    print("    -i|--indir  <DIRECTORY> path to directory of gene PCCs")
    print("    -o|--outdir <DIRECTORY> path to new directory to store gene mutual ranks")
    print("    -k|--skip   <file> path to file containing list of genes to exclude from mutual rank calculation")
    print("    -c|--minpcc    <float> minimum pearson's correlation coefficient between gene pairs (default = -1)\n\n")
    
    
# Read in command line arguments
try:
    options, remainder = getopt.getopt(sys.argv[1:], 'hi:o:k:c:', ['help', 'indir=', 'outdir=', 'skip=', 'minpcc='])
except getopt.GetoptError as err:
    # print help information and exit:
    print(err)  # will print something like "option -a not recognized"
    usage()
    sys.exit(2)

min_pcc = -1
for opt, arg in options:
    if opt in ('-h', '--help'):
        usage()
        sys.exit()
    elif opt in ('-o', '--outdir'):
        outdir = arg
    elif opt in ('-i', '--indir'):
        indir = arg
    elif opt in ('-k', '--skip'):
        skipfile = arg
    elif opt in ('-c', '--minpcc'):
        min_pcc = float(arg)
command = " ".join(sys.argv)

try:
    indir
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

skipList = []
try:
    skipfile
except NameError:
    print('Skipping 0 genes...')
else:
    fi = open(skipfile)
    for gene in fi:
        gene = gene.rstrip()
        skipList.append(gene)
    fi.close()
    print('Skipping', len(skipList), 'genes...')


mrDict = {}
bigDict = {}
idlookup = {}
genelookup = {}
int_count = 0
files = glob.glob(indir + '/*')

for infile in files:
    #print(infile)
    gene = infile.split('/')[-1]
    #print(gene)
    
    int_count = int_count + 1 
    idlookup[gene] = int_count
    genelookup[int_count] = gene


for infile in files:
    #print(infile)
    gene1 = infile.split('/')[-1]
    print("\tReading", gene1)
    id1 = idlookup[gene1]
    count = 0
    if gene1 in skipList:
        #print('\tskipping', gene1)
        continue


    if id1 not in bigDict:
        bigDict[id1] = {}

    
    fi = open(infile)
    for line in fi:
        [gene2, pcc] = line.rstrip().split('\t')
        if gene1 == gene2:
            #print('\tskipping', gene1,gene2)
            continue
        if gene2 in skipList:
            #print('\tskipping', gene2)
            continue

        id2 = idlookup[gene2]
        count += 1 
        bigDict[id1][id2] = count
        
        if id2 not in bigDict:
            continue
        elif id1 in bigDict[id2]:
            runMR(id1,id2)
            del bigDict[id2][id1]
            del bigDict[id1][id2]
            #print('pair found, calculated MR and removed from dictionary')
            continue 

    fi.close()
print('\nREADING COMPLETE!\n\nCalculating MR scores now...\n')

for gid in mrDict:
    gene = genelookup[gid]
    
    pccDict = {}
    pccfile = indir + '/' + gene
    fi = open(pccfile)
    for line in fi:
        [gene2, pcc] = line.rstrip().split('\t')
        pccDict[gene2] = np.round(float(pcc), 3)
    fi.close()
    
    outfile = outdir + '/' + gene
    fo = open(outfile, 'w')
    
    sorted_d = sorted(mrDict[gid].items(), key=operator.itemgetter(1), reverse = False)
    #print(sorted_d)
    for id2mr in sorted_d:
        id2 = id2mr[0]
        mr = id2mr[1]
        gene2 = genelookup[id2]
        line = gene2 + '\t' + str(mr) + '\t' + str(pccDict[gene2]) + '\n'
        fo.write(line)
    
    fo.close()


