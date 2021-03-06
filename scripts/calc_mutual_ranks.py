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

    mra[id1][id2] = mr
    mra[id2][id1] = mr
    
    
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

files = glob.glob(indir + '/*')
try:
    skipfile
except NameError:
    print('Skipping 0 genes...')
else:
    fi = open(skipfile)
    skipcount = 0
    for gene in fi:
        gene = gene.rstrip()
        infile = indir + '/' + gene
        if infile in files:
            files.remove(infile)
            skipcount += 1
    print('Skipping', skipcount, 'genes...')

    fi.close()

mrDict = {}
bigDict = {}
idlookup = {}
genelookup = {}
int_count = 0

numgenes = len(files)
mra = np.zeros(shape=(numgenes,numgenes))

for infile in files:
    #print(infile)
    gene = infile.split('/')[-1]
    #print(gene)
    
    idlookup[gene] = int_count
    genelookup[int_count] = gene
    
    int_count = int_count + 1 


for infile in files:
    #print(infile)
    gene1 = infile.split('/')[-1]
    print("\tReading", gene1)
    id1 = idlookup[gene1]
    count = 0

    if id1 not in bigDict:
        bigDict[id1] = {}

    
    fi = open(infile)
    for line in fi:
        [gene2, pcc] = line.rstrip().split('\t')
        if gene1 == gene2:
            #print('\tskipping', gene1,gene2)
            continue
        if gene2 not in idlookup:
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


for gid1 in genelookup:
    gene = genelookup[gid1]
    
    pccDict = {}
    mrDict = {}
    pccfile = indir + '/' + gene
    fi = open(pccfile)
    for line in fi:
        [gene2, pcc] = line.rstrip().split('\t')
        pccDict[gene2] = np.round(float(pcc), 3)
    fi.close()

    outfile = outdir + '/' + gene
    fo = open(outfile, 'w')

    for gid2 in genelookup:
        if gid1 == gid2:
            continue
        mrDict[gid2] = mra[gid1,gid2]
    
    sorted_d = sorted(mrDict.items(), key=operator.itemgetter(1), reverse = False)
    for gidmr in sorted_d:
        [gid2, mr] = gidmr
        gene2 = genelookup[int(gid2)]
        line = gene2 + '\t' + str(np.round(mr,3)) + '\t' + str(pccDict[gene2]) + '\n'
        fo.write(line)
    
    fo.close()


print('\nMutual ranks and PCCs written to output:', outdir, '\nInput:', indir, 'is now redundant, recommend deleting...\n')

