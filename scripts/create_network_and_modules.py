#!/usr/bin/env python
import math
import time
import os
import os.path
import sys
import getopt
import glob
import subprocess
import numpy as np

def runCMD(cmd):
    val = subprocess.Popen(cmd, shell=True).wait()
    if val == 0:
        pass
    else:
        print ('command failed')
        print (cmd)
        sys.exit(1)


def usage():
    print('\nUsage: '+sys.argv[0]+' -i <mr directory> -c <cluster one jar> -d <decay> [-o <output base name> -w <min edge weight> -p <max pvalue> -q <min quality>')
    print("    -i|--indir       <DIRECTORY> path to directory of gene mutual ranks")
    print("    -c|--clusterone  <FILE> path to clusterone jar file")
    print("    -d|--decay       <INTEGER> decay constant to use for calculating network edge weights")
    print("\n    OPTIONAL:")
    print("    -o|--out         <NAME> base name for all output files (default = same as input directory")
    print("    -w|--minweight   <FLOAT> minimum network edge weight (default = 0.01)")
    print("    -p|--maxpval     <FLOAT> retained modules must have p value less than or equal to this value (default = 1 to retain all modules)")
    print("    -q|--minqual     <FLOAT> retained modules must have weight greater than or equal to this value (default = 0 to retain all modules)")
    print("    -t|--threads     <INTEGER> number of forks to create when multithreading using os.fork (default = 1)")
    print("    -x|--prefix      <NAME> prefix to append to all module names; helpful when comparing between different species/filtering options")


# Read in command line arguments
try:
    options, remainder = getopt.getopt(sys.argv[1:], 'hi:o:c:d:p:q:w:t:x:', ['help', 'indir=', 'out=', 'clusterone=', 'minweight=', 'decay=', 'maxpval=', 'minqual=', 'threads=', 'prefix='])
except getopt.GetoptError as err:
    # print help information and exit:
    print(err)  # will print something like "option -a not recognized"
    usage()
    sys.exit(2)

threads = 1
maxpval = 1
minqual = 0
minweight = 0.01
prefix = ''

for opt, arg in options:
    if opt in ('-h', '--help'):
        usage()
        sys.exit()
    elif opt in ('-i', '--indir'):
        indir = arg
    elif opt in ('-c', '--clusterone'):
        clusterone = arg
    elif opt in ('-d', '--decay'):
        decay = int(arg)
    elif opt in ('-o', '--out'):
        basename = arg
    elif opt in ('-w', '--minweight'):
        minweight = float(arg)
    elif opt in ('-p', '--maxpval'):
        maxpval = float(arg)
    elif opt in ('-q', '--minqual'):
        minqual = float(arg)
    elif opt in ('-t', '--threads'):
        threads = int(arg)
    elif opt in ('-x', '--prefix'):
        prefix = arg
command = " ".join(sys.argv)

try:
    indir
    clusterone
    decay
except NameError:
    usage()
    sys.exit(2)
else:
    if decay < 1 or decay > 999:
        print('--decay must be an integer from 1 to 999')
        usage()
        sys.exit(2)
    print("\nCOMMAND:", command, '\n')

t = time.time()
tmpdir = '/tmp/' + 'ABC' +str(t)

if os.path.exists(tmpdir) == False:
    cmd = 'mkdir ' + tmpdir
    runCMD(cmd)
    print('Creating output directory:', tmpdir, '...\n')
else:
    print('Writing to existing output directory:', tmpdir, '...\n')

try:
    basename
except NameError:
    basename = indir.split('/')[-1]

abcfile = basename + '_' + f"{decay:03d}" + '.abc'
csvfile = basename + '_' + f"{decay:03d}" + '.modules.csv'
outfile = basename + '_' + f"{decay:03d}" + '.modules.txt'

print('Writing network abc file to', abcfile)
print('Writing clusterone csv file to', csvfile)
print('Writing parsed clusterone tab-delimited file to', outfile)

files = glob.glob(indir + '/*')
subDict = {}
num_genes = len(files)
subcount = math.ceil(num_genes/threads)

gct = 0
tct = 1
for file in files:
    gct += 1 
    
    if gct == 1:
        subDict[tct] = []
    
    if gct > subcount:
        gct = 1
        tct += 1
        subDict[tct] = []
        
    subDict[tct].append(file)

sum = 0


for tct in subDict:
    sublen = len(list(subDict[tct]))
    sum = sum + sublen
if num_genes != sum:
    print('ERROR:')
    print('total genes:', num_genes)
    print('total split genes:', sum)
    sys.exit(2)


def fork_abc(tct):
    fileList = subDict[tct]
    tmpfile = tmpdir + '/' + str(tct) + '.abc'
    fo = open(tmpfile, 'w')
    for infile in fileList:
        #print(infile)
        gene1 = infile.split('/')[-1]
    
        fi = open(infile)
        for line in fi:
            [gene2, mr, pcc] = line.rstrip().split('\t')
            weight = np.exp(-1 * (float(mr) - 1) / decay)
        
            if weight < minweight:
                #print('skipping', gene1,gene2,weight)
                break
            else:
                fo.write(gene1 + '\t' + gene2 + '\t' + str(weight) + '\n')
        fi.close()
    fo.close()

print('\nParsing', indir, '...')
for i in range(threads):
    try:
        pid = os.fork()
    except OSError:
        sys.stderr.write("Could not create a child process\n")
        continue
    
    if pid == 0:
        tct = i + 1
        geneList = list(subDict[tct])
        print("  thread", tct, ": Parsing", len(geneList), 'genes...') 
        fork_abc(tct)
        exit()

for i in range(threads):
    finished = os.waitpid(0, 0)
print("\tFinished!\n\n") 


command = 'cat ' + tmpdir + '/*.abc > ' + abcfile
print('\nCombining temporary files:\n', command, '\n')
runCMD(command)

command = 'java -jar ' + clusterone + ' ' + abcfile + ' --output-format csv > ' + csvfile
print('\nRunning clusterone:\n', command, '\n')
runCMD(command)


fi = open(csvfile)
fo = open(outfile, 'w')

print('\nParsing', csvfile, '...')
for line in fi:
    line = line.rstrip().split(',')
    if line[0] == 'Cluster':
        continue
    
    mod = line[0]
    mod = prefix + 'N'+ f"{decay:03d}" + 'M' + f"{int(mod):05d}"
    #print(mod)

    qual = line[5]
    pval = line[6]
    if float(qual) < minqual:
        continue
    if float(pval) > maxpval:
        continue

    genes = line[7]
    genes = genes.replace('"', '')
    fo.write(mod + '\t' + qual + '\t' + pval + '\t' + genes + '\n')
    
fi.close()
fo.close()
print('\tfinished!\n')

