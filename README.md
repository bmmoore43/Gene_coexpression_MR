# mr2mods: Mutual Ranks to Modules workflow
## [Jen Wisecaver](https://www.wisecaverlab.com)
#### 2019-09-25
Set of scripts to identify co-expressed gene sets (i.e., modules) in gene co-expression networks. 

Transforming Pearson’s correlation coefficents (PCCs) into Mutual Ranks (MRs) — first described by [Obayashi & Kinoshita](https://www.ncbi.nlm.nih.gov/pubmed/19767600) — is a good idea if you want to compare between different datasets and/or functional categories ([Obayashi et al, 2018](https://www.ncbi.nlm.nih.gov/pubmed/29216398); [Liesecke et al, 2018](https://www.ncbi.nlm.nih.gov/pubmed/30022075)). However, the transformation requires significant memory and disk space to compute. Moreover, MRs range from 1 to n-1 (where n is the number of genes in the genome), which does not translate well to network edge weights. We’ve implemented a series of R and Python scripts for creating MRs and gene modules directly from a directory of Kallisto gene abundance estimates. The code as been developed to run on multiple threads when possible, which significantly decreasing the total runtime. Finally, the pipeline applies exponential decay functions to transform MRs into edge weights and call co-expressed gene sets (i.e. modules) using the program [ClusterONE](https://www.paccanarolab.org/cluster-one/). 

See [Wisecaver et al. 2017 Plant Cell](https://www.ncbi.nlm.nih.gov/pubmed/28408660) | [PDF](https://static1.squarespace.com/static/59c96b9a51a584c476f1f6f1/t/59dc24f4a9db09b4a109ae77/1507599609212/Plant+Cell+2017+Wisecaver.pdf)

## Steps
An in depth [tutorial](https://github.rcac.purdue.edu/jwisecav/coexp-pipe/blob/master/tutorial/mutual_ranks_to_modules.ipynb) is provided and includes information on files types, run time, memory requirements, etc. The TL:DR steps are provided below. 

1. **Create a gene expression matrix** 
```
Rscript ../scripts/transform_counts.r -a example/abundance_files.txt -l example/transcripts2genes.txt -s example/sample_conditions.txt -t tag -o example/gene_counts 
```


2. **Calculate PCC and MR for all gene pairs**
```
python ../scripts/calc_mr_pcc.py -i example/gene_counts_normalized.matrix -o example/gene_counts_normalized_mr -t 20
```


3. **Run clusterONE and call co-expressed gene modules** 
```
python ../scripts/create_network_and_modules.py -i example/gene_counts_normalized_mr -c ../scripts/cluster_one-1.0.jar -d 5
```
