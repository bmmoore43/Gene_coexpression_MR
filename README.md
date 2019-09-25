# Mutual Ranks and Modules
Set of scripts to identify co-expressed gene sets (i.e., modules) in gene co-expression networks. 

See [Wisecaver et al. 2017 Plant Cell](https://www.ncbi.nlm.nih.gov/pubmed/28408660) | [PDF](https://static1.squarespace.com/static/59c96b9a51a584c476f1f6f1/t/59dc24f4a9db09b4a109ae77/1507599609212/Plant+Cell+2017+Wisecaver.pdf)

## Steps
See the [tutorial](https://github.rcac.purdue.edu/jwisecav/coexp-pipe/blob/master/tutorial/mutual_ranks_to_modules.ipynb) for more information on files types, runtime, memory requirements, etc.  

1. **Create a gene expression matrix.** 
```
Rscript ../scripts/transform_counts.r -a example/abundance_files.txt -l example/transcripts2genes.txt -s example/sample_conditions.txt -t tag -o example/gene_counts 
```


2. **Calculate Pearson's correlation coefficient (PCC) Mutual Rank (MR)for all gene pairs.**
```
python ../scripts/calc_mr_pcc.py -i example/gene_counts_normalized.matrix -o example/gene_counts_normalized_mr -t 20
```


3. **Run clusterONE and call co-expressed gene modules**. 
```
python ../scripts/create_network_and_modules.py -i example/gene_counts_normalized_mr -c ../scripts/cluster_one-1.0.jar -d 5
```
