# Mutual Ranks and Modules
Set of scripts to identify co-expressed gene sets (i.e., modules) in gene co-expression networks. 

See [Wisecaver et al. 2017 Plant Cell](https://www.ncbi.nlm.nih.gov/pubmed/28408660) | [PDF](https://static1.squarespace.com/static/59c96b9a51a584c476f1f6f1/t/59dc24f4a9db09b4a109ae77/1507599609212/Plant+Cell+2017+Wisecaver.pdf)

## Steps
1. **Filter gene expression matrix.** If necessary, remove genes that are not expressed using the prefilter_matrix.pl perl script. Requires [Math::Round](https://metacpan.org/pod/Math::Round) which you may need to install.
```
perl scripts/prefilter_matrix.pl -i example/example_matrix.txt -o example/example_matrix_filtered.txt 
```

2. **Transform the filtered raw counts into variance stabilized abundance estimates.** Requires several R libraries, which you may need to install. Modified from this [DESeq2 vignette](http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#pre-filtering-the-dataset)
```
Rscript scripts/transform_counts.R example/example_matrix_filtered.txt example/example_conditions.txt
```
3. **Calculate Pearson's correlation (PCC) for all gene pairs.** *But first, decide which transformed matrix (vst or rlog) to use. Look at the resulting PDFs and decide which matrix is best for your dataset.* This step can be multithreaded using -t .
```
perl scripts/calc_pearsons_correlation.pl -i rlog_transformed.matrix -o rlog_pcc 
```

4. **Transform PCCs into Mutual Ranks (MRs) and MRs into edge weights.** MRs are transformed to network edge weights using the geometric decay function <a href="https://www.codecogs.com/eqnedit.php?latex=\fn_phv&space;e^{-(MR-1/x)}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\fn_phv&space;e^{-(MR-1/x)}" title="e^{-(MR-1/x)}" /></a>. The output contains results from five different decay functions with x set to 5, 10, 25, 50, and 100, respectively. Only edges greater than set thresholds will be included in the output file. The user can specify different PCC and edge weight thresholds using -c and -w, respectively. 
```
perl scripts/calc_mutual_rank.pl -i rlog_pcc -o rlog_mutual_ranks.txt
```

5. **Run clusterONE and call co-expressed gene modules**. The user must download [ClusterONE](http://www.paccanarolab.org/cluster-one/) and include the path to the ClusterONE jar file using -c . The user must also specify which decay function to use to call modules (either 5, 10, 25, 50, or 100). Optionally, the user can specify P-value and Quality score cutoffs to exclude low scoring modules from the final output. 
```
perl scripts/calc_clusterone_modules.pl -i rlog_mutual_ranks.txt -c $CLUSTERONE -d 5 -p 0.1 -q 0.1
```
