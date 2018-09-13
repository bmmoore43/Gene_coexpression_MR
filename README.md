# coexp-pipe
Steps to identify co-expressed gene sets (i.e. modules) in gene co-expression networks 

1. Filter your gene expression matrix if necessary to remove genes that are not expressed using the prefilter_matrix.pl perl script
```
module load bioinfo perl
perl scripts/prefilter_matrix.pl -i example/example_matrix.txt -o example/example_matrix_filtered.txt 
```

2. Transform the filtered raw counts into variance stabilized abundance estimates
```
module load R
Rscript scripts/transform_counts.R example/example_matrix_filtered.txt example/example_conditions.txt
```

3. Calculate Pearson's correlation for all gene pairs
```
perl scripts/calc_pearsons_correlation.pl -i rlog_transformed.matrix -o rlog_pcc
```

4. Transform Pearson's correlations into Mutual Ranks -> Transform Mutual Ranks into edge weights
```
perl scripts/calc_mutual_rank.pl -i rlog_pcc -o rlog_mutual_ranks.txt
```

5. Run clusterONE and call coexpressed modules
```
perl scripts/calc_clusterone_modules.pl -i rlog_mutual_ranks.txt -c $CLUSTERONE -d 5 -p 0.1 -q 0.1
```
