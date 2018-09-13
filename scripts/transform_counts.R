library("DESeq2")
library("dplyr")
library("ggplot2")
library('hexbin')
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")

## http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#pre-filtering-the-dataset
## Usage: module load R
## Usage: Rscript transform_counts.R counts.matrix conditions.txt

args <- commandArgs(trailingOnly = TRUE)
countsfile <- args[1]  
conditionsfile <- args[2]  

# header rows should have one less column than data rows
countmatrix <- read.table(countsfile, sep="\t", header=TRUE)
coldata <- read.table(conditionsfile, sep="\t", header=TRUE)


## TRANSFORM COUNTS
dds = DESeqDataSetFromMatrix(countData = countmatrix, colData = coldata, design = ~ condition)
vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)

write.table(assay(vsd), "vst_transformed.matrix", sep="\t", quote = FALSE)
write.table(assay(rld), "rlog_transformed.matrix", sep="\t", quote = FALSE)

print('Done writing Tables')
## PLOT SCATTERPLOT COMPARING TRANSFORMED COUNTS
dds <- estimateSizeFactors(dds)
df <- bind_rows(as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>% mutate(transformation = "log2(x + 1)"), as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"), as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))
colnames(df)[1:2] <- c("x", "y")  
pdf(width = 9, height = 3, file = "transformed_counts_comparison.pdf")
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) + coord_fixed() + facet_grid( . ~ transformation)  
dev.off( )

## CALCULATE AND PLOT SAMPLE DISTANCES
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$condition )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf(width = 9, height = 9, file = "sample_distances_vst.pdf", onefile=FALSE)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col = colors)
dev.off( )

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$condition )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf(width = 9, height = 9, file = "sample_distances_rlog.pdf", onefile=FALSE)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col = colors)
dev.off( )

poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$condition )
colnames(samplePoisDistMatrix) <- NULL
pdf(width = 9, height = 9, file = "sample_distances_poisson.pdf", onefile=FALSE )
pheatmap(samplePoisDistMatrix, clustering_distance_rows = poisd$dd, clustering_distance_cols = poisd$dd, col = colors)
dev.off( )
         
          
## CALCULATE SAMPLE PCA PLOTS
pdf(width = 9, height = 9, file = "pca_vst.pdf")
plotPCA(vsd, intgroup = c("condition"))
dev.off( )

pdf(width = 9, height = 9, file = "pca_rlog.pdf")
plotPCA(rld, intgroup = c("condition"))
dev.off( )

