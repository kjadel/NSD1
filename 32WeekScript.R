#install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

#install DESeq2
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")


#install gplots and ggplot2
install.packages(c("ggplot2", "scales", "viridis"))
install.packages("RColorBrewer")
install.packages("gplots")

# load libraries
library(DESeq2)
library(airway)
library(ggplot2)
library(scales)
library(viridis)
library(RColorBrewer)
library(gplots)

sample_info <- as.data.frame(colData(zero_weeks))
sample_info <- sample_info[,c(2,3)]
sample_info$dex <- gsub('trt', 'treated', sample_info$dex)
sample_info$dex <- gsub('untrt', 'untreated', sample_info$dex)
names(sample_info) <- c('cellLine', 'dexamethasone')
write.table(sample_info, file = "sample_info.csv", sep = ',', col.names = T, row.names = T, quote = F)

# script to perform differential gene expression analysis using DESeq2 package
# setwd("~/Desktop/demo/DESeq2_tutorial/data")



zero_weeks <- read.csv("est_counts_genes_kallisto_LizzySet5.csv")
head(zero_weeks)
metaData<- read.csv('32WeeksSampleInfo.csv', header = TRUE, sep = ",")
write.table(metaData, file = "32WeeksSample_Info.csv", sep = ',', col.names = T, row.names = FALSE, quote = F)

countsData <- read.csv('est_counts_genes_kallisto_LizzySet5.csv', header = TRUE, sep = ",")
head(countsData)
write.table(countsData, file = "32WeeksCounts_data.csv", sep = ',', col.names = T, row.names = T, quote = F)


# Step 1: preparing count data ----------------

# read in counts data
counts_data <- read.csv('est_counts_genes_kallisto_LizzySet5.csv')
head(counts_data)


# read in sample info
colData <- read.csv('32WeeksSample_Info.csv')
head(colData)
colnames(counts_data)
# making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data) %in% rownames(colData))


# are they in the same order?
all(colnames(counts_data) == rownames(colData))


# Step 2: construct a DESeqDataSet object ----------

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                       colData = colData,
                       design = ~ Condition)

dds


# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
keep
dds <- dds[keep,]

dds

# set the factor level
dds$Condition <- relevel(dds$Condition, ref = "WT")


# NOTE: collapse technical replicates


# Step 3: Run DESeq ----------------------
dds <- DESeq(dds)
res <- results(dds)

res

# Step 4: Explore Results ----------------

summary(res)

res[order(res$padj),]

res1 <- results(dds, alpha=0.01)
summary(res1)

res2 <- results(dds, alpha=0.01, lfcThreshold=1)
summary(res2)

# contrasts
resultsNames(dds)

# MA plot
plotMA(res2)

#Exporting results to CSV files
write.csv(as.data.frame(res), 
          file="32Weekcondition_treated_results.csv")

resFilt <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
write.csv(resFilt, file="32WeekDE_results_filtered.csv")


# Step 5: Additional ways to visualize the results ----------------------


# 5.1. To Plot the dispersions:
plotDispEsts(dds)


# 5.2. creating a PCA (Principal Components Analyis) plot:
rld <- rlog(dds)
plotPCA(rld,intgroup = "Condition")

# 5.3. Make a better MA-plot

# Coerce to a data frame
deseq2ResDF <- as.data.frame(res)

# Examine this data frame
head(deseq2ResDF)

# Set a boolean column for significance
deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < .1, "Significant", NA)

# Plot the results similar to DEseq2
ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="tomato1", size=2) + labs(x="mean of normalized counts", y="log fold change") + scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") + theme_bw()

# Better plot
ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=padj)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") + labs(x="mean of normalized counts", y="log fold change") + scale_colour_viridis(direction=-1, trans='sqrt') + theme_bw() + geom_density_2d(colour="black", size=2)

#5.4. Viewing normalized counts for a single geneID
plotCounts(dds,gene = "ENSG00000000003",intgroup = "dexamethasone")

#5.5 Make a heatmap of our results

# Transform count data using the variance stablilizing transform
deseq2VST <- vst(dds)

# Convert the DESeq transformed object to a data frame
deseq2VST <- assay(deseq2VST)
deseq2VST <- as.data.frame(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)
head(deseq2VST)

# Keep only the significantly differentiated genes where the fold-change was at least 1
sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= .1 & abs(deseq2ResDF$log2FoldChange) > 2,])
deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]

# Convert the VST counts to long format for ggplot2
install.packages("reshape2")
library(reshape2)

# First compare wide vs long version
deseq2VST_wide <- deseq2VST
deseq2VST_long <- melt(deseq2VST, id.vars=c("Gene"))

head(deseq2VST_wide)
head(deseq2VST_long)

# Now overwrite our original data frame with the long format
deseq2VST <- melt(deseq2VST, id.vars=c("Gene"))

# Make a heatmap
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap

#5.6 Clustering
# Convert the significant genes back to a matrix for clustering
deseq2VSTMatrix <- dcast(deseq2VST, Gene ~ variable)
rownames(deseq2VSTMatrix) <- deseq2VSTMatrix$Gene
deseq2VSTMatrix$Gene <- NULL

# Compute a distance calculation on both dimensions of the matrix
distanceGene <- dist(deseq2VSTMatrix)
distanceSample <- dist(t(deseq2VSTMatrix))

# Cluster based on the distance calculations
clusterGene <- hclust(distanceGene, method="average")
clusterSample <- hclust(distanceSample, method="average")

# Construct a dendogram for samples
install.packages("ggdendro")
library(ggdendro)
sampleModel <- as.dendrogram(clusterSample)
sampleDendrogramData <- segment(dendro_data(sampleModel, type = "rectangle"))
sampleDendrogram <- ggplot(sampleDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + theme_dendro()

# Re-factor samples for ggplot2
deseq2VST$variable <- factor(deseq2VST$variable, levels=clusterSample$labels[clusterSample$order])

# Construct the heatmap. note that at this point we have only clustered the samples NOT the genes
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap

# Combine the dendrogram and the heatmap
install.packages("gridExtra")
library(gridExtra)
grid.arrange(sampleDendrogram, heatmap, ncol=1, heights=c(1,5))

#fixing the alignment
install.packages("gtable")
library(gtable)
library(grid)

# Modify the ggplot objects
sampleDendrogram_1 <- sampleDendrogram + scale_x_continuous(expand=c(.0085, .0085)) + scale_y_continuous(expand=c(0, 0))
heatmap_1 <- heatmap + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0))

# Convert both grid based objects to grobs
sampleDendrogramGrob <- ggplotGrob(sampleDendrogram_1)
heatmapGrob <- ggplotGrob(heatmap_1)

# Check the widths of each grob
sampleDendrogramGrob$widths
heatmapGrob$widths

# Add in the missing columns
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[7], 6)
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[8], 7)

# Make sure every width between the two grobs is the same
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)

# Arrange the grobs into a plot
finalGrob <- arrangeGrob(sampleDendrogramGrob, heatmapGrob, ncol=1, heights=c(2,5))

# Draw the plot
grid.draw(finalGrob)

#add a plot between the dendrogram and the heatmap showing the condition

# Re-order the sample data to match the clustering we did
colData_v2 <- read.csv('sample_info_v2.csv')
head(colData_v2)
colData_v2$Run <- factor(colData_v2$Run, levels=clusterSample$labels[clusterSample$order])

# Construct a plot to show the clinical data
colours <- c("#743B8B", "#8B743B", "#8B3B52")
test_condition <- ggplot(colData_v2, aes(x=Run, y=1, fill=dexamethasone)) + geom_tile() + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0)) + scale_fill_manual(name="Dexa", values=colours) + theme_void()

# Convert the clinical plot to a grob
test_condition_Grob <- ggplotGrob(test_condition)

# Make sure every width between all grobs is the same
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths, test_condition_Grob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)
test_condition_Grob$widths <- as.list(maxWidth)

# Arrange and output the final plot
finalGrob <- arrangeGrob(sampleDendrogramGrob, test_condition_Grob, heatmapGrob, ncol=1, heights=c(2,1,5))
grid.draw(finalGrob)

# Step 6: Pathway Analysis ----------------------
# install gage
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("gage","GO.db","AnnotationDbi","org.Hs.eg.db"), version = "3.16")
library(gage)

# extract the results from the deseq2 data
treated_vs_untreated_DE <- results(dds)

# set up kegg database
kg.hsa <- kegg.gsets(species="hsa")
kegg.sigmet.gs <- kg.hsa$kg.sets[kg.hsa$sigmet.idx]
kegg.dise.gs <- kg.hsa$kg.sets[kg.hsa$dise.idx]

# set up go database
go.hs <- go.gsets(species="human")
go.bp.gs <- go.hs$go.sets[go.hs$go.subs$BP]
go.mf.gs <- go.hs$go.sets[go.hs$go.subs$MF]
go.cc.gs <- go.hs$go.sets[go.hs$go.subs$CC]

# load in libraries to annotate data
library(AnnotationDbi)
library(org.Hs.eg.db)

# annotate the deseq2 results with additional gene identifiers
treated_vs_untreated_DE$symbol <- mapIds(org.Hs.eg.db, keys=row.names(treated_vs_untreated_DE), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
treated_vs_untreated_DE$entrez <- mapIds(org.Hs.eg.db, keys=row.names(treated_vs_untreated_DE), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
treated_vs_untreated_DE$name <- mapIds(org.Hs.eg.db, keys=row.names(treated_vs_untreated_DE), column="GENENAME", keytype="ENSEMBL", multiVals="first")

write.csv(as.data.frame(treated_vs_untreated_DE), 
          file="Pathway_analysis_genes_annotated.csv")

# grab the log fold changes for everything
treated_vs_untreated_DE.fc <- treated_vs_untreated_DE$log2FoldChange
names(treated_vs_untreated_DE.fc) <- treated_vs_untreated_DE$entrez

# Run enrichment analysis on all log fc
fc.kegg.sigmet.p <- gage(treated_vs_untreated_DE.fc, gsets = kegg.sigmet.gs)
fc.kegg.dise.p <- gage(treated_vs_untreated_DE.fc, gsets = kegg.dise.gs)
fc.go.bp.p <- gage(treated_vs_untreated_DE.fc, gsets = go.bp.gs)
fc.go.mf.p <- gage(treated_vs_untreated_DE.fc, gsets = go.mf.gs)
fc.go.cc.p <- gage(treated_vs_untreated_DE.fc, gsets = go.cc.gs)

# covert the kegg results to data frames
fc.kegg.sigmet.p.up <- as.data.frame(fc.kegg.sigmet.p$greater)
fc.kegg.dise.p.up <- as.data.frame(fc.kegg.dise.p$greater)

fc.kegg.sigmet.p.down <- as.data.frame(fc.kegg.sigmet.p$less)
fc.kegg.dise.p.down <- as.data.frame(fc.kegg.dise.p$less)

write.csv(as.data.frame(fc.kegg.sigmet.p.up), 
          file="KEGG_up.csv")

# convert the go results to data frames
fc.go.bp.p.up <- as.data.frame(fc.go.bp.p$greater)
fc.go.mf.p.up <- as.data.frame(fc.go.mf.p$greater)
fc.go.cc.p.up <- as.data.frame(fc.go.cc.p$greater)

fc.go.bp.p.down <- as.data.frame(fc.go.bp.p$less)
fc.go.mf.p.down <- as.data.frame(fc.go.mf.p$less)
fc.go.cc.p.down <- as.data.frame(fc.go.cc.p$less)

# Pathway visualization using Pathview
# Install pathview from bioconductor
BiocManager::install("pathview")
library(pathview)

# View the hsa04140 pathway from the pathway analysis
fc.kegg.sigmet.p.up[grepl("hsa04140", rownames(fc.kegg.sigmet.p.up), fixed=TRUE),]

# Overlay the expression data onto this pathway
pathview(gene.data=treated_vs_untreated_DE.fc, species="hsa", pathway.id="hsa04140")

# Generate Gene symbols for DE genes for webpage analysis input

resFilt <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]

resFilt$symbol <- mapIds(org.Hs.eg.db, keys=row.names(resFilt), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
resFilt$entrez <- mapIds(org.Hs.eg.db, keys=row.names(resFilt), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
resFilt$name <- mapIds(org.Hs.eg.db, keys=row.names(resFilt), column="GENENAME", keytype="ENSEMBL", multiVals="first")
write.csv(resFilt, file="DE_results_filtered_gene_annotated.csv")




