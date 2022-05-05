#loading packages
library(DESeq2)
library(tidyverse)
library(pheatmap)

setwd("/Volumes/Turtle 1/Bioinformatics/Analysis.Static.SCFM2/featurecounts_new/0")

#Inputting data
staticdata1 <- read.csv("Combined_no0.csv", sep = ",", quote = "",skip=2)

#setup data in R
names(staticdata1)[1] <- "locus_tag"

sample_info <-  data.frame(condition=c(rep("ATCC-S", 3),
                                       rep("711a", 3),
                                       rep("1302", 3),
                                       rep("ATCC-R", 3),
                                       rep("711b", 2),
                                       rep("634", 3)))

dds <- DESeqDataSetFromMatrix(countData = staticdata1 %>% column_to_rownames("locus_tag"), colData= sample_info, design = ~condition)


# (A) Normalize for differences in sequencing depth (only for normalized reads output)
# This is for getting normalized counts. Sometimes that's nice for plotting or people like to look at them to get a feel for the data.
#dds <- estimateSizeFactors(dds)
#normalized_counts <- counts(dds, normalized=TRUE)
#normalized_counts %>% write.csv("name.csv")

# (B) VST method of normalization
#DESeq2 offers two different methods to perform a more rigorous analysis:
  #•	rlog — a regularised log, and
  #•	vst — a variance stabilising transformation.
#You’d generally use either of these for downstream analysis, not the above (A)
vst <- varianceStabilizingTransformation(dds, blind = F)

# or (if no conditions/replicates are specified)
#vst2 <- varianceStabilizingTransformation(dds)
________________________________________________________________
#looking at the data
#To see the distribution of values, can also exported as hist(assay)
#assay <- assay(dds)
#write.csv(assay, file = paste("test", ".VST_counts.csv", sep=""))

# Plot top 100 most variable genes
#vars <- rowVars(assay(vst))
#vars <- sort(vars, decreasing=TRUE)
#vars <- vars[1:100]
#topVars <- assay(vst)[names(vars),]
_________________________________________________________________
# Perform differential expression
dds <- DESeq(dds)
#Output differential expression data
Differential_expression_output <- results(dds, cooksCutoff = FALSE, contrast= c("condition", "ATCC-R", "ATCC-S")) 
Differential_expression_output %>% write.csv("Combined_no0_new.old_atcc_stats.csv", row.names=TRUE)
Differential_expression_output <- results(dds, cooksCutoff = FALSE, contrast= c("condition", "711b", "711a")) 
Differential_expression_output %>% write.csv("Combined_no0_new.old_711_stats.csv", row.names=TRUE)
Differential_expression_output <- results(dds, cooksCutoff = FALSE, contrast= c("condition", "ATCC-R", "711a")) 
Differential_expression_output %>% write.csv("Combined_no0_new.old_Rvs711a_stats.csv", row.names=TRUE)
Differential_expression_output <- results(dds, cooksCutoff = FALSE, contrast= c("condition", "711b", "ATCC-S")) 
Differential_expression_output %>% write.csv("Combined_no0_new.old_711vsS_stats.csv", row.names=TRUE)
Differential_expression_output <- results(dds, cooksCutoff = FALSE, contrast= c("condition", "634", "ATCC-S")) 
Differential_expression_output %>% write.csv("Combined_no0_new.old_634vsS_stats.csv", row.names=TRUE)
Differential_expression_output <- results(dds, cooksCutoff = FALSE, contrast= c("condition", "634", "711a")) 
Differential_expression_output %>% write.csv("Combined_no0_new.old_634vs711a_stats.csv", row.names=TRUE)
Differential_expression_output <- results(dds, cooksCutoff = FALSE, contrast= c("condition", "634", "1302")) 
Differential_expression_output %>% write.csv("Combined_no0_new.old_634vs1302_stats.csv", row.names=TRUE)
Differential_expression_output <- results(dds, cooksCutoff = FALSE, contrast= c("condition", "711b", "1302")) 
Differential_expression_output %>% write.csv("Combined_no0_new.old_711bvs1302_stats.csv", row.names=TRUE)
Differential_expression_output <- results(dds, cooksCutoff = FALSE, contrast= c("condition", "ATCC-R", "1302")) 
Differential_expression_output %>% write.csv("Combined_no0_new.old_Rvs1302_stats.csv", row.names=TRUE)
#plot the data
res <- results(dds)
tiff("SvsR_Diffplot.tiff")
plotMA(res, ylim=c(-7,7))
abline(h=c(-1,1), col="blue")
dev.off()


# Build volcano plot
ggplot(data=res, aes(x=log2FoldChange, y=pvalue)) + geom_point()
#flip it
res$delabel <- NA
#res$rn <- row.names(res) 
res$delabel <- res$rn 

res <- as.data.frame(res)
res <- res[with(res, order(log2FoldChange)),]
res$threshold <- as.factor(abs(res$log2FoldChange) > 1 & res$padj < 0.01)
res <- res[!is.na(res$padj),]
ggplot(data=res, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) + 
  geom_point(alpha=0.4, size=1.75) + xlim(c(min(res$log2FoldChange), 
                                            c(max(res$log2FoldChange)))) +
  ylim(c(min(-log10(res$padj)), max(-log10(res$padj)))) + 
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(axis.text=element_text(size=12, face="bold")) +
  theme(axis.title=element_text(size=14)) +
  theme(legend.title=element_text(size=14)) +
  theme(legend.text=element_text(size=12))
ggsave("volcano_plot.tiff", dpi=200, device = "tiff") 

# Plot heatmap of differentially expressed genes
# padj < 0.1 and Fold Change > 2
mostSig <- res[res$padj < 0.1 & abs(res$log2FoldChange) >= 1,]
mostSig <- rownames(mostSig)
norm_counts <- counts(dds, normalize=T)
sig_counts <- norm_counts[mostSig,]

tiff("heatmap_diff_expressed.tiff", width=500, height = 1200)
heatmap.2(sig_counts, trace='none', col=greenred(56), 
          scale="row", mar=c(10, 10), key=T, keysize=1, cexRow = 0.3,
          cexCol = 2, main = "padj < 0.05 & Fold Change >= 1")
dev.off()

tiff("heatmap_diff_expressed_different_style.tiff", width=500, height = 1200)
pheatmap(sig_counts, scale = "row", fontsize_row=1, fontsize_col = 10,
         annotation_col=annotation_col, cluster_cols=F)
dev.off()
_________________________________________________
## Plotting PCA & Heatmap
#first make DESeqDataSet
#normalize data using with the VST method
data1 <- plotPCA(vst, returnData=TRUE)
figS1 <- ggplot(data=data1) + geom_point(aes(x=PC1, y=PC2, color=condition), size=2, shape=16) 
figS1

#PCA with variance 
percentVar <- round(100 * attr(data1, "percentVar"))
ggplot(data1, aes(x=PC1, y=PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

figS2 <- ggplot(data1) + geom_point(aes(x=PC1, y=PC2, color=condition), size=2, shape=16) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +
figS2  
  
## Heatmap code of normalized data
hm <- pheatmap(assay(vst), cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, fontsize_col=12)
hm

___________________________________________
#volcano plot
ggplot(vst=de, aes(x=log2FoldChange, y=pvalue)) + geom_point()
  
