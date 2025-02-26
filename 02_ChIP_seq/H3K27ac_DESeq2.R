#Project: Lymphoma KDM4 Manuscript
#Objective: identify differential H3K27ac peaks in JIB-04 versus DMSO treated OCI-Ly1 cells
#Author: Mohamad Najia

library(data.table)
library(dplyr)
library(DESeq2)
library(BiocParallel)
library(IHW)
library(annotables)
library(stringr)
register(MulticoreParam(2))


#set up environment
project_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
output_dir <- paste0(project_dir, "/DESeq2/")

#import count matrix
fn <- paste0(project_dir,"/H3K27ac.consensus_peaks.featureCounts.txt")
counts.df <- fread(fn, data.table = FALSE)
bed <- counts.df[,c(2,3,4)]
colnames(bed) <- c("chr", "start", "end")
counts.df <- counts.df[,c("DMSO_H3K27ac_REP1.mLb.clN.sorted.bam", "DMSO_H3K27ac_REP2.mLb.clN.sorted.bam", "JIB04_500nM_H3K27ac_REP1.mLb.clN.sorted.bam", "JIB04_500nM_H3K27ac_REP2.mLb.clN.sorted.bam")]

#set up sample metadata
ChIP.condition <- c("DMSO", "DMSO", "JIB04", "JIB04")
colData <- as.data.frame(ChIP.condition)
row.names(colData) <- colnames(counts.df)

#run DESeq2
ChIP.dds <- DESeqDataSetFromMatrix(countData = counts.df, colData = colData, design = ~ ChIP.condition)
ChIP.dds <- DESeq(ChIP.dds, parallel = TRUE)

vsd <- vst(ChIP.dds, blind=FALSE)
plotPCA(vsd, intgroup=c("ChIP.condition"))

rld <- rlog(ChIP.dds, blind = FALSE)
pcaData <- plotPCA(rld, intgroup = c("ChIP.condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plotPCA(rld, intgroup = c("ChIP.condition"))

#perform pairwise comparisons 
res <- results(ChIP.dds, contrast=c("ChIP.condition", "JIB04", "DMSO"), parallel = TRUE, filterFun = ihw, alpha = 0.05)
resSig <- data.frame(res) #pAdjustMethod = "fdr"
resSig$chr <- bed[,1]
resSig$start <- bed[,2]
resSig$end <- bed[,3]

resSig %>% filter(complete.cases(.)) %>% arrange(padj) %>% filter(padj < 0.05) -> out
out$baseMean <- round(out$baseMean, 1)
out$log2FoldChange <- round(out$log2FoldChange, 1)
out$pvalue <- sprintf("%.3e", out$pvalue)
out$padj <- sprintf("%.3e", out$padj)

#output DESeq table of differential peaks
write.table(out[,c("chr", "start", "end", "baseMean", "log2FoldChange", "pvalue", "padj")],
            file = paste0(output_dir, "H3K27ac_DESeq2_JIB04_500nM_v_DMSO_Padj05.tsv"), 
            row.names = FALSE, 
            col.names = TRUE, 
            sep = "\t", 
            quote = FALSE)

#MA plot of differential H3K27ac peaks
resSig$significant <- !is.na(resSig$padj) & resSig$padj < 0.05
resSig$deg_type <- ""
resSig[resSig$significant & resSig$log2FoldChange > 0, "deg_type"] <- "up"
resSig[resSig$significant & resSig$log2FoldChange < 0, "deg_type"] <- "down"


pdf(file = paste0(output_dir, "H3K27ac_MA_plot_JIB04_500nM_v_DMSO.pdf"), width = 3, height = 3, useDingbats = FALSE)  

ggplot(resSig, aes(log2(baseMean), log2FoldChange)) +
  geom_point(aes(color = deg_type), cex = 0.1) + 
  scale_x_continuous(limits = c(-1,15)) + 
  scale_y_continuous(limits = c(-5,5)) + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black") + 
  scale_color_manual(values = c("grey", "red", "blue")) + 
  xlab("log2(Mean Normalized Counts)") +
  ylab("log2(Fold Change) [JIB-04 / DMSO]") + 
  ggtitle("H3K27ac Peaks") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="none",
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))

dev.off()
