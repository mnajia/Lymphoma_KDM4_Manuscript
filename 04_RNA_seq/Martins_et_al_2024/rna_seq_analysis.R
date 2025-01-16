#Project: Lymphoma KDM4 Manuscript
#Author: Mohamad Najia
#Objective: analyze RNA-Seq from Martins et al 2024 to find core genes altered by knockdown of ZNF587/417

library(DESeq2)
library(IHW)
library(jsonlite)
library(dplyr)
library(data.table)
library(tximport)
library(rhdf5)
library(ggplot2)
library(scales)
library(EnhancedVolcano)
library(rjson)
library(ComplexHeatmap)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(circlize)
library(irlba)
library(matrixStats)
library(M3C)
library(clusterProfiler)
library(enrichplot)


####################################################
#Function Declarations
####################################################

scaleRows <- function(x) {
  rm <- rowMeans(x)
  x <- sweep(x, 1, rm)
  sx <- apply(x, 1, sd)
  x <- sweep(x, 1, sx, "/")
}


####################################################
#Set-up Environment
####################################################

# Initialize variables
projectdir <- dirname(rstudioapi::getSourceEditorContext()$path)
tsv_dir <- paste0(projectdir, "/kallisto_output/")

# Import samplesheet
fn <- paste0(projectdir, "/samples_metadata.txt")
df.samples <- fread(fn, header = TRUE, data.table = FALSE)
rownames(df.samples) <- df.samples$srr

#df.samples <- df.samples[c(1:8),]
df.samples <- df.samples[c(9:16),]

# import hg19 USCS genomeStudio gene symbol to transcript id mappings
g2t_map <- read.table(file = "/Volumes/blainey_lab/Mo/genomes/hg19.annot.cdna.gene.symbol.transcript.map", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(g2t_map) <- c("transcript_id", "gene_symbol")

# Import kallisto abundances
tsv_files <- paste0(tsv_dir, df.samples$srr, "/KALLISTO/abundance.tsv")
names(tsv_files) <- df.samples$srr
txi <- tximport(tsv_files, type = "kallisto", tx2gene = g2t_map)



####################################################
#DEG Analysis
####################################################

# Set up DESeq 
sampleTable <- data.frame(condition = df.samples$shRNA)
rownames(sampleTable) <- df.samples$srr

dds <- DESeqDataSetFromTximport(txi, sampleTable, design = ~ condition)
dds <- DESeq(dds)

# Perform PCA
nrow(dds)
keep <- rowSums(counts(dds)) > 1
dds.keep <- dds[keep,]
nrow(dds.keep)
rld <- rlog(dds.keep, blind = FALSE)
pcaData <- plotPCA(rld, intgroup = c("condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x=PC1, y=PC2, fill=condition, color=condition)) + 
  geom_point() + 
  xlab(paste0("PC1: (", percentVar[1], "% Variance Explained)")) +
  ylab(paste0("PC2: (", percentVar[2], "% Variance Explained)")) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank())

# DESeq PCA
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("condition"))

# Extract the DESeq contrast
resSig <- data.frame(results(dds, contrast=c("condition", "ZNF587-ZNF417", "Control"), parallel = TRUE, pAdjustMethod = "fdr", alpha = 0.05))
resSig %>% filter(complete.cases(.)) %>% arrange(-log2FoldChange) %>% filter(padj < 0.05) -> out
out$padj_nlog10 <- -log10(out$padj)
out$gene <- rownames(out)
out$baseMean <- round(out$baseMean, 2)
out$log2FoldChange <- round(out$log2FoldChange, 2)
out$pvalue <-sprintf("%.3e", out$pvalue)
out$padj <-sprintf("%.3e", out$padj)

# Output the differential expression table for statistically significant genes
#write.table(out, file = paste0(projectdir, "/rna_seq_analysis/DESeq/U2932_D6_shZNF587_v_shControl_Padj05.tsv"), 
#            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

write.table(out, file = paste0(projectdir, "/rna_seq_analysis/DESeq/OCI-Ly7_D6_shZNF587_v_shControl_Padj05.tsv"), 
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)



# Output the differential expression table for all genes (significant or not)
resSig %>% filter(complete.cases(.)) %>% arrange(-log2FoldChange) -> out.complete
out.complete$padj_nlog10 <- -log10(out.complete$padj)
out.complete$gene <- rownames(out.complete)

#write.table(out.complete, file = paste0(projectdir, "/rna_seq_analysis/DESeq/U2932_D6_shZNF587_v_shControl.tsv"), 
#            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

write.table(out.complete, file = paste0(projectdir, "/rna_seq_analysis/DESeq/OCI-Ly7_D6_shZNF587_v_shControl.tsv"), 
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)



####################################################
#Determine conserved DEGs across time/cell lines
####################################################

# Import OCI-Ly7 D6 DEGs
fn <- paste0(projectdir, "/rna_seq_analysis/DESeq/OCI-Ly7_D6_shZNF587_v_shControl_Padj05.tsv")
df.ly7 <- fread(fn, data.table = FALSE)

# Import U2932 D6 DEGs
fn <- paste0(projectdir, "/rna_seq_analysis/DESeq/U2932_D6_shZNF587_v_shControl_Padj05.tsv")
df.u2932 <- fread(fn, data.table = FALSE)

# Find conserved up-regulated genes in both cell lines
upgenes <- intersect(df.ly7[(df.ly7$log2FoldChange > 0.75) & (df.ly7$padj < 0.05), "gene"], df.u2932[(df.u2932$log2FoldChange > 0.75) & (df.u2932$padj < 0.05), "gene"])
write.table(paste0(c("shZNF587_UP", "Martins_et_al_2024", upgenes), collapse = "\t"), 
            file = paste0(projectdir, "/rna_seq_analysis/shZNF587_D6_conserved_upregulated_genes.gmt"), 
            quote = FALSE,
            row.names = FALSE, 
            col.names = FALSE,
            sep = "\t")

downgenes <- intersect(df.ly7[(df.ly7$log2FoldChange < 0) & (df.ly7$padj < 0.01), "gene"], df.u2932[(df.u2932$log2FoldChange < 0) & (df.u2932$padj < 0.01), "gene"])
write.table(paste0(c("shZNF587_DOWN", "Martins_et_al_2024", downgenes), collapse = "\t"), 
            file = paste0(projectdir, "/rna_seq_analysis/shZNF587_D6_conserved_downregulated_genes.gmt"), 
            quote = FALSE,
            row.names = FALSE, 
            col.names = FALSE,
            sep = "\t")




