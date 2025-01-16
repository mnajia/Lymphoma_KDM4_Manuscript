#Project: Lymphoma KDM4 Manuscript
#Author: Mohamad Najia
#Objective: analyze RNA-seq of JIB-04 treated OCI-Ly1 cells and determine differentially expressed genes

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
library(msigdbr)
library(org.Hs.eg.db)
library(AnnotationDbi)


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

# Declare color maps
pal_rna <- colorRampPalette(c("#352A86","#343DAE","#0262E0","#1389D2",
                              "#2DB7A3","#A5BE6A","#F8BA43","#F6DA23","#F8FA0D"))(100)
pal_atac <- colorRampPalette(c('#3361A5', '#248AF3', '#14B3FF', 
                               '#88CEEF', '#C1D5DC', '#EAD397', 
                               '#FDB31A','#E42A2A', '#A31D1D'))(100)

# Initialize variables
projectdir <- dirname(rstudioapi::getSourceEditorContext()$path)
tsv_dir <- paste0(projectdir, "/kallisto_output_OCI-Ly1/")

# Import samplesheet
fn <- paste0(projectdir, "/OCI-Ly1_samples.txt")
df.samples <- fread(fn, header = FALSE, data.table = FALSE)
colnames(df.samples) <- c("sample_name", "fq1", "fq2")
rownames(df.samples) <- df.samples$sample_name
df.samples$condition <- c( rep("DMSO", 3), rep("JIB04_150nM", 3), rep("JIB04_500nM", 3) )

# Import hg19 USCS genomeStudio gene symbol to transcript id mappings
g2t_map <- read.table(file = paste0(projectdir, "/hg19.annot.cdna.gene.symbol.transcript.map"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(g2t_map) <- c("transcript_id", "gene_symbol")

# Import kallisto abundances
tsv_files <- paste0(tsv_dir, df.samples$sample_name, "/KALLISTO/abundance.tsv")
names(tsv_files) <- df.samples$sample_name
txi <- tximport(tsv_files, type = "kallisto", tx2gene = g2t_map)

# Save kallisto expression matrices 
write.table(txi$abundance,
            file = paste0(projectdir, "/OCI-Ly1_kallisto_rna_seq_abundance_TPM_matrix.tsv"), 
            row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

write.table(txi$counts,
            file = paste0(projectdir, "/OCI-Ly1_kallisto_rna_seq_counts_matrix.tsv"), 
            row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)



####################################################
#DEG Analysis
####################################################

# Set up DESeq 
sampleTable <- data.frame(condition = df.samples$condition)
rownames(sampleTable) <- df.samples$sample_name

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

pdf(file = paste(projectdir, "/OCI-Ly1_pca_rlog_trandformed_RNA_seq_data.pdf", sep = "/"), width = 4.5, height = 3, useDingbats = FALSE)  
ggplot(pcaData, aes(x=PC1, y=PC2, fill=condition, color=condition)) + 
  geom_point() + 
  xlab(paste0("PC1: (", percentVar[1], "% Variance Explained)")) +
  ylab(paste0("PC2: (", percentVar[2], "% Variance Explained)")) +
  ggtitle("PCA on RNA-Seq of\nJIB-04 treated OCI-Ly1 cells") + 
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank())
dev.off()

# DESeq PCA
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("condition"))

# Extract comparisons for differential expression analysis
matched <- unique(df.samples$condition)
comp <- expand.grid(matched, matched, stringsAsFactors = FALSE)
comp <- subset(comp, Var1!=Var2)
comp$cond <- "condition"
comp <- comp %>% 
  mutate(key = paste0(pmin(Var1, Var2), "_" ,pmax(Var1, Var2), sep = "")) %>%
  dplyr::distinct(key, .keep_all=TRUE)

degs <- list()
for (x in 1:dim(comp)[1]) {
  print(x)
  # extract the DESeq contrast
  resSig <- data.frame(results(dds, contrast=c(comp[x,3], comp[x,1], comp[x,2]), parallel = TRUE, pAdjustMethod = "fdr", alpha = 0.05))
  resSig %>% filter(complete.cases(.)) %>% arrange(-log2FoldChange) %>% filter(padj < 0.05) -> out
  out$padj_nlog10 <- -log10(out$padj)
  out$gene <- rownames(out)

  # round
  out$baseMean <- round(out$baseMean, 2)
  out$log2FoldChange <- round(out$log2FoldChange, 2)
  out$pvalue <-sprintf("%.3e", out$pvalue)
  out$padj <-sprintf("%.3e", out$padj)
  
  # output the differential expression table for statistically significant genes
  write.table(out, file = paste0(projectdir, "/DESeq/OCI-Ly1_", comp[x,1], "_v_", comp[x,2], "_Padj05.tsv"), 
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  
  # output the differential expression table for all genes (significant or not)
  resSig %>% filter(complete.cases(.)) %>% arrange(-log2FoldChange) -> out.complete
  out.complete$padj_nlog10 <- -log10(out.complete$padj)
  out.complete$gene <- rownames(out.complete)
  
  write.table(out.complete, file = paste0(projectdir, "/DESeq/OCI-Ly1_", comp[x,1], "_v_", comp[x,2], ".tsv"), 
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  
  # save the set of DEGs
  degs[[comp[x,4]]] <- out
}



####################################################
#DEG Analysis and Interpretation
####################################################

# Set up MSigDB
m_df <- msigdbr(species = "Homo sapiens")
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)

# 500 nM JIB-04 versus DMSO comparison
currset <- degs[[2]]

currset <- currset %>% mutate(ENTREZID = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "SYMBOL") %>% unname())
geneList <- setNames(currset$log2FoldChange,as.numeric(unlist(currset$ENTREZID)) )

em <- GSEA(geneList, TERM2GENE = m_t2g)
df.gsea <- em@result

core_enrichment <- df.gsea["HALLMARK_P53_PATHWAY", "core_enrichment"] %>% str_split(pattern = "/")
leading_edge_genes <- mapIds(org.Hs.eg.db, keys = core_enrichment[[1]], column = "SYMBOL", keytype = "ENTREZID") %>% unname()

core_enrichment <- df.gsea["HALLMARK_INFLAMMATORY_RESPONSE", "core_enrichment"] %>% str_split(pattern = "/")
leading_edge_genes <- mapIds(org.Hs.eg.db, keys = core_enrichment[[1]], column = "SYMBOL", keytype = "ENTREZID") %>% unname()

core_enrichment <- df.gsea["HALLMARK_TNFA_SIGNALING_VIA_NFKB", "core_enrichment"] %>% str_split(pattern = "/")
leading_edge_genes <- mapIds(org.Hs.eg.db, keys = core_enrichment[[1]], column = "SYMBOL", keytype = "ENTREZID") %>% unname()

core_enrichment <- df.gsea["HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "core_enrichment"] %>% str_split(pattern = "/")
leading_edge_genes <- mapIds(org.Hs.eg.db, keys = core_enrichment[[1]], column = "SYMBOL", keytype = "ENTREZID") %>% unname()

core_enrichment <- df.gsea["HALLMARK_UV_RESPONSE_UP", "core_enrichment"] %>% str_split(pattern = "/")
leading_edge_genes <- mapIds(org.Hs.eg.db, keys = core_enrichment[[1]], column = "SYMBOL", keytype = "ENTREZID") %>% unname()

core_enrichment <- df.gsea["HALLMARK_INTERFERON_GAMMA_RESPONSE", "core_enrichment"] %>% str_split(pattern = "/")
leading_edge_genes <- mapIds(org.Hs.eg.db, keys = core_enrichment[[1]], column = "SYMBOL", keytype = "ENTREZID") %>% unname()



terms <- c("HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_P53_PATHWAY", 
  "HALLMARK_E2F_TARGETS", 
  "HALLMARK_KRAS_SIGNALING_UP", 
  "HALLMARK_HYPOXIA", 
  "HALLMARK_MYC_TARGETS_V2", 
  "HALLMARK_G2M_CHECKPOINT", 
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING", 
  "HALLMARK_UV_RESPONSE_UP", 
  "HALLMARK_IL2_STAT5_SIGNALING", 
  "HALLMARK_IL6_JAK_STAT3_SIGNALING", 
  "HALLMARK_APOPTOSIS", 
  "HALLMARK_UV_RESPONSE_DN", 
  "HALLMARK_INTERFERON_GAMMA_RESPONSE", 
  "HALLMARK_INFLAMMATORY_RESPONSE", 
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE")


df.plot <- df.gsea[terms, c("ID", "NES")]
df.plot$direction <- "up"
df.plot[df.plot$NES < 0, "direction"] <- "down"
df.plot$name <- c("MYC targets v1",
                  "TNFa signaling via NFkB",
                  "Oxidative phosphorylation",
                  "P53 pathway",
                  "E2F targets", 
                  "KRAS signaling up",
                  "Hypoxia", 
                  "MYC targets v2", 
                  "G2M checkpoint",
                  "Wnt beta-catenin signaling", 
                  "UV response up",
                  "IL2-STAT5 signaling", 
                  "IL6-JAK-STAT3 signaling",
                  "Apoptosis",
                  "UV response down",
                  "IFN-g response",
                  "Inflammatory response",
                  "Unfolded protein response")
df.plot <- df.plot[order(df.plot$NES, decreasing = FALSE),]
df.plot$name <- factor(df.plot$name, levels = df.plot$name)

pdf(file = paste0(projectdir, "/OCI-Ly1_GSEA_significant_Hallmark_gene_sets_JIB04_500nM_v_DMSO.pdf"), width = 4, height = 4)  

ggplot(df.plot, aes(x = NES, y = name, fill = direction)) + 
  geom_bar(stat = "identity") +
  xlab("NES") + 
  ylab("Hallmark Gene Set (FDR < 0.05)") +
  ggtitle("GSEA on JIB-04 treated OCI-Ly1 cells") + 
  theme_bw() +
  theme(legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour="black"), 
        axis.text.y = element_text(colour="black"))

dev.off()







