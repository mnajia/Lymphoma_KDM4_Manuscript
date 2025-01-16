#Project: Lymphoma KDM4 Manuscript
#Author: Mohamad Najia
#Objective: plot differentially expressed genes conserved between JIB-04 treated OCI-Ly1 and TMD8 cells

library(dplyr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(enrichplot)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(stringr)


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
#Set up environment
####################################################

# Declare colormaps
pal_rna <- colorRampPalette(c("#352A86","#343DAE","#0262E0","#1389D2",
                              "#2DB7A3","#A5BE6A","#F8BA43","#F6DA23","#F8FA0D"))(100)
pal_atac <- colorRampPalette(c('#3361A5', '#248AF3', '#14B3FF', 
                               '#88CEEF', '#C1D5DC', '#EAD397', 
                               '#FDB31A','#E42A2A', '#A31D1D'))(100)

# Initialize variables
projectdir <- dirname(rstudioapi::getSourceEditorContext()$path)



####################################################
#Plot concordance of differentially expressed genes between TMB8 and OCI-Ly1 lines
####################################################

# Import DESeq results
fn <- paste0(projectdir, "/DESeq/OCI-Ly1_JIB04_500nM_v_DMSO.tsv")
df.ly1 <- fread(fn, data.table = FALSE)
rownames(df.ly1) <- df.ly1$gene

fn <- paste0(projectdir, "/DESeq/TMD8_JIB04_500nM_v_DMSO.tsv")
df.tmd8 <- fread(fn, data.table = FALSE)
rownames(df.tmd8) <- df.tmd8$gene


core_genes <- intersect(df.ly1$gene, df.tmd8$gene)
df.ly1 <- df.ly1[core_genes,]
df.tmd8 <- df.tmd8[core_genes,]

df.diff <- data.frame(gene = df.ly1$gene,
                      log2FC_ly1 = df.ly1$log2FoldChange,
                      padj_ly1 = df.ly1$padj,
                      log2FC_tmd8 = df.tmd8$log2FoldChange,
                      padj_tmd8 = df.tmd8$padj)
df.diff$significant <- "no"
df.diff[ (df.diff$padj_ly1 < 0.1) & (df.diff$padj_tmd8 < 0.1), "significant"] <- "yes"
df.diff$highlight <- ""

highlight_genes <- c("CXCL10", "ZNF587", "ZNF587B", "ZNF417", "HLA-A", "HLA-B", "HLA-C", "IRF7", "CDKN1A", "NFKB1", "IL10RA", "IRF9", "ULBP1")
df.diff[df.diff$gene %in% highlight_genes, "highlight"] <- df.diff[df.diff$gene %in% highlight_genes, "gene"]

ggplot(df.diff, aes(x=log2FC_ly1, y=log2FC_tmd8, color=significant, label=highlight)) + 
  geom_point(size=1) + 
  xlab("OCI-Ly1 - log2(fold change)") + 
  ylab("TMD8 - log2(fold change)") + 
  ggtitle("DEG concordance across DLBCL cells") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size=12),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")) +
  geom_text_repel(
    aes(label = highlight),
    size = 3,
    min.segment.length = 0, 
    seed = 42, 
    box.padding = 0.5,
    max.overlaps = Inf,
    arrow = arrow(length = unit(0.010, "npc")),
    nudge_x = .15,
    nudge_y = .5,
    color = "grey50"
  )




####################################################
#Plot expression heatmap of leading edge genes
####################################################

fn <- paste0(projectdir, "/OCI-Ly1_kallisto_rna_seq_abundance_TPM_matrix.tsv")
df.ly1 <- fread(fn, data.table = FALSE)
rownames(df.ly1) <- df.ly1$V1
df.ly1 <- df.ly1[,grepl("DMSO|500", colnames(df.ly1))]


fn <- paste0(projectdir, "/TMD8_kallisto_rna_seq_abundance_TPM_matrix.tsv")
df.tmd8 <- fread(fn, data.table = FALSE)
rownames(df.tmd8) <- df.tmd8$V1
df.tmd8 <- df.tmd8[,grepl("DMSO|500", colnames(df.tmd8))]


highlight_genes <- c("ZNF587", "ZNF587B", 
                     "IKZF1", "IKZF3", "IKZF4", "SYK", 
                     "E2F5", "FYN", "TP53", "TP53I13", "ATF3", "TP63", "CDKN1A", "CEBPA", "CEBPB", "TXNIP", 
                     "CXCL11", "IRF1","IRF7", "IRF9", "IL6ST", "CALCRL", "JUNB", "FOS", "STAT5A", "ISG20", "ICOSLG",  "CCR7", "IFIT2",
                     "CXCL10", "ISG15", "NFKBIA","CD69", "CD274",  "ICAM1", "IFITM1", "IFI35", "IFI27", "TNFAIP3",  "STAT1"
                      )

df.plot.ly1 <- df.ly1[highlight_genes, ]
df.plot.tmd8 <- df.tmd8[highlight_genes, ]

color_maps <- c("DMSO" = "#4777b2", "JIB04 500 nM" = "#a3062a")
column_ha <- HeatmapAnnotation(condition = c(rep("DMSO",3), rep("JIB04 500 nM", 3)), 
                               col = list(condition = color_maps), 
                               annotation_name_side = "left",
                               border = TRUE)

hm1 <- Heatmap(scaleRows( log2(df.plot.ly1+1) ), 
              column_title = "OCI-Ly1", 
              show_row_names = TRUE, 
              cluster_columns = FALSE, 
              cluster_rows = FALSE,
              show_column_names = FALSE, 
              show_column_dend = FALSE,
              top_annotation = column_ha, 
              col=pal_atac,
              heatmap_legend_param = list(direction = "vertical", border=TRUE),
              border = TRUE,
              name = "Expression Z-score")

hm2 <- Heatmap(scaleRows( log2(df.plot.tmd8+1) ), 
               column_title = "TMD8", 
               show_row_names = TRUE, 
               cluster_columns = FALSE, 
               cluster_rows = FALSE,
               show_column_names = FALSE, 
               show_column_dend = FALSE,
               top_annotation = column_ha, 
               col=pal_atac,
               heatmap_legend_param = list(direction = "vertical", border=TRUE),
               border = TRUE,
               name = "Expression Z-score")

pdf(file = paste0(projectdir, "/OCI-Ly1_TMD8_DEGs_expression_heatmap.pdf"), width = 5, height = 6, useDingbats = FALSE)
draw(hm1+hm2, merge_legend = TRUE)
dev.off()



####################################################
#Plot expression of cell cycle genes
####################################################

# Load cell cycle genes from Seurat
library(Seurat)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cell_cycle_genes <- c(g2m.genes, s.genes)
cell_cycle_phase <- c( rep("G2M", length(g2m.genes)), rep("S", length(s.genes)) )


df.plot.ly1 <- df.ly1[cell_cycle_genes, ] 
df.plot.tmd8 <- df.tmd8[cell_cycle_genes, ] 

color_maps <- c("DMSO" = "#4777b2", "JIB04 500 nM" = "#a3062a")
column_ha <- HeatmapAnnotation(condition = c(rep("DMSO",3),  rep("JIB04 500 nM", 3)), 
                               col = list(condition = color_maps), 
                               annotation_name_side = "left",
                               border = TRUE)

hm1 <- Heatmap(scaleRows( log2(df.plot.ly1+1) ), 
               column_title = "OCI-Ly1", 
               show_row_names = TRUE, 
               cluster_columns = FALSE, 
               cluster_rows = TRUE,
               row_split = cell_cycle_phase,
               row_names_gp = gpar(fontsize = 8),
               show_column_names = FALSE, 
               show_column_dend = FALSE,
               top_annotation = column_ha, 
               col=pal_atac,
               heatmap_legend_param = list(direction = "vertical", border=TRUE),
               border = TRUE,
               name = "Expression Z-score")

hm2 <- Heatmap(scaleRows( log2(df.plot.tmd8+1) ), 
               column_title = "TMD8", 
               show_row_names = TRUE, 
               cluster_columns = FALSE, 
               cluster_rows = TRUE,
               row_split = cell_cycle_phase,
               row_names_gp = gpar(fontsize = 8),
               show_column_names = FALSE, 
               show_column_dend = FALSE,
               top_annotation = column_ha, 
               col=pal_atac,
               heatmap_legend_param = list(direction = "vertical", border=TRUE),
               border = TRUE,
               name = "Expression Z-score")

pdf(file = paste0(projectdir, "/cell_cycle_genes_heatmap.pdf"), width = 5, height = 10, useDingbats = FALSE)
draw(hm1+hm2, merge_legend = TRUE)
dev.off()



####################################################
#Plot expression of zinc finger genes
####################################################

fn <- paste0(projectdir, "/DESeq/OCI-Ly1_JIB04_500nM_v_DMSO_Padj05.tsv")
df.ly1 <- fread(fn, data.table = FALSE)
df.ly1 <- df.ly1[!is.na(df.ly1$log2FoldChange) & !is.na(df.ly1$padj),]
df.ly1$rank <- 1:dim(df.ly1)[1]

fn <- paste0(projectdir, "/DESeq/TMD8_JIB04_500nM_v_DMSO_Padj05.tsv")
df.tmd8 <- fread(fn, data.table = FALSE)
df.tmd8 <- df.tmd8[!is.na(df.tmd8$log2FoldChange) & !is.na(df.tmd8$padj),]
df.tmd8$rank <- 1:dim(df.tmd8)[1]

znfgenes <- intersect(df.ly1[grepl("ZNF", df.ly1$gene), "gene"], df.tmd8[grepl("ZNF", df.tmd8$gene), "gene"])

df.ly1$genelabels <- df.ly1$gene %in% znfgenes
df.ly1$sig <- ""
df.ly1[(df.ly1$log2FoldChange > 0), "sig"] <- "up"
df.ly1[(df.ly1$log2FoldChange < 0), "sig"] <- "down"

gg <- ggplot(df.ly1) + 
  geom_point(aes(rank, log2FoldChange, col=sig)) + 
  xlab("Rank") + 
  ylab("log2(Fold Change)") + 
  ggtitle("JIB-04 versus DMSO treated OCI-Ly1 cells\nDifferentially Expressed ZNF Genes (FDR < 0.05)") + 
  theme_bw() + 
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 10),
    legend.position = "none",
    strip.text = element_text(colour = "black")) + 
  geom_text_repel(aes(rank, log2FoldChange), 
                  label = ifelse(df.ly1$genelabels == TRUE, as.character(df.ly1$gene),""), 
                  box.padding = unit(0.3, "lines"),
                  max.overlaps = Inf) 

pdf(file = paste0(projectdir, "/OCI-Ly1_JIB04_500nM_v_DMSO_volcano_conserved_znf_genes.pdf"), width = 5, height = 5, useDingbats = FALSE)  
gg
dev.off()



df.tmd8$genelabels <- df.tmd8$gene %in% znfgenes
df.tmd8$sig <- ""
df.tmd8[(df.tmd8$log2FoldChange > 0), "sig"] <- "up"
df.tmd8[(df.tmd8$log2FoldChange < 0), "sig"] <- "down"

gg <- ggplot(df.tmd8) + 
  geom_point(aes(rank, log2FoldChange, col=sig)) + 
  xlab("Rank") + 
  ylab("log2(Fold Change)") + 
  ggtitle("JIB-04 versus DMSO treated TMD8 cells\nDifferentially Expressed ZNF Genes (FDR < 0.05)") + 
  theme_bw() + 
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 10),
    legend.position = "none",
    strip.text = element_text(colour = "black")) + 
  geom_text_repel(aes(rank, log2FoldChange), 
                  label = ifelse(df.tmd8$genelabels == TRUE, as.character(df.tmd8$gene),""), 
                  box.padding = unit(0.3, "lines"),
                  max.overlaps = Inf) 

pdf(file = paste0(projectdir, "/TMD8_JIB04_500nM_v_DMSO_volcano_conserved_znf_genes.pdf"), width = 5, height = 5, useDingbats = FALSE)  
gg
dev.off()



####################################################
#Plot log2FC of differentially expressed znf genes between TMD8 and OCI-Ly1 cells
####################################################

fn <- paste0(projectdir, "/DESeq/OCI-Ly1_JIB04_500nM_v_DMSO_Padj05.tsv")
df.ly1 <- fread(fn, data.table = FALSE)
rownames(df.ly1) <- df.ly1$gene

fn <- paste0(projectdir, "/DESeq/TMD8_JIB04_500nM_v_DMSO_Padj05.tsv")
df.tmd8 <- fread(fn, data.table = FALSE)
rownames(df.tmd8) <- df.tmd8$gene

df <- data.frame(znf_genes = znfgenes,
                 ly1_log2FC = df.ly1[znfgenes, "log2FoldChange"],
                 tmd8_log2FC = df.tmd8[znfgenes, "log2FoldChange"])

gg <- ggplot(df) + 
  geom_point(aes(ly1_log2FC, tmd8_log2FC)) + 
  geom_vline(xintercept = 0,  linetype = "dashed", color = "red") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
  scale_x_continuous(limits = c(-4,4)) + 
  scale_y_continuous(limits = c(-2,2)) + 
  xlab("OCI-Ly1 log2(Fold Change)") + 
  ylab("TMD8 log2(Fold Change)") + 
  ggtitle("Differentially Expressed ZNF Genes (FDR < 0.05)\nJIB-04 500 nM versus DMSO") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 10),
    legend.position = "none",
    strip.text = element_text(colour = "black")) + 
  geom_text_repel(aes(ly1_log2FC, tmd8_log2FC), 
                  label = ifelse(df$znf_genes %in% c("ZNF134","ZNF160", "ZNF17", "ZNF543", "ZNF544", "ZNF547", "ZNF586", "ZNF587", "ZNF587B",
                                                     "ZNF671", "ZNF749", "ZNF761", "ZNF765", "ZNF773", "ZNF776", "ZNF805", "ZNF814"), as.character(df$znf_genes),""), 
                  box.padding = unit(0.3, "lines"),
                  max.overlaps = Inf) 

pdf(file = paste0(projectdir, "/OCI-Ly1_TMD8_JIB04_500nM_v_DMSO_log2FC_conserved_znf_genes.pdf"), width = 6, height = 6, useDingbats = FALSE)  
gg
dev.off()



####################################################
#Perform GSEA on differentially expressed genes resulting from knockdown of ZNF587
####################################################

# Import gene sets created from Martins et al 2024
shZNF587_up <- read.gmt(gmtfile = paste0(projectdir, "/Martins_et_al_2024/rna_seq_analysis/shZNF587_D6_conserved_upregulated_genes.gmt")) 
shZNF587_down <- read.gmt(gmtfile = paste0(projectdir, "/Martins_et_al_2024/rna_seq_analysis/shZNF587_D6_conserved_downregulated_genes.gmt")) 

# OCI-Ly1 500 nM dose
fn <- paste0(projectdir, "/DESeq/OCI-Ly1_JIB04_500nM_v_DMSO.tsv")
df.deg <- fread(fn, data.table = FALSE)

geneList <- setNames(df.deg$log2FoldChange, df.deg$gene )
gseares <- GSEA(geneList, TERM2GENE = shZNF587_up, pAdjustMethod = "none", pvalueCutoff = 1)

pdf(file = paste0(projectdir, "/OCI-Ly1_JIB04_500nM_versus_DMSO_GSEA_shZNF587_genes_up_NES_1.67.pdf"), width = 3.5, height = 3.5)  
gseaplot2(gseares, 
          geneSetID = 1, 
          title = "JIB-04 500 nM versus DMSO treated OCI-Ly1 cells\nUp-regulated genes from ZNF587 knockdown (Martins et al 2024)", 
          color = "black", pvalue_table = TRUE,
          rel_heights = c(1,0.1,.4), base_size = 11)
dev.off()

core_enrichment_Ly1 <- gseares["shZNF587_UP", "core_enrichment"] %>% str_split(pattern = "/") %>% unlist()


# OCI-Ly1 150 nM dose
fn <- paste0(projectdir, "/DESeq/OCI-Ly1_JIB04_150nM_v_DMSO.tsv")
df.deg <- fread(fn, data.table = FALSE)

geneList <- setNames(df.deg$log2FoldChange, df.deg$gene )
gseares <- GSEA(geneList, TERM2GENE = shZNF587_up, pAdjustMethod = "none", pvalueCutoff = 1)

pdf(file = paste0(projectdir, "/OCI-Ly1_JIB04_150nM_versus_DMSO_GSEA_shZNF587_genes_up_NES_1.43.pdf"), width = 3.5, height = 3.5)  
gseaplot2(gseares, 
          geneSetID = 1, 
          title = "JIB-04 150 nM versus DMSO treated OCI-Ly1 cells\nUp-regulated genes from ZNF587 knockdown (Martins et al 2024)", 
          color = "black", pvalue_table = TRUE,
          rel_heights = c(1,0.1,.4), base_size = 11)
dev.off()


# TMD8 500 nM dose
fn <- paste0(projectdir, "/DESeq/TMD8_JIB04_500nM_v_DMSO.tsv")
df.deg <- fread(fn, data.table = FALSE)

geneList <- setNames(df.deg$log2FoldChange, df.deg$gene )
gseares <- GSEA(geneList, TERM2GENE = shZNF587_up, pAdjustMethod = "none", pvalueCutoff = 1)

pdf(file = paste0(projectdir, "/TMD8_JIB04_500nM_versus_DMSO_GSEA_shZNF587_genes_up_NES_1.23.pdf"), width = 3.5, height = 3.5)  
gseaplot2(gseares, 
          geneSetID = 1, 
          title = "JIB-04 500 nM versus DMSO treated TMD8 cells\nUp-regulated genes from ZNF587 knockdown (Martins et al 2024)", 
          color = "black", pvalue_table = TRUE,
          rel_heights = c(1,0.1,.4), base_size = 11)
dev.off()

core_enrichment_tmd8 <- gseares["shZNF587_UP", "core_enrichment"] %>% str_split(pattern = "/") %>% unlist()
intersect(core_enrichment_Ly1, core_enrichment_tmd8)

# TMD8 150 nM dose
fn <- paste0(projectdir, "/DESeq/TMD8_JIB04_150nM_v_DMSO.tsv")
df.deg <- fread(fn, data.table = FALSE)

geneList <- setNames(df.deg$log2FoldChange, df.deg$gene )
gseares <- GSEA(geneList, TERM2GENE = shZNF587_up, pAdjustMethod = "none", pvalueCutoff = 1)

pdf(file = paste0(projectdir, "/TMD8_JIB04_150nM_versus_DMSO_GSEA_shZNF587_genes_up_NES_1.33.pdf"), width = 3.5, height = 3.5)  
gseaplot2(gseares, 
          geneSetID = 1, 
          title = "JIB-04 150 nM versus DMSO treated TMD8 cells\nUp-regulated genes from ZNF587 knockdown (Martins et al 2024)", 
          color = "black", pvalue_table = TRUE,
          rel_heights = c(1,0.1,.4), base_size = 11)
dev.off()


