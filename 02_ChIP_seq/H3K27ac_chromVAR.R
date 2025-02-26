#Project: Lymphoma KDM4 Manuscript
#Objective: perform chromVAR TF analysis on H3K27ac ChIP-Seq peaks to determine TF motif archetypes with significant variability across conditions
#Author: Mohamad Najia

library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(irlba)
library(annotables)
library(BuenColors)
library(chromVAR)
library(ChrAccR)
library(SummarizedExperiment)
library(BiocParallel)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(circlize)
library(ComplexHeatmap)
library(GenomicRanges)
library(ggplot2)
library(ggrepel)
register(MulticoreParam(2))
set.seed(14651)


############ Functions ############ 
scaleRows <- function(x) {
  rm <- rowMeans(x)
  x <- sweep(x, 1, rm)
  sx <- apply(x, 1, sd)
  x <- sweep(x, 1, sx, "/")
}
###################################


### Initialize and Import Data ###

#declare colormaps
pal_rna <- colorRampPalette(c("#352A86","#343DAE","#0262E0","#1389D2",
                              "#2DB7A3","#A5BE6A","#F8BA43","#F6DA23","#F8FA0D"))(100)
pal_atac <- colorRampPalette(c('#3361A5', '#248AF3', '#14B3FF', 
                               '#88CEEF', '#C1D5DC', '#EAD397', 
                               '#FDB31A','#E42A2A', '#A31D1D'))(100)

#set up environment
project_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
output_dir <- paste0(project_dir, "/chromVAR/")

#import count matrix
fn <- paste0(project_dir, "/H3K27ac.consensus_peaks.featureCounts.txt")
counts.df <- fread(fn, data.table = FALSE)

df.peaks <- counts.df[,c(2,3,4)]
colnames(df.peaks) <- c("chr", "start", "end")

df.countMatrix <- counts.df[,c("DMSO_H3K27ac_REP1.mLb.clN.sorted.bam", "DMSO_H3K27ac_REP2.mLb.clN.sorted.bam", "JIB04_500nM_H3K27ac_REP1.mLb.clN.sorted.bam", "JIB04_500nM_H3K27ac_REP2.mLb.clN.sorted.bam")]
colnames(df.countMatrix) <- c("DMSO_H3K27ac_rep1", "DMSO_H3K27ac_rep2", "JIB04_500nM_H3K27ac_rep1", "JIB04_500nM_H3K27ac_rep2")



### Apply chromVAR to all peaks ###

#filter peaks to include only standard chromosomes 
inds <- df.peaks$chr %in% paste0("chr", c(1:22, "X", "Y"))
gr.peaks <- makeGRangesFromDataFrame(df.peaks[inds,])
df.countMatrix <- df.countMatrix[inds,]

#sort the peaks
inds <- order(gr.peaks)
gr.peaks <- gr.peaks[inds]
df.countMatrix <- df.countMatrix[inds,]

#chromVAR setup
SE <- SummarizedExperiment(
  rowRanges =  gr.peaks,
  colData = data.frame(Sample = colnames(df.countMatrix), gRNA = c("DMSO", "DMSO", "JIB04", "JIB04")),
  assays = list(counts = as.matrix(df.countMatrix))
)

SE <- filterPeaks(SE)
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg38)


#match motifs to peaks using Vierstra motif archetype database
varchetypes <- readRDS(paste0(project_dir, "/Vierstra_Archetype_Motifs_v2.1.rds"))
motif_ix <- matchMotifs(varchetypes, SE, genome = BSgenome.Hsapiens.UCSC.hg38)

#calculate motif deviations
dev <- computeDeviations(object = SE, annotations = motif_ix)

#get the chromVAR deviations and z-scores
cvdevs <- deviations(dev)
cvzscores <- deviationScores(dev)

#calculate variability across samples
variabilityAll <- computeVariability(dev)
data.frame(variabilityAll) %>% arrange(desc(variability))-> motifVariability
motifVariability$rank <- 1:dim(variabilityAll)[1]

#plot TF variability 
pdf(paste0(output_dir, "chromVAR_overall_vierstra_archetype_motif_variability.pdf"), width = 5, height = 5)
ggplot(motifVariability, aes(x = rank, y = variability, label = name)) + 
  geom_point(size = 0.75) +
  pretty_plot(fontsize = 8) + L_border() + 
  labs(x = "Motif Rank", y = "Accessibility Variability") +
  theme_bw() + 
  geom_label_repel(data = subset(motifVariability, variability > 4),
                   box.padding   = 1, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   max.overlaps = Inf)
dev.off()


#Vierstra archetype motifs: identify motif archetypes with significant variability 
p_cutoff <- 0.05
var_cutoff <- 3 #6

variabilityTop <- filter(variabilityAll, p_value_adj < p_cutoff, variability > var_cutoff)
variabilityTop <- variabilityTop[order(variabilityTop$variability, decreasing = TRUE),]

significantArchetypes <- variabilityTop %>% rownames()
df.chromvar.zscores <- cvzscores[significantArchetypes,]

anames <- str_split_fixed(significantArchetypes, pattern = "\\|", n = 3)
anamesformat <- paste0("(", anames[,1], ") ", anames[,2])

rownames(df.chromvar.zscores) <- anamesformat
x <- variabilityAll[significantArchetypes,"variability"]
names(x) <- anamesformat

#output motif archetypes with significant variability
variabilitySign <- filter(variabilityAll, p_value_adj < p_cutoff)
plotVariability(variabilitySign)
saveRDS(variabilitySign, file = paste0(output_dir, "H3K27ac_chromVAR_vierstra_archetype_motifs_padj0.05_variability.rds"))
write.table(variabilitySign, 
            file = paste0(output_dir, "H3K27ac_chromVAR_vierstra_archetype_motifs_padj0.05_variability.txt"), 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


#plot z-scores of motif archetypes with significant variability across samples 
color_maps <- c("DMSO" = "#a3062a", "JIB04" = "#4777b2")

column_ha <- HeatmapAnnotation(gRNA = c("DMSO", "DMSO", "JIB04", "JIB04"), 
                               col = list(gRNA = color_maps), 
                               annotation_name_side = "left",
                               border = TRUE)

row_ha = rowAnnotation(Variability = anno_barplot(x))

ht1 <- Heatmap(df.chromvar.zscores, 
               col = colorRamp2(c(-8, -6, -4, -2, 0, 2, 4, 6, 8), c('#3361A5', '#248AF3', '#14B3FF', '#88CEEF', '#C1D5DC', '#EAD397', '#FDB31A','#E42A2A', '#A31D1D')), 
               top_annotation = column_ha, 
               right_annotation = row_ha,
               show_row_names = TRUE, 
               show_column_names = FALSE, 
               cluster_rows = FALSE, 
               cluster_columns = FALSE, 
               show_column_dend = FALSE,
               row_km = 2,
               row_title = NULL,
               column_title = "H3K27ac Peaks",
               #row_names_gp = gpar(fontsize = 10),
               #column_names_gp = gpar(fontsize = 10),
               #row_names_side = "left", 
               #row_dend_side = "right",
               #rect_gp = gpar(col = "black", lwd = 0.2),
               heatmap_legend_param = list(direction = "horizontal", border = TRUE),
               name = "chromVAR Z-score", 
               border = TRUE)

pdf(paste0(output_dir, "H3K27ac_chromVAR_vierstra_archetype_motifs_padj0.05_variability3_zscore_heatmap_consensus_peaks.pdf"), width = 5.5, height = 7)
draw(ht1, merge_legend = TRUE)
dev.off()


#plot chromVAR scores for H3K27ac and H3K4me3 together
color_maps <- c("DMSO" = "#a3062a", "JIB04" = "#4777b2")

column_ha <- HeatmapAnnotation(gRNA = c("DMSO", "DMSO", "JIB04", "JIB04"), 
                               col = list(gRNA = color_maps), 
                               annotation_name_side = "left",
                               border = TRUE)


h3k4me3.variability <- readRDS(file = paste0(output_dir, "H3K4me3_chromVAR_vierstra_archetype_motifs_variability.rds"))
h3k4me3variabilitySign <- filter(h3k4me3.variability, p_value_adj < p_cutoff)
overlapArchetypes <- intersect(rownames(variabilityTop), rownames(h3k4me3variabilitySign))


df.chromvar.zscores <- cvzscores[overlapArchetypes,]
anames <- str_split_fixed(overlapArchetypes, pattern = "\\|", n = 3)
anamesformat <- paste0("(", anames[,1], ") ", anames[,2])
rownames(df.chromvar.zscores) <- anamesformat
x <- variabilityAll[overlapArchetypes,"variability"]
names(x) <- anamesformat

row_ha = rowAnnotation(Variability = anno_barplot(x))

ht1 <- Heatmap(df.chromvar.zscores, 
               col = colorRamp2(c(-8, -6, -4, -2, 0, 2, 4, 6, 8), c('#3361A5', '#248AF3', '#14B3FF', '#88CEEF', '#C1D5DC', '#EAD397', '#FDB31A','#E42A2A', '#A31D1D')), 
               top_annotation = column_ha, 
               right_annotation = row_ha,
               show_row_names = TRUE, 
               show_column_names = FALSE, 
               cluster_rows = FALSE, 
               cluster_columns = FALSE, 
               show_column_dend = FALSE,
               row_km = 2,
               row_title = NULL,
               column_title = "H3K27ac Peaks",
               #row_names_gp = gpar(fontsize = 10),
               #column_names_gp = gpar(fontsize = 10),
               #row_names_side = "left", 
               #row_dend_side = "right",
               #rect_gp = gpar(col = "black", lwd = 0.2),
               heatmap_legend_param = list(direction = "horizontal", border = TRUE),
               name = "chromVAR Z-score", 
               border = TRUE)


h3k4me3.cvzscores <- readRDS(file = paste0(output_dir, "H3K4me3_chromVAR_vierstra_archetype_motifs_cvzscores.rds"))
df.chromvar.zscores <- h3k4me3.cvzscores[overlapArchetypes,]

anames <- str_split_fixed(overlapArchetypes, pattern = "\\|", n = 3)
anamesformat <- paste0("(", anames[,1], ") ", anames[,2])
rownames(df.chromvar.zscores) <- anamesformat
y <- h3k4me3variabilitySign[overlapArchetypes,"variability"]
names(y) <- anamesformat

row_ha2 = rowAnnotation(Variability = anno_barplot(y))

ht2 <- Heatmap(df.chromvar.zscores, 
               col = colorRamp2(c(-8, -6, -4, -2, 0, 2, 4, 6, 8), c('#3361A5', '#248AF3', '#14B3FF', '#88CEEF', '#C1D5DC', '#EAD397', '#FDB31A','#E42A2A', '#A31D1D')), 
               top_annotation = column_ha, 
               right_annotation = row_ha2,
               show_row_names = TRUE, 
               show_column_names = FALSE, 
               cluster_rows = FALSE, 
               cluster_columns = FALSE, 
               show_column_dend = FALSE,
               #row_km = 2,
               row_title = NULL,
               column_title = "H3K4me3 Peaks",
               #row_names_gp = gpar(fontsize = 10),
               #column_names_gp = gpar(fontsize = 10),
               #row_names_side = "left", 
               #row_dend_side = "right",
               #rect_gp = gpar(col = "black", lwd = 0.2),
               heatmap_legend_param = list(direction = "horizontal", border = TRUE),
               name = "chromVAR Z-score", 
               border = TRUE)

pdf(paste0(output_dir, "H3K27ac_H3K4me3_chromVAR_vierstra_archetype_motifs_padj0.05_variability3_zscore_heatmap_consensus_peaks.pdf"), width = 6.5, height = 7)
draw(ht1+ht2, merge_legend = TRUE)
dev.off()


