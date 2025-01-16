#Project: Lymphoma KDM4 Manuscript
#Author: Mohamad Najia
#Objective: analyze TEtranscripts quantification of transposable elements from RNA-Seq on JIB-04 treated DLBCL cells

library(dplyr)
library(data.table)
library(ggplot2)
library(scales)
library(stringr)


####################################################
#Set-up Environment
####################################################

# Initialize variables
projectdir <- dirname(rstudioapi::getSourceEditorContext()$path)

# Declare colormaps
cols <- c("Significant" = "red", "Non-significant" = "darkgrey", "Up-regulated" = "#00B2FF", "Down-regulated" = "#00B2FF")



####################################################
#TEtranscripts Analysis
####################################################

## JIB04 150 nM versus DMSO comparisons

# OCI-Ly1 cells
fn <- paste0(projectdir, "/TEtranscripts_output_OCI-Ly1/JIB04_150nM_v_DMSO_TE_analysis.txt")
df.te <- fread(fn, data.table = FALSE)

df.te$nlogFDR <- -1*log10(df.te$padj)
df.te$DEG <- "Non-significant"
df.te[((df.te$log2FoldChange > 1) & (df.te$padj > 0.05)), "DEG"] <- "Up-regulated"
df.te[((df.te$log2FoldChange < -1) & (df.te$padj > 0.05)), "DEG"] <- "Down-regulated"
df.te[((abs(df.te$log2FoldChange) > 1) & (df.te$padj < 0.05)), "DEG"] <- "Significant"

gg1 <- ggplot(df.te, aes(x = log2FoldChange, y = nlogFDR, color = DEG, labels = V1)) +
  scale_colour_manual(values = cols) +
  geom_point(size = 2.5, alpha = 1, na.rm = T) +
  theme_bw(base_size = 14) + 
  xlim(-10, 10) + 
  ylim(0, 20) + 
  ggtitle("OCI-Ly1 JIB04 150 nM v. DMSO TE Expression") + 
  geom_hline(yintercept = -1*log10(0.05), colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = 1, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = -1, colour="#990000", linetype="dashed") +
  xlab("log2(Fold Change)") + 
  ylab("-log10(FDR)") +
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 14)) + 
  geom_text_repel(data=subset(df.te, DEG == "Significant"),
                  aes(log2FoldChange,nlogFDR,label=V1))

# TMD8 cells
fn <- paste0(projectdir, "/TEtranscripts_output_TMD8/JIB04_150nM_v_DMSO_TE_analysis.txt")
df.te <- fread(fn, data.table = FALSE)

df.te$nlogFDR <- -1*log10(df.te$padj)
df.te$DEG <- "Non-significant"
df.te[((df.te$log2FoldChange > 1) & (df.te$padj > 0.05)), "DEG"] <- "Up-regulated"
df.te[((df.te$log2FoldChange < -1) & (df.te$padj > 0.05)), "DEG"] <- "Down-regulated"
df.te[((abs(df.te$log2FoldChange) > 1) & (df.te$padj < 0.05)), "DEG"] <- "Significant"

gg2 <- ggplot(df.te, aes(x = log2FoldChange, y = nlogFDR, color = DEG, labels = V1)) +
  scale_colour_manual(values = cols) +
  geom_point(size = 2.5, alpha = 1, na.rm = T) +
  theme_bw(base_size = 14) + 
  xlim(-10, 10) + 
  ylim(0, 20) + 
  ggtitle("TMD8 JIB04 150 nM v. DMSO TE Expression") + 
  geom_hline(yintercept = -1*log10(0.05), colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = 1, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = -1, colour="#990000", linetype="dashed") +
  xlab("log2(Fold Change)") + 
  ylab("-log10(FDR)") +
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 14)) + 
  geom_text_repel(data=subset(df.te, DEG == "Significant"),
                  aes(log2FoldChange,nlogFDR,label=V1))

pdf(file = paste0(projectdir, "/JIB04_150nM_v_DMSO_TEtranscripts_volcano.pdf"), width = 8, height = 4)  
cowplot::plot_grid(gg1, gg2)
dev.off()



## JIB04 500 nM versus DMSO comparisons

# OCI-Ly1 cells
fn <- paste0(projectdir, "/TEtranscripts_output_OCI-Ly1/JIB04_500nM_v_DMSO_TE_analysis.txt")
df.te <- fread(fn, data.table = FALSE)

df.te$nlogFDR <- -1*log10(df.te$padj)
df.te$DEG <- "Non-significant"
df.te[((df.te$log2FoldChange > 1) & (df.te$padj > 0.05)), "DEG"] <- "Up-regulated"
df.te[((df.te$log2FoldChange < -1) & (df.te$padj > 0.05)), "DEG"] <- "Down-regulated"
df.te[((abs(df.te$log2FoldChange) > 1) & (df.te$padj < 0.05)), "DEG"] <- "Significant"

gg1 <- ggplot(df.te, aes(x = log2FoldChange, y = nlogFDR, color = DEG, labels = V1)) +
  scale_colour_manual(values = cols) +
  geom_point(size = 2.5, alpha = 1, na.rm = T) +
  theme_bw(base_size = 14) + 
  xlim(-10, 10) + 
  ylim(0, 300) + 
  ggtitle("OCI-Ly1 JIB04 500 nM v. DMSO TE Expression") + 
  geom_hline(yintercept = -1*log10(0.05), colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = 1, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = -1, colour="#990000", linetype="dashed") +
  xlab("log2(Fold Change)") + 
  ylab("-log10(FDR)") +
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 14)) + 
  geom_text_repel(data=subset(df.te, DEG == "Significant"),
                  aes(log2FoldChange,nlogFDR,label=V1))

# Quantify significant (FDR < 0.05) differential TE classes
df.te$class <- str_split_fixed(df.te$V1, pattern = ":", n=3)[,3]
df.te$class <- do.call(rbind, lapply(df.te$class, function(x) gsub("\\?$", "", x)))
df.te[df.te$log2FoldChange > 0.75, "direction"] <- "up"
df.te[df.te$log2FoldChange < 0.75, "direction"] <- "down"

tmp <- df.te %>% filter(padj < 0.05 & abs(log2FoldChange) > 0.75) 
df.degs <- table(tmp$direction, tmp$class) %>% as.data.frame()
df.degs$Var1 <- factor(df.degs$Var1, levels = c("up", "down"))

p1 <- ggplot(df.degs, aes(x = as.factor(Var2), 
                    y = Freq * ((-1)^(Var1 == "down")), 
                    fill = Var1)) + 
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, colour="black", linetype="solid") + 
  scale_x_discrete(name = "") +
  scale_y_continuous(limits = c(-25,20), name = "Number of TE families") +
  scale_fill_manual(values = c("blue", "red")) + 
  ggtitle("OCI-Ly1 JIB04 500 nM v DMSO") + 
  theme_bw() +
  theme(legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour="black", angle = 30, vjust = 1, hjust=1), 
        axis.text.y = element_text(colour="black"))


# TMD8 cells
fn <- paste0(projectdir, "/TEtranscripts_output_TMD8/JIB04_500nM_v_DMSO_TE_analysis.txt")
df.te <- fread(fn, data.table = FALSE)

df.te$nlogFDR <- -1*log10(df.te$padj)
df.te$DEG <- "Non-significant"
df.te[((df.te$log2FoldChange > 1) & (df.te$padj > 0.05)), "DEG"] <- "Up-regulated"
df.te[((df.te$log2FoldChange < -1) & (df.te$padj > 0.05)), "DEG"] <- "Down-regulated"
df.te[((abs(df.te$log2FoldChange) > 1) & (df.te$padj < 0.05)), "DEG"] <- "Significant"

gg2 <- ggplot(df.te, aes(x = log2FoldChange, y = nlogFDR, color = DEG, labels = V1)) +
  scale_colour_manual(values = cols) +
  geom_point(size = 2.5, alpha = 1, na.rm = T) +
  theme_bw(base_size = 14) + 
  xlim(-10, 10) + 
  ylim(0, 300) + 
  ggtitle("TMD8 JIB04 500 nM v. DMSO TE Expression") + 
  geom_hline(yintercept = -1*log10(0.05), colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = 1, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = -1, colour="#990000", linetype="dashed") +
  xlab("log2(Fold Change)") + 
  ylab("-log10(FDR)") +
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 14)) + 
  geom_text_repel(data=subset(df.te, DEG == "Significant"),
                  aes(log2FoldChange,nlogFDR,label=V1),
                  max.overlaps = 20)

pdf(file = paste0(projectdir, "/JIB04_500nM_v_DMSO_TEtranscripts_volcano.pdf"), width = 8, height = 4)  
cowplot::plot_grid(gg1, gg2)
dev.off()

# Quantify significant (FDR < 0.05) differential TE classes
df.te$class <- str_split_fixed(df.te$V1, pattern = ":", n=3)[,3]
df.te$class <- do.call(rbind, lapply(df.te$class, function(x) gsub("\\?$", "", x)))
df.te[df.te$log2FoldChange > 0.75, "direction"] <- "up"
df.te[df.te$log2FoldChange < 0.75, "direction"] <- "down"

tmp <- df.te %>% filter(padj < 0.05 & abs(log2FoldChange) > 0.75) 
df.degs <- table(tmp$direction, tmp$class) %>% as.data.frame()
df.degs <- df.degs[-c(7,8),]
df.degs$Var1 <- factor(df.degs$Var1, levels = c("up", "down"))

p2 <- ggplot(df.degs, aes(x = as.factor(Var2), 
                          y = Freq * ((-1)^(Var1 == "down")), 
                          fill = Var1)) + 
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, colour="black", linetype="solid") + 
  scale_x_discrete(name = "") +
  scale_y_continuous(limits = c(-25,20), name = "Number of TE families") +
  scale_fill_manual(values = c("blue", "red")) + 
  ggtitle("TMD8 JIB04 500 nM v DMSO") + 
  theme_bw() +
  theme(legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour="black", angle = 30, vjust = 1, hjust=1), 
        axis.text.y = element_text(colour="black"))

pdf(file = paste0(projectdir, "/JIB04_500nM_v_DMSO_differential_TE_quant_by_class.pdf"), width = 8, height = 3.5)  
cowplot::plot_grid(p1, p2)
dev.off()


