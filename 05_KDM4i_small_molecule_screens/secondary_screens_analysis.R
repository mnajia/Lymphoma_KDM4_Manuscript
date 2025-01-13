#Project: Lymphoma KDM4 Manuscript
#Objective: Analyze secondary ICCB KDM4 small molecule validation screens on histone peptide and nucleosome substrates
#Mohamad Najia

library(data.table)
library(dplyr)
library(ggplot2)
library(ggExtra)
library(ggVennDiagram)
library(ComplexHeatmap)


#################################################### 
#Secondary histone peptide and nucleosome valdiation screens
####################################################

## Set up environment
project_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
fn <- paste0(project_dir, "/secondary_histone_peptide_nucleosome_screens_data.tsv")
df.processed <- fread(fn, data.table = FALSE)

## Process the data for plotting
df.processed[df.processed$Sample_Type == "CPI-455 Pos Control", "Hit_Type"] <- "Positive Control"
df.processed[df.processed$Sample_Type == "DMSO Neg Control", "Hit_Type"] <- "Negative Control"
df.processed[df.processed$Hit_Type == "", "Hit_Type"] <- "No hit: Peptide & Nucleosome"
df.processed$Hit_Type <- factor(df.processed$Hit_Type, levels = c("Negative Control", "Positive Control", "No hit: Peptide & Nucleosome", "Nucleosome Only", "Peptide Only", "Peptide & Nucleosome"))

## Plot secondary nucleosome and peptide screens
pdf(paste0(project_dir, "/secondary_screens_peptide_v_nucleosome.pdf"), width = 7, height = 4, useDingbats = FALSE)

ggplot(df.processed, aes(x=Nuc_Avg_Zscore, y=Peptide_Avg_Zscore)) + 
  geom_point(aes(color = Hit_Type, alpha = Hit_Type), stroke = NA, size = 2) + 
  scale_alpha_manual(values = c(0.2,0.5,0.5,0.5,0.5,0.75)) + 
  geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.25, color = "black") + 
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.25, color = "black") + 
  xlab("Average Z-score, Nucleosome Screen") +
  ylab("Average Z-score, Peptide Screen") + 
  ggtitle("Secondary Histone Peptide & Nucleosome Screens") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))

dev.off()

## Plot replicate reproducibility of the secondary histone peptide screen
pearson_cor <- cor(df.processed$Peptide_Lum_Rep1, df.processed$Peptide_Lum_Rep2, method = "pearson")

pdf(paste0(project_dir, "/secondary_histone_peptide_screen_replicate_reproducibility_lum.pdf"), width = 6, height = 3.5, useDingbats = FALSE)

ggplot(df.processed, aes(x=Peptide_Lum_Rep1, y=Peptide_Lum_Rep2)) + 
  geom_point(aes(color = Hit_Type, alpha = Hit_Type), stroke = NA) + 
  scale_alpha_manual(values = c(0.2,0.5,0.5,0.5,0.5,0.75)) + 
  scale_x_continuous(limits = c(0, 2e5)) + 
  scale_y_continuous(limits = c(0, 2e5)) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") + 
  xlab("Luminescence Replicate 1") +
  ylab("Luminescence Replicate 2") + 
  ggtitle("Secondary Peptide Screen Replicate Reproducibility") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) + 
  annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 1, label = paste0("r = ", round(pearson_cor, digits = 3))) 

dev.off()


pearson_cor <- cor(df.processed$Peptide_Zscore_Rep1, df.processed$Peptide_Zscore_Rep2, method = "pearson")

pdf(paste0(project_dir, "/secondary_histone_peptide_screen_replicate_reproducibility_zscore.pdf"), width = 6, height = 3.5, useDingbats = FALSE)

ggplot(df.processed, aes(x=Peptide_Zscore_Rep1, y=Peptide_Zscore_Rep2)) + 
  geom_point(aes(color = Hit_Type, alpha = Hit_Type), stroke = NA) + 
  scale_alpha_manual(values = c(0.2,0.5,0.5,0.5,0.5,0.75)) + 
  scale_x_continuous(limits = c(-25, 10)) + 
  scale_y_continuous(limits = c(-25, 10)) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") + 
  xlab("Z-score Replicate 1") +
  ylab("Z-score Replicate 2") + 
  ggtitle("Secondary Peptide Screen Replicate Reproducibility") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) + 
  annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 1, label = paste0("r = ", round(pearson_cor, digits = 3))) 

dev.off()

## Plot replicate reproducibility of the secondary nucleosome screen 
pearson_cor <- cor(df.processed$Nuc_Lum_Rep1, df.processed$Nuc_Lum_Rep2, method = "pearson")

pdf(paste0(project_dir, "/secondary_nucleosome_screen_replicate_reproducibility_lum.pdf"), width = 6, height = 3.5, useDingbats = FALSE)

ggplot(df.processed, aes(x=Nuc_Lum_Rep1, y=Nuc_Lum_Rep2)) + 
  geom_point(aes(color = Hit_Type, alpha = Hit_Type), stroke = NA) + 
  scale_alpha_manual(values = c(0.2,0.5,0.5,0.5,0.5,0.75)) + 
  scale_x_continuous(limits = c(0, 4e5)) + 
  scale_y_continuous(limits = c(0, 4e5)) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") + 
  xlab("Luminescence Replicate 1") +
  ylab("Luminescence Replicate 2") + 
  ggtitle("Secondary Nucleosome Screen Replicate Reproducibility") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) + 
  annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 1, label = paste0("r = ", round(pearson_cor, digits = 3))) 

dev.off()


pearson_cor <- cor(df.processed$Nuc_Zscore_Rep1, df.processed$Nuc_Zscore_Rep2, method = "pearson")

pdf(paste0(project_dir, "/secondary_nucleosome_screen_replicate_reproducibility_zscore.pdf"), width = 6, height = 3.5, useDingbats = FALSE)

ggplot(df.processed, aes(x=Nuc_Zscore_Rep1, y=Nuc_Zscore_Rep2)) + 
  geom_point(aes(color = Hit_Type, alpha = Hit_Type), stroke = NA) + 
  scale_alpha_manual(values = c(0.2,0.5,0.5,0.5,0.5,0.75)) + 
  scale_x_continuous(limits = c(-20, 30)) + 
  scale_y_continuous(limits = c(-20, 30)) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") + 
  xlab("Z-score Replicate 1") +
  ylab("Z-score Replicate 2") + 
  ggtitle("Secondary Nucleosome Screen Replicate Reproducibility") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) + 
  annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 1, label = paste0("r = ", round(pearson_cor, digits = 3))) 

dev.off()



#################################################### 
#Secondary AlphaLISA TruHits screens
####################################################
#if the TruHits luminescence is greater than 2 standard deviations of the mean of the negative controls then it is a TruHit and not a false positive that infers with the AlphaLISA assay

## Import data
fn <- paste0(project_dir, "/secondary_alphalisa_truhits_screen_data.tsv")
df.hits <- fread(fn, data.table = FALSE)

## Plot overlaps of hits across all secondary screens
x <- list(Peptide = df.processed[df.processed$Peptide_Screen_Hit == "Y", "Vendor_ID"], 
          Nucleosome = df.processed[df.processed$Nuc_Screen_Hit == "Y", "Vendor_ID"],
          TruHits = df.hits[df.hits$Hit_both == "Y", "Vendor_ID"])

pdf(paste0(project_dir, "/secondary_screens_venn_diagram.pdf"), width = 4, height = 4, useDingbats = FALSE)

ggVennDiagram(x, color = "black", lwd = 0.5, lty = 1) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")

dev.off()

## Process data for plotting
df.hits <- df.hits[df.hits$Sample_Type %in% c("DMSO Neg Control", "Primary Screen Hit"),]
df.hits[df.hits$Sample_Type == "DMSO Neg Control", "Hit_both"] <- "Y"
df.hits[df.hits$Sample_Type == "DMSO Neg Control", "Secondary_Screens_Hit_Type"] <- "Negative Control"
df.hits <- df.hits[df.hits$Secondary_Screens_Hit_Type %in% c("Negative Control", "Nucleosome Only", "Peptide Only", "Peptide & Nucleosome"),]

pdf(paste0(project_dir, "/secondary_alphalisa_truhits_screen.pdf"), width = 6, height = 4, useDingbats = FALSE)

ggplot(df.hits, aes(x=TruHit_Zscore_Rep1, y=TruHit_Zscore_Rep2)) + 
  geom_point(aes(color = Secondary_Screens_Hit_Type, alpha = Hit_both), stroke = NA, size = 3) + 
  scale_alpha_manual(values = c(0.15,1)) + 
  scale_x_continuous(limits = c(-20, 10)) + 
  scale_y_continuous(limits = c(-20, 10)) + 
  geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.25, color = "black") + 
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.25, color = "black") + 
  geom_hline(yintercept = -2, linetype = "dashed", linewidth = 0.25, color = "black") + 
  geom_vline(xintercept = -2, linetype = "dashed", linewidth = 0.25, color = "black") + 
  xlab("Z-score, Replicate 1") +
  ylab("Z-score, Replicate 2") + 
  ggtitle("AlphaLISA TruHits Screen") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) 

dev.off()

## Plot secondary TruHit screen replicate reproducibility
pearson_cor <- cor(df.hits$TruHit_Lum_Rep1, df.hits$TruHit_Lum_Rep2, method = "pearson")

pdf(paste0(project_dir, "/secondary_alphalisa_truhits_screen_replicate_reproducibility_lum.pdf"), width = 6, height = 3.5, useDingbats = FALSE)

ggplot(df.hits, aes(x=TruHit_Lum_Rep1, y=TruHit_Lum_Rep2)) + 
  geom_point(aes(color = Secondary_Screens_Hit_Type, alpha = Secondary_Screens_Hit_Type), stroke = NA) + 
  scale_alpha_manual(values = c(0.2,0.5,0.5,0.5,0.5,0.75)) + 
  scale_x_continuous(limits = c(0, 1.2e6)) + 
  scale_y_continuous(limits = c(0, 1.2e6)) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") + 
  xlab("Luminescence Replicate 1") +
  ylab("Luminescence Replicate 2") + 
  ggtitle("Secondary AlphaLISA TruHits Screen Replicate Reproducibility") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) + 
  annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 1, label = paste0("r = ", round(pearson_cor, digits = 3))) 

dev.off()



