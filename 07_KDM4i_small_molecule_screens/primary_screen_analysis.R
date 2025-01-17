#Project: Lymphoma KDM4 Manuscript
#Objective: Analyze primary ICCB KDM4 small molecule screen on histone peptide substrates
#Mohamad Najia

library(data.table)
library(dplyr)
library(ggplot2)
library(ggExtra)


#################################################### 
#Primary histone peptide screen analysis
####################################################

## Set up environment
project_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
fn <- paste0(project_dir, "/primary_histone_peptide_screen_data.tsv")
df.processed <- fread(fn, data.table = FALSE)

## Process the data for plotting
df.processed$index <- 1:dim(df.processed)[1]
df.processed$classification <- df.processed$Sample_Type

df.processed[df.processed$classification == "Unknown", "classification"] <- ""
df.processed[df.processed$classification == "DMSO Neg Control", "classification"] <- "Negative Control"
df.processed[df.processed$classification == "CPI-455 Pos Control", "classification"] <- "Positive Control"
df.processed[df.processed$Hit_Both == "Y", "classification"] <- "Hit"
df.processed$classification <- factor(df.processed$classification, levels = c("", "Positive Control", "Hit", "Negative Control"))


## Plot primary screen hits
pdf(paste0(project_dir, "/primary_histone_peptide_screen.pdf"), width = 6, height = 3.5, useDingbats = FALSE)

ggplot(df.processed %>% arrange(classification), aes(x=index, y=Average_Zscore)) + 
  geom_point(aes(color = classification, alpha = classification), stroke=NA) + 
  scale_x_continuous(breaks = c(1e4, 2e4, 3e4, 4e4, 5e4)) + 
  scale_color_manual(values = c("grey", "red", "blue", "black")) + #c("blue", "black", "red", "grey")
  scale_alpha_manual(values = c(0.15,0.25,0.4,0.5)) + 
  geom_hline(yintercept = -3, linetype = "dashed", color = "black") + 
  xlab("Small Molecule ID") +
  ylab("Average Z-score") + 
  ggtitle("Primary Peptide Screen") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))

dev.off()


## Plot Z-factor for histone peptide assay
df.controls <- df.processed %>% filter(Sample_Type %in% c("DMSO Neg Control", "CPI-455 Pos Control"), Assay_Plate %in% c("3967", "3968", "3969", "3970", "3971", "3972") )
df.controls$Sample_Type <- factor(df.controls$Sample_Type, levels = c("DMSO Neg Control", "CPI-455 Pos Control"))
positives <- df.controls[df.controls$Sample_Type == "CPI-455 Pos Control", "Crosstalk-Corrected_Luminescence_Rep1"]
negatives <- df.controls[df.controls$Sample_Type == "DMSO Neg Control", "Crosstalk-Corrected_Luminescence_Rep1"]

zfactor <- 1 - ( ( 3*(sd(positives) + sd(negatives)) ) / abs(mean(positives) - mean(negatives)) )

pdf(paste0(project_dir, "/primary_histone_peptide_screen_Z_factor.pdf"), width = 3, height = 2, useDingbats = FALSE)

ggplot(df.controls, aes(x = `Crosstalk-Corrected_Luminescence_Rep1`, fill = Sample_Type)) + 
  geom_histogram( bins = 20) + 
  xlab("Luminescence") + 
  ylab("") + 
  ggtitle("Z-factor Peptide Assay") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, label = paste0("Z' = ", round(zfactor, digits = 3))) 

dev.off()


## Plot primary screen replicate reproducibility
pearson_cor <- cor(df.processed$`Crosstalk-Corrected_Luminescence_Rep1`, df.processed$`Crosstalk-Corrected_Luminescence_Rep2`, method = "pearson")

pdf(paste0(project_dir, "/primary_histone_peptide_screen_replicate_reproducibility_lum.pdf"), width = 5, height = 3.5, useDingbats = FALSE)

ggplot(df.processed %>% arrange(classification), aes(x=`Crosstalk-Corrected_Luminescence_Rep1`, y=`Crosstalk-Corrected_Luminescence_Rep2`)) + 
  geom_point(aes(color = classification, alpha = classification)) + 
  scale_color_manual(values = c("grey", "red", "blue", "black")) + #c("blue", "black", "red", "grey")
  scale_alpha_manual(values = c(0.15,0.15,0.15,0.15)) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") + 
  xlab("Luminescence Replicate 1") +
  ylab("Luminescence Replicate 2") + 
  ggtitle("Primary Peptide Screen Replicate Reproducibility") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 1, label = paste0("r = ", round(pearson_cor, digits = 3))) 

dev.off()


pearson_cor <- cor(df.processed$Zscore_Rep1, df.processed$Zscore_Rep2, method = "pearson")

pdf(paste0(project_dir, "/primary_histone_peptide_screen_replicate_reproducibility_zscore.pdf"), width = 5, height = 3.5, useDingbats = FALSE)

ggplot(df.processed %>% arrange(classification), aes(x=Zscore_Rep1, y=Zscore_Rep2)) + 
  geom_point(aes(color = classification, alpha = classification)) + 
  scale_color_manual(values = c("grey", "red", "blue", "black")) + #c("blue", "black", "red", "grey")
  scale_alpha_manual(values = c(0.15,0.15,0.15,0.15)) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") + 
  xlab("Z-score Replicate 1") +
  ylab("Z-score Replicate 2") + 
  ggtitle("Primary Peptide Screen Replicate Reproducibility") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) + 
  annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 1, label = paste0("r = ", round(pearson_cor, digits = 3))) 

dev.off()



