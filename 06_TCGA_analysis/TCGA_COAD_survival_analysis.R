#Project: Lymphoma KDM4 Manuscript
#Objective: Perform survival analysis on COAD cancer patients from TCGA
#Mohamad Najia

library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RTCGA)
library(RTCGA.clinical)
library(survival)
library(survminer)


#############################
#Functions
#############################

scaleRows <- function(x) {
  rm <- rowMeans(x)
  x <- sweep(x, 1, rm)
  sx <- apply(x, 1, sd)
  x <- sweep(x, 1, sx, "/")
}



#############################
#Initialize environment 
#############################

# Initialize variables
working_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
coad_dir <- paste0(working_dir, "/TCGA-COAD/")

# Import TCGA colorectal cancer RNA-Seq data 
#Publication source: https://www.nature.com/articles/nature11252
#Download source: https://gdc.cancer.gov/about-data/publications/coadread_2012
fn <- paste0(coad_dir, "crc_270_gene_rpkm_datac.txt")
df.rnaseq <- fread(fn, data.table = FALSE)

rownames(df.rnaseq) <- df.rnaseq$Gene
df.rnaseq$Gene <- NULL
colnames(df.rnaseq) <- substr(colnames(df.rnaseq), 1, 12)

# Import TCGA colorectal cancer clinical data
survInfo <- survivalTCGA(RTCGA.clinical::COAD.clinical, extract.cols = "admin.disease_code")
commonIDs <- intersect(survInfo$bcr_patient_barcode, colnames(df.rnaseq))
survInfo.subset <- survInfo %>% filter(bcr_patient_barcode %in% commonIDs)



#############################
#Survival Analysis 
#############################

# Determine normalized, combined expression of KDM4A/C
df.rnaseq.subset <- df.rnaseq[grepl("KDM4A|KDM4C", rownames(df.rnaseq)), commonIDs]
df.z <- scaleRows(df.rnaseq.subset)
df.strata <- data.frame(avgz = colMeans(df.z))

# Stratify patients by KDM4A/C expression
survInfo.subset$avgz <- df.strata[survInfo.subset$bcr_patient_barcode, "avgz"]

quantilecutoffs <- c(1/3, 2/3)
quants <- quantile(survInfo.subset$avgz, quantilecutoffs)

survInfo.subset$expression <- ""
survInfo.subset[survInfo.subset$avgz <= quants[1], "expression"] <- "Low"
survInfo.subset[survInfo.subset$avgz >= quants[2], "expression"] <- "High"

kmTCGA(survInfo.subset, explanatory.names = "expression", pval = TRUE)

sfit <- kmTCGA(survInfo.subset %>% filter(expression %in% c("Low", "High")), 
       explanatory.names = "expression", 
       pval = TRUE,
       conf.int = FALSE,
       return.survfit = TRUE
       )

# Save Kaplan-Meier plot
pdf(paste0(working_dir, "/TCGA_COAD_KM_curve_KDM4AC_expression.pdf"), width = 3.25, height = 4, useDingbats = FALSE)

ggsurvplot(sfit$survfit,
           legend.title = "KDM4A/C Expression",
           legend.labs = c("Top 33%", "Bottom 33%"),
           pval = TRUE,
           conf.int = FALSE,
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = theme_cleantable(),
           palette = c("#E7B800", "#2E9FDF"),
           ggtheme = theme_bw(),
           font.tickslab = c("black")
) 

dev.off()



