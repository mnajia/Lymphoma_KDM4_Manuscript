#Project: Lymphoma KDM4 Manuscript
#Objective: Analyze Cayman Epigenetics Library small molecule screen from OCI-Ly1 and HBL-1 cells
#Mohamad Najia

library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(ComplexHeatmap)


#################################################### 
#Analysis Outline
####################################################
#plot distribution of DMSO controls (total and per plate)
#calculate log2 fold changes of the drugs and controls relative to DMSO
#determine statistical significance relative to the DMSO CTG distribution
#make volcano plots for each dose and each line
#plot log2 fold changes for OCI-Ly1 versus HBL-1 for each dose 
#plot the % inhibition for the hits across doses and in both cell lines as a heatmap



#################################################### 
#Functions Declarations
####################################################

processScreen <- function(df.screen) {
  
  #isolate DMSO controls and small molecules 
  df.dmso <- df.screen %>% filter(Molecule == "DMSO")
  df.sm <- df.screen %>% filter(Molecule != "DMSO")
  
  #iterate across small molecules in the screen
  resfinal <- lapply(split(df.sm, df.sm$Molecule), function(x) {
    
    #iterate across each dose of the small molecule 
    res <- lapply(split(x, x$Dose), function(y) {
      
      #statistical test to determine if there is a difference in the means of small molecule versus DMSO 
      ans <- t.test(y$`CTG Signal`, df.dmso$`CTG Signal`, paired = FALSE, alternative = "two.sided")
      #ans <- wilcox.test(y$`CTG Signal`, df.dmso$`CTG Signal`, paired = FALSE, alternative = "two.sided")
      pval <- ans$p.value
      
      #determine log2 fold change
      log2fc <- log2( mean(y$`CTG Signal`) / mean(df.dmso$`CTG Signal`) )
      
      #determine proportion of DMSO
      pinhib <- mean(y$`CTG Signal`) / mean(df.dmso$`CTG Signal`)
      
      return( c(y$`Cell Type`[1], y$Day[1], y$Molecule[1], y$Dose[1], log2fc, pinhib, pval) )
    })
    
    df <- do.call(rbind,res) %>% as.data.frame()
    
    return(df)
  })
  
  df <- do.call(rbind,resfinal)
  colnames(df) <- c("cell_type", "day", "molecule", "dose", "log2FC", "pinhib", "pval")
  rownames(df) <- 1:dim(df)[1]
  df$dose <- df$dose %>% as.numeric()
  df$log2FC <- df$log2FC %>% as.numeric()
  df$pinhib <- df$pinhib %>% as.numeric()
  df$pval <- df$pval %>% as.numeric()
  
  return(df)
}


data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}



#################################################### 
#Small Molecule Screen Analysis
####################################################

## Set up environment
working_dir <- dirname(rstudioapi::getSourceEditorContext()$path)


## OCI-Ly1 Screen

#import OCI-Ly1 screen data 
fn <- paste0(working_dir, "/OCI-Ly1_compiled_screen_data.txt")
df.ly1 <- fread(fn, data.table = FALSE)

#check the DMSO distribution
df.dmso <- df.ly1 %>% filter(Molecule == "DMSO")

ggplot(df.dmso, aes(x=`CTG Signal`)) + 
  geom_histogram(aes(y=after_stat(density)), bins=50) + 
  geom_vline(aes(xintercept=mean(`CTG Signal`)),
             color="blue", linetype="dashed", linewidth=1) +
  geom_vline(aes(xintercept=300e3),
             color="red", linetype="dashed", linewidth=1) +
  geom_density(alpha=0.2, fill="#FF6666") 


#check outlier values 
table(df.dmso[df.dmso$`CTG Signal` < 300e3, "Plate"])
#...the outliers are equally distributed across all plates and edge wells 

df.dmso.filtered <- df.dmso[df.dmso$`CTG Signal` > 300e3,]

ggplot(df.dmso.filtered, aes(x=`CTG Signal`)) + 
  geom_histogram(aes(y=after_stat(density)), bins=50) + 
  geom_vline(aes(xintercept=mean(`CTG Signal`)),
             color="blue", linetype="dashed", linewidth=1) +
  geom_vline(aes(xintercept=median(`CTG Signal`)),
             color="red", linetype="dashed", linewidth=1) +
  geom_density(alpha=0.2, fill="#FF6666") 


df.ly1.filtered <- df.ly1[ !(df.ly1$Molecule == "DMSO" & df.ly1$`CTG Signal` < 300e3), ]

#determine screen hits 
#df.screen.ly1 <- processScreen(df.ly1.filtered)
df.screen.ly1 <- processScreen(df.ly1)
df.screen.ly1$padj <- p.adjust(df.screen.ly1$pval, method="BH")
df.screen.ly1$log10padj <- -1 * log10(df.screen.ly1$padj)

pvalcutoff <- -1*log10(0.0001)
log2fccutoff <- 1
df.screen.ly1$hit <- "Not a hit"
inds <- df.screen.ly1$log10padj > pvalcutoff & df.screen.ly1$log2FC < -1*abs(log2fccutoff)
df.screen.ly1[inds, "hit"] <- "Hit"
inds <- df.screen.ly1$molecule == "Unused"
df.screen.ly1[inds, "hit"] <- "Control"

cell_type <- "OCI-Ly1"
currdose <- c(1000,500,100)
gglist <- list()

for (i in 1:length(currdose)) {
  gg <- ggplot(df.screen.ly1 %>% filter(dose == currdose[i]), aes(x = log2FC, y = log10padj, color = hit)) + 
    geom_hline(yintercept=pvalcutoff, linetype="dashed", color = "black") +
    geom_vline(xintercept=abs(log2fccutoff), linetype="dashed", color = "black") + 
    geom_vline(xintercept=-1*abs(log2fccutoff), linetype="dashed", color = "black") + 
    geom_point() + 
    scale_color_manual(values=c("blue", "red", "grey")) + 
    labs(x = "log2(Fold Change)", y = "-log10(P-adjust)") + 
    xlim(-10,10) + #ylim(-3,3) + 
    ggtitle(paste0(cell_type, " Screen ", currdose[i], "nM")) + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.x=element_text(colour="black"),
          axis.text.y=element_text(colour="black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.title = element_blank()
    ) +
    geom_label_repel( 
      data=df.screen.ly1 %>% filter(hit == "Hit" & dose == currdose[i]), 
      aes(label=molecule),
      #nudge_x = .1, nudge_y = .1,
      label.size = NA,
      fill = alpha(c("white"),0.1),
      box.padding   = 0.35, 
      point.padding = 0.5,
      max.overlaps = 100,
      segment.color = "black"
      #position = position_dodge(0.8)
    ) 
  
  gglist[[i]] <- gg
}

pdf(paste0(working_dir, "/", cell_type, "_volcano_plots_per_dose_labeled.pdf"), width = 14, height = 4, useDingbats = FALSE)
cowplot::plot_grid(gglist[[1]], gglist[[2]], gglist[[3]], ncol=3)
dev.off()

ggly1 <- gglist 

#save results
write.table(df.screen.ly1, file = paste0(working_dir, "/OCI-Ly1_processed_screen_results.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



## HBL-1 Screen

#import HBL-1 screen data 
fn <- paste0(working_dir, "/HBL-1_compiled_screen_data.txt")
df.hbl1 <- fread(fn, data.table = FALSE)

#check the DMSO distribution
df.dmso <- df.hbl1 %>% filter(Molecule == "DMSO")

ggplot(df.dmso, aes(x=`CTG Signal`)) + 
  geom_histogram(aes(y=after_stat(density)), bins=50) + 
  geom_vline(aes(xintercept=mean(`CTG Signal`)),
             color="blue", linetype="dashed", linewidth=1) +
  geom_vline(aes(xintercept=median(`CTG Signal`)),
             color="red", linetype="dashed", linewidth=1) +
  geom_density(alpha=0.2, fill="#FF6666") 

#determine screen hits 
df.screen.hbl1 <- processScreen(df.hbl1)
df.screen.hbl1$padj <- p.adjust(df.screen.hbl1$pval, method="BH")
df.screen.hbl1$log10padj <- -1 * log10(df.screen.hbl1$padj)

pvalcutoff <- -1*log10(0.0001)
log2fccutoff <- 1
df.screen.hbl1$hit <- "Not a hit"
inds <- df.screen.hbl1$log10padj > pvalcutoff & df.screen.hbl1$log2FC < -1*abs(log2fccutoff)
df.screen.hbl1[inds, "hit"] <- "Hit"
inds <- df.screen.hbl1$molecule == "Unused"
df.screen.hbl1[inds, "hit"] <- "Control"

cell_type <- "HBL-1"
currdose <- c(1000,500,100) 
gglist <- list()

for (i in 1:length(currdose)) {
  gg <- ggplot(df.screen.hbl1 %>% filter(dose == currdose[i]), aes(x = log2FC, y = log10padj, color = hit)) + 
    geom_hline(yintercept=pvalcutoff, linetype="dashed", color = "black") +
    geom_vline(xintercept=abs(log2fccutoff), linetype="dashed", color = "black") + 
    geom_vline(xintercept=-1*abs(log2fccutoff), linetype="dashed", color = "black") + 
    geom_point() + 
    scale_color_manual(values=c("blue", "red", "grey")) + 
    labs(x = "log2(Fold Change)", y = "-log10(P-adjust)") + 
    xlim(-10,10) + #ylim(-3,3) + 
    ggtitle(paste0(cell_type, " Screen ", currdose[i], "nM")) + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.x=element_text(colour="black"),
          axis.text.y=element_text(colour="black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.title = element_blank()
    ) +
    geom_label_repel( 
      data=df.screen.hbl1 %>% filter(hit == "Hit" & dose == currdose[i]), 
      aes(label=molecule),
      #nudge_x = .1, nudge_y = .1,
      label.size = NA,
      fill = alpha(c("white"),0.1),
      box.padding   = 0.35, 
      point.padding = 0.5,
      max.overlaps = 100,
      segment.color = "black"
      #position = position_dodge(0.8)
    ) 
  
  gglist[[i]] <- gg
}

pdf(paste0(working_dir, "/", cell_type, "_volcano_plots_per_dose_labeled.pdf"), width = 14, height = 4, useDingbats = FALSE)
cowplot::plot_grid(gglist[[1]], gglist[[2]], gglist[[3]], ncol=3)
dev.off()

gghbl1 <- gglist

#save results
write.table(df.screen.hbl1, file = paste0(working_dir, "/HBL-1_processed_screen_results.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



## Joint Screen Analysis

#Compare OCI-Ly1 and HBL-1 screens 
colnames(df.screen.hbl1) <- c("cell_type", "day", "molecule", "dose", "HBL1_log2FC", "HBL1_pinhib", "HBL1_pval", "HBL1_padj", "HBL1_log10padj", "hit")
colnames(df.screen.ly1) <- c("cell_type", "day", "molecule", "dose", "Ly1_log2FC", "Ly1_pinhib", "Ly1_pval", "Ly1_padj", "Ly1_log10padj", "hit")
df.screens <- cbind(df.screen.ly1[,2:9], df.screen.hbl1[,5:9])

df.screens$highlight <- "Not a hit"
inds <- df.screen.hbl1$hit == "Hit" & df.screen.ly1$hit == "Hit"
df.screens[inds, "highlight"] <- "Hit: Both"

inds <- df.screen.hbl1$hit == "Hit" & df.screen.ly1$hit == "Not a hit"
df.screens[inds, "highlight"] <- "Hit: HBL-1 only"

inds <- df.screen.hbl1$hit == "Not a hit" & df.screen.ly1$hit == "Hit"
df.screens[inds, "highlight"] <- "Hit: OCI-Ly1 only"

inds <- df.screens$molecule == "Unused"
df.screens[inds, "highlight"] <- "Control"
df.screens$highlight <- factor(df.screens$highlight, levels = c("Not a hit", "Control", "Hit: HBL-1 only", "Hit: OCI-Ly1 only", "Hit: Both"))

currdose <- 1000 #"1 uM"
#currdose <- 500 #"500 nM"
#currdose <- 100 #"100 nM"

pdf(paste0(working_dir, "/D5_Cayman_epigenetics_screen_log2FC_", currdose, "nM.pdf"), width = 5.5, height = 4, useDingbats = FALSE)

ggplot(df.screens %>% filter(dose == currdose), aes(x = Ly1_log2FC, y = HBL1_log2FC, color = highlight)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  geom_vline(xintercept=0, linetype="dashed", color = "black") + 
  geom_point() + 
  scale_color_manual(values=c("grey", "black", "blue", "magenta", "red")) + 
  labs(x = "OCI-Ly1, log2(Fold Change)", y = "HBL-1, log2(Fold Change)") + 
  xlim(-12,2) + ylim(-12,2) + 
  ggtitle(paste0("D5 ", currdose, "nM Cayman Epigenetics Screen")) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank()
  ) +
  geom_label_repel( 
    data=df.screens %>% filter(highlight == "Hit: Both" & dose == currdose), 
    aes(label=molecule),
    label.size = NA,
    fill = alpha(c("white"),0.1),
    box.padding   = 0.35, 
    point.padding = 0.5,
    max.overlaps = 100,
    segment.color = "black"
  ) 

dev.off()

write.table(df.screens, file = paste0(working_dir, "/processed_screen_results.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



#Plot hits at each dose in each screen
inds <- df.screens$dose == 1000 & df.screens$highlight != "Not a hit"
df <- df.screens[inds,c(2,3,5,10,14)]
rownames(df) <- df$molecule
molgroup <- df$highlight
molecules <- df$molecule
df <- df[-c(1,2,5)]
colnames(df) <- c("Ly1_1uM", "HBL1_1uM")

inds <- df.screens$molecule %in% molecules & df.screens$dose == 500
df$Ly1_500nM <- df.screens[inds,5]
df$HBL1_500nM <- df.screens[inds,10]

inds <- df.screens$molecule %in% molecules & df.screens$dose == 100
df$Ly1_100nM <- df.screens[inds,5]
df$HBL1_100nM <- df.screens[inds,10]

df <- df[,c("Ly1_1uM", "Ly1_500nM", "Ly1_100nM", "HBL1_1uM", "HBL1_500nM", "HBL1_100nM")]


column_ha <- HeatmapAnnotation(Dose = factor(rep(c("1000 nM", "500 nM", "100 nM"),2), levels = c("1000 nM", "500 nM", "100 nM")),
                               col = list(Dose = c("1000 nM" = "#b6bb69", "500 nM" = "#7ebb8d", "100 nM" = "#70b69e")),
                               border = TRUE)

row_anno <- rowAnnotation(Class = c("SAH Hydrolase inhibitor",
                                    "Other", 
                                    "BET inhibitor", 
                                    "HDAC inhibitor", 
                                    "EZH2 inhibitor", 
                                    "EZH2 inhibitor", 
                                    "HDAC inhibitor", 
                                    "Other", 
                                    "Other", 
                                    "BET inhibitor", 
                                    "HDAC inhibitor", 
                                    "HDAC inhibitor", 
                                    "H3K9 HMT inhibitor", 
                                    "HDAC inhibitor", 
                                    "BET inhibitor", 
                                    "HDAC inhibitor", 
                                    "DNMT inhibitor", 
                                    "Topoisomerase inhibitor", 
                                    "Other", 
                                    "BET inhibitor", 
                                    "BET inhibitor", 
                                    "HDAC inhibitor", 
                                    "Jumonji Histone Demethylase inhibitor", 
                                    "HDAC inhibitor", 
                                    "Kinase inhibitor", 
                                    "HDAC inhibitor", 
                                    "HDAC inhibitor", 
                                    "BET inhibitor", 
                                    "HDAC inhibitor", 
                                    "HDAC inhibitor", 
                                    "Kinase inhibitor", 
                                    "HDAC inhibitor", 
                                    "HDAC inhibitor", 
                                    "HDAC inhibitor", 
                                    "HDAC inhibitor", 
                                    "DNMT inhibitor", 
                                    "p53 activator", 
                                    "HDAC inhibitor", 
                                    NA), 
                          col = list(Class = setNames(RColorBrewer::brewer.pal(name = "Paired", n = 12), list("BET inhibitor", "DNMT inhibitor", "EZH2 inhibitor", "H3K9 HMT inhibitor", "HDAC inhibitor", "Jumonji Histone Demethylase inhibitor", "Kinase inhibitor", "Other", "p53 activator", "SAH Hydrolase inhibitor", "Topoisomerase inhibitor", "NA"))),
                          border = TRUE)

hm <- Heatmap(df,
        col = pals::warmcool(100), #circlize::colorRamp2(c(0,1), c("white", "blue")),
        cluster_columns = FALSE, 
        cluster_rows = TRUE,
        show_row_dend = FALSE,
        show_column_names = FALSE,
        top_annotation = column_ha,
        right_annotation = row_anno, 
        row_split = molgroup,
        column_split = factor( c(rep("OCI-Ly1",3), rep("HBL-1", 3)) ),
        border = TRUE,
        rect_gp = gpar(col = "black", lwd = 1),
        name = "CTG Signal (Relative to DMSO)", 
        heatmap_legend_param = list(direction = "horizontal", border = TRUE)
        )

draw(hm, merge_legend = TRUE, heatmap_legend_side = "right")

hmrows <- row_order(hm)
new_hmrows <- c(hmrows$Control, rev(hmrows$`Hit: Both`), rev(hmrows$`Hit: HBL-1 only`), rev(hmrows$`Hit: OCI-Ly1 only`))
molgroup <- factor(molgroup, levels = c("Control", "Hit: Both", "Hit: HBL-1 only", "Hit: OCI-Ly1 only"))


hm2 <- Heatmap(df,
        col = circlize::colorRamp2(c(0,0.5,1), c("#B40426", "#DEDCDB", "#3B4CC0")),
        cluster_columns = FALSE, 
        cluster_rows = FALSE,
        show_row_dend = FALSE,
        show_column_names = FALSE,
        top_annotation = column_ha,
        right_annotation = row_anno,
        row_order = new_hmrows,
        row_split = molgroup,
        column_split = factor( c(rep("OCI-Ly1",3), rep("HBL-1", 3)) ),
        border = TRUE,
        rect_gp = gpar(col = "black", lwd = 1),
        name = "CTG Signal (Relative to DMSO)", 
        heatmap_legend_param = list(direction = "horizontal", border = TRUE)
)

pdf(paste0(working_dir, "/D5_Cayman_epigenetics_screen_hits_dose_response_drug_annotation.pdf"), width = 7, height = 7, useDingbats = FALSE)
draw(hm2, merge_legend = TRUE, heatmap_legend_side = "right")
dev.off()



#Plot dose response of positive controls
pc.sm <- c("GSK343", "UNC1999", "3-Deazaneplanocin A", "I-BET762", "I-BET151", "(+)-JQ1", "Panobinostat", "Chaetocin", "JIB-04")

df.pc <- df.ly1 %>% filter(Molecule %in% pc.sm)
df.dmso <- df.ly1 %>% filter(Molecule == "DMSO")
inds <- sample(1:dim(df.dmso)[1], 12, replace=FALSE)
df.pc$NormCTG <- df.pc$`CTG Signal` / median(df.dmso$`CTG Signal`) #mean(df.dmso[inds, "CTG Signal"])

df2 <- data_summary(df.pc, varname="NormCTG", 
                    groupnames=c("Dose", "Molecule"))
df2$Molecule <- factor(df2$Molecule, levels = pc.sm)

ggly1 <- ggplot(df2, aes(x=Dose, y=NormCTG)) +
  geom_line()+
  geom_pointrange(aes(ymin=NormCTG-sd, ymax=NormCTG+sd)) +
  scale_x_continuous(limits = c(0,1000), breaks = c(0,250,500,750,1000)) + 
  scale_y_continuous(limits = c(0,1.4), breaks = c(0,0.4,0.8,1.2)) + 
  xlab("") + 
  ylab("Relative CTG Signal") + 
  facet_wrap(~Molecule, nrow = 1) + 
  theme_linedraw(base_size = 12) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        strip.background = element_blank(),
        strip.text = element_text(colour = "black"),
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 12), 
        axis.text.x=element_blank(),
        axis.title.x=element_blank())


df.pc <- df.hbl1 %>% filter(Molecule %in% pc.sm)
df.dmso <- df.hbl1 %>% filter(Molecule == "DMSO")
inds <- sample(1:dim(df.dmso)[1], 12, replace=FALSE)
df.pc$NormCTG <- df.pc$`CTG Signal` / median(df.dmso$`CTG Signal`) #mean(df.dmso[inds, "CTG Signal"])

df2 <- data_summary(df.pc, varname="NormCTG", 
                    groupnames=c("Dose", "Molecule"))
df2$Molecule <- factor(df2$Molecule, levels = pc.sm)

gghbl1 <- ggplot(df2, aes(x=Dose, y=NormCTG)) +
  geom_line()+
  geom_pointrange(aes(ymin=NormCTG-sd, ymax=NormCTG+sd)) +
  scale_x_continuous(limits = c(0,1000), breaks = c(0,250,500,750,1000)) + 
  scale_y_continuous(limits = c(0,1.4), breaks = c(0,0.4,0.8,1.2)) + 
  xlab("Dose (nM)") + 
  ylab("Relative CTG Signal") + 
  facet_wrap(~Molecule, nrow = 1) + 
  theme_linedraw(base_size = 12) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 12), 
        axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1))

pdf(paste0(working_dir, "/D5_Cayman_epigenetics_screen_hits_dose_response_select_molecules.pdf"), width = 11, height = 3.5, useDingbats = FALSE)
cowplot::plot_grid(ggly1, gghbl1, nrow = 2)
dev.off()



