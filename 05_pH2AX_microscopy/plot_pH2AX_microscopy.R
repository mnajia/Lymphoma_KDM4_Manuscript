#Project: Lymphoma KDM4 Manuscript
#Objective: analyze microscopy of pH2A.X from OCI-Ly1 and HBL-1 cells treated with 100 nM QC6352 for 48 hours
#Mohamad Najia

library(data.table)
library(dplyr)
library(ggplot2)
library(ggExtra)
library(scales) 


# Import data
project_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

fn <- paste0(project_dir, "/QC6352_100nM_48hrs_pH2AX_microscopy.csv")
df <- fread(fn, data.table = FALSE)

df <- filter(df, WELL %in% c("C09", "C10", "C11", "C12"))
df$cell_type <- "OCI-Ly1"
df[df$WELL %in% c("C11", "C12"), "cell_type"] <- "HBL-1"
df$condition <- ""
df[df$WELL %in% c("C09", "C11"), "condition"] <- "DMSO"
df[df$WELL %in% c("C10", "C12"), "condition"] <- "QC6352 (100 nM)"


# Perform filtering to remove putative debris

#check distribution of DAPI intensity 
ggplot(df, aes(x = DAPI_IntDen, fill = WELL)) + 
  geom_density(alpha = 0.7) + 
  geom_vline(xintercept = 2e5) + 
  xlab("") + 
  ylab("") + 
  ggtitle("") + 
  facet_grid(~cell_type) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))

#filter based on DAPI signal
df.hbl <- filter(df, cell_type == "HBL-1" & DAPI_IntDen > 2e5)
df.ly1 <- filter(df, cell_type == "OCI-Ly1" & DAPI_IntDen > 5e5)
df <- rbind(df.hbl, df.ly1)


# Plot pH2A.X distributions
gg <- ggplot(df, aes(x = H2AX_speckles_IntDen, fill = condition)) + # H2AX_speckles_RawIntDen
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  xlab("pH2A.X Speckles Intensity per Cell") + 
  ylab("") + 
  ggtitle("") + 
  scale_x_continuous(limits = c(0,40)) + 
  facet_grid(~cell_type) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))

pdf(paste0(project_dir, "/D2_QC6352_100nM_pH2AX_speckles_intensity_density.pdf"), width = 7, height = 3, useDingbats = FALSE)
gg
dev.off()


ggplot(df, aes(x = H2AX_speckles_Mean, fill = condition)) + 
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  xlab("pH2A.X mean number of speckles per cell") + 
  ylab("") + 
  ggtitle("") + 
  facet_grid(~cell_type) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))



