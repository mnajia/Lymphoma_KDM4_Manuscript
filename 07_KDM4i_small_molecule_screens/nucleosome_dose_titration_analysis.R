#Project: Lymphoma KDM4 Manuscript
#Objective: Analyze dose titration data on nucleosome substrates of the top hits from histone peptide and nucleosome secondary screens
#Mohamad Najia

library(data.table)
library(dplyr)
library(ggplot2)
library(ggExtra)
library(scales) 


#################################################### 
#Analysis of small molecule dose titration data on nucleosomes
####################################################

## Set up environment
project_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
fn <- paste0(project_dir, "/nucleosome_dose_titration_prioritized_hits_10122022.txt")
df <- fread(fn, data.table = FALSE)
df <- df %>% filter(! Molecule_Name %in% c("DMSO", "CPI-455"))

## Plot all small molecules
data_summary <- df %>%
  group_by(Concentration, Molecule_Name) %>%
  summarise(
    MeanValue = mean(Signal_Relative_to_DMSO),
    SD = sd(Signal_Relative_to_DMSO)
  )

pdf(file = paste0(project_dir, "/nucleosome_dose_titration_prioritized_hits_10122022_all_molecules.pdf"), width = 8, height = 8, useDingbats = FALSE)

ggplot() +
  geom_line(data = data_summary, aes(x = Concentration, y = MeanValue, group = Molecule_Name)) + 
  geom_errorbar(data = data_summary, aes(x = Concentration, y = MeanValue, ymin = MeanValue-SD, ymax = MeanValue+SD, width = 0.2)) + 
  geom_point(data = df, aes(x = Concentration, y = Signal_Relative_to_DMSO, group = Molecule_Name),
             position = position_dodge(0.9), color = "blue", size = 3, alpha = 0.25) +
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x)) + 
  labs(x = "Small Molecule Dose (uM)", y = "Signal Relative to DMSO") + 
  ggtitle("Nucleosome H3K9me3 Demethylation Activity") + 
  facet_wrap(~Molecule_Name, ncol = 5, nrow = 4) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black"),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.ticks = element_line(color = "black"),
        panel.grid.minor = element_blank())

dev.off()


## Plot promising small molecules with IC50s less than or equal to 10 uM
#df2 <- df %>% filter(Source_Well %in% c("G05", "J03", "J05", "M03"))
df2 <- df %>% filter(Molecule_Name %in% c("CB67730375", "CB40551911", "CB43433036", "CB32128362"))
df2$Molecule_Name <- factor(df2$Molecule_Name, levels = c("CB67730375", "CB40551911", "CB43433036", "CB32128362"))

data_summary <- df2 %>%
  group_by(Concentration, Molecule_Name) %>%
  summarise(
    MeanValue = mean(Signal_Relative_to_DMSO),
    SD = sd(Signal_Relative_to_DMSO)
  )

pdf(file = paste0(project_dir, "/nucleosome_dose_titration_prioritized_hits_10122022_promising_molecules.pdf"), width = 8.5, height = 3, useDingbats = FALSE)

ggplot() +
  geom_line(data = data_summary, aes(x = Concentration, y = MeanValue, group = Molecule_Name)) + 
  geom_errorbar(data = data_summary, aes(x = Concentration, y = MeanValue, ymin = MeanValue-SD, ymax = MeanValue+SD, width = 0.2)) + 
  geom_point(data = df2, aes(x = Concentration, y = Signal_Relative_to_DMSO, group = Molecule_Name),
             position = position_dodge(0.9), color = "blue", size = 3, alpha = 0.25) +
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x)) + 
  scale_y_continuous(limits = c(0,1.25), breaks = c(0,0.25,0.5,0.75,1,1.25)) + 
  xlab("Small Molecule Dose (uM)") + 
  ylab("Signal Relative to DMSO") + 
  ggtitle("Nucleosome H3K9me3 Demethylation Activity") + 
  facet_wrap(~Molecule_Name, ncol = 4, nrow = 1) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black"),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.ticks = element_line(color = "black"),
        panel.grid.minor = element_blank())

dev.off()


## Plot the remaining small molecules for the supplement
df2 <- df[!(df$Molecule_Name %in% c("CB67730375", "CB40551911", "CB43433036", "CB32128362")),]

data_summary <- df2 %>%
  group_by(Concentration, Molecule_Name) %>%
  summarise(
    MeanValue = mean(Signal_Relative_to_DMSO),
    SD = sd(Signal_Relative_to_DMSO)
  )

pdf(file = paste0(project_dir, "/nucleosome_dose_titration_prioritized_hits_10122022_remaining_molecules.pdf"), width = 10, height = 3.5, useDingbats = FALSE)

ggplot() +
  geom_line(data = data_summary, aes(x = Concentration, y = MeanValue, group = Molecule_Name)) +
  geom_errorbar(data = data_summary, aes(x = Concentration, y = MeanValue, ymin = MeanValue-SD, ymax = MeanValue+SD, width = 0.2)) + 
  geom_point(data = df2, aes(x = Concentration, y = Signal_Relative_to_DMSO, group = Molecule_Name),
             position = position_dodge(0.9), color = "blue", size = 3, alpha = 0.25) +
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x)) + 
  #scale_y_continuous(limits = c(0,1.25), breaks = c(0,0.25,0.5,0.75,1,1.25)) + 
  xlab("Small Molecule Dose (uM)") + 
  ylab("Signal Relative to DMSO") + 
  ggtitle("Nucleosome H3K9me3 Demethylation Activity") + 
  facet_wrap(~Molecule_Name, ncol = 8, nrow = 2) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black"),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.ticks = element_line(color = "black"),
        panel.grid.minor = element_blank())

dev.off()



#################################################### 
#Analysis of dose titration of ChemBridge biosimilar molecules on nucleosomes
####################################################

## Import data
fn <- paste0(project_dir, "/nucleosome_dose_titration_chembridge_biosimilars_11092023.txt")
df <- fread(fn, data.table = FALSE)
df <- df %>% filter(! Molecule_Name %in% c("DMSO", "CPI-455"))
df <- df[df$Exclude != "Y",]

## Plot all small molecules
data_summary <- df %>%
  group_by(Concentration, Molecule_Name) %>%
  summarise(
    MeanValue = mean(Signal_Relative_to_DMSO),
    SD = sd(Signal_Relative_to_DMSO)
  )

pdf(file = paste0(project_dir, "/nucleosome_dose_titration_chembridge_biosimilars_11092023_all_molecules.pdf"), width = 8, height = 12, useDingbats = FALSE)

ggplot() +
  geom_line(data = data_summary, aes(x = Concentration, y = MeanValue, group = Molecule_Name)) + 
  geom_errorbar(data = data_summary, aes(x = Concentration, y = MeanValue, ymin = MeanValue-SD, ymax = MeanValue+SD, width = 0.2)) + 
  geom_point(data = df, aes(x = Concentration, y = Signal_Relative_to_DMSO, group = Molecule_Name),
             position = position_dodge(0.9), color = "blue", size = 3, alpha = 0.25) +
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x)) + 
  #labels = trans_format("log10", math_format(10^.x))) +
  labs(x = "Small Molecule Dose (uM)", y = "Signal Relative to DMSO") + 
  ggtitle("Nucleosome H3K9me3 Demethylation Activity") + 
  facet_wrap(~Molecule_Name, ncol = 7, nrow = 10) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black"),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.ticks = element_line(color = "black"),
        panel.grid.minor = element_blank())

dev.off()

## Plot promising small molecules
#df2 <- df %>% filter(Source_Well %in% c("C07", "C11", "G05",
#                                        "G09", "H07", "I11",
#                                        "J05"))
df2 <- df %>% filter(Molecule_Name %in% c("CB50272401", "CB91488421", "CB70948541",
                                        "CB90889252", "CB69140372", "CB66138813",
                                        "CB61634775"))

data_summary <- df2 %>%
  group_by(Concentration, Molecule_Name) %>%
  summarise(
    MeanValue = mean(Signal_Relative_to_DMSO),
    SD = sd(Signal_Relative_to_DMSO)
  )

pdf(file = paste0(project_dir, "/nucleosome_dose_titration_chembridge_biosimilars_11092023_promising_molecules.pdf"), width = 10, height = 2, useDingbats = FALSE)

ggplot() +
  geom_line(data = data_summary, aes(x = Concentration, y = MeanValue, group = Molecule_Name)) + 
  geom_errorbar(data = data_summary, aes(x = Concentration, y = MeanValue, ymin = MeanValue-SD, ymax = MeanValue+SD, width = 0.2)) + 
  geom_point(data = df2, aes(x = Concentration, y = Signal_Relative_to_DMSO, group = Molecule_Name),
             position = position_dodge(0.9), color = "blue", size = 3, alpha = 0.25) +
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x)) + 
  scale_y_continuous(limits = c(0,1.25), breaks = c(0,0.25,0.5,0.75,1,1.25)) + 
  xlab("Small Molecule Dose (uM)") + 
  ylab("Signal Relative to DMSO") + 
  ggtitle("Nucleosome H3K9me3 Demethylation Activity") + 
  facet_wrap(~Molecule_Name, ncol = 8, nrow = 1) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black"),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.ticks = element_line(color = "black"),
        panel.grid.minor = element_blank())

dev.off()



#################################################### 
#Analysis of dose titration of Enamine biosimilar molecules on nucleosomes
####################################################

## Import data
fn <- paste0(project_dir, "/nucleosome_dose_titration_enamine_biosimilars_11082023.txt")
df <- fread(fn, data.table = FALSE)
df <- df %>% filter(! Molecule_Name %in% c("DMSO", "CPI-455"))

## Plot all small molecules
data_summary <- df %>%
  group_by(Concentration, Molecule_Name) %>%
  summarise(
    MeanValue = mean(Signal_Relative_to_DMSO),
    SD = sd(Signal_Relative_to_DMSO)
  )

pdf(file = paste0(project_dir, "/nucleosome_dose_titration_enamine_biosimilars_11092023_all_molecules.pdf"), width = 8, height = 12, useDingbats = FALSE)

ggplot() +
  geom_line(data = data_summary, aes(x = Concentration, y = MeanValue, group = Molecule_Name)) + 
  geom_errorbar(data = data_summary, aes(x = Concentration, y = MeanValue, ymin = MeanValue-SD, ymax = MeanValue+SD, width = 0.2)) + 
  geom_point(data = df, aes(x = Concentration, y = Signal_Relative_to_DMSO, group = Molecule_Name),
             position = position_dodge(0.9), color = "blue", size = 3, alpha = 0.25) +
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x)) + 
  #labels = trans_format("log10", math_format(10^.x))) +
  labs(x = "Small Molecule Dose (uM)", y = "Signal Relative to DMSO") + 
  ggtitle("Nucleosome H3K9me3 Demethylation Activity") + 
  facet_wrap(~Molecule_Name, ncol = 7, nrow = 10) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black"),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.ticks = element_line(color = "black"),
        panel.grid.minor = element_blank())

dev.off()

## Plot promising small molecules
#df2 <- df %>% filter(Source_Well %in% c("C09", "F09", "G05", "M07", "N07")) 
df2 <- df %>% filter(Molecule_Name %in% c("Z1256159460", "Z363389670", "Z2372639318", "Z1143402100", "Z969007916")) 

data_summary <- df2 %>%
  group_by(Concentration, Molecule_Name) %>%
  summarise(
    MeanValue = mean(Signal_Relative_to_DMSO),
    SD = sd(Signal_Relative_to_DMSO)
  )

pdf(file = paste0(project_dir, "/nucleosome_dose_titration_enamine_biosimilars_11092023_promising_molecules.pdf"), width = 6, height = 2, useDingbats = FALSE)

ggplot() +
  geom_line(data = data_summary, aes(x = Concentration, y = MeanValue, group = Molecule_Name)) + 
  geom_errorbar(data = data_summary, aes(x = Concentration, y = MeanValue, ymin = MeanValue-SD, ymax = MeanValue+SD, width = 0.2)) + 
  geom_point(data = df2, aes(x = Concentration, y = Signal_Relative_to_DMSO, group = Molecule_Name),
             position = position_dodge(0.9), color = "blue", size = 3, alpha = 0.25) +
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x)) + 
  scale_y_continuous(limits = c(0,1.25), breaks = c(0,0.25,0.5,0.75,1,1.25)) + 
  xlab("Small Molecule Dose (uM)") + 
  ylab("Signal Relative to DMSO") + 
  ggtitle("Nucleosome H3K9me3 Demethylation Activity") + 
  facet_wrap(~Molecule_Name, ncol = 5, nrow = 2) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black"),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.ticks = element_line(color = "black"),
        panel.grid.minor = element_blank())

dev.off()







