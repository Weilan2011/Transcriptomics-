#install.packages("patchwork")
#install.packages("pheatmap")
#install.packages("gridExtra")
#install.packages("tidyverse")
#install.packages("egg")
#install.packages("git2r")
#install.packages("ggpubr")
#install.packages("gridExtra")
#install.packages("gridExtra")
#install.packages("gridExtra")
#install.packages("gridExtra")
#install.packages("gridExtra")
#install.packages("gridExtra")
#install.packages("gridExtra")

library(pheatmap)
library(gridExtra)
library(tidyverse)
library(egg)
library(grid)
library(git2r)
library(ggpubr)
library(DESeq2)
library(ggplot2)
library(readxl)
library(rstatix)
library(dplyr)
library(Hmisc)
library(GGally)
library(ggcorrplot)
library(cowplot)
library(RColorBrewer)
library(tidyverse)
library(tidyr)
library(lattice)
library(nlme)
library(lsmeans)
library(vegan)
library(mgcv)
library(EnhancedVolcano)
library(patchwork)

#### setup themes for figures ###
mytheme1 <- theme(panel.background = element_rect(fill = "white", color="black"), 
                  axis.title = element_text(face="bold",size=12, color = "black"),
                  axis.text.y  = element_text(face="bold",size=12, color = "black"),
                  axis.text.x  = element_text(face="bold",size=12, color = "black"),
                  panel.grid.minor.y=element_blank(), 
                  panel.grid.minor.x=element_blank(),
                  panel.grid.major.y=element_blank(),
                  panel.grid.major.x=element_blank(),
                  strip.background = element_blank(),
                  strip.text.x = element_text(face="bold",size=14, color = "black"),
                  legend.position = "right",
                  legend.box = "horizontal",
                  legend.title= element_blank(),
                  legend.text =element_text(size=14),
                  legend.key = element_blank())

mytheme2 <- theme(panel.background = element_rect(fill = "white", color="black"), 
                  axis.title = element_text(face="bold",size=16, color = "black"),
                  axis.text.y  = element_text(face="bold",size=12, color = "black"),
                  axis.text.x  = element_text(face="bold",size=12, color = "black"),
                  panel.grid.minor.y=element_blank(), 
                  panel.grid.minor.x=element_blank(),
                  panel.grid.major.y=element_blank(),
                  panel.grid.major.x=element_blank(),
                  strip.background = element_blank(),
                  strip.text.x = element_text(face="bold",size=10, color = "black"),
                  legend.position = "right",
                  legend.box = "vertical",
                  legend.title= element_text(face="bold",size=10, color = "black"),
                  legend.text =element_text(size=10),
                  legend.key = element_blank())

setwd("/home/weilan/ENCODE/Subset_regtools_out")
subset_combined_data <- read.csv("STAR_HISAT2_combined_data2.csv", header = TRUE, sep = ",")
subset_combined_data2 <- read.csv("TEST_combined_data.csv", header = TRUE, sep = ",")

Junctions_per_reads_data <- subset_combined_data[-1,]

str(Junctions_per_reads_data)
str(subset_combined_data2)

# Convert factors
Junctions_per_reads_data$Junction.Type  <- factor(Junctions_per_reads_data$Junction.Type , 
                                                  levels = c('DA','N','D','A','NDA'))
subset_combined_data2$Junction_Type  <- factor(subset_combined_data2$Junction_Type , 
                                                  levels = c('DA','N','D','A','NDA'))

Junctions_per_reads_data$Known_junction <- factor(Junctions_per_reads_data$Known_junction, 
                                                  levels = c('0','1'))
subset_combined_data2$Known_junction <- factor(subset_combined_data2$Known_junction, 
                                                  levels = c('0','1'))

str(subset_combined_data2)


scatter_plot1 <- ggplot(Junctions_per_reads_data, 
                        aes(x=K562_HDGF_rep1_hisat2, 
                            y=K562_HDGF_rep1_STAR, color=Junction.Type)) + 
  geom_point(size=2) + 
  scale_x_log10() +
  scale_y_log10() +
  labs(title= "K562_HDGF", y="STAR", x="HISAT2")+ 
  mytheme2

scatter_plot2 <- ggplot(Junctions_per_reads_data, 
                        aes(x=K562_NT_rep1_hisat2, 
                            y=K562_NT_rep1_STAR, color=Junction.Type)) +
  geom_point(size=2) + 
  scale_x_log10() +
  scale_y_log10() +
  labs(title= "K562_NT", y="STAR", x="HISAT2")+ 
  mytheme2

scatter_plot3 <- ggplot(Junctions_per_reads_data, 
                        aes(x=HepG2_NUP35_rep1_hisat2, 
                            y=HepG2_NUP35_rep1_STAR, color=Junction.Type)) +
  geom_point(size=2) + 
  scale_x_log10() +
  scale_y_log10() +
  labs(title= "HepG2_NUP35", y="STAR", x="HISAT2")+ mytheme2

scatter_plot4 <- ggplot(Junctions_per_reads_data, aes(x=HepG2_NT_rep1_hisat2, 
                                                      y=HepG2_NT_rep1_STAR, color=Junction.Type)) +
  geom_point(size=2) + 
  scale_x_log10() +
  scale_y_log10() +
  labs(title= "HepG2_NT", y="STAR", x="HISAT2")+ mytheme2

scatter_plot5 <- ggplot(subset_combined_data2, aes(x=RBM_8162_hisat2, 
                                                      y=RBM_8162_STAR, 
                                                      color=Junction_Type)) +
  geom_point(size=2) + 
  scale_x_log10() +
  scale_y_log10() +
  labs(title= "RBM_8162", y="STAR", x="HISAT2")+ mytheme2
scatter_plot5


# Combine the scatter plots into one
combined_plot <- scatter_plot1 + scatter_plot2 + scatter_plot3 + scatter_plot4 
combined_plot 
