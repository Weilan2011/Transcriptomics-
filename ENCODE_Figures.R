# install.packages("UpSetR")

library(UpSetR)
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
                  axis.title = element_text(size=14, color = "black"),
                  axis.text.y  = element_text(size=12, color = "black"),
                  axis.text.x  = element_text(size=12, color = "black"),
                  panel.grid.major.y=element_blank(), 
                  panel.grid.major.x=element_blank(),
                  panel.grid.minor.y=element_line(color = "lightgrey"),
                  panel.grid.minor.x=element_line(color = "lightgrey"),
                  axis.line = element_line(color = "black"),
                  strip.background = element_blank(),
                  strip.text.x = element_text(face="bold",size=10, color = "black"),
                  legend.position = "right",
                  legend.box = "vertical",
                  legend.title= element_text(face="bold",size=10, color = "black"),
                  legend.text =element_text(size=10),
                  legend.key = element_blank())

setwd("/home/weilan/ENCODE/4_ALS/")

ALS_rMATS_Summary_data <- read.csv("ENCODE_ALS_Summary.csv", header = TRUE, sep = ",")
str(ALS_rMATS_Summary_data)

# Convert factors
ALS_rMATS_Summary_data$As_Type  <- factor(ALS_rMATS_Summary_data$As_Type, 
                                             levels = c( 'A3SS','A5SS', 'MXE', 'RI','SE'))

ALS_rMATS_Summary_data$Target  <- factor(ALS_rMATS_Summary_data$Target, 
                                              levels = c('EWSR1', 'FUS', 'HNRNPA1','HNRNPA2B1', 
                                                         'MATR3', 'TAF15', 'TARDBP','TIA1'))

ALS_rMATS_Summary_data$Cell_line<- factor(ALS_rMATS_Summary_data$Cell_line, 
                                              levels = c('HepG2','K562'))

str(ALS_rMATS_Summary_data)

AS_counts_Significant.Events <-  ALS_rMATS_Summary_data %>% 
  ggplot()+
  geom_col(aes(x=As_Type, y=Significant.Events, fill=Cell_line), position = "dodge")+
  mytheme2+
  scale_fill_manual(values = c("#285D92", "orange"))  +
  facet_wrap(~Target, nrow=4)+
  coord_flip()+
  labs(y="# of Significant events 
       (FDR<0.05, |IncLevelDifference|>0.2)", x="Alternative Splicing Events") 
AS_counts_Significant.Events


ALS_rMATS_ercent_data <- read.csv("ENCODE_ALS_Summary_Percent.csv", header = TRUE, sep = ",")
str(ALS_rMATS_ercent_data)

# Convert factors
ALS_rMATS_ercent_data$As_Type  <- factor(ALS_rMATS_ercent_data$As_Type, 
                                          levels = c('SE','RI', 'MXE','A5SS','A3SS'))

ALS_rMATS_ercent_data$Cell_Target  <- factor(ALS_rMATS_ercent_data$Cell_Target, 
                                              levels = c('HepG2_EWSR1', 'K562_EWSR1',
                                                         'HepG2_FUS', 'K562_FUS',
                                                         'HepG2_HNRNPA1','K562_HNRNPA1',
                                                         'HepG2_HNRNPA2B1','K562_HNRNPA2B1', 
                                                         'HepG2_MATR3', 'K562_MATR3',
                                                         'HepG2_TAF15', 'K562_TAF15',
                                                         'HepG2_TARDBP','K562_TARDBP',
                                                         'HepG2_TIA1','K562_TIA1'))

ALS_rMATS_ercent_data$AS_Category  <- factor(ALS_rMATS_ercent_data$AS_Category, 
                                           levels = c('IJC Percent','SJC Percent'))

str(ALS_rMATS_ercent_data)

AS_counts_event_percent <-  ggplot(data = ALS_rMATS_ercent_data)+
  geom_bar(aes(x= As_Type, y=Group_2, fill=AS_Category), stat = 'identity')+
  mytheme2+
  scale_fill_manual(values = c("dimgray","lightslateblue"))  +
  facet_wrap(~Cell_Target, nrow=4)+
  coord_flip()+
  labs(y="% of inclusion & exlcusion events ", x="Alternative Splicing Events") 
AS_counts_event_percent


######### HepG2_K562 Alternative Splicing events overlaps ########
# install.packages("patchwork")
library(ggplot2)
library(dplyr)
library(patchwork)
library(purrr)
library(ggrepel)

# Define the file paths
files <- c(
  "ENCODE_ALS_SE_Summary.csv",
  "ENCODE_ALS_RI_Summary.csv",
  "ENCODE_ALS_MXE_Summary.csv",
  "ENCODE_ALS_A5SS_Summary.csv",
  "ENCODE_ALS_A3SS_Summary.csv"
)

# Define the targets
targets <- c("EWSR1","FUS","HNRNPA1","HNRNPA2B1", "MATR3", "TAF15", "TARDBP", "TIA1")

# Process each file
for (file in files) {
  plot_list <- list()  # Initialize an empty list to store the plots
  df <- read.csv(file)  # Load the dataset
  
  # Extract the type of splicing event from the file name for the plot title
  splicing_event_type <- gsub(".*_ALS_(.*)_Summary\\.csv", "\\1", basename(file))
  
  for (target in targets) {
    # Check if the columns exist and convert to numeric if they are not
    hepG2_col <- paste0("HepG2_", target)
    k562_col <- paste0("K562_", target)
    
    if (!(hepG2_col %in% names(df)) || !(k562_col %in% names(df))) {
      message(paste("Skipping: Columns for", target, "not found in", file))
      next  # Skip this iteration if columns do not exist
    }
    
    df[[hepG2_col]] <- as.numeric(as.character(df[[hepG2_col]]))
    df[[k562_col]] <- as.numeric(as.character(df[[k562_col]]))
    
    # Add columns to indicate significance
    df$HepG2_significant <- abs(df[[hepG2_col]]) > 0.2
    df$K562_significant <- abs(df[[k562_col]]) > 0.2
    
    # Categorize events based on significance
    df$Category <- case_when(
      df$HepG2_significant & df$K562_significant ~ "Both",
      df$HepG2_significant & !df$K562_significant ~ "HepG2 Only",
      !df$HepG2_significant & df$K562_significant ~ "K562 Only"
    )
    
    # Count the number of events for each category
    counts <- df %>% group_by(Category) %>% summarise(n = n(), .groups = 'drop')
    
    # Generate scatter plot
    p <- ggplot(df, aes(x = .data[[hepG2_col]], y = .data[[k562_col]])) +
      geom_point(aes(color = Category), alpha = 0.5) +
      scale_color_manual(values = c("Both" = "red", "HepG2 Only" = "#285D92", "K562 Only" = "orange", "neither" = "grey")) +
      labs(title = paste(splicing_event_type, target),
           x = paste("HepG2", "IncLevelDifference"),
           y = paste("K562", "IncLevelDifference")) +
      geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
      geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
      mytheme2 # Use a simple theme
    
    # Add annotations for the counts in the upper left corner
    annotations_y <- max(df[[k562_col]], na.rm = TRUE) * c(1.1, 0.98, 0.86)  # Adjust these based on your data range
    for(i in 1:nrow(counts)) {
      p <- p + annotate("text", x = min(df[[hepG2_col]], na.rm = TRUE), y = annotations_y[i], 
                        label = paste(counts$Category[i], counts$n[i]), hjust = 0, vjust = 1, size = 3)
    }
    
    plot_list[[target]] <- p  # Add the plot to the list
  }
  
  # Combine the plots using patchwork OUTSIDE the targets loop but INSIDE the files loop
  combined_plot <- purrr::reduce(plot_list, `+`) + plot_layout(ncol = 2) # Arrange in 2 columns
  
  # Print the combined plot for the current AS type
  print(combined_plot) # Display the plot
  
  # Optionally, pause in non-interactive environments
  if (!interactive()) Sys.sleep(time = 2) # Pause for 2 seconds; adjust time as needed
  
  # Save the combined plot to a file
  plot_filename <- paste0("ALS_rMATS_", splicing_event_type, "_combined_Significance_Plot.png")
  ggsave(plot_filename, combined_plot, width = 20, height = 20)  # Adjust size as needed
}

###### Heatmap ######
ECODE_ALS_AS_GENE_heatmap <- read.csv("ENCODE_ALS_RBP_Gene_AS_heatmap.csv", 
                                             header = TRUE, row.names = 1) %>% filter(Total_SE<8)

ECODE_ALS_AS_GENE_heatmap_Binary <- read.csv("ENCODE_ALS_RBP_Gene_AS_heatmap_Binary.csv", 
                                      header = TRUE, row.names = 1) %>% filter(Total_SE>8)
str(ECODE_ALS_AS_GENE_heatmap)
str(ECODE_ALS_AS_GENE_heatmap_Binary)

ECODE_ALS_SE_HepG2 <- ECODE_ALS_AS_GENE_heatmap[16:23]
ECODE_ALS_SE_K562 <- ECODE_ALS_AS_GENE_heatmap[24:31]
ECODE_ALS_SE <- ECODE_ALS_AS_GENE_heatmap[16:31]
ECODE_ALS_SE_Binary <- ECODE_ALS_AS_GENE_heatmap_Binary[16:31]

ECODE_ALS_RI_HepG2 <- ECODE_ALS_AS_GENE_heatmap[32:39]
ECODE_ALS_RI_K562 <- ECODE_ALS_AS_GENE_heatmap[40:47]
ECODE_ALS_RI <- ECODE_ALS_AS_GENE_heatmap[32:47]
ECODE_ALS_RI_Binary <- ECODE_ALS_AS_GENE_heatmap_Binary[32:47]

ECODE_ALS_MXE_HepG2 <- ECODE_ALS_AS_GENE_heatmap[48:55]
ECODE_ALS_MXE_K562 <- ECODE_ALS_AS_GENE_heatmap[56:63]
ECODE_ALS_MXE <- ECODE_ALS_AS_GENE_heatmap[48:63]
ECODE_ALS_MXE <- ECODE_ALS_AS_GENE_heatmap[48:63]

ECODE_ALS_A3SS_HepG2 <- ECODE_ALS_AS_GENE_heatmap[64:71]
ECODE_ALS_A3SS_K562 <- ECODE_ALS_AS_GENE_heatmap[72:79]
ECODE_ALS_A3SS <- ECODE_ALS_AS_GENE_heatmap[64:79]

ECODE_ALS_A5SS_HepG2 <- ECODE_ALS_AS_GENE_heatmap[80:87]
ECODE_ALS_A5SS_K562 <- ECODE_ALS_AS_GENE_heatmap[88:95]
ECODE_ALS_A5SS <- ECODE_ALS_AS_GENE_heatmap[80:95]

str(ECODE_ALS_SE_HepG2)
str(ECODE_ALS_RI_HepG2)
str(ECODE_ALS_MXE_HepG2)
str(ECODE_ALS_A3SS_HepG2)
str(ECODE_ALS_A5SS_HepG2)

str(ECODE_ALS_SE_K562)
str(ECODE_ALS_RI_K562)
str(ECODE_ALS_MXE_K562)
str(ECODE_ALS_A3SS_K562)
str(ECODE_ALS_A5SS_K562)

ECODE_ALS_SE_HepG2 <- as.matrix(ECODE_ALS_SE_HepG2)
ECODE_ALS_RI_HepG2 <- as.matrix(ECODE_ALS_RI_HepG2)
ECODE_ALS_MXE_HepG2 <- as.matrix(ECODE_ALS_MXE_HepG2)
ECODE_ALS_A3SS_HepG2 <- as.matrix(ECODE_ALS_A3SS_HepG2)
ECODE_ALS_A5SS_HepG2 <- as.matrix(ECODE_ALS_A5SS_HepG2)

ECODE_ALS_SE_K562 <- as.matrix(ECODE_ALS_SE_K562)
ECODE_ALS_RI_K562 <- as.matrix(ECODE_ALS_RI_K562)
ECODE_ALS_MXE_K562 <- as.matrix(ECODE_ALS_MXE_K562)
ECODE_ALS_A3SS_K562 <- as.matrix(ECODE_ALS_A3SS_K562)
ECODE_ALS_A5SS_K562 <- as.matrix(ECODE_ALS_A5SS_K562)


pheatmap::pheatmap(ECODE_ALS_SE_HepG2, scale = "none",  
                   cluster_cols = T, cluster_rows =T,
                   color=colorRampPalette(c("navy", "white", "red"))(50),
                   cutree_cols=4, cutree_rows=6)

pheatmap::pheatmap(ECODE_ALS_SE_K562, scale = "none",  
                   cluster_cols = T, cluster_rows =T,
                   color=colorRampPalette(c("navy", "white", "red"))(50),
                   cutree_cols=4, cutree_rows=6)

pheatmap::pheatmap(ECODE_ALS_RI, scale = "none",  
                   cluster_cols = T, cluster_rows =T,
                   color=colorRampPalette(c("navy", "white", "red"))(50),
                   cutree_cols=4, cutree_rows=6)

pheatmap::pheatmap(ECODE_ALS_SE, scale = "none",  
                   cluster_cols = T, cluster_rows =T,
                   color=colorRampPalette(c("navy", "white", "red"))(50),
                   cutree_cols=4, cutree_rows=6)

pheatmap::pheatmap(ECODE_ALS_SE_Binary, scale = "none",  
                   cluster_cols = T, cluster_rows =T,
                   color=colorRampPalette(c("navy", "white", "red"))(50),
                   cutree_cols=4, cutree_rows=6)


#### ALS RBP CrypSplice ####
setwd("/home/weilan/ENCODE/5_CrypSplice/ALS_consolidated/Gene_symbols_udated/")
ECODE_ALS_CrypSplice <- read.csv("aggregated_significant_juncIDs_with_annotation.csv", 
                                 header = TRUE) 

# Remove "JS_diff_" from column names
names(ECODE_ALS_CrypSplice) <- gsub("JS_diff_", "", names(ECODE_ALS_CrypSplice))

library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)

str(ECODE_ALS_CrypSplice)

# Adjusted function to include nrow parameter
create_combined_scatter_plot_with_significance <- function(ECODE_ALS_CrypSplice, nrow=4) {
  gene_targets <- unique(gsub("^(K562_|HepG2_)", "", names(ECODE_ALS_CrypSplice)[grep("^(K562_|HepG2_)", names(ECODE_ALS_CrypSplice))]))
  combined_plot_data <- data.frame()
  
  for (target in gene_targets) {
    k562_col <- paste("K562", target, sep="_")
    hepg2_col <- paste("HepG2", target, sep="_")
    
    plot_data <- ECODE_ALS_CrypSplice %>%
      select(juncID, gene_name, annotation, all_of(k562_col), all_of(hepg2_col)) %>%
      mutate(
        K562 = as.numeric(as.character(!!sym(k562_col))),
        HepG2 = as.numeric(as.character(!!sym(hepg2_col))),
        HepG2_significant = abs(HepG2) > 0.1,
        K562_significant = abs(K562) > 0.1,
        Category = case_when(
          HepG2_significant & K562_significant ~ "Both",
          HepG2_significant & !K562_significant ~ "HepG2 Only",
          !HepG2_significant & K562_significant ~ "K562 Only",
          TRUE ~ "Not Significant"
        ),
        GeneTarget = target
      ) %>%
      select(juncID, gene_name, annotation, K562, HepG2, Category, GeneTarget)
    
    combined_plot_data <- bind_rows(combined_plot_data, plot_data)
  }
  
  # Using nrow parameter in facet_wrap
  p <- ggplot(combined_plot_data, aes(x=K562, y=HepG2, color=Category)) +
    geom_point(alpha=1, size=3) +
    scale_color_manual(values = c("Both" = "red", "HepG2 Only" = "#285D92", "K562 Only" = "orange", "Not Significant" = "grey")) +
    facet_wrap(~GeneTarget, scales = "free", nrow = nrow) +
    labs(x="K562", y="HepG2") +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
    mytheme2 # Ensure mytheme2 is defined or use theme_minimal()
  
  print(p)
}

# Create a long-format DataFrame that lists each juncID with the associated gene targets where it's significant
# Determine significance for each gene target within each cell line
ECODE_ALS_CrypSplice_long <- pivot_longer(ECODE_ALS_CrypSplice, 
                                          cols = starts_with(c("K562", "HepG2")), 
                                          names_to = "Condition_GeneTarget", 
                                          values_to = "Value")

# Determine significance and extract cell line and gene target info
significance_threshold <- 0.1
ECODE_ALS_CrypSplice_long <- ECODE_ALS_CrypSplice_long %>%
  mutate(Significant = abs(Value) > significance_threshold,
         CellLine = ifelse(str_detect(Condition_GeneTarget, "^K562"), "K562", 
                           ifelse(str_detect(Condition_GeneTarget, "^HepG2"), "HepG2", NA)),
         GeneTarget = str_remove(Condition_GeneTarget, "^K562_|^HepG2_")) %>%
  filter(Significant)

# Group by juncID and CellLine, and concatenate gene targets
gene_target_lists <- ECODE_ALS_CrypSplice_long %>%
  group_by(juncID, CellLine) %>%
  summarise(GeneTargets = toString(unique(GeneTarget)), .groups = 'drop')

# Spread concatenated gene targets into separate columns for K562 and HepG2
gene_targets_spread <- gene_target_lists %>%
  pivot_wider(names_from = CellLine, values_from = GeneTargets, values_fill = list(GeneTargets = ""))

# Count distinct gene targets affecting each juncID and combine with gene target names
counts_with_names <- ECODE_ALS_CrypSplice_long %>%
  group_by(juncID) %>%
  summarise(K562_Count = sum(CellLine == "K562", na.rm = TRUE),
            HepG2_Count = sum(CellLine == "HepG2", na.rm = TRUE),
            Both_Count = sum(CellLine %in% c("K562", "HepG2"), na.rm = TRUE)) %>%
  left_join(gene_targets_spread, by = "juncID") %>%
  ungroup()

# Optionally, merge these counts and names with the original data for context
final_details <- left_join(ECODE_ALS_CrypSplice %>% select(juncID, gene_name, annotation) %>% distinct(),
                           counts_with_names, by = "juncID")

# Filter for juncIDs with Both_Count >= 2
final_details_filtered <- final_details %>%
  filter(Both_Count >= 2)

# View the filtered resulting table
head(final_details_filtered)

# Save the filtered final table
write.csv(final_details_filtered, "juncID_gene_target_impact_details_filtered.csv", row.names = FALSE)