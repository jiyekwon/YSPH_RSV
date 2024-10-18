install.packages("patchwork")

# Reseq plate heatmaps ----------------------------------------------------


library(tidyverse)
library(reshape2)
library(viridis)

data <- read.table("/vast/palmer/pi/grubaugh/datasets/pgcoe/RSV/241009_reseq_plate1/summary/RSVB_ampdepth.txt", 
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Ensure depth is numeric
data$depth <- as.numeric(data$depth)
data$name <- factor(data$name, levels = paste("Amplicon", 1:50))
data$sample <- factor(data$sample, levels = rev(unique(data$sample)))
data$depth[data$depth == 0] <- 0.0001
str(data)


# Create the plot
plate1_reseq <- ggplot(data, aes(x = name, y = sample, fill = pmin(depth, 500))) +
  geom_tile() +
  scale_fill_gradient2(
    colors = c("white", "#ffcccc", "#ff9999", "#ff6666", "#ff3333"),
    values = c(0, 0.1, 0.25, 0.5, 1),  
    scale_x_continuous(limits = c(0,500)),  
    breaks = c(0, 100, 250, 500),  
    labels = c("0", "100", "250", "500+"),
    name = "Coverage\nDepth"
  ) +
  # Highlight negative values in red
  geom_tile(data = subset(data, depth < 0), 
            aes(x = name, y = sample), 
            fill = "red") +
  # Customize theme
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    panel.grid.major = element_blank(),  # Remove grid
    panel.grid.minor = element_blank(),  # Remove grid
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("Coverage Depth Against RSVB, Plate1 ReSeq") +
  # Add labels for negative values
  geom_text(data = subset(data, depth < 0),
            aes(label = "Error"),
            color = "white",
            size = 2)
print(plate1_reseq)

# old plate heatmap -------------------------------------------------------
library(dplyr)
library(ggplot2)

file_path <- "/vast/palmer/pi/grubaugh/datasets/pgcoe/RSV/240613_RSV001-006/summary/RSVB_ampdepth.txt"
data94 <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
print(str(data94))

amplicon_data_94 <- data94 %>%
  filter(grepl("^Yale-RSV-00[0-9]{1,2}$", sample)) %>%  
  filter(as.numeric(gsub("Yale-RSV-00", "", sample)) <= 94)
print(head(amplicon_data_94))
str(amplicon_data_94)
amplicon_data_94$depth <- as.numeric(amplicon_data_94$depth)
amplicon_data_94$name <- factor(amplicon_data_94$name, levels = paste("Amplicon", 1:50))
amplicon_data_94$sample <- factor(amplicon_data_94$sample, levels = rev(unique(amplicon_data_94$sample)))
data94$depth[data94$depth == 0] <- 0.0001

plate1_orig <- ggplot(amplicon_data_94, aes(x = name, y = sample, fill = pmin(depth, 500))) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("white", "#ffcccc", "#ff9999", "#ff6666", "#ff3333"),
    values = c(0, 0.1, 0.25, 0.5, 1),
    scale_x_continuous(limits = c(0,500)),
    breaks = c(0, 100, 250, 500),
    labels = c("0", "100", "250", "500+"),
    name = "Coverage\nDepth"
  ) +
  geom_tile(data = subset(amplicon_data_94, depth < 0), 
            aes(x = name, y = sample), 
            fill = "red") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 3),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("Coverage Depth Against RSVB, Plate1 Original") +
  geom_text(data = subset(amplicon_data_94, depth < 0),
            aes(label = "Error"),
            color = "white",
            size = 2)
print(plate1_orig)

# patchwork ---------------------------------------------------------------
plate1_total <- plate1_orig / plate1_reseq
ggsave("plate1_comparison.png", plot = plate1_total, width = 10, height = 8, dpi = 300)
