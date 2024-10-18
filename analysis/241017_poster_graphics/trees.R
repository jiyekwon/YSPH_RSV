# install everything ------------------------------------------------------
install.packages("BiocManager")
install.packages("rlang")
install.packages("ggtree")
install.packages("tidytree")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("ape")
install.packages("jsonlite")
install.packages("tidyr")
install.packages("crayon")
install.packages("tidyverse")

library(tidyr)
library(BiocManager)
library(rlang)
library(ggtree)
library(ggplot2)
library(tidytree)
library(dplyr)
library(ape)
library(jsonlite)
library(crayon)
library(tidyverse)
# pull meta and trees -----------------------------------------------------
treeA <- read.tree("RSVA_tree.nwk")
treeB <- read.tree("RSVB_tree.nwk")
json_dataA <- fromJSON("RSVA_auspice.json")
json_dataB <- fromJSON("RSVB_auspice.json")
metadata_file_gisaid <- read.table("rsv_metadata_gisaid.txt", header = TRUE, stringsAsFactors = FALSE)
metadata_file_yale <- read.table("rsv_metadata.txt", header = TRUE, stringsAsFactors = FALSE)
metadata_file_gisaid <- metadata_file_gisaid %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d"))
metadata_separated_yale <- metadata_file_yale %>%
  separate(location, into = c("continent", "country", "region"), sep = " / ", extra = "merge", fill = "right")
metadata_separated_yale <- metadata_separated_yale %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d"))
metadata_total <- bind_rows(metadata_file_gisaid, metadata_separated_yale)

# plot trees --------------------------------------------------------------
treeA <- read.tree("RSVA_tree.nwk")
tip_dataA <- data.frame(
  label = treeA$tip.label,
  continent = sample(c("North America", "South America", "Europe", "Asia", "Africa", "Oceania"), 
                     size = length(treeA$tip.label), replace = TRUE),
  stringsAsFactors = FALSE
)
tip_dataA$category <- ifelse(grepl("^Yale-", tip_data$label), "Yale", tip_data$continent)
color_palette <- c(
  "North America" = "#d6604d",
  "South America" = "#f4a582",
  "Europe" = "#9970ab",
  "Asia" = "#2166ac",
  "Africa" = "#b8e186",
  "Oceania" = "#053061",
  "Yale" = "#ff000d"
)
p <- ggtree(treeA) %<+% tip_dataA +
  geom_point(aes(color = category), size = 2, alpha = 0.8, na.rm = TRUE) +
  scale_color_manual(values = color_palette, na.translate = FALSE) +
  theme(legend.position = "left") +
  labs(color = "Category")
p <- p + geom_point(data = subset(p$data, !isTip & is.na(category)), 
                    color = "transparent", size = 0)

treeB <- read.tree("RSVB_tree.nwk")
tip_dataB <- data.frame(
  label = treeB$tip.label,
  continent = sample(c("North America", "South America", "Europe", "Asia", "Africa", "Oceania"), 
                     size = length(treeB$tip.label), replace = TRUE),
  stringsAsFactors = FALSE
)
tip_dataB$category <- ifelse(grepl("^Yale-", tip_dataB$label), "Yale", tip_dataB$continent)
pB <- ggtree(treeB) %<+% tip_dataB +
  geom_point(aes(color = category), size = 2, alpha = 0.8, na.rm = TRUE) +
  scale_color_manual(values = color_palette, na.translate = FALSE) +
  theme(legend.position = "left") +
  labs(color = "Category")
pB <- pB + geom_point(data = subset(p$data, !isTip & is.na(category)), 
                    color = "transparent", size = 0)
print(pB)
print(p)
ggsave("RSVA_tree_R.png", p, width = 12, height = 8, dpi = 300)
ggsave("RSVA_tree_R.png", p, width = 12, height = 8, dpi = 300)
