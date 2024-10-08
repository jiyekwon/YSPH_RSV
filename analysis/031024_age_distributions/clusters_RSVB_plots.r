#################################original plots################
install.packages("tidyverse")
library(readxl)
clusters_B <- read_xlsx("Clusters_B.xlsx")
library(dplyr)

convert_to_years <- function(age) {
  # Extract the numeric part and the unit
  num <- as.numeric(gsub("[^0-9]", "", age))  # Extract the number
  unit <- gsub("[0-9]", "", age)              # Extract the unit
  
  # Convert based on the unit
  if (grepl("dys", unit)) {
    return(num / 365)  # Convert days to years
  } else if (grepl("mos", unit)) {
    return(num / 12)   # Convert months to years
  } else if (grepl("yrs", unit)) {
    return(num)        # Leave years as is
  } else {
    return(NA)         # Handle any unexpected cases
  }
}

clusters_B <- clusters_B %>%
  mutate(age_in_years = sapply(Age, convert_to_years))

library(ggplot2)
library(gridExtra)
library(dplyr)

p1 <- ggplot(clusters_B, aes(x = Expansion, y = age_in_years, fill = expansion_colors)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Age Distribution by Expansion (Box Plot)",
       y = "Age (years)",
       x = "Expansion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(ggridges)
p5 <- ggplot(clusters_B, aes(x = age_in_years, y = Expansion, fill = Expansion)) +
  geom_density_ridges(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Age Distribution by Expansion (Ridge Plot)",
       x = "Age (years)",
       y = "Expansion")

# Alternatively, save plots to PDF
pdf("age_distribution_plots.pdf", width = 12, height = 15)
grid.arrange(p1), ncol = 2)
dev.off()

#############################################################################
# Load required libraries
library(ggplot2)
library(gridExtra)
library(dplyr)
library(ggridges)
library(lubridate)

# Convert date string to Date object
clusters_B <- clusters_B %>%
  mutate(date = as.Date(Date, format = "%m/%d/%Y"))

# Define custom colors for each expansion
# Replace these with your actual expansion names and desired colors
expansion_colors <- c(
  "A" = "red",
  "B" = "blue",
  "C" = "green",
  "D" = "purple",
  "E" = "orange"
)

# First set of plots: Age Distribution
# 1. Box Plot
p1 <- ggplot(clusters_B, aes(x = Expansion, y = age_in_years, fill = Expansion)) +
  geom_boxplot() +
  scale_fill_manual(values = expansion_colors) +
  theme_minimal() +
  labs(title = "Age Distribution by Expansion (Box Plot)",
       y = "Age (years)",
       x = "Expansion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Second set of plots: Temporal Analysis
# 1. Age vs Time Scatter Plot
t1 <- ggplot(clusters_B, aes(x = date, y = age_in_years, color = Expansion)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = expansion_colors) +
  theme_minimal() +
  labs(title = "Age Distribution Over Time",
       x = "Date",
       y = "Age (years)")
