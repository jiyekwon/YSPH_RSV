# Load required libraries
library(tidyverse)
library(lubridate)

# Assuming your data is in a dataframe called 'df'
# If it's not, you'll need to read it in first

# Data preparation
df <- df %>%
  mutate(
    date_collected = ymd(date_collected),
    season = case_when(
      month(date_collected) %in% c(10, 11) ~ "Early (Oct/Nov)",
      month(date_collected) %in% c(12, 1) ~ "Middle (Dec/Jan)",
      month(date_collected) %in% c(2, 3) ~ "Late (Feb/Mar)",
      TRUE ~ "Other"
    )
  )

# Function to create plot for each RSV type
create_rsv_plot <- function(data, rsv_type) {
  ggplot(data %>% filter(str_detect(combined_column, rsv_type)), 
         aes(x = age_group, fill = age_group)) +
    geom_bar() +
    facet_wrap(~ season, scales = "free_y") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(title = paste("RSV", rsv_type, "Cases by Age Group and Season"),
         x = "Age Group",
         y = "Number of Cases") +
    scale_fill_brewer(palette = "Set3")
}

# Create plots
plot_a <- create_rsv_plot(df, "A")
plot_b <- create_rsv_plot(df, "B")
plot_ab <- create_rsv_plot(df, "A/B")

# Display plots
print(plot_a)
print(plot_b)
print(plot_ab)

# If you want to save the plots, you can use ggsave:
# ggsave("rsv_a_plot.png", plot_a, width = 12, height = 8)
# ggsave("rsv_b_plot.png", plot_b, width = 12, height = 8)
# ggsave("rsv_ab_plot.png", plot_ab, width = 12, height = 8)





####by genotype #####

install.packages("lubridate")
library(lubridate)
install.packages("stringr")
library(stringr)


df3 <- read.csv("final_calls.csv")
df2 <- read_xlsx("Accurate_names_tube_IDS.xlsx")

names(df3)[names(df3) == "sample"] <- "SeqID"
df <- merge(df3, df2, by = "SeqID")
unmatched_in_merged_df <- merged_df[!merged_df$SeqID %in% df2$SeqID, ]
unmatched_in_df2 <- df2[!df2$SeqID %in% df$SeqID, ]
print(unmatched_in_df3)
print(unmatched_in_df2)


# Data preparation
df <- merged_df %>%
  mutate(
    date_collected = ymd(Collected),  # Ensure Collected is in the proper date format
    season = case_when(
      month(date_collected) %in% c(10, 11) ~ "Early (Oct/Nov)",
      month(date_collected) %in% c(12, 1) ~ "Middle (Dec/Jan)",
      month(date_collected) %in% c(2, 3) ~ "Late (Feb/Mar)",
      TRUE ~ "Other"
    ),
    season = factor(season, levels = c("Early (Oct/Nov)", "Middle (Dec/Jan)", "Late (Feb/Mar)", "Other")),
    rsv_type = case_when(
      type == "RSVA" ~ "RSV A",
      type == "RSVB" ~ "RSV B",
      type == "RSVA/B" ~ "RSV A/B",
      TRUE ~ "Unknown"
    ),
    age_group = case_when(
      str_detect(Age, "mos") ~ "<1 year",  # Handle ages in months
      str_detect(Age, "yr") ~ paste(str_extract(Age, "\\d+"), "year"),
      str_detect(Age, "yrs") ~ paste(str_extract(Age, "\\d+"), "years"),
      TRUE ~ Age
    ),
    age_group = case_when(
      age_group == "<1 year" ~ "<1 year",  # Keep the <1 year group intact
      age_group == "1 year" ~ "1-4 years",
      as.numeric(str_extract(age_group, "\\d+")) < 1 ~ "<1 year",  # Retain this check for safety
      as.numeric(str_extract(age_group, "\\d+")) <= 4 ~ "1-4 years",
      as.numeric(str_extract(age_group, "\\d+")) <= 17 ~ "5-17 years",
      as.numeric(str_extract(age_group, "\\d+")) <= 29 ~ "18-29 years",
      as.numeric(str_extract(age_group, "\\d+")) <= 49 ~ "30-49 years",
      as.numeric(str_extract(age_group, "\\d+")) <= 64 ~ "50-64 years",
      TRUE ~ "65+ years"
    ),
    age_group = factor(age_group, levels = c("<1 year", "1-4 years", "5-17 years", "18-29 years", "30-49 years", "50-64 years", "65+ years"))
  ) %>%
  filter(season != "Other")

# Define custom color palette
season_colors <- c("Early (Oct/Nov)" = "salmon", "Middle (Dec/Jan)" = "lightgreen", "Late (Feb/Mar)" = "deepskyblue", "Other" = "gray")

combined_plot_AB_A_B <- ggplot(df, aes(x = age_group, fill = season)) +
  geom_bar(position = "dodge") +
  facet_grid(rsv_type ~ season, scales = "free_y") +  # Added scales = "free_y"
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(color = "black", face = "bold"),
    panel.spacing.y = unit(1, "lines")
  ) +
  labs(
    title = "RSV Cases by Type, Age Group, and Season",
    x = "Age Group",
    y = "Number of Cases",
    fill = "Season"
  ) +
  scale_fill_manual(values = season_colors) +
  scale_x_discrete(drop = FALSE)

print(combined_plot_AB_A_B)
ggsave("combined_rsv_plot.png", combined_plot, width = 15, height = 10)

summary_table <- df %>%
  group_by(rsv_type, season, age_group) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = season, values_from = count, values_fill = 0) %>%
  arrange(rsv_type, age_group)

print(kable(summary_table, caption = "Summary of RSV Cases by Type, Season, and Age Group"))
write_xlsx(summary_table, "rsv_summary_table.xlsx")

##stack plots###

summary_df <- df %>%
  group_by(age_group, season) %>%
  summarise(total_cases = n(), .groups = 'drop')  # Count occurrences

totalseasonplot_combined <- ggplot(summary_df, aes(x = age_group, y = total_cases, fill = season)) +
  geom_bar(stat = "identity", position = "dodge") +  # Use stat = "identity" for total counts
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(color = "black", face = "bold"),
    panel.spacing.y = unit(1, "lines")
  ) +
  labs(
    title = "Total RSV Cases by Age Group and Season",
    x = "Age Group",
    y = "Total Number of Cases",
    fill = "Season"
  ) +
  scale_fill_manual(values = season_colors) +
  scale_x_discrete(drop = FALSE)

print(totalseasonplot_combined)
ggsave("totalseasonplot_combined.png", combined_plot, width = 15, height = 10)


totalseasonplot <- ggplot(summary_df, aes(x = age_group, y = total_cases, fill = season)) +
  geom_bar(stat = "identity", position = "dodge") +  # Use stat = "identity" for total counts
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(color = "black", face = "bold"),
    panel.spacing.y = unit(1, "lines")
  ) +
  labs(
    title = "Total RSV Cases by Age Group and Season",
    x = "Age Group",
    y = "Total Number of Cases",
    fill = "Season"
  ) +
  scale_fill_manual(values = season_colors) +
  scale_x_discrete(drop = FALSE) +
  facet_wrap(~season)  # Facet by season

print(totalseasonplot)

install.packages("patchwork")
library(patchwork)

combinedplots <- (totalseasonplot / combined_plot_AB_A_B)
ggsave(filename = "RSV_Season_24.png", plot = combinedplots, width = 12, height = 8, dpi = 300)
