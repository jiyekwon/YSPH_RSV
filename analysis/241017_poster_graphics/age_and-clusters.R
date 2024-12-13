# adding in meta sections -------------------------------------------------

install.packages("readxl")
install.packages("tidyverse")
library(readxl)
library(tidyverse)

allmeta <- read.table("data/rsv_metadata.txt",sep="\t",header=T)
cluster_list <- read_xlsx("FINAL_clusters_B_Post_April_1_23.xlsx")
call_list <- read.table("/vast/palmer/pi/grubaugh/datasets/pgcoe/RSV/240613_RSV001-006/summary/final_calls.txt", header = TRUE, sep = "\t", na.strings = "", fill = TRUE)
call_list <- call_list %>%
  mutate(merged_call = coalesce(call, mashcall))
call_list <- call_list %>% select(-call, -mashcall)
call_list <- call_list %>% 
  rename(
    name = sample
  )
cluster_list <- cluster_list %>%
  rename(
    name = SeqID
  )

allmeta <- merge(allmeta, cluster_list, by = "name", all = TRUE)
allmeta <- merge(allmeta, call_list, by = "name", all = TRUE)
allmeta <- allmeta %>% 
  mutate(Clade_Name = ifelse(is.na(Clade_Name), "no_clade", Clade_Name))

write.table(allmeta, file = "full_rsv_metadata.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
allmeta <- allmeta %>%
  filter(!is.na(merged_call), merged_call != "Other") %>%
  mutate(
    season = case_when(
      month(date) %in% c(10, 11) ~ "Early (Oct/Nov)",
      month(date) %in% c(12, 1) ~ "Middle (Dec/Jan)",
      month(date) %in% c(2, 3) ~ "Late (Feb/Mar)",
      TRUE ~ "Other"  # This should be filtered out later
    ),
    season = factor(season, levels = c("Early (Oct/Nov)", "Middle (Dec/Jan)", "Late (Feb/Mar)")),
    agecat = factor(agecat, levels = c("Infants", "Preschool", "Children", "Adults", "Geriatrics"))
  ) %>%
  filter(season != "Other")

# cluster graphing bullshit -----------------------------------------------

age_order <- c("Geriatrics", "Adults", "Children", "Preschool", "Infants")
allmeta$date <- as.Date(allmeta$date, format = "%Y-%m-%d")

colour_codes <- c("Infants" = "#f1eef6",
                  "Preschool" = "#bdc9e1",
                  "Children" = "#74a9cf",
                  "Adults" = "#2b8cbe",
                    "Geriatrics" = "#045a8d")
season_colors <- c("Early (Oct/Nov)" = "salmon", "Middle (Dec/Jan)" = "lightgreen", "Late (Feb/Mar)" = "deepskyblue", "Other" = "gray")

# plotting clusters -------------------------------------------------------
#all vs none, by agecat
cladegrouping <- ggplot(allmeta %>% filter(!is.na(agecat)), 
                        aes(x = ifelse(Clade_Name == "no_clade", "No Clade", "Has Clade"), 
                            fill = factor(agecat, levels = age_order))) +  
  geom_bar(position = "fill") + 
  scale_fill_manual(values = colour_codes) + 
  theme_minimal() +
  labs(title = "Proportion of Age Categories by Clade Status",
       x = "Clade Category",
       y = "Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#clades across time, by agecat
cladetime <- ggplot() + 
  geom_dotplot(data = allmeta %>% filter(Clade_Name != "no_clade" & !is.na(agecat)),
               aes(x = Clade_Name, y = as.Date(date), fill = factor(agecat, levels = age_order)),
               binaxis = "y", 
               stackdir = "center", 
               alpha = 1, 
               binwidth = 1, 
               stackgroups = TRUE) + 
  coord_flip() +  scale_fill_manual(values = colour_codes) +
  scale_y_date(date_breaks = "2 weeks", date_labels = "%b %d, %Y") +
  labs(x = "Clade Name", 
       y = "Sampling Date", 
       title = "Stacked Dotplot of Clades by Sampling Date and Age Category") +
  theme_minimal()

# staked histogram across season ages ------------------------------------------
create_rsv_plot <- function(allmeta, rsv_type) {
  ggplot(data %>% filter(str_detect(merged_call, rsv_type)), 
         aes(x = agecat, fill = agecat)) +
    geom_bar() +
    facet_wrap(~ season, scales = "free_y") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(title = paste("RSV", rsv_type, "Cases by Age Category and Season"),
         x = "Age Category",
         y = "Number of Cases") +
    scale_fill_manual(values = colour_codes)
}

plot_a <- create_rsv_plot(allmeta, "A")
plot_b <- create_rsv_plot(df, "B")
plot_ab <- create_rsv_plot(df, "A/B")
print(plot_a)
print(plot_b)
print(plot_ab)


combined_plot_AB_A_B <- ggplot(allmeta, aes(x = agecat, fill = agecat)) +
  geom_bar(position = "dodge") +
  facet_grid(merged_call ~ season, scales = "free_y") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(color = "black", face = "bold"),
    panel.spacing.y = unit(1, "lines")
  ) +
  labs(
    title = "RSV Cases by Merged Call, Age Group, and Season",
    x = "Age Group",
    y = "Number of Cases",
    fill = "Season"
  ) +
  scale_fill_manual(values = colour_codes) +
  scale_x_discrete(drop = FALSE)

# Print the plot
print(combined_plot_AB_A_B)
ggsave("combined_rsv_plot.png", combined_plot_AB_A_B, width = 15, height = 10)



