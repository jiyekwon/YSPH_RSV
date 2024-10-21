# adding in meta sections -------------------------------------------------

# install.packages("readxl")
# install.packages("tidyverse")
# install.packages("svglite")
library(readxl)
library(tidyverse)
library(dplyr)
library(svglite)


allmeta <- read.table("data/rsv_metadata.txt",sep="\t",header=T)
cluster_list <- read_xlsx("FINAL_clusters_B_Post_April_1_23.xlsx")
call_list <- read.table("data/RSVAB_final_calls.txt", header = TRUE, sep = "\t", na.strings = "", fill = TRUE)
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

allmeta$month <- month(as.Date(allmeta$date))
allmeta$season_month <- allmeta$month-9
allmeta$season_month[allmeta$season_month < 0] <- allmeta$month[allmeta$season_month < 0]+3

allmeta$agecat <- factor(allmeta$agecat,
                         levels=rev(c("Infants","Preschool","Children","Adults","Geriatrics")),
                         ordered=T)

allmeta$month <- factor(allmeta$month, 
                        levels = c(10, 11, 12, 1, 2, 3), 
                        labels = c("October", "November", "December", "January", "February", "March"))

# cluster graphing bullshit -----------------------------------------------

age_order <- c("Geriatrics", "Adults", "Children", "Preschool", "Infants")
allmeta$date <- as.Date(allmeta$date, format = "%Y-%m-%d")

colour_codes <- c("Infants" = "#023858",
                  "Preschool" = "#045a8d",
                  "Children" = "#0570b0",
                  "Adults" = "#2b8cbe",
                    "Geriatrics" = "#74a9cf")

season_colors <- c("Early (Oct/Nov)" = "salmon", "Middle (Dec/Jan)" = "lightgreen", "Late (Feb/Mar)" = "deepskyblue", "Other" = "gray")

# plotting clusters -------------------------------------------------------
#all vs none, by agecat
allmeta$Clade_Name

cladegrouping <- ggplot(allmeta %>% filter(!is.na(agecat)), 
                        aes(x = ifelse(Clade_Name == "no_clade", "No Clade", "Has Clade"), 
                            fill = factor(agecat, levels = age_order))) +  # Set levels for agecat
  geom_bar(position = "fill") + 
  scale_fill_manual(values = colour_codes) + 
  theme_minimal() +
  labs(x = "",
       y = "Proportion",
       fill = "Ages") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_blank(),
        legend.position = "none")
print(cladegrouping)

ggsave("cladegrouping.svg", cladegrouping, device = "svg", width = 5, height = 7, dpi = 400)


#clades across time, by agecat
cladetime <- ggplot() + 
  geom_dotplot(data = allmeta %>% filter(Clade_Name != "no_clade" & !is.na(agecat)),
               aes(x = Clade_Name, y = as.Date(date), fill = factor(agecat, levels = age_order)),
               binaxis = "y", 
               stackdir = "center", 
               alpha = 1, 
               binwidth = 1, 
               stackgroups = TRUE, 
               show.legend = FALSE) +  # Remove the legend
  coord_flip() + 
  scale_fill_manual(values = colour_codes) +
  scale_y_date(date_breaks = "2 weeks", date_labels = "%b %d, %Y") +
  labs(x = "Clade Name", 
       y = "Sampling Date", 
       title = "Stacked Dotplot of Clades by Sampling Date and Age Category") +
  theme_minimal()

print(cladetime)

ggsave("cladetime.png", cladetime, width = 15, height = 10)

# staked histogram across season genotype by ages ------------------------------------------
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



rsv_type <- "A"
agedist_across_season_plot <- ggplot(subset(allmeta,merged_call != "RSVA/B"), 
       aes(x = month, fill = agecat)) +
  geom_bar(position="fill") +
  facet_grid(merged_call ~ ., scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        #legend.position = "none"
        ) +
  labs(title = paste("RSV", "Cases by Age Category and Season"),
       x = "",
       y = "Proportion of cases") +
  scale_fill_manual(values = colour_codes,name="")
agedist_across_season_plot
ggsave("combined_RSVAB_age_dists_across_season.png",dpi=400,units="mm",width=300,height=250)
ggsave("combined_RSVAB_age_dists_across_season.svg",dpi=400,units="mm",width=300,height=250)


combined_plot_AB_A_B <- ggplot(allmeta, aes(x = merged_call, fill = agecat)) +
  geom_bar(position = "stack") +  # Stack the bars by age category
  facet_grid(. ~ season, scales = "free_y") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(color = "black", face = "bold"),
    panel.spacing.y = unit(1, "lines"),
    legend.position = "none"
  ) +
  labs(
    x = "Type",
    y = "",
    fill = "Age Group"
  ) +
  scale_fill_manual(values = colour_codes) +
  scale_x_discrete(drop = FALSE)  
print(combined_plot_AB_A_B)
ggsave("combined_rsv_plot.png", combined_plot_AB_A_B, width = 15, height = 10)

#just A and B not AB poor souls

combined_plot_A_B <- ggplot(allmeta %>% filter(merged_call %in% c("RSVA", "RSVB")), 
                            aes(x = merged_call, fill = factor(agecat, levels = age_order))) + 
  geom_bar(position = "stack") +  # Stack the bars by age category
  facet_grid(. ~ season, scales = "free_y") +  # Facet by season
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(color = "black", face = "bold"),
    panel.spacing.y = unit(1, "lines"),
    legend.position = "none"  # Position legend on the right
  ) +
  labs(
    x = "Type",
    y = "Number of Cases",
    fill = "Age Group"
  ) +
  scale_fill_manual(values = colour_codes) +  
  scale_x_discrete(drop = FALSE)
print(combined_plot_A_B)

ggsave("combined_A_B.svg", combined_plot_A_B, device = "svg", width = 6, height = 10, dpi = 400)

# age graphs -----------------------------------------------------------------

agehist <- ggplot(allmeta,aes(x=floor(age),fill=agecat)) +
  geom_bar(stat="count")  +
  ylab("n")+xlab("age (yrs)") +
  theme(axis.title.y.left = element_text(angle=90,vjust=0.5),
        axis.title.x=element_text(angle = 0),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        legend.position = c(1, 0),
        legend.justification = c(1, 0)) +
  labs(
    x = "Age",
    y = "Number of Cases",
  ) +
  scale_fill_manual(values = colour_codes)
print(agehist)

catdist <- ggplot(allmeta, aes(x = agecat, fill = agecat)) + 
  geom_bar() +
  scale_y_continuous(name = "Number of Cases", sec.axis = sec_axis(trans = ~./nrow(allmeta), name = "%")) +
  theme(axis.title.y.left = element_text(angle = 90, vjust = 0.5),
        axis.title.y.right = element_text(angle = 0, vjust = 0.5),
        axis.title.x = element_text(angle = 0),  
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray93"),  # Major grid lines
        panel.grid.minor = element_line(color = "lightgrey", linetype = "dotted"),  # Minor grid lines (optional)
        legend.position = "none") +
  labs(
    x = "Age",  
    y = "Number of Cases"
  ) +
  scale_fill_manual(values = colour_codes)

print(catdist)

ages_by_month_n <- ggplot(allmeta, aes(x = month, fill = agecat)) + 
  geom_bar(position = "stack") +
  scale_y_continuous(name = "Number of cases", 
                     sec.axis = sec_axis(trans = ~./nrow(allmeta), name = "Proportion")) +
  theme(
    axis.title.y.left = element_text(angle = 90, vjust = 0.5),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 10, b = 40, l = 10),
    panel.background = element_blank(),
    plot.background = element_blank()
  ) +
  xlab("Month") +
  scale_fill_manual(values = colour_codes) +
  ggtitle("Number of cases by age dempgraphic across transmission months in Connecticut, 2024")


print(ages_by_month_n)
ggsave("ages_by_month_N.svg", ages_by_month_n, device = "svg", width = 12, height = 6, dpi = 400)


# formatting --------------------------------------------------------------
(agehist / (cladegrouping | combined_plot_AB_A_B))

poster_graphic <- catdist + combined_plot_AB_A_B + cladegrouping + plot_layout(ncol = 3, widths = c(3,3,1))
print(poster_graphic)
ggsave("age_distribution_plot.png", poster_graphic, width = 25, height = 10)

install.packages("cowplot")
library(cowplot)

combined_plot <- ggdraw(agehist) + draw_plot(combined_plot_AB_A_B, x = 0.5, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(cladegrouping, x = 0.10, y = 0.5, width = 0.45, height = 0.5)
print(combined_plot)


library("patchwork")
layout="AAB
AAC"
agedist_across_season_plot + guide_area() + cladegrouping + plot_layout(design=layout,guides="collect")
ggsave("combined_RSVAB_age_season_inclades.png",dpi=400,units="mm",width=400,height=300)
ggsave("combined_RSVAB_age_season_inclades.svg",dpi=400,units="mm",width=400,height=300)
