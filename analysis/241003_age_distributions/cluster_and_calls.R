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


# cluster graphing bullshit -----------------------------------------------

age_order <- c("Geriatrics", "Adults", "Children", "Preschool", "Infants")
allmeta$date <- as.Date(allmeta$date, format = "%Y-%m-%d")

# plotting clusters -------------------------------------------------------
#all vs none, by agecat
cladegrouping <- ggplot(allmeta, 
                        aes(x = ifelse(Clade_Name == "no_clade", "No Clade", "Has Clade"), 
                            fill = factor(agecat, levels = age_order))) +  
  geom_bar(position = "fill") +  
  theme_minimal() +
  labs(title = "Proportion of Age Categories by Clade Status",
       x = "Clade Category",
       y = "Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#clades across time, by agecat
cladetime <- ggplot() + 
  geom_dotplot(data = allmeta %>% filter(Clade_Name != "no_clade" & !is.na(agecat)),
               aes(x = Clade_Name, y = as.Date(date), fill = agecat),
               binaxis = "y", 
               stackdir = "center", 
               alpha = 1, 
               binwidth = 1, 
               stackgroups = TRUE) + 
  coord_flip() +  scale_fill_brewer(palette = "Set2", name = "Age Category") +
  labs(x = "Clade Name", 
       y = "Sampling Date", 
       title = "Stacked Dotplot of Clades by Sampling Date and Age Category") +
  theme_minimal()
