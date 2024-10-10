# Load metadata
data <- read.delim("rsv_metadata.txt", header = TRUE, stringsAsFactors = FALSE)

data <- read.csv('rsv_metadata.txt')

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(lubridate)

data$date <- as.Date(data$date, format="%Y-%m-%d")

# Example data
data <- data.frame(
  agestr = c("120 days", "6 months", "2 years", "3 months", "500 days"),
  stringsAsFactors = FALSE
)

# Function to convert ages to days
convert_to_days <- function(age_str) {
  # Extract numeric value
  numeric_value <- as.numeric(gsub("[^0-9]", "", age_str))
  
  # Convert based on unit
  if (grepl("day", age_str, ignore.case = TRUE) || grepl("days", age_str, ignore.case = TRUE)) {
    return(numeric_value)
  } else if (grepl("month", age_str, ignore.case = TRUE) || grepl("mos", age_str, ignore.case = TRUE)) {
    return(numeric_value * 30) # Approximate conversion
  } else if (grepl("year", age_str, ignore.case = TRUE) || grepl("yrs", age_str, ignore.case = TRUE)) {
    return(numeric_value * 365) # Approximate conversion
  } else {
    return(NA) # Handle unknown units
  }
}
data$AgeInDays <- sapply(data$agestr, function(x) convert_to_days(x))
data$Age_in_Years <- data$AgeInDays / 365

print(data)

# Define age groups
data <- data %>%
  mutate(AgeGroup = case_when(
    Age <= 5 ~ "0-5",
    Age <= 10 ~ "6-10",
    Age <= 18 ~ "11-18",
    Age <= 65 ~ "19-65",
    TRUE ~ "65+"
  ))

# Assume all samples are positive
data$Status <- "Positive"

ggplot(data, aes(x = date, fill = AgeGroup)) +
  geom_histogram(binwidth = 7, position = "stack") +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  labs(
    title = "Transmission of Virus from Young Children to Adults",
    x = "Collection Date",
    y = "Number of Positive Cases",
    fill = "Age Group"
  ) +
  theme_minimal()

# Save the plot
ggsave("virus_transmission_plot.png")
