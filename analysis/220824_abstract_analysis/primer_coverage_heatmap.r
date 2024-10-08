install.packages(c("ggplot2", "reshape2"))
library(ggplot2)
library(reshape2)

read_depth_files <- function(folder_path, genome_type) {
  # Adjust the pattern to match either RSVA or RSVB specifically
  files <- list.files(folder_path, pattern = paste0("\\.", genome_type, "\\.depth\\.txt$"), full.names = TRUE)
  
  # Check if any files were found
  if (length(files) == 0) {
    warning(paste("No files found for genome type", genome_type))
    return(NULL)
  }
  
  # Read and combine all files
  combined_data <- lapply(files, function(file) {
    # Read the file with explicit column names (assuming Position and Depth)
    data <- read.table(file, header = FALSE, col.names = c("Genome", "Position", "Depth"))
    if (nrow(data) == 0) {
      warning(paste("File", file, "is empty or could not be read"))
    } else {
      # Extract sample name from file name
      data$Sample <- gsub(paste0("\\.", genome_type, "\\.depth\\.txt$"), "", basename(file))
    }
    return(data)
  })
  
  # Combine the data for all samples
  combined_data <- do.call(rbind, combined_data)
  return(combined_data)
}

folder_path <- "/gpfs/gibbs/project/grubaugh/emd224/PGCOE/Analysis/Primer_Efficency/test_primers/depthfiles"

reshape_for_heatmap <- function(data) {
  # Ensure the 'Depth' column exists
  if (!"Depth" %in% colnames(data)) stop("Column 'Depth' not found in data")
  # Reshape data to wide format, using Position as the row identifier and Sample as columns
  data_wide <- reshape2::dcast(data, Position ~ Sample, value.var = "Depth", fill = 0)
  return(data_wide)
}

rsva_wide <- reshape_for_heatmap(rsva_data)
rsvb_wide <- reshape_for_heatmap(rsvb_data)

highlight_regions_rsva <- data.frame(
  xmin = c(4673, 5648),  # Example coordinates for RSVA
  xmax = c(5595, 7550),
  ymin = -Inf,
  ymax = Inf
)

highlight_regions_rsvb <- data.frame(
  xmin = c(4675, 5653),  # Example coordinates for RSVB
  xmax = c(5600, 7552),
  ymin = -Inf,
  ymax = Inf
)

plot_heatmap <- function(data_wide, genome_type) {
  melted_data <- reshape2::melt(data_wide, id.vars = "Position")
  
  ggplot(melted_data, aes(x = Position, y = variable, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue") +
    labs(title = paste("Depth Heatmap -", genome_type), x = "Position", y = "Sample", fill = "Depth") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}
plot_heatmap <- function(data_wide, genome_type) {
  # Melt the data
  melted_data <- reshape2::melt(data_wide, id.vars = "Position")
  
  # Custom color scale with specific breaks for easier interpretation
  ggplot(melted_data, aes(x = Position, y = variable, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = c("white", "yellow", "orange", "red", "blue"),
      values = scales::rescale(c(0, 10, 100, 1000, max(melted_data$value))),
      breaks = c(0, 10, 100, 1000),
      labels = c("0X", "10X", "100X", "1000X+"),
      limits = c(0, 1000)
    ) +
    labs(title = paste("Depth Heatmap -", genome_type), x = "Position", y = "Sample", fill = "Depth") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


plot_heatmap <- function(data_wide, genome_type) {
  melted_data <- reshape2::melt(data_wide, id.vars = "Position")
  p <- ggplot(melted_data, aes(x = Position, y = variable, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = c("white", "lightblue", "yellow", "orange", "red", "purple", "black"), 
      values = scales::rescale(c(0, 1, 10, 100, 1000, 10000, 100000)),
      breaks = c(0, 10, 100, 1000, 10000, 100000),  # Define the breaks (6 breaks)
      labels = c("0", "1-10", "11-100", "101-1000", "1001-10000", "10001-100000"),  
      na.value = "white"
    ) +
    labs(title = paste("Depth Heatmap -", genome_type), x = "Position", y = "Sample", fill = "Depth") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.key.height = unit(2, "cm"),
          legend.key.width = unit(0.5, "cm"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.spacing.y = unit(0.5, "cm")) +
    guides(fill = guide_legend(
      keywidth = unit(1.5, "cm"),  # Width of legend keys
      keyheight = unit(1.5, "cm"), # Height of legend keys
      title.position = "top",      # Position of the legend title
      title.hjust = 0.5,           # Center the legend title
      label.position = "bottom",   # Position of labels relative to keys
      label.hjust = 0.5            # Center the labels
    ))
  p + geom_rect(data = highlight_regions, 
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = "red", alpha = 0.2, color = "black")
}


plot_heatmap <- function(data_wide, genome_type, highlight_regions = NULL) {
  # Melt the data
  melted_data <- reshape2::melt(data_wide, id.vars = "Position")
  
  # Define color breaks and labels
  breaks <- c(0, 1, 10, 100, 1000, 10000, 100000)
  labels <- c("0", "1-10", "11-100", "101-1000", "1001-10000", "10001-100000", "100000+")
  
  # Plot the heatmap
  p <- ggplot(melted_data, aes(x = Position, y = variable, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = c("white", "lightblue", "yellow", "orange", "red", "purple", "black"),
      values = scales::rescale(breaks),
      breaks = breaks,
      labels = labels,
      na.value = "white"
    ) +
    labs(title = paste("Depth Heatmap -", genome_type), x = "Position", y = "Sample", fill = "Depth") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.key.height = unit(2, "cm"),
          legend.key.width = unit(0.5, "cm"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.spacing.y = unit(0.5, "cm")) +
    guides(fill = guide_legend(
      keywidth = unit(1.5, "cm"),
      keyheight = unit(1.5, "cm"),
      title.position = "top",
      title.hjust = 0.5,
      label.position = "bottom",
      label.hjust = 0.5
    ))
  
  if (!is.null(highlight_regions) && nrow(highlight_regions) > 0) {
    p <- p + geom_rect(data = highlight_regions, 
                       aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                       fill = "transparent", color = "black", linewidth = 1.5, linetype = "solid",
                       inherit.aes = FALSE)
  }
  return(p)
}




rsva_data <- read_depth_files(folder_path, "RSVA")
rsvb_data <- read_depth_files(folder_path, "RSVB")

print(head(rsva_data))
print(head(rsvb_data))

plot_heatmap(rsva_wide, "RSVA", highlight_regions_rsva)
plot_heatmap(rsvb_wide, "RSVB", highlight_regions_rsvb)
