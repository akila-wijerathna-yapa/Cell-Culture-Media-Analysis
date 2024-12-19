################################################################################

# Read the Sample_Names.txt file
sample_names <- read.csv("Data/Sample_Names.txt", header = FALSE, stringsAsFactors = FALSE)

# Assign column names for clarity
colnames(sample_names) <- c("Original", "GroupName_Replicate")

# Step 1: Extract Group Names
sample_names$Group <- sub("_\\d+$", "", sample_names$GroupName_Replicate)  # Remove replicate numbers

# Step 2: Get the Number of Groups and Samples
num_groups <- length(unique(sample_names$Group))  # Unique group names
num_samples <- nrow(sample_names)                # Total number of samples

# Print the results
print(paste("Number of groups:", num_groups))
print(paste("Number of samples:", num_samples))

# Optional: View groups and their sample counts
group_counts <- table(sample_names$Group)
print(group_counts)



################################################################################

# Check for any negative values across numeric columns
any_negative <- any(df_normalized  %>% 
                      select(where(is.numeric)) %>% 
                      unlist() < 0, na.rm = TRUE)

# Print result
if (any_negative) {
  print("There are negative values in the dataframe.")
} else {
  print("No negative values found in the dataframe.")
}

################################################################################


# Define a color palette with 5 distinct colors
color_palette <- c( "#A6CEE3" , "#1F78B4" , "#33A02C" , "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00")

# c( "#A6CEE3" , "#1F78B4" , "#B2DF8A" , "#33A02C" , "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00" , "#CAB2D6" ,"darkgreen")

# Manually define 10 intervals (breaks)
breaks <- seq(min(heatmap_data), max(heatmap_data), length.out = 8)

# Use custom breaks in pheatmap
pheatmap(
  heatmap_data,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = color_palette,  # 10 colors
  breaks = breaks,        # Ensure values are evenly mapped to colors
  main = "Heatmap with 10 Distinct Colors and Custom Breaks",
  show_rownames = FALSE,
  show_colnames = TRUE
)

