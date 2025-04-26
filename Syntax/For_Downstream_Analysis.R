## For upstream analysis syntax refer to "Untargted_POS_NEG.R" file

df_combined <- read.csv("Tables/df_combined_POS_NEG.csv") 

# Remove the 'X' column
df_combined <- df_combined %>% select(-X)

str(df_combined)

### Important

# This file has normalized data.
# For PCA and Heatmap scale the data as below
# After you find metabolites of interest (From Top50, UpSetPlots etc...), for other statistical analysis use "df_combined_POS_NEG.csv" file (only normalized, without scaling); 


###############################################################################

## Scaling

# Select numeric data (exclude metadata)

df_scaled_numeric <- df_combined %>%
  select(-Metabolite.name, -Adduct.type, -Average.Rt.min.) %>% 
  scale() %>%
  as.data.frame()

# Combine the scaled numeric data with metadata
df_combined_scaled <- cbind(
  df_combined %>% select(Metabolite.name, Adduct.type, Average.Rt.min.),
  df_scaled_numeric
)

# Check the structure of the scaled data
str(df_combined_scaled)

# Check scaling results
summary(df_combined_scaled %>% select(-Metabolite.name, -Adduct.type, -Average.Rt.min.))

################################################################################
## PCA

# Prepare data for PCA (exclude metadata)
df_combined_scaled_pca_data <- df_combined_scaled %>%
  select(-Metabolite.name, -Adduct.type, -Average.Rt.min.) %>%
  t() %>%
  as.data.frame()

# Assign rownames to the PCA data
rownames(df_combined_scaled_pca_data) <- colnames(df_pos_normalized %>% select(-Metabolite.name,-Adduct.type, -Average.Rt.min.))

# Perform PCA on scaled data
df_combined_scaled_pca_data_result <- prcomp(df_combined_scaled_pca_data, center = FALSE, scale. = FALSE)

# Create a data frame for plotting PCA results
df_combined_scaled_pca_data <- as.data.frame(df_combined_scaled_pca_data_result$x[, 1:2])  # Extract PC1 and PC2
df_combined_scaled_pca_data$Group <- sapply(rownames(df_combined_scaled_pca_data), function(x) strsplit(x, "_")[[1]][1])
df_combined_scaled_pca_data$Sample <- rownames(df_combined_scaled_pca_data)

# Predefined solid colors for 9 groups
combined_solid_colors <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
  "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22"
)

# Dynamically map colors to groups
combined_group_colors <- setNames(solid_colors, unique(df_combined_scaled_pca_data$Group))

# Step 1: Calculate variance explained by each PC
combined_explained_variance <- df_combined_scaled_pca_data_result$sdev^2 / sum(df_combined_scaled_pca_data_result$sdev^2)
combined_pc1_var <- round(explained_variance[1] * 100, 2)  # PC1 percentage
combined_pc2_var <- round(explained_variance[2] * 100, 2)  # PC2 percentage

# Step 2: Create the PCA plot with PC values in axis labels

combined_pca_plot <- ggplot(df_combined_scaled_pca_data, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, max.overlaps = 20, box.padding = 0.5, point.padding = 0.5) +
  scale_color_manual(values = combined_group_colors) +
  labs(
    title = "PCA Plot with Percent Variance Explained",
    x = paste0("Principal Component 1 (", pc1_var, "%)"),
    y = paste0("Principal Component 2 (", pc2_var, "%)"),
    color = "Group"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  )

# Print the PCA plot
print(combined_pca_plot)

# Save the PCA plot as PDF and JPEG
ggsave("Combined_PCA_Plot_Custom_Colors.pdf", plot = combined_pca_plot, device = "pdf", width = 10, height = 8, dpi = 300)
ggsave("Combined_PCA_Plot_Custom_Colors.jpeg", plot = combined_pca_plot, device = "jpeg", width = 10, height = 8, dpi = 300)


################################################################################
## Heat Map

# Extract numeric data for heatmap
combined_heatmap_data <- df_combined_scaled %>% select(-Metabolite.name,-Adduct.type, -Average.Rt.min.) %>% as.matrix()

color_palette <- colorRampPalette(c( "#A6CEE3" , "#33A02C", "#E31A1C"))(50)

pheatmap(
  combined_heatmap_data,
  cluster_rows = TRUE,    # Perform hierarchical clustering on rows (metabolites).
  cluster_cols = TRUE,    # Perform hierarchical clustering on columns (samples).
  color = color_palette,  # Apply the defined color palette.
  main = "Heatmap of Scaled Data",  # Title of the heatmap.
  show_rownames = FALSE,  # Optional: Hide row names (metabolites) to improve visualization.
  show_colnames = TRUE    # Optional: Show column names (sample names).
)


pdf_file <- "Heatmap_Scaled_Data.pdf"
jpeg_file <- "Heatmap_Scaled_Data.jpeg"

# Save as PDF
pdf(file = pdf_file, width = 10, height = 8)  # Adjust width and height for desired size
pheatmap(
  combined_heatmap_data,
  cluster_rows = TRUE,    # Perform hierarchical clustering on rows (metabolites).
  cluster_cols = TRUE,    # Perform hierarchical clustering on columns (samples).
  color = color_palette,  # Apply the defined color palette.
  main = "Heatmap of Scaled Data",  # Title of the heatmap.
  show_rownames = FALSE,  # Hide row names (metabolites).
  show_colnames = TRUE    # Show column names (sample names).
)
dev.off()  # Close the PDF device

# Save as JPEG
jpeg(file = jpeg_file, width = 1200, height = 1000, res = 150)  # High resolution
pheatmap(
  combined_heatmap_data,
  cluster_rows = TRUE,    # Perform hierarchical clustering on rows (metabolites).
  cluster_cols = TRUE,    # Perform hierarchical clustering on columns (samples).
  color = color_palette,  # Apply the defined color palette.
  main = "Heatmap of Scaled Data",  # Title of the heatmap.
  show_rownames = FALSE,  # Hide row names (metabolites).
  show_colnames = TRUE    # Show column names (sample names).
)
dev.off()  # Close the JPEG device

# Print confirmation
cat("Heatmap saved as", pdf_file, "and", jpeg_file, "in the working directory.\n")

################################################################################
# Summarize multiple copies of metabolites coming from different libraries

# Step 1: Create a new column 'Metabolite.name.mod' by extracting text before the first ';'
df_combined_scaled_v2 <- df_combined_scaled %>%
  mutate(Metabolite.name.mod = str_extract(Metabolite.name, "^[^;]+"))

# Step 2: Group by 'Metabolite.name.mod' and calculate the average for all sample columns
df_combined_scaled_v2_averaged <- df_combined_scaled_v2 %>%
  group_by(Metabolite.name.mod) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  ungroup()


## Heat Map

# Extract numeric data for heatmap
heatmap_data_v2 <- df_combined_scaled_v2_averaged %>% select(-Metabolite.name.mod, -Average.Rt.min.) %>% as.matrix()

color_palette <- colorRampPalette(c( "#A6CEE3" , "#33A02C", "#E31A1C"))(2)

pdf_file <- "Heatmap_combined_scaled_v2_averaged.pdf"
jpeg_file <- "Heatmap_combined_scaled_v2_averaged.jpeg"

# Save as PDF
pdf(file = pdf_file, width = 10, height = 8)  # Adjust width and height for desired size
pheatmap(
  heatmap_data_v2,
  cluster_rows = TRUE,    
  cluster_cols = TRUE,    
  color = color_palette,  
  main = "Heatmap of Scaled Data",  
  show_rownames = FALSE,  
  show_colnames = TRUE)

dev.off()  

# Save as JPEG
jpeg(file = jpeg_file, width = 1200, height = 1000, res = 150)  # High resolution
pheatmap(
  heatmap_data_v2,
  cluster_rows = TRUE,    
  cluster_cols = TRUE,    
  color = color_palette,  
  main = "Heatmap of Scaled Data",  
  show_rownames = FALSE,  
  show_colnames = TRUE)

dev.off()  

# Print confirmation
cat("Heatmap saved as", pdf_file, "and", jpeg_file, "in the working directory.\n")


################################################################################
# Find Top metabolites

# Summarize multiple copies of metabolites coming from different libraries

# Step 1: Create a new column 'Metabolite.name.mod' by extracting text before the first ';'
df_combined_v2 <- df_combined %>%
  mutate(Metabolite.name.mod = str_extract(Metabolite.name, "^[^;]+"))

# Step 2: Group by 'Metabolite.name.mod' and calculate the average for all sample columns
df_combined_v2_averaged <- df_combined_v2 %>%
  group_by(Metabolite.name.mod) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  ungroup()


# Calculate the Average Intensity Across All Samples

# Exclude metadata columns
metabolite_data <- df_combined_v2_averaged %>% 
  select(-Metabolite.name.mod, -Average.Rt.min.)


df_combined_v2_averaged <- df_combined_v2_averaged %>% 
  mutate(Median_Intensity = apply(metabolite_data, 1, median, na.rm = TRUE))

# Sort and select the top metabolites by median intensity
highest_metabolites <- df_combined_v2_averaged %>%
  arrange(desc(Median_Intensity)) %>%
  select(Metabolite.name.mod, Median_Intensity)

# Print the top metabolites by median intensity
print(highest_metabolites)


write.csv(df_combined_v2_averaged, "df_combined_v2_averaged.csv")


################################################################################
#  find the top metabolites for each group

df_combined_v3 <- df_combined %>%
  mutate(Metabolite.name.mod = str_extract(Metabolite.name, "^[^;]+"))

# Create a new DataFrame by dropping the specified columns and calculating the mean of rows grouped by "Metabolite.name.mod"
df_combined_v3_Met_Avg <- df_combined_v3 %>%
  select(-Metabolite.name, -Adduct.type, -Average.Rt.min.) %>% # Drop the specified columns
  group_by(Metabolite.name.mod) %>% # Group by "Metabolite.name.mod"
  summarise(across(everything(), mean, na.rm = TRUE)) %>% # Calculate the mean for each sample column
  ungroup() # Remove grouping to get a clean DataFrame

# View the resulting DataFrame
head(df_combined_v3_Met_Avg)

# Create a new data frame with the mean of each group
df_combined_v3_sample_mean <- df_combined_v3_Met_Avg

# Extract group names from sample columns
sample_columns <- colnames(df_combined_v3_sample_mean)[-1] # Exclude the "Metabolite.name.mod" column
group_names <- sapply(sample_columns, function(x) sub("_\\d+$", "", x)) # Remove replicate numbers (_n) to get group names

# Calculate the mean for each group and add it as a new column
for (group in unique(group_names)) {
  # Get the columns corresponding to the current group
  group_columns <- sample_columns[grepl(paste0("^", group, "_"), sample_columns)]
  
  # Calculate the mean of the group columns and add a new column to the data frame
  df_combined_v3_sample_mean[[group]] <- rowMeans(
    df_combined_v3_sample_mean[group_columns],
    na.rm = TRUE
  )
}

# View the resulting data frame
head(df_combined_v3_sample_mean)


################################################################################
## Find the top 25 metabolites for each group


# 1. Find the top 25 metabolites for each group
group_columns <- colnames(df_combined_v3_sample_mean)[30:38] # Extract group column names

top_metabolites <- lapply(group_columns, function(group) {
  df_combined_v3_sample_mean %>%
    arrange(desc(.data[[group]])) %>% # Sort by the group column in descending order
    slice_head(n = 25) %>% # Get the top 25 rows
    pull(Metabolite.name.mod) # Extract the metabolite names
})

# Name the list elements by group
names(top_metabolites) <- group_columns

# 2. Create a unique list of all metabolites across all groups
top_metabolite_list <- unique(unlist(top_metabolites))

# 3. Create a new dataframe for these metabolites
# Filter the original sample data for the selected metabolites
df_top_metabolites_replicates <- df_combined_v3_Met_Avg %>%
  filter(Metabolite.name.mod %in% top_metabolite_list)

# View the resulting dataframes
print("Top metabolites for each group:")
print(top_metabolites)

print("Final dataframe showcasing replicates behavior:")
head(df_top_metabolites_replicates)


################################################################################



# Define the group column names
group_columns <- colnames(df_combined_v3_sample_mean)[30:38]

# Create a list of top 25 metabolites for each group
top_metabolites <- lapply(group_columns, function(group) {
  df_combined_v3_sample_mean %>%
    arrange(desc(.data[[group]])) %>% # Sort by the group column in descending order
    slice_head(n = 50) %>% # Get the top 25 rows
    pull(Metabolite.name.mod) # Extract the metabolite names
})

# Name the list elements by group
names(top_metabolites) <- group_columns

# Combine the list into a presence/absence dataframe
top_metabolites_df <- tibble(Metabolite.name.mod = unique(unlist(top_metabolites)))

# Create presence/absence columns for each group
for (group in group_columns) {
  top_metabolites_df[[group]] <- ifelse(
    top_metabolites_df$Metabolite.name.mod %in% top_metabolites[[group]],
    1,
    0
  )
}

# UpSet Plot using ComplexUpset
ComplexUpset::upset(
  top_metabolites_df,
  group_columns,
  name = "Top Metabolites by Group",
  base_annotations = list(
    'Intersection Size' = intersection_size(
      counts = TRUE,
      text = list(size = 4)
    )
  ),
  set_sizes = upset_set_size(
    geom = geom_bar(
      width = 0.5, # Removed aes(fill = ...) to fix the issue
      fill = "steelblue" # You can set a static fill color if needed
    )
  ),
  width_ratio = 0.3
)

################################################################################
# Find metabolites in each intersection
intersection_list <- top_metabolites_df %>%
  pivot_longer(cols = group_columns, names_to = "Group", values_to = "Presence") %>%
  group_by(Metabolite.name.mod) %>%
  summarise(Intersection = paste(Group[Presence == 1], collapse = ", ")) %>%
  arrange(Intersection)

# View the intersections
print(intersection_list)

# Find metabolites shared between specific groups
shared_metabolites <- intersection_list %>%
  filter(str_detect(Intersection, "HyPea-7404") & str_detect(Intersection, "HyPep-4601N"))

# Find metabolites only in HyPea-7404
only_in_HyPea_7404 <- top_metabolites_df %>%
  filter(`HyPea-7404` == 1 & rowSums(select(., -Metabolite.name.mod, -`HyPea-7404`)) == 0)

# View the result
print(only_in_HyPea_7404)

# Find metabolites only in HyPep-7504
only_in_HyPep_7504 <- top_metabolites_df %>%
  filter(`HyPep-7504` == 1 & rowSums(select(., -Metabolite.name.mod, -`HyPep-7504`)) == 0)

# View the result
print(only_in_HyPep_7504)

# Find metabolites shared between HyPea-7404, HyPep-4601N, and HY-YEST-503
specific_intersection <- intersection_list %>%
  filter(Intersection == "HyPea-7404, HyPep-4601N, HY-YEST-503")

# View specific intersection
print(specific_intersection)

################################################################################
# After you find metabolites of interest (From Top50, UpSetPlots etc...), for other statistical analysis use "df_combined_POS_NEG.csv" file (only normalized, without scaling); 


