#---- Libraries Loading ----

library(tidyverse)
library(pheatmap)
library(vsn)
library(patchwork)
library(naniar)
library(simputation)
library(missForest)



#---- Data Loading ----

df_pos <- read.csv("Data/POS_Height_1_2024_09_24_05_19_14.csv", stringsAsFactors = FALSE)


#---- Run Quality Check ----

df_run <- df_pos %>%
  filter(Metabolite.name != "Unknown",
         !grepl("no MS2:", Metabolite.name),
         MS.MS.matched != "FALSE")

colnames(df_run)


# Creating a new dataframe with the specified columns

df_run_selected <- df_run %>%
  select(Average.Rt.min., Metabolite.name, 
         starts_with("X240830_1137_hydrolysates_blank"),
         starts_with("X240830_1137_hydrolysates_pool"),
         starts_with("sample") & ends_with("_sb.cn.pos_dda1"))

# Verify the selected column names

colnames(df_run_selected)


#Rename Samples
# Create a lookup table for renaming columns

sample_names <- read.csv("Data/Sample_Names.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(sample_names) <- c("Original", "Renamed")  # Assign meaningful column names

# Create a lookup table for renaming

rename_lookup <- setNames(sample_names$Renamed, sample_names$Original)


# Rename columns in df_run_selected using the lookup table

df_run_selected <- df_run_selected %>%
  rename_with(~ rename_lookup[gsub("_sb.cn.pos_dda1", "", .x)], .cols = contains("sample"))

# Explicitly rename blank and pool columns, if needed

df_run_selected <- df_run_selected %>%
  rename(
    Blank_1 = X240830_1137_hydrolysates_blank_sb.cn.pos_dda1,
    Blank_2 = X240830_1137_hydrolysates_blank_sb.cn.pos_dda2,
    Blank_3 = X240830_1137_hydrolysates_blank_sb.cn.pos_dda3,
    Blank_4 = X240830_1137_hydrolysates_blank_sb.cn.pos_dda4,
    Pool_1 = X240830_1137_hydrolysates_pool_sb.cn.pos_dda1,
    Pool_2 = X240830_1137_hydrolysates_pool_sb.cn.pos_dda2
  )

# Trim any leading or trailing spaces in column names

colnames(df_run_selected) <- trimws(colnames(df_run_selected))

# Verify renamed column names

colnames(df_run_selected)
write.csv(df_run_selected, "df_run_selected.csv" ) 



# Converting selected columns to long format

df_run_selected_long <- df_run_selected %>%
  pivot_longer(
    cols = -c(Average.Rt.min., Metabolite.name), # Exclude metadata columns
    names_to = "Samples",
    values_to = "Area"
  )


### Raw Data distribution ----

# Create the box plot

ggplot(df_run_selected_long, aes(x = Samples, y = Area)) +
  geom_boxplot() +
  geom_jitter(width=0.15)+
  theme(axis.text.x = element_text(angle = 90))+  
  labs(x = "Samples", y = "Area", title = "Box Plot of Area by Sample")

# Create the box plot with log-transformed data

ggplot(df_run_selected_long, aes(x = Samples, y = log(Area))) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.5) +  # Added alpha for better visualization of points
  theme(axis.text.x = element_text(angle = 90))+  
  labs(x = "Samples", y = "Log10(Area)", title = "Box Plot of Log-Transformed Area by Sample")


### PCA ----
# 1. Prepare data for PCA
# Select only the columns representing samples and exclude metadata columns

df_pca_data <- df_run_selected %>%
  select(-Average.Rt.min., -Metabolite.name) %>%
  t() %>%
  as.data.frame()

# Standardize the data (mean=0, variance=1) for PCA

df_pca_data_scaled <- scale(df_pca_data)

# 2. Perform PCA

pca_result <- prcomp(df_pca_data_scaled, center = TRUE, scale. = TRUE)

# 3. Extract sample names and group by ignoring replicate numbers

sample_names <- rownames(df_pca_data)
group_names <- str_extract(sample_names, "^[^_]+")  # Extract the part before "_"

# 4. Create a data frame with PCA results and group labels

pca_data <- data.frame(Sample = sample_names, 
                       Group = group_names,
                       PC1 = pca_result$x[,1], 
                       PC2 = pca_result$x[,2])

# 5. Plot the PCA results with colors based on sample groups

ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 3) +
  labs(title = "PCA Plot of Samples by Group", x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal() +
  theme(legend.position = "right")



################################
### Missing values ----
################################


#### Missing vlues % ----

# Calculate missing value percentage for each sample
missing_percentage <- df_run_selected %>%
  select(-Average.Rt.min., -Metabolite.name) %>%  # Exclude metadata columns
  summarise(across(everything(), ~ mean(is.na(.) | . == 0) * 100)) %>%  # Calculate percentage of missing values
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Missing_Percentage")  # Convert to long format

# Plot missing value percentage as a bar chart

ggplot(missing_percentage, aes(x = Sample, y = Missing_Percentage)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Missing Value Percentage by Sample",
       x = "Sample Names",
       y = "Missing Value Percentage (%)")


################ ----
# Step 1: Identify sample columns and separate group names from replicate numbers

# Define the sample columns (exclude metadata columns like "Average.Rt.min." and "Metabolite.name")

sample_mv_columns <- names(df_run_selected)[!(names(df_run_selected) %in% c("Average.Rt.min.", "Metabolite.name", "Blank_1", "Blank_2", "Blank_3", "Blank_4", "Pool_1", "Pool_2"))]

# Convert the data to long format to make it easier to work with groups and replicates

df_mv_handle_long <- df_run_selected %>%
  pivot_longer(cols = all_of(sample_mv_columns), 
               names_to = "Sample", 
               values_to = "Area") %>%
  mutate(Group = sub("_\\d+$", "", Sample),       # Extract the base sample group name (everything before the last underscore and digit)
         Replicate = as.numeric(sub(".*_", "", Sample))  # Extract replicate number (the last number after the underscore)
  )

# Step 2: Count the non-missing values in each group for each metabolite

# Group by metabolite and group, and count non-missing values in each group

df_group_counts <- df_mv_handle_long %>%
  group_by(Metabolite.name, Group) %>%
  summarise(NonMissingCount = sum(!is.na(Area)), .groups = "drop")

# Step 3: Filter metabolites that meet the threshold (at least 2 non-missing values in any group)

# Find metabolites with at least 2 replicates detected in any group

metabolites_to_keep <- df_group_counts %>%
  filter(NonMissingCount >= 2) %>%
  pull(Metabolite.name) %>%
  unique()  # Unique list of metabolites that meet the criteria

# Step 4: Filter the original df_run_selected based on the metabolites_to_keep list

df_filtered <- df_run_selected %>%
  filter(Metabolite.name %in% metabolites_to_keep)


# Step 5: Count the number of identified metabolites (non-missing values) per sample
# Select sample columns and count non-NA values per column

metabolite_counts <- df_filtered %>%
  select(-Average.Rt.min., -Metabolite.name) %>%
  summarise(across(everything(), ~ sum(!is.na(.)))) %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Count")

# Step 6: Plot a bar plot for the number of identified metabolites per sample

ggplot(metabolite_counts, aes(x = Sample, y = Count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Number of Identified Metabolites per Sample", x = "Sample", y = "Count of Identified Metabolites")



# missing value visualization

# Step 1: Create a binary matrix where 1 represents missing values (0 in this case) and 0 represents present values

missing_matrix <- df_filtered %>%
  select(-Average.Rt.min., -Metabolite.name) %>%
  mutate(across(everything(), ~ ifelse(. == 0, 1, 0))) %>%  # Mark 0 as missing (1), others as present (0)
  as.matrix()

# Double-check the binary matrix for any missing values represented as 1
table(missing_matrix)


# Define color and breaks explicitly for binary data
color_palette <- c("white", "black")
breaks <- c(-0.5, 0.5, 1.5)

# Plot the heatmap
pheatmap(missing_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = color_palette,
         breaks = breaks,
         show_rownames = FALSE,
         show_colnames = TRUE,
         main = "Heatmap of Missing Values Across Metabolites and Samples")


# Missing value identification using library(naniar) 

# Convert 0 values to NA in the dataset temporarily for missing data analysis

df_temp <- df_filtered %>%
  mutate(across(everything(), ~ ifelse(. == 0, NA, .)))

vis_miss(df_temp, cluster = TRUE)  # Visualize clustered missing values

# Upset plot for missing data combinations
gg_miss_upset(df_temp)



# Group-wise missingness summary
missing_summary <- df_temp %>%
  pivot_longer(cols = -c(Average.Rt.min., Metabolite.name), names_to = "Sample", values_to = "Value") %>%
  mutate(Type = case_when(
    grepl("Blank", Sample) ~ "Blank",
    grepl("Pool", Sample) ~ "Pool",
    TRUE ~ "Sample"
  )) %>%
  group_by(Type, Sample) %>%
  summarise(Missing_Percentage = mean(is.na(Value)) * 100, .groups = "drop")

# Plot missingness percentages by groups
ggplot(missing_summary, aes(x = Sample, y = Missing_Percentage, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Missingness by Group and Sample",
       x = "Sample Names",
       y = "Missing Value Percentage (%)")



# Count the occurrences of each type of "missing" value in df_run_selected
missing_summary <- df_run_selected %>%
  summarise(
    Total_Zero = sum(across(where(is.numeric), ~ sum(. == 0, na.rm = TRUE))),
    Total_NA = sum(is.na(.)),
    Total_NaN = sum(across(where(is.numeric), ~ sum(is.nan(.))))
  )

# Display the summary of missing values
print(missing_summary)

# "0 " indicates a true biological absence (e.g., a metabolite is genuinely undetectable in a sample), these values should remain as 0.
# Keeping 0 preserves the biological meaning of "no signal" or "no presence."
# Hence downstream analysis should consider 0 as a meaningful observation rather than missing data.


###########################################################
### Let's try this withour Missing value imputation ----
# DDA; these are random observation; not everything was quantified
##########################################################

# If blanks have high intensities for certain metabolites, subtracting them first avoids propagating those effects through normalization.
# When using scaling techniques (e.g., Pareto scaling or mean centering) for multivariate analysis (e.g., PCA), it is usually better to subtract blanks first.



# Step 1: Calculate Average Blank Values

# Calculate the average of blank intensities
df_with_blank_avg <- df_filtered %>%
  rowwise() %>%
  mutate(Average_Blank = mean(c_across(starts_with("Blank_")), na.rm = TRUE))


# Step 2: Subtract Blank Values from Samples
# Define the sample columns (exclude metadata and blank columns)
sample_columns <- setdiff(names(df_with_blank_avg), 
                          c("Metabolite.name", "Average.Rt.min.", 
                            "Blank_1", "Blank_2", "Blank_3", "Blank_4", "Pool_1", "Pool_2", "Average_Blank"))

## Subtract blanks and replace negatives with NA or 0
## Negative values can be replaced with NA to treat them as missing data for downstream analyses or with zero if they signify an absence of signal rather than an artifact.

# # Subtract blanks and replace negatives with NA
# df_blank_corrected <- df_with_blank_avg %>%
#   mutate(across(all_of(sample_columns), ~ ifelse((. - Average_Blank) < 0, NA, . - Average_Blank))) %>%
#   select(-starts_with("Blank_"), -starts_with("Pool_"), -Average_Blank)  # Drop Blank, Pool, and Average_Blank columns

# Subtract blanks and replace negatives with zero
df_blank_corrected <- df_with_blank_avg %>%
  mutate(across(all_of(sample_columns), ~ ifelse((. - Average_Blank) < 0, 0, . - Average_Blank))) %>%
  select(-starts_with("Blank_"), -starts_with("Pool_"), -Average_Blank)  # Drop Blank, Pool, and Average_Blank columns

################################
# Check for any negative values across numeric columns
any_negative <- any(df_blank_corrected %>% 
                      select(where(is.numeric)) %>% 
                      unlist() < 0, na.rm = TRUE)

# Print result
if (any_negative) {
  print("There are negative values in the dataframe.")
} else {
  print("No negative values found in the dataframe.")
}
#####################################

# Identify metabolites with >20% zero values in samples
metabolites_to_keep <- df_blank_corrected %>%
  select(-Metabolite.name, -Average.Rt.min.) %>%
  summarise(across(everything(), ~ mean(. == 0, na.rm = TRUE))) %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Zero_Percentage") %>%
  filter(Zero_Percentage <= 0.20) %>%
  pull(Sample)

# Filter the dataset to keep only those metabolites
df_filtered_80 <- df_blank_corrected %>%
  select(Metabolite.name, Average.Rt.min., all_of(metabolites_to_keep))


####################################
## Visualize missing values
# Prepare the dataset for heatmap
# Remove metadata columns (Metabolite.name, Average.Rt.min.)
heatmap_data <- df_filtered_80 %>%
  select(-Metabolite.name, -Average.Rt.min.)

# Create a binary-like matrix for coloring: 0 stays as 0, all other values become 1
binary_heatmap_data <- heatmap_data %>%
  mutate(across(everything(), ~ ifelse(. == 0, 0, 1))) %>%
  as.matrix()

# Define color palette and breaks for binary data
color_palette <- c("black", "lightblue") # Black for 0, light blue for all other values
breaks <- c(-0.5, 0.5, 1.5) # Ensure two distinct bins: 0 and non-zero values

# Plot the heatmap
pheatmap(
  binary_heatmap_data,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = color_palette,
  breaks = breaks,
  main = "Two-Color Heatmap of Zero Patterns in df_filtered_80",
  show_rownames = FALSE,
  show_colnames = TRUE
)



################################################

# Many feature extraction software programs, such as MZmine 3, often generate tables with missing values denoted as ‘NA’, ‘NaN’ or 0. This means that for several m/z and RT traces in a given sample, there may not be a peak detected and therefore no value is available. However, many statistical approaches, such as PCA, require numerical values for each observation. Hence, these features with missing values need to be removed or imputed. In this section, we handle the zero values in our blank-removed feature quantification table. https://www.nature.com/articles/s41596-024-01046-3

# To handle log transformation when dataset contains zeros, we can use one of the following approaches to avoid undefined values (since log ⁡ ( 0 ) log(0) is undefined):

# Approach 1: Add a Small Constant Before Log Transformation
# Add a small constant (𝑘) to all values, ensuring no zeros remain.
# k is typically chosen as a fraction of the smallest non-zero value in the dataset to minimize distortion.

# # Identify the smallest non-zero value in the dataset
# min_nonzero <- df_filtered_80 %>%
#   select(-Metabolite.name, -Average.Rt.min.) %>%
#   unlist() %>%
#   .[. > 0] %>%
#   min(na.rm = TRUE)
# 
# # Add a small constant (e.g., 10% of the smallest non-zero value)
# k <- min_nonzero * 0.1
# 
# # Apply log transformation with the constant
# df_log_transformed <- df_blank_corrected %>%
#   mutate(across(where(is.numeric), ~ log(. + k)))
# 
# # 2. Histogram of Log-Transformed Values
# 
# # Convert data to long format for easier visualization
# df_long_log <- df_log_transformed %>%
#   pivot_longer(cols = -c(Metabolite.name, Average.Rt.min.), names_to = "Sample", values_to = "Log_Area")
# 
# # Plot histogram of log-transformed values
# ggplot(df_long_log, aes(x = Log_Area)) +
#   geom_histogram(bins = 50, fill = "skyblue", color = "black") +
#   labs(title = "Distribution of Log-Transformed Metabolite Areas",
#        x = "Log(Area + k)", y = "Frequency") +
#   theme_minimal()



# Extract numeric data (excluding metadata columns)
data_matrix <- as.matrix(df_filtered_80 %>% select(-Metabolite.name, -Average.Rt.min.))

# Apply VSN normalization
vsn_fit <- vsn2(data_matrix)

# Extract normalized data
normalized_data <- as.data.frame(exprs(vsn_fit))

# Add metadata columns back
df_normalized <- cbind(
  df_filtered_80 %>% select(Metabolite.name, Average.Rt.min.),
  normalized_data
)

# Extract normalized data for imputation
numeric_data <- df_normalized %>% select(-Metabolite.name, -Average.Rt.min.)

# Perform imputation using missForest
imputed_result <- missForest(as.data.frame(numeric_data))

# Extract the imputed data
imputed_data <- as.data.frame(imputed_result$ximp)

# Combine imputed data with metadata
df_imputed <- cbind(
  df_normalized %>% select(Metabolite.name, Average.Rt.min.),
  imputed_data
)

################################
# visualization
# Prepare data for boxplot / Log2
# Prepare data for plotting (before normalization)
df_long_before <- df_filtered_80 %>%
  pivot_longer(cols = -c(Metabolite.name, Average.Rt.min.), names_to = "Sample", values_to = "Intensity") %>%
  mutate(Stage = "Before Normalization", Log2_Intensity = log2(Intensity + 1))  # Add log2 transformation

# Prepare data for plotting (after normalization)
df_long_after_norm <- df_normalized %>%negati
  pivot_longer(cols = -c(Metabolite.name, Average.Rt.min.), names_to = "Sample", values_to = "Intensity") %>%
  mutate(Stage = "After Normalization", Log2_Intensity = log2(Intensity + 1))  # Add log2 transformation

# Prepare data for plotting (after imputation)
df_long_after_impute <- df_imputed %>%
  pivot_longer(cols = -c(Metabolite.name, Average.Rt.min.), names_to = "Sample", values_to = "Intensity") %>%
  mutate(Stage = "After Imputation", Log2_Intensity = log2(Intensity + 1))  # Add log2 transformation



# Boxplot for before normalization
plot_before <- ggplot(df_long_before, aes(x = Sample, y = Log2_Intensity)) +
  geom_boxplot(outlier.alpha = 0.5, fill = "lightblue") +
  labs(title = "Before Normalization", x = "Samples", y = "Log2(Intensity)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)  # Rotate x-axis labels
  )

# Boxplot for after normalization
plot_after_norm <- ggplot(df_long_after_norm, aes(x = Sample, y = Log2_Intensity)) +
  geom_boxplot(outlier.alpha = 0.5, fill = "lightgreen") +
  labs(title = "After Normalization", x = "Samples", y = "Log2(Intensity)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)  # Rotate x-axis labels
  )

# Boxplot for after imputation
plot_after_impute <- ggplot(df_long_after_impute, aes(x = Sample, y = Log2_Intensity)) +
  geom_boxplot(outlier.alpha = 0.5, fill = "lightcoral") +
  labs(title = "After Imputation", x = "Samples", y = "Log2(Intensity)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)  # Rotate x-axis labels
  )

# Combine the three plots using patchwork
combined_plot <- plot_before + plot_after_norm + plot_after_impute +
  plot_layout(ncol = 1) & theme(legend.position = "none")

# Display the combined plot
print(combined_plot)


##############################
## PCA

# Prepare numeric data for PCA (exclude metadata)
pca_data <- df_imputed %>%
  select(-Metabolite.name, -Average.Rt.min.) %>%
  t() %>%
  as.data.frame()

# Extract group names by splitting sample names at "_"
rownames(pca_data) <- colnames(df_imputed %>% select(-Metabolite.name, -Average.Rt.min.))
pca_data$Group <- sapply(rownames(pca_data), function(x) strsplit(x, "_")[[1]][1])

# Standardize the data
pca_data_scaled <- scale(pca_data[, -ncol(pca_data)])  # Exclude 'Group' column

# Perform PCA
pca_result <- prcomp(pca_data_scaled, center = TRUE, scale. = TRUE)

# Create a data frame for plotting PCA results
pca_plot_data <- as.data.frame(pca_result$x[, 1:2])  # Extract PC1 and PC2
pca_plot_data$Group <- pca_data$Group  # Add group information
pca_plot_data$Sample <- rownames(pca_data)  # Add sample names


library(viridis)

# Dynamically create a larger palette using viridis
num_groups <- length(unique(pca_plot_data$Group))  # Count unique groups
darker_colors <- viridis(num_groups, option = "D")  # Generate dark colors using the "D" viridis option

ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text(vjust = -1, size = 3) +
  scale_color_manual(values = darker_colors) +  # Use dynamically generated darker colors
  labs(
    title = "PCA Plot with Grouped Samples (Dynamic Dark Colors)",
    x = "Principal Component 1",
    y = "Principal Component 2",
    color = "Group"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )


#######################################################
# Summarize multiple copies of metabolites coming from different libraries

# Step 1: Create a new column 'Metabolite.name.mod' by extracting text before the first ';'
df_imputed <- df_imputed %>%
  mutate(Metabolite.name.mod = str_extract(Metabolite.name, "^[^;]+"))

# Step 2: Group by 'Metabolite.name.mod' and calculate the average for all sample columns
df_averaged <- df_imputed %>%
  group_by(Metabolite.name.mod) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  ungroup()

# Step 2: Group by 'Metabolite.name.mod' and calculate the average for all sample columns
df_averaged <- df_imputed %>%
  group_by(Metabolite.name.mod) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE))) %>%
  ungroup()




# Find Top metabolites


# Step 1: Find the maximum value across all sample columns for each metabolite
df_highest <- df_averaged %>%
  rowwise() %>%
  mutate(Max_Intensity = max(c_across(where(is.numeric)), na.rm = TRUE)) %>%
  ungroup()

# Step 2: Arrange metabolites by the maximum intensity value in descending order
df_highest_sorted <- df_highest %>%
  arrange(desc(Max_Intensity))

# Step 3: Select the top metabolites (e.g., top 10)
top_metabolites <- df_highest_sorted %>%
  select(Metabolite.name.mod, Max_Intensity) %>%
  head(10)

# Print the top metabolites
print(top_metabolites)


# Exclude metadata columns and calculate row-wise maximum intensity
metabolite_max <- df_imputed %>%
  rowwise() %>%
  mutate(Max_Intensity = max(c_across(-c(Metabolite.name, Average.Rt.min.)), na.rm = TRUE)) %>%
  ungroup()

# Select the top 10 metabolites with the highest intensities
top_metabolites <- metabolite_max %>%
  arrange(desc(Max_Intensity)) %>%
  select(Metabolite.name, Max_Intensity) %>%
  head(10)

# Print the top metabolites
print(top_metabolites)

# Optional: Visualize the top metabolites
ggplot(top_metabolites, aes(x = reorder(Metabolite.name, -Max_Intensity), y = Max_Intensity)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Top 10 Metabolites by Maximum Intensity",
       x = "Metabolite Name",
       y = "Max Intensity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))














