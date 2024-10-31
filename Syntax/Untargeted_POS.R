#---- Libraries Loading ----

library(tidyverse)
library(pheatmap)
library(vsn)
library(patchwork)


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
         X240830_1137_hydrolysates_blank_sb.cn.pos_dda1,
         X240830_1137_hydrolysates_blank_sb.cn.pos_dda2,
         X240830_1137_hydrolysates_blank_sb.cn.pos_dda3,
         X240830_1137_hydrolysates_blank_sb.cn.pos_dda4,
         X240830_1137_hydrolysates_pool_sb.cn.pos_dda1,
         X240830_1137_hydrolysates_pool_sb.cn.pos_dda2,
         sample1_sb.cn.pos_dda1, sample10_sb.cn.pos_dda1, sample11_sb.cn.pos_dda1,
         sample12_sb.cn.pos_dda1, sample13_sb.cn.pos_dda1, sample14_sb.cn.pos_dda1,
         sample15_sb.cn.pos_dda1, sample16_sb.cn.pos_dda1, sample17_sb.cn.pos_dda1,
         sample18_sb.cn.pos_dda1, sample19_sb.cn.pos_dda1, sample2_sb.cn.pos_dda1,
         sample20_sb.cn.pos_dda1, sample21_sb.cn.pos_dda1, sample22_sb.cn.pos_dda1,
         sample23_sb.cn.pos_dda1, sample24_sb.cn.pos_dda1, sample25_sb.cn.pos_dda1,
         sample26_sb.cn.pos_dda1, sample27_sb.cn.pos_dda1, sample28_sb.cn.pos_dda1,
         sample3_sb.cn.pos_dda1, sample4_sb.cn.pos_dda1, sample5_sb.cn.pos_dda1,
         sample6_sb.cn.pos_dda1, sample7_sb.cn.pos_dda1, sample8_sb.cn.pos_dda1,
         sample9_sb.cn.pos_dda1)

colnames(df_run_selected)

#Rename Samples
# Create a lookup table for renaming columns
products <- c(
  "sample14" = " HyPep-YE_1",
  "sample15" = " HyPep-YE_2",
  "sample16" = " HyPep-YE_3",
  "sample1" = " HyPea-7404_1",
  "sample2" = " HyPea-7404_2",
  "sample3" = " HyPea-7404_3",
  "sample17" = " HY-YEST-412_1",
  "sample18" = " HY-YEST-412_2",
  "sample19" = " HY-YEST-412_3",
  "sample26" = " HY-YEST-555_1",
  "sample27" = " HY-YEST-555_2",
  "sample28" = " HY-YEST-555_3",
  "sample20" = " HY-YEST-466_1",
  "sample21" = " HY-YEST-466_2",
  "sample22" = " HY-YEST-466_3",
  "sample23" = " HY-YEST-503_1",
  "sample24" = " HY-YEST-503_2",
  "sample25" = " HY-YEST-503_3",
  "sample8" = " HyPep-4601N_1",
  "sample9" = " HyPep-4601N_2",
  "sample10" = " HyPep-4601N_3",
  "sample4" = " HyPep-1510_1",
  "sample5" = " HyPep-1510_2",
  "sample6" = " HyPep-1510_3",
  "sample7" = " HyPep-1510_4",
  "sample11" = " HyPep-7504_1",
  "sample12" = " HyPep-7504_2",
  "sample13" = " HyPep-7504_3")



# Function to ensure unique names by appending suffixes to duplicates
make_unique_names <- function(names_vec) {
  make.unique(names_vec, sep = "_")
}

# Renaming columns using gsub to extract the relevant part and match with the lookup table
df_run_selected <- df_run_selected %>%
  rename_with(~ make_unique_names(products[gsub("_sb.cn.pos_dda1", "", .x)]), 
              .cols = contains("sample"))

# Rename the "blank" and "pool" columns explicitly
df_run_selected<- df_run_selected %>%
  rename(
    Blank_1 = X240830_1137_hydrolysates_blank_sb.cn.pos_dda1,
    Blank_2 = X240830_1137_hydrolysates_blank_sb.cn.pos_dda2,
    Blank_3 = X240830_1137_hydrolysates_blank_sb.cn.pos_dda3,
    Blank_4 = X240830_1137_hydrolysates_blank_sb.cn.pos_dda4,
    Pool_1 = X240830_1137_hydrolysates_pool_sb.cn.pos_dda1,
    Pool_2 = X240830_1137_hydrolysates_pool_sb.cn.pos_dda2
  )

colnames(df_run_selected)

# Remove leading and trailing spaces from column names
colnames(df_run_selected) <- trimws(colnames(df_run_selected))


# Converting selected columns to long format
df_run_selected_long <- df_run_selected %>%
  pivot_longer(cols = c("Blank_1", 
                        "Blank_2", 
                        "Blank_3", 
                        "Blank_4", 
                        "Pool_1", 
                        "Pool_2", 
                        "HyPea-7404_1", 
                        "HyPep-4601N_3", 
                        "HyPep-7504_1", 
                        "HyPep-7504_2", 
                        "HyPep-7504_3", 
                        "HyPep-YE_1", 
                        "HyPep-YE_2", 
                        "HyPep-YE_3", 
                        "HY-YEST-412_1", 
                        "HY-YEST-412_2", 
                        "HY-YEST-412_3", 
                        "HyPea-7404_2", 
                        "HY-YEST-466_1", 
                        "HY-YEST-466_2", 
                        "HY-YEST-466_3", 
                        "HY-YEST-503_1", 
                        "HY-YEST-503_2", 
                        "HY-YEST-503_3", 
                        "HY-YEST-555_1", 
                        "HY-YEST-555_2", 
                        "HY-YEST-555_3", 
                        "HyPea-7404_3", 
                        "HyPep-1510_1", 
                        "HyPep-1510_2", 
                        "HyPep-1510_3", 
                        "HyPep-1510_4", 
                        "HyPep-4601N_1", 
                        "HyPep-4601N_2"), 
                        names_to = "Samples", values_to = "Area")


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
### Missing value handling ----
################################

# Step 1: Identify sample columns and separate group names from replicate numbers

# Define the sample columns (exclude metadata columns like "Average.Rt.min." and "Metabolite.name")
sample_mv_columns <- names(df_run_selected)[!(names(df_run_selected) %in% c("Average.Rt.min.", "Metabolite.name", "Blank_1", "Blank_2", "Blank_3", "Blank_4", "Pool_1", "Pool_2"))]

# Convert the data to long format to make it easier to work with groups and replicates
df_mv_handle_long <- df_run_selected %>%
  pivot_longer(cols = all_of(sample_mv_columns), 
               names_to = "Sample", 
               values_to = "Area") %>%
  mutate(
    Group = sub("_\\d+$", "", Sample),       # Extract the base sample group name (everything before the last underscore and digit)
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





# Install and load naniar if not installed
install.packages("naniar")
library(naniar)

# Convert 0 values to NA in the dataset temporarily for missing data analysis
df_temp <- df_filtered %>%
  mutate(across(everything(), ~ ifelse(. == 0, NA, .)))

# Check for missing values
library(naniar)
vis_miss(df_temp, cluster = TRUE)  # Visualize clustered missing values

# Upset plot for missing data combinations
gg_miss_upset(df_temp)



# Exclude non-numeric columns for pivoting
df_temp_long <- df_temp %>%
  select(where(is.numeric)) %>%  # Keep only numeric columns
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value")

# Plot histograms for all variables in one facet
ggplot(df_temp_long, aes(x = Value, fill = is.na(Value))) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("blue", "red"), labels = c("Non-Missing", "Missing")) +
  facet_wrap(~ Variable, scales = "free_x") +
  labs(title = "Histograms of Missing vs Non-Missing Values Across Variables",
       x = "Values",
       y = "Frequency",
       fill = "Data Status") +
  theme_minimal()


##################################
### Missing value imputation ----
#################################

### This data looks like missing not at random (MNAR). Data imputation ----

library(simputation)

# Step 1: Convert 0s to NAs in the dataset (if 0 represents missing values)
df_temp <- df_filtered %>%
  mutate(across(where(is.numeric), ~ ifelse(. == 0, NA, .)))  # Target only numeric columns

# Step 2: Confirm the presence of NAs after conversion
converted_na_count <- sum(is.na(df_temp))
print(paste("NA count after conversion:", converted_na_count))

# Step 3: Reapply k-NN imputation specifically to only the columns with remaining NA values
df_filtered_imputed <- impute_knn(df_temp, formula = ~ ., k = 10)

# Step 4: Verify that imputation has been completed
remaining_NA_count <- sum(is.na(df_filtered_imputed))  # Should be 0 if imputation is complete
print(paste("Remaining NA count after imputation:", remaining_NA_count))


# Verify the result
sum(is.na(df_filtered_imputed))

# Replace both 0s and NAs with 1 in numeric columns only
df_filtered_imputed <- df_filtered %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.) | . == 0, 1, .)))

# Verify if all NAs and 0s are replaced
sum(is.na(df_filtered_imputed))  # Should be 0 if imputation is complete


# Verify the result
sum(is.na(df_filtered_imputed))  # Should be 0 if imputation is complete



## now do donesstream after checking these imputation plots????

# Prepare data as before
df_temp_long <- df_temp %>%
  pivot_longer(cols = starts_with("Blank") | starts_with("Hy") | starts_with("Pool"), 
               names_to = "Sample", values_to = "Area") %>%
  mutate(Log2_Area = log2(Area + 1),  # Add 1 to avoid log of zero
         Source = "Before Imputation")

df_filtered_imputed_long <- df_filtered_imputed %>%
  pivot_longer(cols = starts_with("Blank") | starts_with("Hy") | starts_with("Pool"), 
               names_to = "Sample", values_to = "Area") %>%
  mutate(Log2_Area = log2(Area + 1),
         Source = "After Imputation")

# Combine both data frames for plotting
df_combined <- bind_rows(df_temp_long, df_filtered_imputed_long)

# Plot with faceting
ggplot(df_combined, aes(x = Log2_Area)) +
  geom_density(color = "blue", size = 1) +  # Set line color and thickness
  labs(title = "Density Plot of Log2 Area Before and After Imputation",
       x = "Log2(Area)",
       y = "Density") +
  facet_wrap(~ Source, ncol = 1) +  # Separate plots vertically by Source
  theme_minimal() +
  theme(strip.text = element_text(size = 12),  # Adjust facet title font size
        legend.position = "none")



###########################
### Data Normalization ----
###########################

# Load necessary libraries
library(vsn)
library(patchwork)

# Step 1: Prepare data for visualization (before normalization)
df_long_raw <- df_filtered_imputed %>%
  pivot_longer(cols = -c(Average_Rt_min_, Metabolite_name), names_to = "Sample", values_to = "Area")

# Step 2: Plot boxplot of raw data using log2 transformation
plot_raw <- ggplot(df_long_raw, aes(x = log2(Area + 1), y = Sample)) +  # log2 transformation for visualization
  geom_boxplot(fill = "lightblue") +
  theme(axis.text.y = element_text(angle = 0)) +
  labs(title = "Boxplot of Metabolite Areas (Before Normalization)", y = "Sample", x = "Log2(Area + 1)")

# Step 3: Apply VSN normalization on df_filtered_imputed
# Ensure the selection excludes Average_Rt_min_ and Metabolite_name columns
data_matrix <- as.matrix(df_filtered_imputed %>% select(-Average_Rt_min_, -Metabolite_name))
vsn_normalized <- vsn::vsn2(data_matrix)

# Extract normalized data
normalized_data <- as.data.frame(exprs(vsn_normalized))
normalized_data$Metabolite_name <- df_filtered_imputed$Metabolite_name  # Add back Metabolite names for plotting

# Step 4: Convert normalized data to long format for plotting
df_long_normalized <- normalized_data %>%
  pivot_longer(cols = -Metabolite_name, names_to = "Sample", values_to = "Area")

# Plot boxplot of normalized data
plot_normalized <- ggplot(df_long_normalized, aes(x = Area, y = Sample)) +
  geom_boxplot(fill = "lightgreen") +
  theme(axis.text.y = element_text(angle = 0)) +
  labs(title = "Boxplot of Metabolite Areas (After VSN Normalization)", y = "Sample", x = "VSN Normalized Area")

# Step 5: Combine plots with patchwork (top: raw data, bottom: normalized data)
combined_plot <- plot_raw / plot_normalized  # Stack the plots vertically
print(combined_plot)



####################################
#### plot and see imputed data -----
####################################


# Step 1: Select sample columns (excluding metadata) for PCA
df_pca_data <- df_filtered_imputed %>%
  select(-Average_Rt_min_, -Metabolite_name) %>%
  t() %>%
  as.data.frame()

# Standardize the data (mean=0, variance=1) for PCA
df_pca_data_scaled <- scale(df_pca_data)

# Step 2: Perform PCA
pca_result <- prcomp(df_pca_data_scaled, center = TRUE, scale. = TRUE)


# Step 3: Extract sample names and group by keeping the full prefix for "HY_YEST" samples
sample_names <- rownames(df_pca_data)
group_names <- ifelse(grepl("^HY_YEST", sample_names), 
                      str_extract(sample_names, "^[^_]+_[^_]+"),  # Extract "HY_YEST" as one unit
                      str_extract(sample_names, "^[^_]+"))       # Extract first part for others

# Step 4: Create a data frame with PCA results and group labels
pca_data <- data.frame(Sample = sample_names, 
                       Group = group_names,
                       PC1 = pca_result$x[,1], 
                       PC2 = pca_result$x[,2])

# Step 5: Plot the PCA results with color by group and labeled sample names
ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 3) +
  geom_text(aes(label = Sample), hjust = 1.1, vjust = 1.1, size = 2.5, check_overlap = TRUE) +  # Adjust text size and placement
  labs(title = "PCA Plot of Samples by Group", x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal() +
  theme(legend.position = "right")




### Boxplot with Log2 Transformation


# Step 6: Prepare the imputed data for visualization (before normalization)
df_long_raw <- df_filtered_imputed %>%
  pivot_longer(cols = -c(Average_Rt_min_, Metabolite_name), names_to = "Sample", values_to = "Area")

# Step 7: Plot boxplot of raw data using log2 transformation
plot_raw <- ggplot(df_long_raw, aes(x = log2(Area + 1), y = Sample)) +  # log2 transformation for visualization
  geom_boxplot(fill = "lightblue") +
  theme(axis.text.y = element_text(angle = 0)) +
  labs(title = "Boxplot of Metabolite Areas (Before Normalization)", y = "Sample", x = "Log2(Area + 1)")

# Step 8: Apply VSN normalization on df_filtered_imputed
data_matrix <- as.matrix(df_filtered_imputed %>% select(-Average_Rt_min_, -Metabolite_name))
vsn_normalized <- vsn::vsn2(data_matrix)

# Extract normalized data and convert it for plotting
normalized_data <- as.data.frame(exprs(vsn_normalized))
normalized_data$Metabolite_name <- df_filtered_imputed$Metabolite_name
df_long_normalized <- normalized_data %>%
  pivot_longer(cols = -Metabolite_name, names_to = "Sample", values_to = "Area")

# Step 9: Plot boxplot of normalized data
plot_normalized <- ggplot(df_long_normalized, aes(x = Area, y = Sample)) +
  geom_boxplot(fill = "lightgreen") +
  theme(axis.text.y = element_text(angle = 0)) +
  labs(title = "Boxplot of Metabolite Areas (After VSN Normalization)", y = "Sample", x = "VSN Normalized Area")

# Step 10: Combine plots with patchwork (top: raw data, bottom: normalized data)
combined_plot <- plot_raw / plot_normalized  # Stack the plots vertically
print(combined_plot)


### Box plot ----

# Convert normalized data to long format for ggplot2
df_long_normalized <- normalized_data %>%
  pivot_longer(cols = -Metabolite_name, names_to = "Sample", values_to = "Area")

# Plot the box plot with log-transformed data
ggplot(df_long_normalized, aes(x = Sample, y = log2(Area))) +
  geom_boxplot(fill = "lightgreen") +
  geom_jitter(width = 0.15, alpha = 0.5, color = "blue") +  # Optional: Color jitter for better contrast
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Samples", y = "Log2(Area)", title = "Box Plot of Log-Transformed Area by Sample (Normalized Data)")



##########################################
### Filter out potential contaminants ----
##########################################


# Step 1. Define blank columns

# Define blank columns
blank_columns <- c("Blank_1", "Blank_2", "Blank_3", "Blank_4")

# Calculate Q3 thresholds for each blank sample in normalized_data
blank_thresholds <- normalized_data %>%
  select(all_of(blank_columns)) %>%
  summarise(across(everything(), ~ quantile(.x, 0.75, na.rm = TRUE)))

# Step2. Transform Data to Long Format and Flag Points Above Q3 in Blanks

# Transform normalized_data to long format and flag values above Q3
df_normalized_data_long <- normalized_data %>%
  pivot_longer(cols = -Metabolite_name, names_to = "Samples", values_to = "Area") %>%
  mutate(flag = case_when(
    Samples == "Blank_1" & Area > blank_thresholds$Blank_1 ~ "Above Q3",
    Samples == "Blank_2" & Area > blank_thresholds$Blank_2 ~ "Above Q3",
    Samples == "Blank_3" & Area > blank_thresholds$Blank_3 ~ "Above Q3",
    Samples == "Blank_4" & Area > blank_thresholds$Blank_4 ~ "Above Q3",
    TRUE ~ "Below Q3"
  ))



# Step 3: Plot Box Plot with Highlighted Points in Long Format

ggplot(df_normalized_data_long, aes(x = Samples, y = log10(Area))) +
geom_boxplot(outlier.shape = NA) +  
  geom_jitter(aes(color = flag), width = 0.15, alpha = 0.5) +  
  scale_color_manual(values = c("Above Q3" = "red", "Below Q3" = "black")) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Samples", y = "Log10(Area)", title = "Box Plot of Log10(Area) by Sample with Flagged Metabolites")


# Step 4. Filter for Metabolites Below Q3 in All Blanks
# Identify metabolites below Q3 in all blanks
below_Q3_metabolites <- normalized_data %>%
  rowwise() %>%
  filter(all(c_across(all_of(blank_columns)) <= blank_thresholds)) %>%
  pull(Metabolite_name)

# Filter data to keep only below Q3 metabolites
df_below_Q3 <- normalized_data %>%
  filter(Metabolite_name %in% below_Q3_metabolites)


# Step 5. Create New PCA with df_below_Q3
# Prepare data for PCA by selecting sample columns
df_pca_data <- df_below_Q3 %>%
  select(-Metabolite_name) %>%
  t() %>%
  as.data.frame()

# Standardize and perform PCA
df_pca_data_scaled <- scale(df_pca_data)
pca_result <- prcomp(df_pca_data_scaled, center = TRUE, scale. = TRUE)

# Extract sample names and groups
sample_names <- rownames(df_pca_data)
group_names <- ifelse(grepl("^HY_YEST", sample_names), 
                      str_extract(sample_names, "^[^_]+_[^_]+"), 
                      str_extract(sample_names, "^[^_]+"))

# Create PCA plot with color based on group
pca_data <- data.frame(Sample = sample_names, Group = group_names, PC1 = pca_result$x[,1], PC2 = pca_result$x[,2])
ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 3) +
  labs(title = "PCA Plot of Samples by Group (Below Q3 Metabolites)", x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal() +
  theme(legend.position = "right")


# Step

# Step

# Step



# # Step 1: Define blank columns
# blank_columns <- c("Blank_1", "Blank_2", "Blank_3", "Blank_4")
# 
# # Step 2: Calculate Q3 thresholds for each blank sample
# blank_thresholds <- df_run_selected %>%
#   select(all_of(blank_columns)) %>%
#   summarise(across(everything(), ~ quantile(.x, 0.75, na.rm = TRUE)))
# 
# # Step 3: Transform the data to long format and flag specific data points above Q3 in blanks
# df_long <- df_run_selected %>%
#   pivot_longer(cols = -c(Average.Rt.min., Metabolite.name),
#                names_to = "Samples", values_to = "Area") %>%
#   mutate(flag = case_when(
#     Samples == "Blank_1" & Area > blank_thresholds$Blank_1 ~ "Above Q3",
#     Samples == "Blank_2" & Area > blank_thresholds$Blank_2 ~ "Above Q3",
#     Samples == "Blank_3" & Area > blank_thresholds$Blank_3 ~ "Above Q3",
#     Samples == "Blank_4" & Area > blank_thresholds$Blank_4 ~ "Above Q3",
#     TRUE ~ "Below Q3"
#   ))
# 
# # Step 4: Plot the box plot with log10 scale, highlighting only specific points above Q3 in red
# ggplot(df_long, aes(x = Samples, y = log10(Area))) +
#   geom_boxplot(outlier.shape = NA) +  # Single box plot per sample
#   geom_jitter(aes(color = flag), width = 0.15, alpha = 0.5) +  # Highlight only Above Q3 points in red
#   scale_color_manual(values = c("Above Q3" = "red", "Below Q3" = "black")) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   labs(x = "Samples", y = "Log10(Area)", title = "Box Plot of Log10(Area) by Sample with Flagged Metabolites")
# 
# 
# 
# 
# # Q3 red in Blank and corresponding metabolites in blue
# 
# 
# # Step 1: Define blank columns
# blank_columns <- c("Blank_1", "Blank_2", "Blank_3", "Blank_4")
# 
# # Step 2: Calculate Q3 thresholds for each blank sample
# blank_thresholds <- df_run_selected %>%
#   select(all_of(blank_columns)) %>%
#   summarise(across(everything(), ~ quantile(.x, 0.75, na.rm = TRUE)))
# 
# # Step 3: Flag specific data points above Q3 in blanks, and identify metabolites for blue highlighting
# df_long <- df_run_selected %>%
#   pivot_longer(cols = -c(Average.Rt.min., Metabolite.name),
#                names_to = "Samples", values_to = "Area") %>%
#   mutate(flag = case_when(
#     Samples == "Blank_1" & Area > blank_thresholds$Blank_1 ~ "Above Q3",  # Red in blanks
#     Samples == "Blank_2" & Area > blank_thresholds$Blank_2 ~ "Above Q3",
#     Samples == "Blank_3" & Area > blank_thresholds$Blank_3 ~ "Above Q3",
#     Samples == "Blank_4" & Area > blank_thresholds$Blank_4 ~ "Above Q3",
#     Metabolite.name %in% (df_run_selected %>%
#                             filter((Blank_1 > blank_thresholds$Blank_1) |
#                                      (Blank_2 > blank_thresholds$Blank_2) |
#                                      (Blank_3 > blank_thresholds$Blank_3) |
#                                      (Blank_4 > blank_thresholds$Blank_4)) %>%
#                             pull(Metabolite.name)) ~ "Above Q3 in Blanks",  # Blue in all samples for these metabolites
#     TRUE ~ "Below Q3"  # Default color for others
#   ))
# 
# # Step 4: Plot the box plot with log10 scale, using red for Q3 points in blanks and blue across all samples
# ggplot(df_long, aes(x = Samples, y = log10(Area))) +
#   geom_boxplot(outlier.shape = NA) +  # Single box plot per sample
#   geom_jitter(aes(color = flag), width = 0.15, alpha = 0.5) +  # Highlight with colors based on flags
#   scale_color_manual(values = c("Above Q3" = "red", "Above Q3 in Blanks" = "blue", "Below Q3" = "black")) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   labs(x = "Samples", y = "Log10(Area)", title = "Box Plot of Log10(Area) by Sample with Highlighted Metabolites")
# 
# 
# ## Select Below Q3 for downstream data analysis
# 
# 
# # Step 1: Define blank columns
# blank_columns <- c("Blank_1", "Blank_2", "Blank_3", "Blank_4")
# 
# # Step 2: Calculate Q3 thresholds for each blank sample
# blank_thresholds <- df_run_selected %>%
#   select(all_of(blank_columns)) %>%
#   summarise(across(everything(), ~ quantile(.x, 0.75, na.rm = TRUE)))
# 
# # Step 3: Identify metabolites that are below Q3 in all blank samples
# below_Q3_metabolites <- df_run_selected %>%
#   rowwise() %>%
#   filter(all(c_across(all_of(blank_columns)) <= blank_thresholds)) %>%
#   pull(Metabolite.name)
# 
# # Step 4: Filter the original data frame to keep only "Below Q3" metabolites
# df_below_Q3 <- df_run_selected %>%
#   filter(Metabolite.name %in% below_Q3_metabolites)
# 
# 
# ## Plotting
# 
# # Step 5: Convert the filtered data frame to long format
# df_below_Q3_long <- df_below_Q3 %>%
#   pivot_longer(cols = -c(Average.Rt.min., Metabolite.name),
#                names_to = "Samples", values_to = "Area")
# 
# # Step 6: Plot the long-format data in log10 scale
# ggplot(df_below_Q3_long, aes(x = Samples, y = log10(Area))) +
#   geom_boxplot(outlier.shape = NA) +  # Box plot for each sample without outliers
#   geom_jitter(width = 0.15, alpha = 0.5, color = "blue") +  # Jitter points for visibility
#   theme(axis.text.x = element_text(angle = 90)) +  
#   labs(x = "Samples", y = "Log10(Area)", title = "Box Plot of Log10(Area) by Sample for Below Q3 Metabolites")
# 
# 
# 
# # New PCA for clean data
# 
# # Step 1: Prepare data for PCA by removing non-numeric columns and transposing the data
# df_pca_data <- df_below_Q3 %>%
#   select(-Average.Rt.min., -Metabolite.name) %>%
#   t() %>%
#   as.data.frame()
# 
# # Step 2: Standardize the data for PCA
# df_pca_data_scaled <- scale(df_pca_data)
# 
# # Step 3: Perform PCA
# pca_result <- prcomp(df_pca_data_scaled, center = TRUE, scale. = TRUE)
# 
# # Step 4: Extract sample names and group them by ignoring replicate numbers
# sample_names <- rownames(df_pca_data)
# group_names <- str_extract(sample_names, "^[^_]+")  # Extract everything before the final underscore
# 
# # Step 5: Create a data frame with PCA results and grouping labels
# pca_data <- data.frame(Sample = sample_names, 
#                        Group = group_names,
#                        PC1 = pca_result$x[,1], 
#                        PC2 = pca_result$x[,2])
# 
# # Step 6: Plot the PCA results with colors based on sample groups
# ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
#   geom_point(size = 3) +
#   labs(title = "PCA Plot of Samples by Group (Below Q3 Metabolites)", 
#        x = "Principal Component 1", 
#        y = "Principal Component 2") +
#   theme_minimal() +
#   theme(legend.position = "right")
# 
# 
# 
# 

###############################################
### Clean the data for downstream analysis ----
###############################################

# Step 1: Calculate the average of the blank columns for each metabolite in df_below_Q3
df_with_blank_avg <- df_below_Q3 %>%
  rowwise() %>%
  mutate(Average_Blank = mean(c_across(starts_with("Blank_")), na.rm = TRUE))

# Step 2: Define sample columns, excluding `Pool_1`, `Pool_2`, and all `Blank` columns
sample_columns <- setdiff(names(df_below_Q3), c("Metabolite_name", "Blank_1", "Blank_2", "Blank_3", "Blank_4", "Pool_1", "Pool_2"))

# Step 3: Subtract the average blank value from each sample column to remove background noise
df_clean_to_stat <- df_with_blank_avg %>%
  mutate(across(all_of(sample_columns), ~ . - Average_Blank)) %>%
  select(Metabolite_name, all_of(sample_columns))  # Keep only the relevant columns for downstream analysis

# View the resulting dataframe
str(df_clean_to_stat)


#######################################
## Summarize replicate Metabolite.name`
#######################################

# Step 1: Clean `Metabolite.name` by removing extra information after ";"
df_clean_to_stat_v1 <- df_clean_to_stat %>%
  mutate(Metabolite_name = str_remove(Metabolite_name, ";.*"))

# Step 2: Group by the cleaned `Metabolite.name` and calculate the average for each unique metabolite
df_clean_to_stat_v1 <- df_clean_to_stat_v1 %>%
  group_by(Metabolite_name) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  ungroup()

# Step 3: Reshape to long format for averaging replicates without duplicate column names
df_clean_to_stat_v1_long <- df_clean_to_stat_v1 %>%
  pivot_longer(cols = -Metabolite_name, names_to = "Sample", values_to = "Value")

# Step 4: Extract base sample names by removing replicate suffixes
df_clean_to_stat_v1_long <- df_clean_to_stat_v1_long %>%
  mutate(Base_Sample = str_remove(Sample, "_\\d+$"))

# Step 5: Group by `Metabolite_name` and `Base_Sample`, then calculate the average for each base sample
df_clean_to_stat_v2 <- df_clean_to_stat_v1_long %>%
  group_by(Metabolite_name, Base_Sample) %>%
  summarise(Average_Value = mean(Value, na.rm = TRUE), .groups = "drop")

# Step 6: Pivot back to wide format to get a final data frame with averaged values for each base sample
df_clean_to_stat_v2 <- df_clean_to_stat_v2 %>%
  pivot_wider(names_from = Base_Sample, values_from = Average_Value)

# View the resulting dataframe structure
str(df_clean_to_stat_v2)

# Step 7: Prepare data for heatmap plotting
# Truncate Metabolite.name to the first 30 characters for readability
df_clean_to_stat_v2 <- df_clean_to_stat_v2 %>%
  mutate(Metabolite_name = str_sub(Metabolite_name, 1, 30))

# Convert to matrix format for heatmap
heatmap_matrix <- df_clean_to_stat_v2 %>%
  column_to_rownames(var = "Metabolite_name") %>%
  as.matrix()

# Step 8: Plot the heatmap
pheatmap(heatmap_matrix,
         scale = "row",  # Normalize each row for better visualization
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Heatmap of Metabolite Data (Averaged Across Replicates)",
         fontsize_row = 6)  # Adjust font size for readability




# Save heatmap as a PDF with high resolution
pdf("Heatmap_Metabolite_Data.pdf", width = 12, height = 10)  # Set appropriate width and height
pheatmap(heatmap_matrix,
         scale = "row",  # Normalize each row for better visualization
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Heatmap of Metabolite Data (Averaged Across Replicates)",
         fontsize_row = 6)  # Adjust font size for row labels
dev.off()  # Close the PDF device

# Save heatmap as a high-resolution JPEG
jpeg("Heatmap_Metabolite_Data.jpg", width = 1200, height = 1000, res = 300)  # Set high resolution with res = 300
pheatmap(heatmap_matrix,
         scale = "row",  # Normalize each row for better visualization
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Heatmap of Metabolite Data (Averaged Across Replicates)",
         fontsize_row = 6)  # Adjust font size for row labels
dev.off()  # Close the JPEG device












### Retrieve infromation for  below Q3 metabolites ----

# Step 1: Retrieve the list of metabolite names from df_below_Q3
metabolites_below_Q3 <- df_below_Q3$Metabolite.name

# Step 2: Filter df_pos to only include rows with metabolite names present in df_below_Q3
df_pos_filtered <- df_pos %>%
  filter(Metabolite.name %in% metabolites_below_Q3)

# View the resulting filtered data frame
print(df_pos_filtered)

