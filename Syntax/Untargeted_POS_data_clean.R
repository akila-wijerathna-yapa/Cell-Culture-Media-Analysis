
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


### Data Normalization ----

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
