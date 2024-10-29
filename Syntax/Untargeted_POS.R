#---- Libraries Loading ----

library(tidyverse)
library(pheatmap)


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


##########################################
### Filter out potential contaminants ----
##########################################


# Step 1: Define blank columns
blank_columns <- c("Blank_1", "Blank_2", "Blank_3", "Blank_4")

# Step 2: Calculate Q3 thresholds for each blank sample
blank_thresholds <- df_run_selected %>%
  select(all_of(blank_columns)) %>%
  summarise(across(everything(), ~ quantile(.x, 0.75, na.rm = TRUE)))

# Step 3: Transform the data to long format and flag specific data points above Q3 in blanks
df_long <- df_run_selected %>%
  pivot_longer(cols = -c(Average.Rt.min., Metabolite.name),
               names_to = "Samples", values_to = "Area") %>%
  mutate(flag = case_when(
    Samples == "Blank_1" & Area > blank_thresholds$Blank_1 ~ "Above Q3",
    Samples == "Blank_2" & Area > blank_thresholds$Blank_2 ~ "Above Q3",
    Samples == "Blank_3" & Area > blank_thresholds$Blank_3 ~ "Above Q3",
    Samples == "Blank_4" & Area > blank_thresholds$Blank_4 ~ "Above Q3",
    TRUE ~ "Below Q3"
  ))

# Step 4: Plot the box plot with log10 scale, highlighting only specific points above Q3 in red
ggplot(df_long, aes(x = Samples, y = log10(Area))) +
  geom_boxplot(outlier.shape = NA) +  # Single box plot per sample
  geom_jitter(aes(color = flag), width = 0.15, alpha = 0.5) +  # Highlight only Above Q3 points in red
  scale_color_manual(values = c("Above Q3" = "red", "Below Q3" = "black")) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Samples", y = "Log10(Area)", title = "Box Plot of Log10(Area) by Sample with Flagged Metabolites")




# Q3 red in Blank and corresponding metabolites in blue


# Step 1: Define blank columns
blank_columns <- c("Blank_1", "Blank_2", "Blank_3", "Blank_4")

# Step 2: Calculate Q3 thresholds for each blank sample
blank_thresholds <- df_run_selected %>%
  select(all_of(blank_columns)) %>%
  summarise(across(everything(), ~ quantile(.x, 0.75, na.rm = TRUE)))

# Step 3: Flag specific data points above Q3 in blanks, and identify metabolites for blue highlighting
df_long <- df_run_selected %>%
  pivot_longer(cols = -c(Average.Rt.min., Metabolite.name),
               names_to = "Samples", values_to = "Area") %>%
  mutate(flag = case_when(
    Samples == "Blank_1" & Area > blank_thresholds$Blank_1 ~ "Above Q3",  # Red in blanks
    Samples == "Blank_2" & Area > blank_thresholds$Blank_2 ~ "Above Q3",
    Samples == "Blank_3" & Area > blank_thresholds$Blank_3 ~ "Above Q3",
    Samples == "Blank_4" & Area > blank_thresholds$Blank_4 ~ "Above Q3",
    Metabolite.name %in% (df_run_selected %>%
                            filter((Blank_1 > blank_thresholds$Blank_1) |
                                     (Blank_2 > blank_thresholds$Blank_2) |
                                     (Blank_3 > blank_thresholds$Blank_3) |
                                     (Blank_4 > blank_thresholds$Blank_4)) %>%
                            pull(Metabolite.name)) ~ "Above Q3 in Blanks",  # Blue in all samples for these metabolites
    TRUE ~ "Below Q3"  # Default color for others
  ))

# Step 4: Plot the box plot with log10 scale, using red for Q3 points in blanks and blue across all samples
ggplot(df_long, aes(x = Samples, y = log10(Area))) +
  geom_boxplot(outlier.shape = NA) +  # Single box plot per sample
  geom_jitter(aes(color = flag), width = 0.15, alpha = 0.5) +  # Highlight with colors based on flags
  scale_color_manual(values = c("Above Q3" = "red", "Above Q3 in Blanks" = "blue", "Below Q3" = "black")) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Samples", y = "Log10(Area)", title = "Box Plot of Log10(Area) by Sample with Highlighted Metabolites")


## Select Below Q3 for downstream data analysis


# Step 1: Define blank columns
blank_columns <- c("Blank_1", "Blank_2", "Blank_3", "Blank_4")

# Step 2: Calculate Q3 thresholds for each blank sample
blank_thresholds <- df_run_selected %>%
  select(all_of(blank_columns)) %>%
  summarise(across(everything(), ~ quantile(.x, 0.75, na.rm = TRUE)))

# Step 3: Identify metabolites that are below Q3 in all blank samples
below_Q3_metabolites <- df_run_selected %>%
  rowwise() %>%
  filter(all(c_across(all_of(blank_columns)) <= blank_thresholds)) %>%
  pull(Metabolite.name)

# Step 4: Filter the original data frame to keep only "Below Q3" metabolites
df_below_Q3 <- df_run_selected %>%
  filter(Metabolite.name %in% below_Q3_metabolites)


## Plotting

# Step 5: Convert the filtered data frame to long format
df_below_Q3_long <- df_below_Q3 %>%
  pivot_longer(cols = -c(Average.Rt.min., Metabolite.name),
               names_to = "Samples", values_to = "Area")

# Step 6: Plot the long-format data in log10 scale
ggplot(df_below_Q3_long, aes(x = Samples, y = log10(Area))) +
  geom_boxplot(outlier.shape = NA) +  # Box plot for each sample without outliers
  geom_jitter(width = 0.15, alpha = 0.5, color = "blue") +  # Jitter points for visibility
  theme(axis.text.x = element_text(angle = 90)) +  
  labs(x = "Samples", y = "Log10(Area)", title = "Box Plot of Log10(Area) by Sample for Below Q3 Metabolites")



# New PCA for clean data

# Step 1: Prepare data for PCA by removing non-numeric columns and transposing the data
df_pca_data <- df_below_Q3 %>%
  select(-Average.Rt.min., -Metabolite.name) %>%
  t() %>%
  as.data.frame()

# Step 2: Standardize the data for PCA
df_pca_data_scaled <- scale(df_pca_data)

# Step 3: Perform PCA
pca_result <- prcomp(df_pca_data_scaled, center = TRUE, scale. = TRUE)

# Step 4: Extract sample names and group them by ignoring replicate numbers
sample_names <- rownames(df_pca_data)
group_names <- str_extract(sample_names, "^[^_]+")  # Extract everything before the final underscore

# Step 5: Create a data frame with PCA results and grouping labels
pca_data <- data.frame(Sample = sample_names, 
                       Group = group_names,
                       PC1 = pca_result$x[,1], 
                       PC2 = pca_result$x[,2])

# Step 6: Plot the PCA results with colors based on sample groups
ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 3) +
  labs(title = "PCA Plot of Samples by Group (Below Q3 Metabolites)", 
       x = "Principal Component 1", 
       y = "Principal Component 2") +
  theme_minimal() +
  theme(legend.position = "right")






# Remove Blank sample values from samples ----

# Step 1: Calculate the average of the blank columns for each metabolite
df_with_blank_avg <- df_below_Q3 %>%
  rowwise() %>%
  mutate(Average_Blank = mean(c_across(starts_with("Blank_")), na.rm = TRUE))

# Step 2: Define sample columns, excluding `Pool_1` and `Pool_2`
sample_columns <- setdiff(names(df_below_Q3), c("Average.Rt.min.", "Metabolite.name", "Blank_1", "Blank_2", "Blank_3", "Blank_4", "Pool_1", "Pool_2"))

# Step 3: Subtract the average blank value from each sample column and create the new data frame
df_clean_to_stat <- df_with_blank_avg %>%
  mutate(across(all_of(sample_columns), ~ . - Average_Blank)) %>%
  select(Metabolite.name, all_of(sample_columns))

# View the resulting dataframe
str(df_clean_to_stat)


## Summarize replicate Metabolite.name`  

# Step 1: Clean the `Metabolite.name` by removing extra information after ";"
df_clean_to_stat <- df_clean_to_stat %>%
  mutate(Metabolite.name = sub(";.*", "", Metabolite.name))

# Step 2: Group by the cleaned `Metabolite.name` and calculate the average for each unique metabolite
df_clean_to_stat_v1 <- df_clean_to_stat %>%
  group_by(Metabolite.name) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  ungroup()

# View the resulting dataframe
str(df_clean_to_stat_v1)


# Summarize sample replicates


# Step 1: Reshape to long format to handle replicate averaging without creating duplicate column names
df_long <- df_clean_to_stat_v1 %>%
  pivot_longer(cols = -Metabolite.name, names_to = "Sample", values_to = "Value")

# Step 2: Extract base sample names by removing replicate suffixes
df_long <- df_long %>%
  mutate(Base_Sample = str_remove(Sample, "_\\d+$"))

# Step 3: Group by Metabolite.name and Base_Sample, then calculate the average for each base sample
df_clean_to_stat_v2 <- df_long %>%
  group_by(Metabolite.name, Base_Sample) %>%
  summarise(Average_Value = mean(Value, na.rm = TRUE), .groups = "drop")  # Calculate mean across replicates

# Step 4: Pivot back to wide format to get a final data frame with averaged values for each base sample
df_clean_to_stat_v2 <- df_clean_to_stat_v2 %>%
  pivot_wider(names_from = Base_Sample, values_from = Average_Value)

# View the resulting dataframe
str(df_clean_to_stat_v2)



# Heatmap

# Step 1: Truncate Metabolite.name to the first 30 characters
df_clean_to_stat_v2 <- df_clean_to_stat_v2 %>%
  mutate(Metabolite.name = str_sub(Metabolite.name, 1, 30))

# Step 2: Prepare the data matrix for the heatmap
# Remove the Metabolite.name column and convert the data to a matrix
heatmap_matrix <- df_clean_to_stat_v2 %>%
  column_to_rownames(var = "Metabolite.name") %>%
  as.matrix()

# Step 3a: Save heatmap as a PDF with high resolution
pdf("Heatmap_Metabolite_Data.pdf", width = 12, height = 10)  # Set appropriate width and height
pheatmap(heatmap_matrix,
         scale = "row",  # Normalize each row for better visualization
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Heatmap of Metabolite Data (Averaged Across Replicates)",
         fontsize_row = 6)  # Adjust font size for row labels
dev.off()

# Step 3b: Save heatmap as a high-resolution JPG
jpeg("Heatmap_Metabolite_Data.jpg", width = 1200, height = 1000, res = 300)  # Set high resolution with res = 300
pheatmap(heatmap_matrix,
         scale = "row",  # Normalize each row for better visualization
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Heatmap of Metabolite Data (Averaged Across Replicates)",
         fontsize_row = 6)  # Adjust font size for row labels
dev.off()

