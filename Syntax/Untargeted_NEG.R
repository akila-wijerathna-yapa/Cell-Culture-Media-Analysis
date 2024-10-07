library(tidyverse)
library(ggplot2)

#---- Data Loading ----
df_neg <- read.csv("Data/NEG_Height_1_2024_09_25_03_48_22.csv", stringsAsFactors = FALSE)

#---- Sample QC for running ----

df_run <- df_neg %>%
  filter(Metabolite.name != "Unknown",
         !grepl("no MS2:", Metabolite.name),
         MS.MS.matched != "FALSE")

colnames(df_run)

# Creating a new dataframe with the specified columns

df_run_selected <- df_run %>%
  select(Average.Rt.min., Metabolite.name, 
         X240830_1137_hydrolysates_blank_sb.cn.neg_dda1,
         X240830_1137_hydrolysates_pool_sb.cn.neg_dda1, 
         X240830_1137_hydrolysates_pool_sb.cn.neg_dda2, 
         sample1_sb.cn.neg_dda1, sample10_sb.cn.neg_dda1, sample11_sb.cn.neg_dda1,
         sample12_sb.cn.neg_dda1, sample13_sb.cn.neg_dda1, sample14_sb.cn.neg_dda1, 
         sample15_sb.cn.neg_dda1, sample16_sb.cn.neg_dda1, sample17_sb.cn.neg_dda1,
         sample18_sb.cn.neg_dda1, sample19_sb.cn.neg_dda1, sample2_sb.cn.neg_dda1,
         sample20_sb.cn.neg_dda1, sample21_sb.cn.neg_dda1, sample22_sb.cn.neg_dda1,
         sample23_sb.cn.neg_dda1, sample24_sb.cn.neg_dda1, sample25_sb.cn.neg_dda1,
         sample26_sb.cn.neg_dda1, sample27_sb.cn.neg_dda1, sample28_sb.cn.neg_dda1,
         sample3_sb.cn.neg_dda1, sample4_sb.cn.neg_dda1, sample5_sb.cn.neg_dda1,
         sample6_sb.cn.neg_dda1, sample7_sb.cn.neg_dda1, sample8_sb.cn.neg_dda1,
         sample9_sb.cn.neg_dda1
  )


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
  rename_with(~ make_unique_names(products[gsub("_sb.cn.neg_dda1", "", .x)]), 
              .cols = contains("sample"))

# Rename the "blank" and "pool" columns explicitly
df_run_selected <- df_run_selected %>%
  rename(
    Blank = X240830_1137_hydrolysates_blank_sb.cn.neg_dda1,
    Pool_1 = X240830_1137_hydrolysates_pool_sb.cn.neg_dda1,
    Pool_2 = X240830_1137_hydrolysates_pool_sb.cn.neg_dda2
  )

# Inspect the renamed dataset
head(df_run_selected)

colnames(df_run_selected)

# Remove leading and trailing spaces from column names
colnames(df_run_selected) <- trimws(colnames(df_run_selected))


# Convert df_run_selected to long format
df_run_selected_long <- df_run_selected %>%
  pivot_longer(
    cols = c("Blank", "Pool_1", "Pool_2", "HyPea-7404_1", "HyPep-4601N_3", "HyPep-7504_1", 
             "HyPep-7504_2", "HyPep-7504_3", "HyPep-YE_1", "HyPep-YE_2", "HyPep-YE_3", 
             "HY-YEST-412_1", "HY-YEST-412_2", "HY-YEST-412_3", "HyPea-7404_2", "HY-YEST-466_1", 
             "HY-YEST-466_2", "HY-YEST-466_3", "HY-YEST-503_1", "HY-YEST-503_2", "HY-YEST-503_3", 
             "HY-YEST-555_1", "HY-YEST-555_2", "HY-YEST-555_3", "HyPea-7404_3", "HyPep-1510_1", 
             "HyPep-1510_2", "HyPep-1510_3", "HyPep-1510_4", "HyPep-4601N_1", "HyPep-4601N_2"),
    names_to = "Samples",
    values_to = "Area"
  )


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



#----
# Filter Data

df_clean <- df_neg %>%
  filter(Metabolite.name != "Unknown",
           !grepl("no MS2:", Metabolite.name),
           MS.MS.matched != "FALSE")

# Inspect the cleaned data
head(df_clean)

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
df_clean <- df_clean %>%
  rename_with(~ make_unique_names(products[gsub("_sb.cn.neg_dda1", "", .x)]), 
              .cols = contains("sample"))

# Rename the "blank" and "pool" columns explicitly
df_clean <- df_clean %>%
  rename(
    Blank = X240830_1137_hydrolysates_blank_sb.cn.neg_dda1,
    Pool_1 = X240830_1137_hydrolysates_pool_sb.cn.neg_dda1,
    Pool_2 = X240830_1137_hydrolysates_pool_sb.cn.neg_dda2
  )

# Inspect the renamed dataset
head(df_clean)



# Check unique values in the 'Metabolite.name' column

unique_metabolites <- unique(df_clean$Metabolite)

# Get the count of unique values
num_unique_metabolites <- length(unique_metabolites)

# Print the count
print(num_unique_metabolites)

