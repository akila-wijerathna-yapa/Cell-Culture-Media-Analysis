library(tidyverse)
library(ggplot2)

df_pos <- read.csv("Data/POS_Height_1_2024_09_24_05_19_14.csv", stringsAsFactors = FALSE)



# Filter Data

df_clean <- df_pos %>%
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
  rename_with(~ make_unique_names(products[gsub("_sb.cn.pos_dda1", "", .x)]), 
              .cols = contains("sample"))

# Rename the "blank" and "pool" columns explicitly
df_clean <- df_clean %>%
  rename(
    Blank_1 = X240830_1137_hydrolysates_blank_sb.cn.pos_dda1,
    Blank_2 = X240830_1137_hydrolysates_blank_sb.cn.pos_dda2,
    Blank_3 = X240830_1137_hydrolysates_blank_sb.cn.pos_dda3,
    Blank_4 = X240830_1137_hydrolysates_blank_sb.cn.pos_dda4,
    Pool_1 = X240830_1137_hydrolysates_pool_sb.cn.pos_dda1,
    Pool_2 = X240830_1137_hydrolysates_pool_sb.cn.pos_dda2
  )

# Inspect the renamed dataset
head(df_clean)



# Check unique values in the 'Metabolite.name' column

unique_metabolites <- unique(df_clean$Metabolite)

# Get the count of unique values
num_unique_metabolites <- length(unique_metabolites)

# Print the count
print(num_unique_metabolites)

