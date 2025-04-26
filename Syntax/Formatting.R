# Load necessary packages
library(tidyverse)

# Read the CSV file
df <- read.csv("Data/Metab.csv")

# Convert to wide format
df_wide <- df %>%
  pivot_wider(
    id_cols = c(Metabolite, Category),  # rows defined by Metabolite and Category
    names_from = Hyd,                   # new column names from Hyd
    values_from = avg                   # values come from avg column
  )

# View the wide format table
head(df_wide)
write.csv(df_wide, "Data/Metabolomics_WideFormat.csv", row.names = FALSE)
