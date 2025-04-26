# Load required packages
library(clustvis)
library(pheatmap)
library(ggplot2)
library(tidyverse)

# Read dataset
df_combined <- read.csv("data/Metab_merged.csv", stringsAsFactors = FALSE)

# Set row names to the metabolite column and remove it from the dataframe
rownames(df_combined) <- df_combined$Metabolite
df_combined <- df_combined %>% select(-Metabolite)

# Import data into ClustVis
data_imported <- importData(
  file = "data/Metab_merged.csv",
  sep = ",",
  nbrRowAnnos = 1, # First column is the annotation (Metabolite names)
  nbrColAnnos = 0, # No column annotations
  transpose = FALSE
)

data_processed <- processData(
  data_imported,
  transformation = "ln(x + 1)", # Log transformation
  rowScaling = "uv", # Unit variance scaling
  pcaMethod = "svdImpute", # Singular Value Decomposition with imputation
  maxNaRows = 0.2, # Remove rows with more than 20% missing values
  maxNaCols = 0.2  # Remove columns with more than 20% missing values
)

pca_plot <- generatePCA(
  proc = data_processed,
  pcx = 1, pcy = 2, # Select PC1 and PC2
  colorAnno = NULL, # No color annotations
  showEllipses = TRUE,
  showSampleIds = TRUE
)

# Save the PCA plot
savePCA(pca_plot, file = "PCA_plot.pdf")


heatmap_plot <- generateHeatmap(
  proc = data_processed,
  showImputed = TRUE, # Show imputed values
  clustDistRows = "correlation", # Cluster rows using correlation distance
  clustMethodRows = "average", # Use average linkage for hierarchical clustering
  clustDistCols = "correlation",
  clustMethodCols = "average",
  matrixColorScheme = "RdBu", # Red-Blue color scheme
  showNumbers = FALSE,
  showRownames = TRUE,
  showColnames = TRUE
)

# Save the heatmap plot
saveHeatmap(heatmap_plot, file = "Heatmap_plot.pdf")

