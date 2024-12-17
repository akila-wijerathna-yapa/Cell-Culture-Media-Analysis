# Check for any negative values across numeric columns
any_negative <- any(df_filtered_80  %>% 
                      select(where(is.numeric)) %>% 
                      unlist() < 0, na.rm = TRUE)

# Print result
if (any_negative) {
  print("There are negative values in the dataframe.")
} else {
  print("No negative values found in the dataframe.")
}
