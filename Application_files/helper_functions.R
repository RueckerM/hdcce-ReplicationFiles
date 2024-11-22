# Function to individual specifically demean panel data variables 
# and report back original, individual means and the demeaned data 
#------------------------------------------------------------------
demean_panel_variables <- function(data, id_column, variables) {
  # Calculate individual-specific means
  individual_means <- aggregate(
    data[variables], 
    by = list(id = data[[id_column]]), 
    FUN = mean, 
    na.rm = TRUE
  )
  
  # Merge individual means back to original dataset
  demeaned_data <- merge(data, individual_means, by.x = id_column, by.y = "id", suffixes = c("", "_mean"))
  
  # Subtract individual means from original variables
  for (var in variables) {
    demeaned_data[[paste0(var, "_demeaned")]] <- 
      demeaned_data[[var]] - demeaned_data[[paste0(var, "_mean")]]
  }
  
  return(demeaned_data)
}


# Function to remove digits from the exact zeroes.
#------------------------------------------------------------------
round_custom <- function(x , digits=5) {
  # Check if input is a matrix
  if (!is.matrix(x)) {
    stop("Input must be a matrix")
  }
  
  # Create a matrix to store the rounded results
  # Using the same dimensions as the input matrix
  rounded <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  
  # Round each element of the matrix
  for (i in 1:nrow(x)) {
    for (j in 1:ncol(x)) {
      if (x[i, j] == 0) {
        rounded[i, j] <- "-"
      } else {
        rounded[i, j] <- round(x[i, j], digits = digits)
      }
    }
  }
  
  return(rounded)
}
