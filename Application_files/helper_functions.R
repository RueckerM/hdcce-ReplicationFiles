# Function to remove digits from the exact zeroes.

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
