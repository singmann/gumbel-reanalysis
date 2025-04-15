get_max_rhat <- function(brms_model) {
  # Ensure posterior package is available
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("The 'posterior' package is required but not installed.")
  }
  
  # Extract R-hat values
  rhat_values <- posterior::rhat(brms_model)
  
  # Return the maximum R-hat
  max_rhat <- max(rhat_values, na.rm = TRUE)
  return(max_rhat)
}

