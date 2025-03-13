logit_transform <- function(x, min_val = 0.703, max_val = 0.78) {
  # If min and max aren't provided, use the data min/max with a small buffer
  if(is.null(min_val)) min_val <- min(x, na.rm = TRUE) - 0.001
  if(is.null(max_val)) max_val <- max(x, na.rm = TRUE) + 0.001

  # Scale to [0,1]
  x_scaled <- (x - min_val) / (max_val - min_val)

  # Handle boundary cases to avoid Inf/-Inf
  x_scaled[x_scaled <= 0] <- 0.001
  x_scaled[x_scaled >= 1] <- 0.999

  # Apply logit
  log(x_scaled / (1 - x_scaled))
}

# Inverse logit for back-transformation
inverse_logit_transform <- function(x, min_val = 0.703, max_val = 0.78) {
  p <- 1 / (1 + exp(-x))
  p * (max_val - min_val) + min_val
}

