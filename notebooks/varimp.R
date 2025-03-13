varimp <- function(model) {
  # Check dimensions of model$varcount
  if (length(dim(model$varcount)) == 3) {
    # For 3D array: Convert to list of matrices (one per chain)
    # Then combine into a single matrix
    chains <- dim(model$varcount)[1]
    counts <- do.call(rbind, lapply(1:chains, function(i) model$varcount[i,,]))

    # Preserve variable names if they exist
    if (!is.null(dimnames(model$varcount)[[3]])) {
      colnames(counts) <- dimnames(model$varcount)[[3]]
    }
  } else {
    # For 2D matrix, use as is
    counts <- model$varcount
  }

  # Calculate variable importances
  varimps <- colMeans(counts/rowSums(counts))

  # Create data frame using the names from varimps
  var.df <- data.frame(
    names = names(varimps),
    varimps = as.vector(varimps)
  )

  # Reorder factors for plotting
  var.df <- transform(var.df, names = reorder(names, -varimps))

  # Get relative importances for plotting
  rel <- counts/rowSums(counts)

  # Create and display plot
  p <- rel %>%
    data.frame() %>%
    gather() %>%
    group_by(key) %>%
    summarise(mean = mean(value), sd = sd(value, na.rm = TRUE)) %>%
    transform(Var = reorder(key, mean)) %>%
    ggplot(aes(x = Var, y = mean)) +
    geom_pointrange(aes(y = mean, x = Var, ymin = mean - sd, ymax = mean + sd),
                    color = "#00AFDD") +
    xlab(NULL) +
    ylab("Variable importance") +
    coord_flip() +
    theme_bw()

  print(p)

  # Return variable importance data frame
  return(var.df)
}
