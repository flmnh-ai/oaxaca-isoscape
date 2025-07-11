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

eck4 <- sf::st_crs('PROJCRS["Eckert_IV",
        BASEGEOGCRS["WGS 84",
                    DATUM["unknown",
                          ELLIPSOID["WGS84",6378137,298.257223563,
                                    LENGTHUNIT["metre",1,
                                               ID["EPSG",9001]]]],
                    PRIMEM["Greenwich",0,
                           ANGLEUNIT["Degree",0.0174532925199433]]],
        CONVERSION["unnamed",
                   METHOD["Eckert IV"],
                   PARAMETER["Longitude of natural origin",0,
                             ANGLEUNIT["Degree",0.0174532925199433],
                             ID["EPSG",8802]],
                   PARAMETER["False easting",0,
                             LENGTHUNIT["metre",1],
                             ID["EPSG",8806]],
                   PARAMETER["False northing",0,
                             LENGTHUNIT["metre",1],
                             ID["EPSG",8807]]],
        CS[Cartesian,2],
        AXIS["(E)",east,
             ORDER[1],
             LENGTHUNIT["metre",1,
                        ID["EPSG",9001]]],
        AXIS["(N)",north,
             ORDER[2],
             LENGTHUNIT["metre",1,
                        ID["EPSG",9001]]]]')

convergence <- function(workflow){
  bart_fit_engine <- extract_fit_engine(workflow)
  plot(bart_fit_engine)

  print(t(bart_fit_engine$sigma) |>
    as_tibble() |>
    pivot_longer(everything()) |>
    ggplot(aes(value, color = name)) +
    geom_density())

  print(t(bart_fit_engine$sigma) |>
    as_tibble() |>
    mutate(iter = 1:n()) |>
    pivot_longer(-iter) |>
    ggplot(aes(iter, value, color = name)) +
    geom_line() +
    facet_wrap(~name))

  sigma_mcmc <- t(bart_fit_engine$sigma) |>
    as_tibble() |>
    map(mcmc) |>
    mcmc.list()
  plot(sigma_mcmc)

  # Check Gelman-Rubin statistic
  print(gelman.diag(sigma_mcmc))

  # Effective sample size
  effectiveSize(sigma_mcmc)
}

varimp <- function(model) {
  model <- extract_fit_engine(model)
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


#' Create a prediction map for strontium isotope values
#'
#' @param predictions Stars object containing Sr predictions
#' @param region Optional sf object defining region boundaries
#' @param samples Optional sf object with sample points to overlay
#' @param title Plot title
#' @return ggplot object with the map
plot_isoscape <- function(predictions, region = NULL, samples = NULL, title = "Strontium Isoscape") {
  p <- ggplot()

  if (!is.null(region)) {
    p <- p + geom_sf(data = region, fill = 'black', color = NA)
  }

  p <- p +
    geom_stars(data = predictions |> setNames('Sr')) +
    scale_fill_viridis_c(option = 'magma', na.value = NA, name = "87Sr/86Sr") +
    theme_bw() +
    theme(
      panel.background = element_rect(fill = "lightblue"),
      panel.grid.major = element_line(color = "white", linewidth = 0.5),
      plot.margin = margin(0, 0, 0, 0)
    ) +
    labs(title = title, x = '', y = '') +
    ggspatial::annotation_scale(location = "bl", width_hint = 0.2) +
    ggspatial::annotation_north_arrow(
      location = "bl", which_north = "true",
      pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"),
      style = ggspatial::north_arrow_fancy_orienteering
    )

  if (!is.null(region)) {
    p <- p + geom_sf(data = region, fill = NA, color = 'white', linewidth = 0.5)
  }

  if (!is.null(samples)) {
    p <- p + geom_sf(data = samples, aes(color = Sr), size = 2) +
      scale_color_viridis_c(option = 'magma', name = "Measured 87Sr/86Sr")
  }

  return(p)
}
