utils::globalVariables("population")
#' Dasymetric Mapping
#'
#' Maps X-grid covariates to Y-grid using centroid-based spatial join
#'
#' @param misaligned_data List with gridx and gridy from simulate_misaligned_data
#'
#' @return sf object with Y grid containing mapped X covariates
#' @importFrom sf st_centroid st_within st_drop_geometry
#' @importFrom stats glm as.formula poisson binomial lm confint
#' @export
dasymetric_mapping <- function(misaligned_data) {
  # Extract components
  gridx <- misaligned_data$gridx
  gridy <- misaligned_data$gridy
  
  # Convert grids to centroids for gridx
  gridx_centroids <- st_centroid(gridx)
  
  # Find which y polygon contains each x centroid
  intersections <- st_within(gridx_centroids, gridy)
  
  # Create mapping table
  mapping_table <- data.frame(
    x_id = seq_len(nrow(gridx)),
    y_id = sapply(intersections, function(x) if(length(x) > 0) x[1] else NA)
  )
  
  # Function to aggregate values
  aggregate_to_y <- function(x_values, mapping) {
    y_values <- numeric(nrow(gridy))
    for(i in seq_along(y_values)) {
      x_indices <- mapping$x_id[mapping$y_id == i]
      if(length(x_indices) > 0) {
        y_values[i] <- mean(x_values[x_indices], na.rm = TRUE)
      } else {
        y_values[i] <- NA
      }
    }
    return(y_values)
  }
  
  # Transfer X variables to Y grid
  x_vars <- grep("covariate_x_", names(gridx), value = TRUE)
  result <- gridy
  
  for(var in x_vars) {
    result[[var]] <- aggregate_to_y(gridx[[var]], mapping_table)
  }
  
  return(result)
}

#' Fit Dasymetric Model
#'
#' Fits regression model to dasymetrically mapped data
#'
#' @param mapped_data Output from dasymetric_mapping
#' @param outcome_type Distribution type: 'normal', 'poisson', or 'binomial'
#'
#' @return Data frame with parameter estimates and confidence intervals
#' @export
#' @importFrom sf st_drop_geometry
#' @importFrom stats glm as.formula poisson binomial lm confint
fit_dasymetric_model <- function(mapped_data,outcome_type) {
  # Remove geometry and convert to regular data frame
  mapped_data <- st_drop_geometry(mapped_data)
  
  # Get x and y variable names
  x_vars <- grep("covariate_x_", names(mapped_data), value = TRUE)
  y_vars <- grep("covariate_y_", names(mapped_data), value = TRUE)
  all_vars <- c(x_vars, y_vars)
  
  if (outcome_type=='poisson'){
    
    # Fit Poisson model: y as the count outcome with offset
    formula <- as.formula(paste("y ~", paste(all_vars, collapse = " + ")))
    model <- glm(formula, family = poisson(link = "log"),offset=log(population),data = mapped_data)
  
  } else if (outcome_type=='binomial'){
    
    # Fit binomial model: y as successes, population as trials
    formula <- as.formula(paste("cbind(y, population - y) ~", paste(all_vars, collapse = " + ")))
    model <- glm(formula, family = binomial(link = "logit"), data = mapped_data)
    
  } else if (outcome_type=='normal'){
    
    mapped_data$y_percap<-mapped_data$y/mapped_data$population
    formula <- as.formula(paste("y_percap ~", paste(all_vars, collapse = " + ")))
    model <- lm(formula, data = mapped_data)
  
  } else{
    
      stop('Outcome type not recognized!')
    
    }
  
  # Calculate rate ratios and CIs
  coef_table <- summary(model)$coefficients
  ci_lower <- coef_table[,1] - 1.96*coef_table[,2]
  ci_upper <- coef_table[,1] + 1.96*coef_table[,2]
  
  results <- data.frame(
    beta_hat = coef_table[,1],
    ci_lower = ci_lower, 
    ci_upper = ci_upper
  )
  
  return(results)
}
