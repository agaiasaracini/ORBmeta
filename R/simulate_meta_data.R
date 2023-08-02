
#' @export

#To simulate one data set of a meta analysis with  ORB
simulate_meta_data <- function(n_studies, mu, tau_squared, n_treatment=50, n_control=50, gamma){

  if (!require("dplyr")) {
    # If the package is not loaded, try to load it
    if (!requireNamespace("dplyr", quietly = TRUE)) {
      # If the package is not installed, give an error message
      stop("The 'dplyr' package is not installed.")
    } else {
      # If the package is installed but not loaded, load it
      library("dplyr")
    }
  }



  #Number of studies per meta-analysis
  n_studies <- n_studies

  #True treatment effect and true heterogeneity parameter
  mu <- mu
  tau_squared <- tau_squared

  #Empty data frame
  meta_data <- data.frame(study = 1:n_studies,
                          n_control = numeric(n_studies),
                          n_treatment = numeric(n_studies),
                          y = numeric(n_studies),
                          s = numeric(n_studies),
                          p_value = numeric(n_studies))


  n_treatment <- rep(n_treatment, n_studies)
  n_control <- rep(n_control, n_studies)


  # Generate true treatment effect mean and standard errors for all studies at once
  theta_i <- rnorm(n_studies, mean = mu, sd = sqrt(tau_squared))
  sigma_i <- sqrt(rchisq(n_studies, df = 2*n_control - 2) / ((n_control - 1) * n_control))

  # Generate observed treatment effect for all studies at once
  y_i <- rnorm(n_studies, mean = theta_i, sd = sigma_i)

  # Calculate one-sided p-values for all studies at once
  p_i <- pnorm(y_i / sigma_i, lower.tail = FALSE)
  p_i <- round(p_i, 10)


  # Fill in data frame
  meta_data$n_control <- n_control
  meta_data$n_treatment <- n_treatment
  meta_data$y <- y_i
  meta_data$s <- sigma_i
  meta_data$p_value <- p_i


  #Missing data mechanism
  meta_data_miss <- meta_data


  #probability of being missing, ie 1-prob of selection
  replace_prob <- exp(-4 *as.numeric(meta_data$p_value)^gamma)
  high_index <- runif(n_studies) >= replace_prob
  meta_data_miss$y[high_index] <- "high"
  meta_data_miss$s[high_index] <- "high"


  #among the rest, choose some at random some to be low risk
  if (length(which(meta_data_miss$y == "high")) > 0) { # if there is at least one HR
    non_high_indices <- which(meta_data_miss$y != "high") # indices reported
    n_high_rows <- length(which(meta_data_miss$y == "high")) # number of HR
    n_values_to_replace <- min(sample(0:n_high_rows, 1), length(non_high_indices)) # number of values to replace so that never more low than high
    replace_indices <- ifelse(length(non_high_indices) == 1,
                              non_high_indecies, #if there is only one non-high study, pick that one to be low
                              sample(non_high_indices, n_values_to_replace)) #if there are more than one, pick
    meta_data_miss$y[replace_indices] <- "low"
    meta_data_miss$s[replace_indices] <- "low"

  } else {

    meta_data_miss$y <- meta_data_miss$y
    meta_data_miss$s <- meta_data_miss$s

  }

  #make sure we have at least 3 studies that are reported otheriwse calculationg of ci unstable
  #pls also v unrealistic that there would be a meta analysis conducted on 2 studies!
  non_high_low_rows <- which(!(meta_data_miss$y %in% c("high", "low")))
  n_non_high_low_rows <- length(non_high_low_rows)

  if (n_non_high_low_rows < 3) {
    high_low_rows <- which(meta_data_miss$y %in% c("high", "low"))

    if (n_non_high_low_rows == 0) {
      random_rows <- sample(high_low_rows, 3)
    } else {
      random_rows <- sample(high_low_rows, 3 - n_non_high_low_rows)
    }

    meta_data_miss$y[random_rows] <- meta_data$y[random_rows]
    meta_data_miss$s[random_rows] <- meta_data$s[random_rows]
  }

  meta_data_miss$s_true <- meta_data$s

  return(meta_data_miss)

}


