

#To simulate one data set of a meta analysis with and without ORB
simulate_meta_data <- function(n_studies, mu, tau_squared, min_n=50, max_n=500, seed, gamma){

  if (!require("biostatUZH")) {
    # If the package is not loaded, try to load it
    if (!requireNamespace("biostatUZH", quietly = TRUE)) {
      # If the package is not installed, give an error message
      stop("The 'biostatUZH' package is not installed.")
    } else {
      # If the package is installed but not loaded, load it
      library("biostatUZH")
    }
  }

set.seed(seed)

#Number of studies per meta-analysis
n_studies <- n_studies

#Minimum and maximum sample sizes per treatment arm
min_n <- min_n
max_n <- max_n

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

#Loop over studies and fill in data frame
for (i in 1:n_studies) {

  #Generate sample sizes for control and treatment groups
  n_control <- sample(min_n:max_n, 1)
  n_treatment <- sample(min_n:max_n, 1)

  #Ensure similar sample sizes in control and treatment
  if (abs(n_control - n_treatment) > 50) {
    if (n_control > n_treatment) {
      n_control <- n_treatment + 50
    } else {
      n_treatment <- n_control + 50
    }
  }

  #Generate true treatment effect mean and standard errors
  theta_i <- rnorm(1, mean = mu, sd = sqrt(tau_squared))

  #Standard errors from Chisquare distribution
  n_i <- n_control+n_treatment
  sigma_i <- sqrt(rchisq(1, df= 2*n_i - 2)/((n_i -1)*n_i))

  #Generate observed treatment effect
  y_i <- rnorm(1, mean = theta_i, sd = sigma_i)

  #Calculate one-sided p-value
  p_i <- pnorm(y_i / sigma_i, lower.tail = FALSE)
  p_i <- round(p_i, 10)

  #Fill in data frame
  meta_data[i, "n_control"] <- n_control
  meta_data[i, "n_treatment"] <- n_treatment
  meta_data[i, "y"] <- y_i
  meta_data[i, "s"] <- sigma_i
  meta_data[i, "p_value"] <- p_i
}


#Missing data mechanism:
#if not significant, missing and classified as high risk of bias

#dicht selection: if non-significant>missing
#NB this might not alawys hold bc we need at least 3 studies so we might have to pu back
#some non-significant ones
#meta_data_missing <- meta_data
#meta_data_missing$y <- ifelse(meta_data_missing$p_value >= 0.05,
                             # "high",
                            #  meta_data_missing$y)
#meta_data_missing$s <- ifelse(meta_data_missing$p_value >= 0.05,
                             # "high",
                              #meta_data_missing$s)

#probability of missing as a function of the pvalue
# Replace studies that reported with "high" according to selection prob.
# Calculate the replacement probability
#replace_prob <- 1 - exp(-4 * pnorm(-as.numeric(meta_data$y) / as.numeric(meta_data$s))^gamma)


#probability of being missing, ie 1-prob of selection
replace_prob <- exp(-4 *as.numeric(meta_data$p_value)^gamma)
# Generate a vector of "high" based on the replacement probability
meta_data_missing <- meta_data
for (i in 1:n_studies) {
  if (runif(1) >= replace_prob[i]) {
    meta_data_missing$y[i] <- "high"
    meta_data_missing$s[i] <- "high"
  }
}
#replace_vector <- ifelse(runif(n_studies) < replace_prob, "high", "")
# Replace the values of y and s with "high" where needed

#meta_data_missing$y <- ifelse(replace_vector == "high", "high", meta_data$y)
#meta_data_missing$s <- ifelse(replace_vector == "high", "high", meta_data$s)





#Format p-value & other columns
meta_data$p_value <- formatPval(meta_data$p_value)
meta_data$n_control <- format(meta_data$n_control)
meta_data$n_treatment <- format(meta_data$n_treatment)
meta_data$y <- format(meta_data$y)
meta_data$s <- format(meta_data$s)


#make sure we have at least 3 studies that are reported otheriwse calculationg of ci unstable
#pls also v unrealistic that there would be a meta analysis conducted on 2 studies!
if (length(which(meta_data_missing$y != "high")) < 3) {

  non_high_rows <- which(meta_data_missing$y != "high")
  high_rows <- which(meta_data_missing$y == "high")

  random_rows <- sample(high_rows, 3 - length(non_high_rows))
  meta_data_missing$y[random_rows] <- meta_data$y[random_rows]
  meta_data_missing$s[random_rows] <- meta_data$s[random_rows]

}

colnames(meta_data_missing) <- c("study ID",
                                 "N control",
                                 "N treatment",
                                 "y",
                                 "sd",
                                 "p-value")

return(list(meta_data = meta_data,
            meta_data_missing = meta_data_missing))

}


