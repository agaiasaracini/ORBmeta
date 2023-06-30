
#The feORBbin function which is described in the feORBbin.Rd file
#Implementation of the Copas et al. (2019) method for ORB adjustment
#in meta-anlyses of clinical trials, for the random effects model.
#fe: fixed effect
#ORB: correction
#bin: takes observed cell counts for treatment and control

library(rootSolve)

feORBbin <- function(a, c, n1, n2, outcome, init_param, alpha_ben, alpha_harm, true.SE=NULL) {

  #Indecies where we have the reported outcomes (not high nor low) and the HR outcomes
  rep_index <- which(a != "high" & a != "low")
  HR_index <- which(a == "high")

  # a,c,n1,n2 values for the reported outcomes
  a_rep <- as.numeric(a[rep_index])
  c_rep <- as.numeric(c[rep_index])
  n1_rep <- as.numeric(n1[rep_index])
  n2_rep <- as.numeric(n2[rep_index])

  #CONTINUITY CORRECTION
  if (length(a_rep == 0) | length(c_rep == 0) > 0){

    a_0_index <- c(which(a_rep == 0), which(c_rep == 0)) #WHERE ARE THE ZERO CELL COUNTS
    a_rep[a_0_index] <- a_rep[a_0_index] + 0.5
    c_rep[a_0_index] <- c_rep[a_0_index] + 0.5
    n1_rep[a_0_index] <- n1_rep[a_0_index] + 0.5
    n2_rep[a_0_index] <- n2_rep[a_0_index] + 0.5 #ADD 0.5 TO ALL OF THE CELLS IN THAT STUDY!

  } else {

    a_rep <- a_rep
    c_rep <- c_rep
    n1_rep <- n1_rep
    n2_rep <- n2_rep
  }


  # n1, n2 values for the unreported outcomes, needed for imputation of sigma sq for HR studies
  n_HR <- as.numeric(n1[HR_index]) + as.numeric(n2[HR_index])

  # Unadjusted logRR and sigma_sq of logRR for each study
  logRR <- log((a_rep*n2_rep)/(c_rep*n1_rep))
  sigma_squared <- ((n1_rep - a_rep)/(n1_rep*a_rep)) + ((n2_rep - c_rep)/(n2_rep*c_rep))


  # k estiamtion, based on reported studies
  k <- sum(1/sigma_squared)/sum((n1_rep + n2_rep))

  #Possibility to pass the true SE to the function
  if (!is.null(true.SE)){

    sigma_squared_imputed <- (as.numeric(true.SE)[HR_index])^2
  } else {
    #Imputed variances for the HR studies
    sigma_squared_imputed <- 1/(k*n_HR)

  }

  #Average sigma squared value that is used when we do not adjust for ORB and when we do
  sigma_squared_average_unadjusted <- mean(sigma_squared)
  sigma_squared_average_adjusted <- mean(c(sigma_squared, sigma_squared_imputed))

  # unadjusted likelihood function to be maximized
  f.u <- function(mu, logRR, sigma_squared) {

    -(1/2)*sum((((logRR - mu )^2)/(sigma_squared)))
  }

  # set initial values for mu and tau_squared
  init_param <- init_param

  # maximize unadjusted function optim() L-BFGS-B
  fit.u <- optim(init_param, f.u, logRR = logRR, sigma_squared = sigma_squared,
                 method = "Brent",
                 control = list(fnscale = -1),
                 lower = -10,
                 upper = 10,
                 hessian=TRUE)

  mle.u <- fit.u$par[1] #MLE ESTIMATE


  if(outcome == "benefit"){

    # adjusted likelihood function for harm to be maximized
    f.adj.b <- function(mu, logRR, sigma_squared, sigma_squared_imputed) {

      z_alpha <- qnorm(1 - alpha_ben/2)
      -(1/2)*sum((((logRR - mu )^2)/(sigma_squared))) +

        if (!is.null(sigma_squared_imputed)){
        sum(log(pnorm(z_alpha - mu/sqrt(sigma_squared_imputed))))
               #- pnorm(-z_alpha -mu/sqrt(sigma_squared_imputed))))
        } else {

          0
        }

    }

    fit.adj.b <- optim(init_param, f.adj.b, logRR = logRR, sigma_squared = sigma_squared, sigma_squared_imputed = sigma_squared_imputed,
                       method = "Brent",
                       control = list(fnscale = -1),
                       lower = -10,
                       upper=10,
                       hessian=TRUE)

    mle.b <- fit.adj.b$par[1]

    #Confidence Intervals WALD
    #Wald CI based on Hessian output of optim() function
    a <- alpha_ben
    #Unadjsted
    fisher_info.u <- solve(-fit.u$hessian)
    s.u <- sqrt(diag(fisher_info.u))
    ci.u <- fit.u$par + qnorm(c(a/2, 1-a/2)) * s.u
    #Adjusted
    fisher_info.adj.b <- solve(-fit.adj.b$hessian)
    s.adj.b <- sqrt(diag(fisher_info.adj.b))
    ci.u.adj.b <- fit.adj.b$par + qnorm(c(a/2, 1-a/2)) * s.adj.b

    #Likelihood Ratio Confidence Intervals
    z <- qchisq(1-a, df=1) #3.841

    #Copas et al. 2019
    #f.u.CI <- fucntion(mu){
    # -(f.u(mle.u, logRR, sigma_squared) - f.u(mu, logRR, sigma_squared) + (1/2)*(z^2))^2
    #}

    #Likelihood Ratio Statistics from likelihood theory
    #Unadjusted
    LRstat.u <- function(mu){

      2*(f.u(mle.u, logRR, sigma_squared) - f.u(mu, logRR, sigma_squared))

    }

    #From Held's book solutions likelihood inference
    #From Held's book solutions likelihood inference
    LR.lower.u <- uniroot(Vectorize(function(mu){LRstat.u(mu) - z}),
                          interval=c(-10,mle.u))$root

    LR.upper.u <- uniroot(Vectorize(function(mu){LRstat.u(mu) - z}),
                          interval=c(mle.u, 10))$root


    #Benefit Adjusted
    LRstat.b <- function(mu){

      2*(f.adj.b(mle.b, logRR, sigma_squared, sigma_squared_imputed) - f.adj.b(mu, logRR, sigma_squared, sigma_squared_imputed))

    }

    LR.lower.b <- uniroot(Vectorize(function(mu){LRstat.b(mu) - z}),
                          interval=c(-10,mle.b))$root

    LR.upper.b <- uniroot(Vectorize(function(mu){LRstat.b(mu) - z}),
                          interval=c(mle.b, 10))$root


    p_val <- pnorm(logRR / sigma_squared, lower.tail = FALSE) #one sided pvalue benefit

    return(list(RR_unadjusted = mle.u,
                CI_unadjusted_low = LR.lower.u,
                CI_unadjusted_up = LR.upper.u,

                CI_unadjusted_low_WALD = ci.u[1],
                CI_unadjusted_up_WALD = ci.u[2],

                RR_adjusted_benefit = mle.b,

                CI_adjusted_benefit_low = LR.lower.b,
                CI_adjusted_benefit_up = LR.upper.b,

                CI_adjusted_benefit_low_WALD = ci.u.adj.b[1],
                CI_adjusted_benefit_up_WALD = ci.u.adj.b[2],
                logRR = logRR,
                p_value = p_val


    ))


  } else if (outcome == "harm"){

    # adjusted likelihood function for benefit to be maximized
    f.adj.h <- function(mu, logRR, sigma_squared, sigma_squared_imputed) {

      z_alpha <- qnorm(1 - alpha_harm/2)
      -(1/2)*sum((((logRR - mu )^2)/(sigma_squared))) +

        if (!is.null(sigma_squared_imputed)){
        sum(log(pnorm(mu/sqrt(sigma_squared_imputed))))

        } else {
          0
        }
    }

    fit.adj.h <- optim(init_param, f.adj.h, logRR = logRR, sigma_squared = sigma_squared, sigma_squared_imputed = sigma_squared_imputed,
                       method = "Brent",
                       control = list(fnscale = -1),
                       lower = -10,
                       upper=10,
                       hessian=TRUE)

    mle.h <- fit.adj.h$par[1]

    #CONFIDENCE INTERVALS WALD
    a <- alpha_harm #for harm Copas et al use 99% conf level
    #Unadjusted
    fisher_info.u <- solve(-fit.u$hessian)
    s.u <- sqrt(diag(fisher_info.u))
    ci.u <- fit.u$par + qnorm(c(a/2, 1-a/2)) * s.u
    #Adjusted
    fisher_info.adj.h <- solve(-fit.adj.h$hessian)
    s.adj.h <- sqrt(diag(fisher_info.adj.h))
    ci.u.adj.h <- fit.adj.h$par + qnorm(c(a/2, 1-a/2)) * s.adj.h



    #Likelihood Ratio Confidence Intervals
    z <- qchisq(1-alpha_harm, df=1) #Copas et al use alt

    #Unadjsuted
    LRstat.u <- function(mu){

      2*(f.u(mle.u, logRR, sigma_squared) - f.u(mu, logRR, sigma_squared))

    }

    #From Held's book solutions likelihood inference
    LR.lower.u <- uniroot(Vectorize(function(mu){LRstat.u(mu) - z}),
                          interval=c(-10,mle.u))$root

    LR.upper.u <- uniroot(Vectorize(function(mu){LRstat.u(mu) - z}),
                          interval=c(mle.u, 10))$root

    #harm djusted
    LRstat.h <- function(mu){

      2*(f.adj.h(mle.h, logRR, sigma_squared, sigma_squared_imputed) - f.adj.h(mu, logRR, sigma_squared, sigma_squared_imputed))

    }

    LR.lower.h <- uniroot(Vectorize(function(mu){LRstat.h(mu) - z}),
                          interval=c(-10,mle.u))$root

    LR.upper.h <- uniroot(Vectorize(function(mu){LRstat.h(mu) - z}),
                          interval=c(mle.u, 10))$root

    p_val = pnorm(logRR / sigma_squared, lower.tail = FALSE)

    return(list(RR_unadjusted = mle.u,

                CI_unadjusted_low = LR.lower.u,
                CI_unadjusted_up = LR.upper.u,

                CI_unadjusted_low_WALD = ci.u[1],
                CI_unadjusted_up_WALD = ci.u[2],

                RR_adjusted_harm = mle.h,

                CI_adjusted_harm_low = LR.lower.h,
                CI_adjusted_harm_up = LR.upper.h,

                CI_adjusted_harm_low_WALD = ci.u.adj.h[1],
                CI_adjusted_harm_up_WALD = ci.u.adj.h[2],
                logRR = logRR,
                p_value = p_val


    ))

  } else {

    return("invalid outcome input")
  }


}
