
#The reORBgen function which is described in the reORBgen.Rd file
#Implementation of the Copas et al. (2019) method for ORB adjustment
#in meta-anlyses of clinical trials, for the random effects model.
#re: random effect
#ORB: correction
#gen: generic, because it takes in an observed treatment effect, assumed to be
#     normally distributed, with standard error and sample sizes


#'
library(rootSolve)

#Takes in  y: observed treatment effect y normally distributed
#          s: standard error
#          n1: sample size treatment
#          n2: sample size control
#          init_param: c(mean, heterogeneity) initial guess of true treatment effect
#                      and heterogeneity parameters for initialization of optimiziation
#          alpha: confidence level

#' @export
reORBgen <- function(y, s, n1, n2, outcome, init_param, alpha, true.SE=NULL, LR.CI = TRUE) {

  #Indecies where we have the reported outcomes and the unreported with high risk (HR)
  rep_index <- which(!is.na(y) & y != "high")
  HR_index <- which(y == "high")

  #obs treatment effect, standard error^2, and sample sizes of reported outcomes
  logRR <- as.numeric(y[rep_index])
  sigma_squared <- (as.numeric(s[rep_index]))^2
  n1_rep <- as.numeric(n1[rep_index])
  n2_rep <- as.numeric(n2[rep_index])


  #total sample size of unreported outcomes
  n_HR <- as.numeric(n1[HR_index]) + as.numeric(n2[HR_index])

  #Copas et al. 2014, 2019 method for imputation of missing variances

  #k estimation, based on reported studies
  k <- sum(1/sigma_squared)/sum((n1_rep + n2_rep))


  if (!is.null(true.SE)){

    sigma_squared_imputed <- (as.numeric(true.SE)[HR_index])^2
  } else {
    # imputed variances for the HR studies                          #ok until now
    sigma_squared_imputed <- 1/(k*n_HR)

  }


  #Unadjusted log-likelihood function to be maximized
  f.u <- function(params, logRR, sigma_squared) {
    mu <- params[1]
    tau_squared <- params[2]
    -(1/2)*sum(log(sigma_squared + tau_squared) + ((logRR - mu)^2)/(sigma_squared + tau_squared))
  }

  #Set initial values for mu and tau_squared
  init_params <- init_param

  #Maximize unadjusted function optim()
  fit.u <- optim(init_params, f.u, logRR = logRR, sigma_squared = sigma_squared,
                 method = "Nelder-Mead",
                 control = list(fnscale = -1),
                 hessian=FALSE)

  #Return unadjusted mu and tau_squared
  mle.u <- fit.u$par[1]
  mle.tau <- max(fit.u$par[2],0)



  if(outcome == "benefit"){

    # adjusted likelihood function for beneficial outcome to be maximized
    f.adj.b <- function(params, logRR, sigma_squared, sigma_squared_imputed) {
      mu <- params[1]
      tau_squared <- params[2]
      z_alpha <- qnorm(1 - alpha/2)
      -(1/2)*sum(log(sigma_squared + tau_squared) + ((logRR - mu)^2)/(sigma_squared + tau_squared)) +
        sum(log(pnorm((z_alpha*sqrt(sigma_squared_imputed) - mu)/sqrt(sigma_squared_imputed + tau_squared)) -
                  pnorm((-z_alpha*sqrt(sigma_squared_imputed) - mu)/sqrt(sigma_squared_imputed + tau_squared))))

    }

    #Maximize log-likelihood
    fit.adj.b <- optim(init_params, f.adj.b, logRR = logRR, sigma_squared = sigma_squared, sigma_squared_imputed = sigma_squared_imputed,
                       method = "Nelder-Mead",
                       control = list(fnscale = -1),

                       hessian=FALSE)

    #Return adjusted mu and tau_squared
    mle.b <- fit.adj.b$par[1]
    mle.b.tau <- max(fit.adj.b$par[2],0)

    #WALD CONFIDENCE INTERVALS
    #a <- 0.05 #for harm Copas et al use 99% conf level
    #Unadjusted
    #fisher_info.u <- solve(-fit.u$hessian)
    #s.u <- sqrt(diag(fisher_info.u)[1])
    #ci.u <- fit.u$par[1] + qnorm(c(a/2, 1-a/2)) * s.u
    #Adjusted benefit
    #fisher_info.adj.b <- solve(-fit.adj.b$hessian)
    #s.adj.b <- sqrt(diag(fisher_info.adj.b)[1])
    #ci.u.adj.b <- fit.adj.b$par[1] + qnorm(c(a/2, 1-a/2)) * s.adj.b

    #LIKELIHOOD RATIO CONFIDENCE INTERVALS

    if (LR.CI){

    #Unadjusted
    z <- qchisq(1-alpha, df=1) #3.841

    #Re-write log likelihood with two inputs instead of param = c(mu, tau_squared)
    ll.u <- function(mu, tau_squared, logRR, sigma_squared) {

      -(1/2)*sum(log(sigma_squared + tau_squared) + ((logRR - mu)^2)/(sigma_squared + tau_squared))
    }


    pl.u <- function(mu, logRR, sigma_squared) { #take in vector of mus

      res <- mu

      for (i in seq_along(mu)) { #for all these values of mu
        optimResult <- optim(par = init_param[2],
                             fn = function(tau_squared) ll.u(mu[i], tau_squared, logRR=logRR, sigma_squared = sigma_squared),
                             method = "Nelder-Mead",
                             #lower= 0,
                             control = list(fnscale = -1))

        res[i] <- optimResult$value
      }
      return(res)
    }

    f <- function(mu, logRR, sigma_squared){
      pl.u(mu, logRR=logRR, sigma_squared=sigma_squared) - pl.u(mle.u, logRR=logRR, sigma_squared=sigma_squared) + 1/2*qchisq(0.95, df=1)
    }

    eps <- sqrt(.Machine$double.eps)
    lowerBound.u <- uniroot(f, interval = c(-1 + eps, mle.u), logRR=logRR, sigma_squared=sigma_squared)$root
    upperBound.u <- uniroot(f, interval = c( mle.u, 1 - eps), logRR=logRR, sigma_squared=sigma_squared)$root


    #Adjusted benefit

    #Re-write log likelihood with two inputs instead of param = c(mu, tau_squared)
    ll.b <- function(mu, tau_squared, logRR, sigma_squared, sigma_squared_imputed) {

      z_alpha <- qnorm(1-alpha/2)

      -(1/2)*sum(log(sigma_squared + tau_squared) + ((logRR - mu)^2)/(sigma_squared + tau_squared)) +
        sum(log(pnorm((z_alpha*sqrt(sigma_squared_imputed) - mu)/sqrt(sigma_squared_imputed + tau_squared)) -
                  pnorm((-z_alpha*sqrt(sigma_squared_imputed) - mu)/sqrt(sigma_squared_imputed + tau_squared))))
    }


    pl.b <- function(mu, logRR, sigma_squared, sigma_squared_imputed) { #take in vector of mus

      res <- mu

      for (i in seq_along(mu)) { #for all these values of mu
        optimResult <- optim(par = init_param[2],
                             fn = function(tau_squared) ll.b(mu[i], tau_squared,
                                                             logRR=logRR,
                                                             sigma_squared = sigma_squared,
                                                             sigma_squared_imputed = sigma_squared_imputed),
                             method = "Nelder-Mead",

                             control = list(fnscale = -1))

        res[i] <- optimResult$value
      }
      return(res)
    }

    f.b <- function(mu, logRR, sigma_squared, sigma_squared_imputed){

      pl.b(mu, logRR=logRR, sigma_squared=sigma_squared, sigma_squared_imputed = sigma_squared_imputed) - pl.b(mle.b, logRR=logRR, sigma_squared=sigma_squared, sigma_squared_imputed = sigma_squared_imputed) + 1/2*qchisq(0.95, df=1)
    }

    lowerBound.b <- uniroot(f.b, interval = c(-0.5, mle.b), logRR=logRR,
                                sigma_squared=sigma_squared,
                                sigma_squared_imputed=sigma_squared_imputed)$root
    upperBound.b <- uniroot(f.b, interval = c(mle.b, 0.5), logRR=logRR,
                                sigma_squared=sigma_squared,
                                sigma_squared_imputed = sigma_squared_imputed)$root


    return(list(mu_unadjusted = mle.u,
                LR_mu_unadjusted_low = lowerBound.u,
                LR_mu_unadjusted_up = upperBound.u,


                #CI_unadjusted_low_WALD = ci.u[1],
                #CI_unadjusted_up_WALD = ci.u[2],

                mu_adjusted_benefit = mle.b,
                LR_mu_adjusted_low = lowerBound.b,
                LR_mu_adjusted_up = upperBound.b,

                tau_squared_unadjusted = mle.tau,
                tau_squared_adjusted = mle.b.tau


                #CI_adjusted_benefit_low_WALD = ci.u.adj.b[1],
                #CI_adjusted_benefit_up_WALD = ci.u.adj.b[2]

    ))


    } else {



      return(list(
        mu_unadjusted = mle.u,
        mu_adjusted_benefit = mle.b,
        tau_squared_unadjusted = mle.tau,
        tau_squared_adjusted = mle.b.tau
      ))


    }
  } else if (outcome == "harm"){

    #Adjusted log-likelihood function for harmful outcome to be maximized
    f.adj.h <- function(params, logRR, sigma_squared, sigma_squared_imputed) {
      mu <- params[1]
      tau_squared <- params[2]
      z_alpha <- qnorm(1 - alpha/2)
      -(1/2)*sum(log(sigma_squared + tau_squared) + ((logRR - mu)^2)/(sigma_squared + tau_squared)) +
        sum(log(pnorm(mu/sqrt(sigma_squared_imputed + tau_squared))))
    }

    #Maximize log-likelihoood
    fit.adj.h <- optim(init_params, f.adj.h, logRR = logRR, sigma_squared = sigma_squared, sigma_squared_imputed = sigma_squared_imputed,
                       method = "Nelder-Mead",
                       control = list(fnscale = -1),

                       hessian = FALSE)

    #Adjusted mu and tau_squared estimates
    mle.h <- fit.adj.h$par[1]
    mle.h.tau <- max(fit.adj.h$par[2],0)



    #LIKELIHOOD RATIO CONFIDENCE INTERVALS

    if (LR.CI){


    #Unadjusted
    z <- qchisq(1-alpha, df=1) #3.841

    #Re-write log likelihood with two inputs instead of param = c(mu, tau_squared)
    ll.u <- function(mu, tau_squared, logRR, sigma_squared) {

      -(1/2)*sum(log(sigma_squared + tau_squared) + ((logRR - mu)^2)/(sigma_squared + tau_squared))
    }


    pl.u <- function(mu, logRR, sigma_squared) { #take in vector of mus

      res <- mu

      for (i in seq_along(mu)) { #for all these values of mu
        optimResult <- optim(par = init_param[2],
                             fn = function(tau_squared) ll.u(mu[i], tau_squared, logRR=logRR, sigma_squared = sigma_squared),
                             method = "Nelder-Mead",
                             #lower= 0,
                             control = list(fnscale = -1))

        res[i] <- optimResult$value
      }
      return(res)
    }

    f <- function(mu, logRR, sigma_squared){
      pl.u(mu, logRR=logRR, sigma_squared=sigma_squared) - pl.u(mle.u, logRR=logRR, sigma_squared=sigma_squared) + 1/2*qchisq(0.95, df=1)
    }

    lowerBound.u <- uniroot(f, interval = c(-0.5, mle.u), logRR=logRR, sigma_squared=sigma_squared)$root
    upperBound.u <- uniroot(f, interval = c(mle.u, 0.5), logRR=logRR, sigma_squared=sigma_squared)$root


    ll.h <- function(mu, tau_squared, logRR, sigma_squared, sigma_squared_imputed) {

      z_alpha <- qnorm(1 - alpha/2)
      -(1/2)*sum(log(sigma_squared + tau_squared) + ((logRR - mu)^2)/(sigma_squared + tau_squared)) +
        sum(log(pnorm(mu/sqrt(sigma_squared_imputed + tau_squared))))
    }


    pl.h <- function(mu, logRR, sigma_squared, sigma_squared_imputed) { #take in vector of mus

      res <- mu

      for (i in seq_along(mu)) { #for all these values of mu
        optimResult <- optim(par = init_param[2],
                             fn = function(tau_squared) ll.h(mu[i], tau_squared, logRR=logRR, sigma_squared = sigma_squared, sigma_squared_imputed = sigma_squared_imputed),
                             method = "L-BFGS-B",
                             lower= 0,
                             control = list(fnscale = -1))

        res[i] <- optimResult$value
      }
      return(res)
    }

    f.h <- function(mu, logRR, sigma_squared, sigma_squared_imputed){
      pl.h(mu, logRR=logRR, sigma_squared=sigma_squared, sigma_squared_imputed = sigma_squared_imputed) - pl.h(mle.h, logRR=logRR, sigma_squared=sigma_squared, sigma_squared_imputed = sigma_squared_imputed) + 1/2*qchisq(0.99, df=1)
    }

    lowerBound.h <- uniroot(f.h, interval = c(-0.5, mle.h), logRR=logRR, sigma_squared=sigma_squared, sigma_squared_imputed = sigma_squared_imputed)$root
    upperBound.h <- uniroot(f.h, interval = c(mle.h, 0.5), logRR=logRR, sigma_squared=sigma_squared, sigma_squared_imputed = sigma_squared_imputed)$root





    #WALD CONFIDENCE INTERVALS
    #a <- 0.01 #for harm Copas et al use 99% conf level
    #Unadjusted
    #fisher_info.u <- solve(-fit.u$hessian)
    #s.u <- sqrt(diag(fisher_info.u)[1])
    #ci.u <- fit.u$par[1] + qnorm(c(a/2, 1-a/2)) * s.u
    #Adjusted harm
    #fisher_info.adj.h <- solve(-fit.adj.h$hessian)
    #s.adj.h <- sqrt(diag(fisher_info.adj.h)[1])
    #ci.u.adj.h <- fit.adj.h$par[1] + qnorm(c(a/2, 1-a/2)) * s.adj.h



    return(list(mu_unadjusted = mle.u,
                LR_mu_unadjusted_low_h = lowerBound.u,
                LR_mu_unadjusted_up_h = upperBound.y,


                #CI_unadjusted_low_WALD = exp(ci.u)[1],
                #CI_unadjusted_up_WALD = exp(ci.u)[2],

                mu_adjusted_harm = mle.h,
                LR_mu_adjusted_low_h = lowerBound.h,
                LR_mu_adjusted_up_h = upperBound.h,


                tau_squared_unadjsuted = mle.tau,
                tau_squared_adjusted = mle.h.tau


                #CI_adjusted_harm_low_WALD = exp(ci.u.adj.h)[1],
                #CI_adjusted_harm_up_WALD = exp(ci.u.adj.h)[2]

    ))

    }else {
      return(list(
        mu_unadjusted = mle.u,
        mu_adjusted_harm = mle.h,
        tau_squared_unadjusted = mle.tau,
        tau_squared_adjusted = mle.h.tau
      ))
}


  } else {

    return("invalid outcome input")
  }


}



