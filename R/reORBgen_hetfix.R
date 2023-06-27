

library(rootSolve)
#If tau squared is known, init_param is just one number, an initial guess of mu!
#Not a vector of two entries anymore

reORBgen_hetfix <- function(y, s, n1, n2, outcome, init_param, alpha, true.SE=NULL, fixed.tau.squared) {

  #Indecies where we have the reported outcomes and the unreported with high risk (HR)
  rep_index <- which(!is.na(y) & y != "high")
  HR_index <- which(y == "high")

  #Observed treatment effect, standard error^2, and sample sizes of reported outcomes
  logRR <- as.numeric(y[rep_index])
  sigma_squared <- (as.numeric(s[rep_index]))^2
  n1_rep <- as.numeric(n1[rep_index])
  n2_rep <- as.numeric(n2[rep_index])


  #Total sample size of unreported outcomes
  n_HR <- as.numeric(n1[HR_index]) + as.numeric(n2[HR_index])

  #Copas et al. (2014, 2019) method for imputation of missing variances

  #k estimation, based on reported studies
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


  #Fixed tau squared
  tau_squared <- fixed.tau.squared

  #Likelihood function is a function of mu only
  f.u <- function(mu, logRR, sigma_squared, tau_squared){

    -(1/2)*sum(log(sigma_squared + tau_squared) + ((logRR - mu)^2)/(sigma_squared + tau_squared))
  }

  #Initial parameter is only mu
  init_params <- init_param

  #Optimize only with respect to mu
  fit.u <- optim(init_params, f.u, logRR = logRR, sigma_squared = sigma_squared, tau_squared=tau_squared,
                 method = "Brent",

                 lower = -10,
                 upper = 10,
                 control = list(fnscale = -1),
                 hessian=FALSE)
  mle.u <- fit.u$par



  #Beneficial outcome adjustment for ORB

  if(outcome == "benefit"){


    #Adjusted log-likelihood function for beneficial outcome to be maximized
    f.adj.b <- function(params, logRR, sigma_squared, sigma_squared_imputed, tau_squared) {
      mu <- params
      z_alpha <- qnorm(1 - alpha/2)
      -(1/2)*sum(log(sigma_squared + tau_squared) + ((logRR - mu)^2)/(sigma_squared + tau_squared)) +

        if (length(sigma_squared_imputed)>0){
          sum(log(pnorm((z_alpha*sqrt(sigma_squared_imputed) - mu)/sqrt(sigma_squared_imputed + tau_squared)))) #one-sided
        }else{
          0
        }

    }

    #Set initial values for mu and tau_squared
    init_params <- init_param

    #Maximize log-likelihood
    fit.adj.b <- optim(init_params, f.adj.b, logRR = logRR, sigma_squared = sigma_squared, sigma_squared_imputed = sigma_squared_imputed, tau_squared=tau_squared,
                       method = "Brent",
                       control = list(fnscale = -1),
                       lower = -10,
                       upper = 10,

                       hessian=FALSE)

    #Return adjusted mu and tau_squared
    mle.b <- fit.adj.b$par


    #LIKELIHOOD RATIO CONFIDENCE INTERVALS


    #Unadjusted
    z <- qchisq(1-alpha, df=1) #3.841

    LRstat.u <- function(mu){

      2*(f.u(mle.u, logRR, sigma_squared, tau_squared) - f.u(mu, logRR, sigma_squared, tau_squared))

    }

    #From Held's book solutions likelihood inference
    LR.lower.u <- uniroot(Vectorize(function(mu){LRstat.u(mu) - z}),
                          interval=c(-10,mle.u))$root

    LR.upper.u <- uniroot(Vectorize(function(mu){LRstat.u(mu) - z}),
                          interval=c(mle.u, 10))$root


    #Benefit Adjusted
    LRstat.b <- function(mu){

      2*(f.adj.b(mle.b, logRR, sigma_squared, sigma_squared_imputed, tau_squared) - f.adj.b(mu, logRR, sigma_squared, sigma_squared_imputed, tau_squared))

    }

    LR.lower.b <- uniroot(Vectorize(function(mu){LRstat.b(mu) - z}),
                          interval=c(-10,mle.b))$root

    LR.upper.b <- uniroot(Vectorize(function(mu){LRstat.b(mu) - z}),
                          interval=c(mle.b, 10))$root

    return(list(mu_unadjusted = mle.u,
                LR_mu_unadjusted_low = LR.lower.u,
                LR_mu_unadjusted_up = LR.upper.u,


                #CI_unadjusted_low_WALD = ci.u[1],
                #CI_unadjusted_up_WALD = ci.u[2],

                mu_adjusted_benefit = mle.b,
                LR_mu_adjusted_low = LR.lower.b,
                LR_mu_adjusted_up = LR.upper.b,


                average_sigma_squared_unadjusted = sigma_squared_average_unadjusted,
                sigma_squared_average_adjusted = sigma_squared_average_adjusted


                #CI_adjusted_benefit_low_WALD = ci.u.adj.b[1],
                #CI_adjusted_benefit_up_WALD = ci.u.adj.b[2]

    ))





    #Adjustment for ORB in harmful outcome
  } else if (outcome == "harm"){

    #Adjusted log-likelihood function for harmful outcome to be maximized
    f.adj.h <- function(mu, logRR, sigma_squared, sigma_squared_imputed, tau_squared) {


      z_alpha <- qnorm(1 - alpha/2)
      -(1/2)*sum(log(sigma_squared + tau_squared) + ((logRR - mu)^2)/(sigma_squared + tau_squared)) +
        sum(log(pnorm(mu/sqrt(sigma_squared_imputed + tau_squared))))
    }

    init_params <- init_param

    #Maximize log-likelihoood
    fit.adj.h <- optim(init_params, f.adj.h, logRR = logRR, sigma_squared = sigma_squared, sigma_squared_imputed = sigma_squared_imputed, tau_squared=tau_squared,
                       method = "Brent",
                       lower = -10,
                       upper = 10,
                       control = list(fnscale = -1),

                       hessian = FALSE)

    #Adjusted mu and tau_squared estimates
    mle.h <- fit.adj.h$par


    z <- qchisq(0.99, df=1) #Copas et al use alt

    #Unadjsuted
    LRstat.u <- function(mu){

      2*(f.u(mle.u, logRR, sigma_squared, tau_squared) - f.u(mu, logRR, sigma_squared, tau_squared))

    }

    #From Held's book solutions likelihood inference
    LR.lower.u <- uniroot.all(Vectorize(function(mu){LRstat.u(mu) - z}),
                              interval=c(-10, 10))[1]

    LR.upper.u <- uniroot.all(Vectorize(function(mu){LRstat.u(mu) - z}),
                              interval=c(-10, 10))[2]

    #harm djusted
    LRstat.h <- function(mu){

      2*(f.adj.h(mle.h, logRR, sigma_squared, sigma_squared_imputed, tau_squared) - f.adj.h(mu, logRR, sigma_squared, sigma_squared_imputed, tau_squared))

    }

    LR.lower.h <- uniroot.all(Vectorize(function(mu){LRstat.h(mu) - z}),
                              interval=c(-10,10))[1]

    LR.upper.h <- uniroot.all(Vectorize(function(mu){LRstat.h(mu) - z}),
                              interval=c(-10,10))[2]



    return(list(mu_unadjusted = mle.u,
                LR_mu_unadjusted_low_h = LR.lower.u,
                LR_mu_unadjusted_up_h = LR.upper.u,


                #CI_unadjusted_low_WALD = exp(ci.u)[1],
                #CI_unadjusted_up_WALD = exp(ci.u)[2],

                mu_adjusted_harm = mle.h,
                LR_mu_adjusted_low_h = LR.lower.h,
                LR_mu_adjusted_up_h = LR.upper.h,

                average_sigma_squared_unadjusted = sigma_squared_average_unadjusted,
                sigma_squared_average_adjusted = sigma_squared_average_adjusted


                #CI_adjusted_harm_low_WALD = exp(ci.u.adj.h)[1],
                #CI_adjusted_harm_up_WALD = exp(ci.u.adj.h)[2]

    ))




  } else {

    return("invalid outcome input")
  }


}
