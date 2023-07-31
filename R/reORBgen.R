
#The reORBgen function which is described in the reORBgen.Rd file
#Implementation of the Copas et al. (2019) method for ORB adjustment
#in meta-anlyses of clinical trials, for the random effects model.
#re: random effect
#ORB: correction
#gen: generic, because it takes in an observed treatment effect, assumed to be
#     normally distributed, with standard error and sample sizes



#library(rootSolve)

#Takes in  a: cell counts treatment
#          c: cell counts control
#          n1: sample size treatment
#          n2: sample size control
#          init_param: c(log RR, heterogeneity) initial guess of true treatment effect
#                      and heterogeneity parameters for initialization of optimization
#          alpha: confidence level
#          true.SE=NULL if we know the standard errors
#          LR.CI=FALSE if we want it to calculate the profile likelihood confidence intervals

#' @export

reORBgen <- function(a=NULL, c=NULL,
                     mu1=NULL, mu2=NULL, sd1=NULL, sd2=NULL,
                     y=NULL, s=NULL,
                     n1,
                     n2,
                     outcome,
                     init_param,
                     alpha_ben=NULL,
                     alpha_ben_one.sided = TRUE,
                     alpha_harm=NULL,
                     true.SE=NULL, LR.CI = FALSE,
                     Wald.CI = FALSE,
                     selection.benefit = "Copas.oneside",
                     weight_param=15,
                     opt_method="L-BFGS-B",
                     lower = c(-5, 0.00001),
                     upper = c(5,5),
                     low.risk = FALSE
) {


  if (!require("rootSolve")) {
    # If the package is not loaded, try to load it
    if (!requireNamespace("rootSolve", quietly = TRUE)) {
      # If the package is not installed, give an error message
      stop("The 'rootSolve' package is not installed.")
    } else {
      # If the package is installed but not loaded, load it
      library("rootSolve")
    }
  }


  #Default is to ignore the low risk of bias, but if TRUE then we
  #use all the unreported studies fro adjustment
  #this supposedly can be done for the two smooth selection functions

  #Can take in binary counts and calculate log RR: a, c
  #Means and calculate differences in means: y1, y2, sd1, sd2
  #Directly normally distributed effect measure: y, s


  #Beneficial outcome
  #In the original paper by Copas et al., they use a piece-wise selection function
  #which has P(unreporting)=1 if not-significant, where they consider two tails
  #however, more realistic that, if significant in the "wrong" direction,
  #even less likley to be reported! hence Copas.oneside = default in which we
  #dont consider the siginficance in the "wrong direction"

  #Ideally have
  #Copas.twoside from original paper here always two sided by nature of the function
  #Copas.oneside DEFAULT
  #Neg.exp.piecewise
  #Sigmoid.cont
  #For the first two, L-BFGS-B sometimes doesnt work well so Nelder-Mead needed

  #okay so if alpha_ben_one.sided=TRUE, apart from copas.twosides, which is always two sided by nature,
  #we use the one-sided threshold z=1.64
  #for the one.side copas and for the otehr selectiosn functions
  #depends on the context but for my simulations this makes most sense since im removeing/creatign ORB based on
  #the one sided p value
  sel.ben <- selection.benefit

  #If we have Exponential.piecewise or Sigmoid.continuous, we also need a weight parameter
  w <- weight_param

  if (Wald.CI){
    my.hessian <- TRUE
  } else{
    my.hessian <- FALSE
  }

  method <- opt_method
  lower <- lower
  upper <- upper

  #Significance
  # z_alpha <- qnorm(1-alpha/2)

  if (!is.null(a) & !is.null(c)){
    #Reported outcomes
    #Indecies where we don't have low/high, i.e., reported outcomes
    #If we turn C into numeric, the low/high (unreported) become NA
    #Where do we have low and were do we have high
    Rep_index <- which(!is.na(as.numeric(a)))

    if (low.risk == FALSE){

      HR_index <- which(a == "high")

    } else { ##if low.risk=TRUE, we include all unreported studies in the adjustment
      HR_index <- which(a == "high" | a == "low")
    }

    # a,c,n1,n2, n values for the reported studies
    a_rep <- as.numeric(a[Rep_index])
    c_rep <- as.numeric(c[Rep_index])
    n1_rep <- as.numeric(n1[Rep_index])
    n2_rep <- as.numeric(n2[Rep_index])

    if (length(a_rep == 0) | length(c_rep == 0) > 0){ #Continuity correction

      a_0_index <- c(which(a_rep == 0), which(c_rep == 0)) #Where are the zero cell counts
      a_rep[a_0_index] <- a_rep[a_0_index] + 0.5
      c_rep[a_0_index] <- c_rep[a_0_index] + 0.5
      n1_rep[a_0_index] <- n1_rep[a_0_index] + 0.5
      n2_rep[a_0_index] <- n2_rep[a_0_index] + 0.5 #Add 0.5 to the cell counts with 0

    } else {

      a_rep <- a_rep
      c_rep <- c_rep
      n1_rep <- n1_rep
      n2_rep <- n2_rep
    }


    #How many studies are reported?
    N_rep <- length(Rep_index)

    #Unreported study sizes, we might have info from n1,n2 or just the total
    ntot <- as.numeric(n1) + as.numeric(n2)
    n_HR <- as.numeric(ntot[HR_index])


    #logRR and standard error needed to evaluate if the study is significant or not
    #Outcome is significant if |logRR| > z_alpha*sigma
    logRR <- log((a_rep*n2_rep)/(c_rep*n1_rep))
    s <- sqrt(((n1_rep - a_rep)/(n1_rep*a_rep)) + ((n2_rep - c_rep)/(n2_rep*c_rep)))
    sigma_squared <- s^2



    #Imputed values of sigma squared for the unreported studies
    #K value based on the reported studies
    k <- sum(1/sigma_squared)/sum((n1_rep + n2_rep))

  } else if (!is.null(mu1) & !is.null(mu2) & !is.null(sd1) & !is.null(sd2)){


    Rep_index <- which(!is.na(as.numeric(mu1)))

    if (low.risk == FALSE){

      HR_index <- which(mu1 == "high")

    } else {

      HR_index <- which(mu1 == "high" | mu1 == "low")

    }

    #Unreported study sizes, we might have info from n1,n2 or just the total
    ntot <- as.numeric(n1) +as.numeric(n2)
    n_HR <- as.numeric(ntot[HR_index])

    #mu1,mu2,n1,n2,n values for the reported studies
    mu1_rep <- as.numeric(mu1[Rep_index])
    mu2_rep <- as.numeric(mu2[Rep_index])
    sd1_rep <- as.numeric(sd1[Rep_index])
    sd2_rep <- as.numeric(sd2[Rep_index])
    n1_rep <- as.numeric(n1[Rep_index])
    n2_rep <- as.numeric(n2[Rep_index])

    #How many studies are reported?
    N_rep <- length(Rep_index)

    #Differenece in means
    #Standard errors given
    logRR <- mu1_rep - mu2_rep
    s <- sqrt((as.numeric(sd1_rep)^2)/(as.numeric(n1_rep)) + (as.numeric(sd2_rep)^2)/(as.numeric(n2_rep)))
    sigma_squared <- s^2

    k <- sum(1/sigma_squared)/sum((n1_rep + n2_rep))

  } else if (!is.null(y) & !is.null(s)) {



    #Indecies where we have the reported outcomes and the unreported with high risk (HR)
    Rep_index <- which(!is.na(as.numeric(y)))

    if (low.risk == FALSE){

      HR_index <- which(y == "high")
    } else {

      HR_index <- which(y == "high" | y == "low")
    }

    #Observed treatment effect, standard error^2, and sample sizes of reported outcomes
    logRR <- as.numeric(y[Rep_index])
    sigma_squared <- (as.numeric(s[Rep_index]))^2
    n1_rep <- as.numeric(n1[Rep_index])
    n2_rep <- as.numeric(n2[Rep_index])


    #Total sample size of unreported outcomes
    n_HR <- as.numeric(n1[HR_index]) + as.numeric(n2[HR_index])

    #Copas et al. (2014, 2019) method for imputation of missing variances

    #k estimation, based on reported studies
    k <- sum(1/sigma_squared)/sum((n1_rep + n2_rep))



  } else {

    return("Error: invalid inputs. Input either a and c values or y1,sd1 and y2,sd2 values.")
  }


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

  #Unadjusted log-likelihood function to be maximized
  f.u <- function(params, logRR, sigma_squared) {
    mu <- params[1]
    tau_squared <- params[2]
    -(1/2)*sum(log(sigma_squared + tau_squared) + ((logRR - mu)^2)/(sigma_squared + tau_squared))
  }

  #Set initial values for mu and tau_squared
  init_params <- init_param

  #Maximize unadjusted function optim()
  if (method == "L-BFGS-B"){

    fit.u <- optim(init_params, f.u, logRR = logRR, sigma_squared = sigma_squared,
                   #method = "Nelder-Mead",
                   method=method,
                   lower = lower,
                   upper = upper,
                   control = list(fnscale = -1),
                   hessian=my.hessian)
  } else {

    fit.u <- optim(init_params, f.u, logRR = logRR, sigma_squared = sigma_squared,
                   method = method,
                   control = list(fnscale = -1),
                   hessian=my.hessian)

  }

  #Return unadjusted mu and tau_squared
  mle.u <- fit.u$par[1]
  mle.tau <- max(fit.u$par[2],0)


  #Beneficial outcome adjustment for ORB

  if(outcome == "benefit"){

    z_alpha.copas <- qnorm(1-alpha_ben/2) #when using the copas.twosides it's always two sided!

    if (alpha_ben_one.sided == TRUE){ #for the other selection functions can choose

      z_alpha <- qnorm(1-alpha_ben)
    } else {
      z_alpha <- qnorm(1-alpha_ben/2)
    }

    #Adjusted log-likelihood function for beneficial outcome to be maximized
    f.adj.b <- function(params, logRR, sigma_squared, sigma_squared_imputed) {
      mu <- params[1]
      tau_squared <- params[2]



      #The contribution from the reported studies is always present
      -(1/2)*sum(log(sigma_squared + tau_squared) + ((logRR - mu)^2)/(sigma_squared + tau_squared)) +


        if (length(sigma_squared_imputed) > 0) {

          if (sel.ben == "Copas.oneside"){

            sum(log(pnorm((z_alpha*sqrt(sigma_squared_imputed) - mu)/sqrt(sigma_squared_imputed + tau_squared))))

          } else if (sel.ben == "Copas.twoside") {

            sum(log(pnorm((z_alpha.copas*sqrt(sigma_squared_imputed) - mu)/sqrt(sigma_squared_imputed + tau_squared))
                    - pnorm((-z_alpha.copas*sqrt(sigma_squared_imputed) - mu)/sqrt(sigma_squared_imputed + tau_squared))))

          } else if (sel.ben == "Neg.exp.piecewise") {

            #Probability of unreporting
            piece_wise_neg_exp <- function(y, sigma_squared_imputed, z_alpha) {
              ifelse(y < sqrt(sigma_squared_imputed)*z_alpha,
                     1,
                     exp(-w * (y - z_alpha*sqrt(sigma_squared_imputed))))
            }

            #Integrand
            integrand <- function(y, mu, tau_squared, sigma_squared_imputed, z_alpha) {

              weight <- piece_wise_neg_exp(y, sigma_squared_imputed, z_alpha)
              (1 / sqrt(2 * pi * (sigma_squared_imputed + tau_squared))) * exp(-0.5 * ((y - mu)^2 / (sigma_squared_imputed + tau_squared))) * weight
            }

            #log-lik contribution
            sum(log(sapply(sigma_squared_imputed, function(sigma_sq_imputed) {
              integrate(function(y) integrand(y, mu, tau_squared, sigma_sq_imputed, z_alpha),
                        lower = -10, upper = 10, subdivisions = 200, stop.on.error=FALSE)$value
            })))

          } else if (sel.ben == "Sigmoid.cont") {

            #Probability of unreporting
            cont_sigmoid <- function(y, sigma_squared_imputed) {

              1/(1+exp(w*(y-z_alpha*sqrt(sigma_squared_imputed))))

            }

            #Integrand
            integrand <- function(y, mu, tau_squared, sigma_squared_imputed) {

              weight <- cont_sigmoid(y, sigma_squared_imputed)
              (1 / sqrt(2 * pi * (sigma_squared_imputed + tau_squared))) * exp(-0.5 * ((y - mu)^2 / (sigma_squared_imputed + tau_squared))) * weight
            }

            #log-lik contribution
            sum(log(sapply(sigma_squared_imputed, function(sigma_sq_imputed) {
              integrate(function(y) integrand(y, mu, tau_squared, sigma_sq_imputed),
                        lower = -10, upper = 10, subdivisions = 200, stop.on.error=FALSE)$value
            })))

          } else {

            stop("Error: Invalid weight function specified. ")
          }



        } else {
          0  # Return 0 if sigma_squared_imputed is empty
        }


    }

    #Maximize log-likelihood
    if (method == "L-BFGS-B"){

      fit.adj.b <- optim(init_params, f.adj.b, logRR = logRR, sigma_squared = sigma_squared, sigma_squared_imputed = sigma_squared_imputed,
                         method = method,
                         lower = lower,
                         upper = upper,
                         control = list(fnscale = -1),

                         hessian=my.hessian)
    } else {

      fit.adj.b <- optim(init_params, f.adj.b, logRR = logRR, sigma_squared = sigma_squared, sigma_squared_imputed = sigma_squared_imputed,
                         method = method,
                         control = list(fnscale = -1),
                         hessian=my.hessian)


    }

    #Return adjusted mu and tau_squared
    mle.b <- fit.adj.b$par[1]
    mle.b.tau <- max(fit.adj.b$par[2],0)


    #REML ESTIMATION OF TAU SQUARED UNADJUSTED AND ADJUSTED

    #unadjusted REML likelihood
    lik.REML <- function(tau_squared, sigma_squared, logRR ){

      w_i <- 1/(sigma_squared + tau_squared)
      mu.REML <- sum(w_i*logRR)/sum(w_i)


      (1/2)*sum(log(w_i)) -(1/2)*log(sum(w_i)) - (1/2)*sum(w_i*((logRR - mu.REML)^2))

    }

    init <- init_param[2]

    fit.lik.REML <- optim(init, lik.REML, sigma_squared=sigma_squared, logRR=logRR,
                          method= "Brent",
                          lower = 0, upper=10,
                          control = list(fnscale = -1))

    tau.REML <- max(fit.lik.REML$par,0)


    #Adjusted REML function benefit
    lik.adj.REML <- function(tau_squared, sigma_squared, sigma_squared_imputed, logRR ){

      w_i <- 1/(sigma_squared + tau_squared)
      mu.REML <- sum(w_i*logRR)/sum(w_i)



      (1/2)*sum(log(w_i)) -(1/2)*log(sum(w_i)) - (1/2)*sum(w_i*((logRR - mu.REML)^2)) +

        if (length(sigma_squared_imputed) > 0) {

          if (sel.ben == "Copas.oneside"){

            sum(log(pnorm((z_alpha*sqrt(sigma_squared_imputed) - mu.REML)/sqrt(sigma_squared_imputed + tau_squared))))

          } else if (sel.ben == "Copas.twoside") {

            sum(log(pnorm((z_alpha.copas*sqrt(sigma_squared_imputed) - mu.REML)/sqrt(sigma_squared_imputed + tau_squared))
                    - pnorm((-z_alpha.copas*sqrt(sigma_squared_imputed) - mu.REML)/sqrt(sigma_squared_imputed + tau_squared))))

          } else if (sel.ben == "Neg.exp.piecewise") {

            #Probability of unreporting
            piece_wise_neg_exp <- function(y, sigma_squared_imputed, z_alpha) {
              ifelse(y < sqrt(sigma_squared_imputed)*z_alpha,
                     1,
                     exp(-w * (y - z_alpha*sqrt(sigma_squared_imputed))))
            }

            #Integrand
            integrand <- function(y, mu.REML, tau_squared, sigma_squared_imputed, z_alpha) {

              weight <- piece_wise_neg_exp(y, sigma_squared_imputed, z_alpha)
              (1 / sqrt(2 * pi * (sigma_squared_imputed + tau_squared))) * exp(-0.5 * ((y - mu.REML)^2 / (sigma_squared_imputed + tau_squared))) * weight
            }

            #log-lik contribution
            sum(log(sapply(sigma_squared_imputed, function(sigma_sq_imputed) {
              integrate(function(y) integrand(y, mu.REML, tau_squared, sigma_sq_imputed, z_alpha),
                        lower = -20, upper = 20, subdivisions = 200, stop.on.error=FALSE)$value
            })))

          } else if (sel.ben == "Sigmoid.cont") {

            #Probability of unreporting
            cont_sigmoid <- function(y, sigma_squared_imputed, z_alpha) {
              ifelse(y < sqrt(sigma_squared_imputed)*z_alpha,
                     1/(1+exp(w*(y-z_alpha*sqrt(sigma_squared_imputed)))),
                     1/(1+exp(w*(y-z_alpha*sqrt(sigma_squared_imputed)))))

            }

            #Integrand
            integrand <- function(y, mu.REML, tau_squared, sigma_squared_imputed, z_alpha) {

              weight <- cont_sigmoid(y, sigma_squared_imputed, z_alpha)
              (1 / sqrt(2 * pi * (sigma_squared_imputed + tau_squared))) * exp(-0.5 * ((y - mu.REML)^2 / (sigma_squared_imputed + tau_squared))) * weight
            }

            #log-lik contribution
            sum(log(sapply(sigma_squared_imputed, function(sigma_sq_imputed) {
              integrate(function(y) integrand(y, mu.REML, tau_squared, sigma_sq_imputed, z_alpha),
                        lower = -20, upper = 20, subdivisions = 200, stop.on.error=FALSE)$value
            })))

          } else {

            stop("Error: Invalid weight function specified. ")
          }



        } else {
          0  # Return 0 if sigma_squared_imputed is empty
        }


    }


    fit.lik.adj.REML <- optim(init, lik.adj.REML, sigma_squared=sigma_squared, sigma_squared_imputed= sigma_squared_imputed, logRR=logRR,
                              method= "Brent",
                              lower = 0, upper=10,
                              control = list(fnscale = -1))

    tau.adj.REML <- max(fit.lik.adj.REML$par,0)



    #LIKELIHOOD RATIO CONFIDENCE INTERVALS

    if (LR.CI){

      #Unadjusted
      z <- qchisq(1-alpha_ben, df=1) #3.841

      #Re-write log likelihood with two inputs instead of param = c(mu, tau_squared)
      ll.u <- function(mu, tau_squared, logRR, sigma_squared) {

        -(1/2)*sum(log(sigma_squared + tau_squared) + ((logRR - mu)^2)/(sigma_squared + tau_squared))
      }

      #Profile log-likelihood
      pl.u <- function(mu, logRR, sigma_squared) { #take in vector of mus

        res <- mu

        for (i in seq_along(mu)) { #for all these values of mu
          optimResult <- optim(par = init_param[2],
                               fn = function(tau_squared) ll.u(mu[i], tau_squared, logRR=logRR, sigma_squared = sigma_squared),
                               method = "Brent",
                               lower= 0.0001,
                               upper=10,
                               control = list(fnscale = -1))

          res[i] <- optimResult$value
        }
        return(res)
      }

      f <- function(mu, logRR, sigma_squared){
        pl.u(mu, logRR=logRR, sigma_squared=sigma_squared) - pl.u(mle.u, logRR=logRR, sigma_squared=sigma_squared) + 1/2*qchisq(0.95, df=1)
      }

      eps <- sqrt(.Machine$double.eps)
      lowerBound.u <- uniroot(f, interval = c(-10, mle.u), logRR=logRR, sigma_squared=sigma_squared)$root
      upperBound.u <- uniroot(f, interval = c( mle.u, 10), logRR=logRR, sigma_squared=sigma_squared)$root


      #Adjusted benefit

      #Re-write log likelihood with two inputs instead of param = c(mu, tau_squared)
      ll.b <- function(mu, tau_squared, logRR, sigma_squared, sigma_squared_imputed) {



        #The contribution from the reported studies is always present
        -(1/2)*sum(log(sigma_squared + tau_squared) + ((logRR - mu)^2)/(sigma_squared + tau_squared)) +


          if (length(sigma_squared_imputed) > 0) {

            if (sel.ben == "Copas.oneside"){

              sum(log(pnorm((z_alpha*sqrt(sigma_squared_imputed) - mu)/sqrt(sigma_squared_imputed + tau_squared))))

            } else if (sel.ben == "Copas.twoside") {

              sum(log(pnorm((z_alpha.copas*sqrt(sigma_squared_imputed) - mu)/sqrt(sigma_squared_imputed + tau_squared))
                      - pnorm((-z_alpha.copas*sqrt(sigma_squared_imputed) - mu)/sqrt(sigma_squared_imputed + tau_squared))))

            } else if (sel.ben == "Neg.exp.piecewise") {

              #Probability of unreporting
              piece_wise_neg_exp <- function(y, sigma_squared_imputed, z_alpha) {
                ifelse(y < sqrt(sigma_squared_imputed)*z_alpha,
                       1,
                       exp(-w * (y - z_alpha*sqrt(sigma_squared_imputed))))
              }

              #Integrand
              integrand <- function(y, mu, tau_squared, sigma_squared_imputed, z_alpha) {

                weight <- piece_wise_neg_exp(y, sigma_squared_imputed, z_alpha)
                (1 / sqrt(2 * pi * (sigma_squared_imputed + tau_squared))) * exp(-0.5 * ((y - mu)^2 / (sigma_squared_imputed + tau_squared))) * weight
              }

              #log-lik contribution
              sum(log(sapply(sigma_squared_imputed, function(sigma_sq_imputed) {
                integrate(function(y) integrand(y, mu, tau_squared, sigma_sq_imputed, z_alpha),
                          lower = -20, upper = 20, subdivisions = 200, stop.on.error=FALSE)$value
              })))

            } else if (sel.ben == "Sigmoid.cont") {

              #Probability of unreporting
              cont_sigmoid <- function(y, sigma_squared_imputed, z_alpha) {
                ifelse(y < sqrt(sigma_squared_imputed)*z_alpha,
                       1/(1+exp(w*(y-z_alpha*sqrt(sigma_squared_imputed)))),
                       1/(1+exp(w*(y-z_alpha*sqrt(sigma_squared_imputed)))))

              }

              #Integrand
              integrand <- function(y, mu, tau_squared, sigma_squared_imputed, z_alpha) {

                weight <- cont_sigmoid(y, sigma_squared_imputed, z_alpha)
                (1 / sqrt(2 * pi * (sigma_squared_imputed + tau_squared))) * exp(-0.5 * ((y - mu)^2 / (sigma_squared_imputed + tau_squared))) * weight
              }

              #log-lik contribution
              sum(log(sapply(sigma_squared_imputed, function(sigma_sq_imputed) {
                integrate(function(y) integrand(y, mu, tau_squared, sigma_sq_imputed, z_alpha),
                          lower = -20, upper = 20, subdivisions = 200, stop.on.error=FALSE)$value
              })))

            } else {

              stop("Error: Invalid weight function specified. ")
            }



          } else {
            0  # Return 0 if sigma_squared_imputed is empty
          }



      }


      pl.b <- function(mu, logRR, sigma_squared, sigma_squared_imputed) { #take in vector of mus

        res <- mu

        for (i in seq_along(mu)) { #for all these values of mu
          optimResult <- optim(par = init_param[2],
                               fn = function(tau_squared) ll.b(mu[i], tau_squared,
                                                               logRR=logRR,
                                                               sigma_squared = sigma_squared,
                                                               sigma_squared_imputed = sigma_squared_imputed),
                               method = "Brent",
                               lower=0.00001,
                               upper=10,

                               control = list(fnscale = -1))

          res[i] <- optimResult$value
        }
        return(res)
      }

      f.b <- function(mu, logRR, sigma_squared, sigma_squared_imputed){

        pl.b(mu, logRR=logRR, sigma_squared=sigma_squared, sigma_squared_imputed = sigma_squared_imputed) - pl.b(mle.b, logRR=logRR, sigma_squared=sigma_squared, sigma_squared_imputed = sigma_squared_imputed) + 1/2*qchisq(0.95, df=1)
      }

      lowerBound.b <- uniroot(f.b, interval = c(-010, mle.b), logRR=logRR,
                              sigma_squared=sigma_squared,
                              sigma_squared_imputed=sigma_squared_imputed)$root
      upperBound.b <- uniroot(f.b, interval = c(mle.b, 10), logRR=logRR,
                              sigma_squared=sigma_squared,
                              sigma_squared_imputed = sigma_squared_imputed)$root



      if (Wald.CI){


        #WALD CONFIDENCE INTERVALS
        a <- alpha_ben #for harm Copas et al use 99% conf level
        #Unadjusted
        fisher_info.u <- solve(-fit.u$hessian)
        s.u <- sqrt(diag(fisher_info.u)[1])
        ci.u <- fit.u$par[1] + qnorm(c(a/2, 1-a/2)) * s.u
        #Adjusted benefit
        fisher_info.adj.b <- solve(-fit.adj.b$hessian)
        s.adj.b <- sqrt(diag(fisher_info.adj.b)[1])
        ci.u.adj.b <- fit.adj.b$par[1] + qnorm(c(a/2, 1-a/2)) * s.adj.b



        return(list(mu_unadjusted = mle.u,
                    LR_mu_unadjusted_low = lowerBound.u,
                    LR_mu_unadjusted_up = upperBound.u,


                    CI_unadjusted_low_WALD = ci.u[1],
                    CI_unadjusted_up_WALD = ci.u[2],

                    mu_adjusted_benefit = mle.b,
                    LR_mu_adjusted_low = lowerBound.b,
                    LR_mu_adjusted_up = upperBound.b,

                    tau_squared_unadjusted = mle.tau,
                    tau_squared_adjusted = mle.b.tau,

                    tau_squared_unadjusted_REML = tau.REML,
                    tau_squared_adjusted_REML = tau.adj.REML,

                    average_sigma_squared_unadjusted = sigma_squared_average_unadjusted,
                    sigma_squared_average_adjusted = sigma_squared_average_adjusted,


                    CI_adjusted_benefit_low_WALD = ci.u.adj.b[1],
                    CI_adjusted_benefit_up_WALD = ci.u.adj.b[2]

        ))

      } else {

        return(list(mu_unadjusted = mle.u,
                    LR_mu_unadjusted_low = lowerBound.u,
                    LR_mu_unadjusted_up = upperBound.u,



                    mu_adjusted_benefit = mle.b,
                    LR_mu_adjusted_low = lowerBound.b,
                    LR_mu_adjusted_up = upperBound.b,

                    tau_squared_unadjusted = mle.tau,
                    tau_squared_adjusted = mle.b.tau,

                    tau_squared_unadjusted_REML = tau.REML,
                    tau_squared_adjusted_REML = tau.adj.REML,

                    average_sigma_squared_unadjusted = sigma_squared_average_unadjusted,
                    sigma_squared_average_adjusted = sigma_squared_average_adjusted



        ))


      }


    } else {



      return(list(
        mu_unadjusted = mle.u,
        mu_adjusted_benefit = mle.b,
        tau_squared_unadjusted = mle.tau,
        tau_squared_adjusted = mle.b.tau,
        tau_squared_unadjusted_REML = tau.REML,
        tau_squared_adjusted_REML = tau.adj.REML,
        average_sigma_squared_unadjusted = sigma_squared_average_unadjusted,
        sigma_squared_average_adjusted = sigma_squared_average_adjusted
      ))


    }


    #Adjustment for ORB in harmful outcome
  } else if (outcome == "harm"){

    #Adjusted log-likelihood function for harmful outcome to be maximized
    f.adj.h <- function(params, logRR, sigma_squared, sigma_squared_imputed) {
      mu <- params[1]
      tau_squared <- params[2]

      -(1/2)*sum(log(sigma_squared + tau_squared) + ((logRR - mu)^2)/(sigma_squared + tau_squared)) +

        if (length(sigma_squared_imputed)>0){

          sum(log(pnorm(mu/sqrt(sigma_squared_imputed + tau_squared))))

        } else {

          0
        }
    }

    #Maximize log-likelihoood
    if (method == "L-BFGS-B"){
      fit.adj.h <- optim(init_params, f.adj.h, logRR = logRR, sigma_squared = sigma_squared, sigma_squared_imputed = sigma_squared_imputed,
                         method = method,
                         lower = lower,
                         upper = upper,
                         control = list(fnscale = -1),

                         hessian = my.hessian)
    } else {

      fit.adj.h <- optim(init_params, f.adj.h, logRR = logRR, sigma_squared = sigma_squared, sigma_squared_imputed = sigma_squared_imputed,
                         method = method,
                         control = list(fnscale = -1),

                         hessian = my.hessian)
    }

    #Adjusted mu and tau_squared estimates
    mle.h <- fit.adj.h$par[1]
    mle.h.tau <- max(fit.adj.h$par[2],0)



    #LIKELIHOOD RATIO CONFIDENCE INTERVALS

    if (LR.CI){


      #Unadjusted
      z <- qchisq(1-alpha_harm, df=1) #3.841

      #Re-write log likelihood with two inputs instead of param = c(mu, tau_squared)
      ll.u <- function(mu, tau_squared, logRR, sigma_squared) {

        -(1/2)*sum(log(sigma_squared + tau_squared) + ((logRR - mu)^2)/(sigma_squared + tau_squared))
      }


      pl.u <- function(mu, logRR, sigma_squared) { #take in vector of mus

        res <- mu

        for (i in seq_along(mu)) { #for all these values of mu
          optimResult <- optim(par = init_param[2],
                               fn = function(tau_squared) ll.u(mu[i], tau_squared, logRR=logRR, sigma_squared = sigma_squared),
                               method = "Brent",
                               lower= 0.00001,
                               upper=10,
                               control = list(fnscale = -1))

          res[i] <- optimResult$value
        }
        return(res)
      }

      f <- function(mu, logRR, sigma_squared){
        pl.u(mu, logRR=logRR, sigma_squared=sigma_squared) - pl.u(mle.u, logRR=logRR, sigma_squared=sigma_squared) + 1/2*qchisq(0.95, df=1)
      }

      lowerBound.u <- uniroot(f, interval = c(-10, mle.u), logRR=logRR, sigma_squared=sigma_squared)$root
      upperBound.u <- uniroot(f, interval = c(mle.u, 10), logRR=logRR, sigma_squared=sigma_squared)$root


      ll.h <- function(mu, tau_squared, logRR, sigma_squared, sigma_squared_imputed) {

        z_alpha <- qnorm(1 - alpha_harm/2)
        -(1/2)*sum(log(sigma_squared + tau_squared) + ((logRR - mu)^2)/(sigma_squared + tau_squared)) +
          sum(log(pnorm(mu/sqrt(sigma_squared_imputed + tau_squared))))
      }


      pl.h <- function(mu, logRR, sigma_squared, sigma_squared_imputed) { #take in vector of mus

        res <- mu

        for (i in seq_along(mu)) { #for all these values of mu
          optimResult <- optim(par = init_param[2],
                               fn = function(tau_squared) ll.h(mu[i], tau_squared, logRR=logRR, sigma_squared = sigma_squared, sigma_squared_imputed = sigma_squared_imputed),
                               method = "Brent",
                               lower= 0.00001,
                               upper = 10,
                               control = list(fnscale = -1))

          res[i] <- optimResult$value
        }
        return(res)
      }

      f.h <- function(mu, logRR, sigma_squared, sigma_squared_imputed){
        pl.h(mu, logRR=logRR, sigma_squared=sigma_squared, sigma_squared_imputed = sigma_squared_imputed) - pl.h(mle.h, logRR=logRR, sigma_squared=sigma_squared, sigma_squared_imputed = sigma_squared_imputed) + 1/2*qchisq(0.99, df=1)
      }

      lowerBound.h <- uniroot(f.h, interval = c(-10, mle.h), logRR=logRR, sigma_squared=sigma_squared, sigma_squared_imputed = sigma_squared_imputed)$root
      upperBound.h <- uniroot(f.h, interval = c(mle.h, 10), logRR=logRR, sigma_squared=sigma_squared, sigma_squared_imputed = sigma_squared_imputed)$root





      if (Wald.CI){
        #WALD CONFIDENCE INTERVALS
        a <- alpha_harm #for harm Copas et al use 99% conf level
        #Unadjusted
        fisher_info.u <- solve(-fit.u$hessian)
        s.u <- sqrt(diag(fisher_info.u)[1])
        ci.u <- fit.u$par[1] + qnorm(c(a/2, 1-a/2)) * s.u
        #Adjusted harm
        fisher_info.adj.h <- solve(-fit.adj.h$hessian)
        s.adj.h <- sqrt(diag(fisher_info.adj.h)[1])
        ci.u.adj.h <- fit.adj.h$par[1] + qnorm(c(a/2, 1-a/2)) * s.adj.h



        return(list(mu_unadjusted = mle.u,
                    LR_mu_unadjusted_low = lowerBound.u,
                    LR_mu_unadjusted_up = upperBound.u,


                    CI_unadjusted_low_WALD = ci.u[1],
                    CI_unadjusted_up_WALD = ci.u[2],

                    mu_adjusted_harm = mle.h,
                    LR_mu_adjusted_low = lowerBound.h,
                    LR_mu_adjusted_up = upperBound.h,


                    tau_squared_unadjsuted = mle.tau,
                    tau_squared_adjusted = mle.h.tau,

                    average_sigma_squared_unadjusted = sigma_squared_average_unadjusted,
                    sigma_squared_average_adjusted = sigma_squared_average_adjusted,


                    CI_adjusted_harm_low_WALD = ci.u.adj.h[1],
                    CI_adjusted_harm_up_WALD = ci.u.adj.h[2]

        ))

      } else {


        return(list(mu_unadjusted = mle.u,
                    LR_mu_unadjusted_low = lowerBound.u,
                    LR_mu_unadjusted_up = upperBound.y,

                    mu_adjusted_harm = mle.h,
                    LR_mu_adjusted_low = lowerBound.h,
                    LR_mu_adjusted_up = upperBound.h,


                    tau_squared_unadjsuted = mle.tau,
                    tau_squared_adjusted = mle.h.tau,

                    average_sigma_squared_unadjusted = sigma_squared_average_unadjusted,
                    sigma_squared_average_adjusted = sigma_squared_average_adjusted



        ))

      }

    }else {

      return(list(
        mu_unadjusted = mle.u,
        mu_adjusted_harm = mle.h,
        tau_squared_unadjusted = mle.tau,
        tau_squared_adjusted = mle.h.tau,

        average_sigma_squared_unadjusted = sigma_squared_average_unadjusted,
        sigma_squared_average_adjusted = sigma_squared_average_adjusted
      ))
    }


  } else {

    return("invalid outcome input")
  }


}




