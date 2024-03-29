
#' @export

#If tau squared is known

reORBgen_hetfixed <- function(a=NULL, c=NULL,
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
                              opt_method="Brent",
                              lower = -10,
                              upper = 10,
                              low.risk = FALSE,
                              tau_squared_fixed
) {


  #Tau squared is now a known and fixed parameter and we only optimize with respect to mu
  #So all optimization problems become one dimensional
  #LR intervals instead of PL
  tau_squared <- tau_squared_fixed

  #Default is to ignore the low risk of bias, but if TRUE then we
  #use all the unreported studies fro adjustment
  #this supposedly can be done for the two smooth selection functions

  #Can take in binary counts and calculate log RR: a, c
  #Means and calculate differences in means: y1, y2, sd1, sd2
  #Directly normally distributed effect measure: y, s


  #Beneficial outcome
  #In the original paper by Copas et al. (2019), they use a piece-wise selection function
  #which has P(unreporting)=1 if not-significant, where they consider two tails
  #however, more realistic that, if significant in the "wrong" direction,
  #even less likely to be reported! hence Copas.oneside = default in which we
  #don't consider the significance in the "wrong direction"

  #Selection functions:
  #Copas.twoside from original paper here always two sided by nature of the function
  #Copas.oneside DEFAULT
  #Neg.exp.piecewise
  #Sigmoid.cont
  #Brent optimization because one dimensional

  #If alpha_ben_one.sided=TRUE, apart from copas.twosides, which is always two sided by nature,
  #we use the one-sided threshold z=1.64
  #Preferred to used alpha_ben_one.sided=FALSE for all selection functions, so z=1.96

  sel.ben <- selection.benefit

  #If we have Neg.exp.piecewise or Sigmoid.cont, we also need a weight parameter
  w <- weight_param

  if (Wald.CI){
    my.hessian <- TRUE
  } else{
    my.hessian <- FALSE
  }

  method <- opt_method
  lower <- lower
  upper <- upper



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

      a_0_both <- which(a_rep == 0.5 & c_rep == 0.5) #If both zero, remove the study!

      if (length(a_0_both)>0){

        a_rep <- a_rep[-a_0_both]
        c_rep <- c_rep[-a_0_both]
        n1_rep <- n1_rep[-a_0_both]
        n2_rep <- n2_rep[-a_0_both]
      } else {

        a_rep <- a_rep
        c_rep <- c_rep
        n1_rep <- n1_rep
        n2_rep <- n2_rep

      }

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
  f.u <- function(mu, logRR, sigma_squared, tau_squared) {

    -(1/2)*sum(log(sigma_squared + tau_squared) + ((logRR - mu)^2)/(sigma_squared + tau_squared))
  }

  #Set initial values for mu and tau_squared
  init_params <- init_param

  #Maximize unadjusted function optim()

  fit.u <- optim(init_params, f.u, logRR = logRR, sigma_squared = sigma_squared, tau_squared = tau_squared,

                 method=method,
                 lower = lower,
                 upper = upper,
                 control = list(fnscale = -1),
                 hessian=my.hessian)


  #Return unadjusted mu and tau_squared
  mle.u <- fit.u$par


  #Beneficial outcome adjustment for ORB

  if(outcome == "benefit"){

    z_alpha.copas <- qnorm(1-alpha_ben/2) #when using the copas.twosides it's always two sided!

    if (alpha_ben_one.sided == TRUE){ #for the other selection functions can choose

      z_alpha <- qnorm(1-alpha_ben)
    } else {
      z_alpha <- qnorm(1-alpha_ben/2)
    }

    #Adjusted log-likelihood function for beneficial outcome to be maximized
    f.adj.b <- function(mu, logRR, sigma_squared, sigma_squared_imputed, tau_squared) {

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

    fit.adj.b <- optim(init_params, f.adj.b, logRR = logRR, sigma_squared = sigma_squared, sigma_squared_imputed = sigma_squared_imputed, tau_squared=tau_squared,
                       method = method,
                       lower = lower,
                       upper = upper,
                       control = list(fnscale = -1),

                       hessian=my.hessian)


    #Return adjusted mu and tau_squared
    mle.b <- fit.adj.b$par



    #LIKELIHOOD RATIO CONFIDENCE INTERVALS

    if (LR.CI){

      #Unadjusted
      z <- qchisq(1-alpha_ben, df=1) #3.841

      LRstat.u <- function(mu){

        2*(f.u(mle.u, logRR, sigma_squared, tau_squared) - f.u(mu, logRR, sigma_squared, tau_squared))

      }

      #From prof Held's book likelihood inference
      LR.lower.u <- uniroot(Vectorize(function(mu){LRstat.u(mu) - z}),
                            interval=c(-10,mle.u))$root

      LR.upper.u <- uniroot(Vectorize(function(mu){LRstat.u(mu) - z}),
                            interval=c(mle.u, 10))$root

      #Adjusted benefit

      #Re-write log likelihood with two inputs instead of param = c(mu, tau_squared)
      LRstat.b <- function(mu){

        2*(f.adj.b(mle.b, logRR, sigma_squared, sigma_squared_imputed, tau_squared) - f.adj.b(mu, logRR, sigma_squared, sigma_squared_imputed, tau_squared))

      }

      LR.lower.b <- uniroot(Vectorize(function(mu){LRstat.b(mu) - z}),
                            interval=c(-10,mle.b))$root

      LR.upper.b <- uniroot(Vectorize(function(mu){LRstat.b(mu) - z}),
                            interval=c(mle.b, 10))$root

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
                    LR_mu_unadjusted_low = LR.lower.u,
                    LR_mu_unadjusted_up = LR.upper.u,


                    CI_unadjusted_low_WALD = ci.u[1],
                    CI_unadjusted_up_WALD = ci.u[2],

                    mu_adjusted_benefit = mle.b,
                    LR_mu_adjusted_low = LR.lower.b,
                    LR_mu_adjusted_up = LR.upper.b,


                    average_sigma_squared_unadjusted = sigma_squared_average_unadjusted,
                    sigma_squared_average_adjusted = sigma_squared_average_adjusted,


                    CI_adjusted_benefit_low_WALD = ci.u.adj.b[1],
                    CI_adjusted_benefit_up_WALD = ci.u.adj.b[2]

        ))

      } else {

        return(list(mu_unadjusted = mle.u,
                    LR_mu_unadjusted_low = LR.lower.u,
                    LR_mu_unadjusted_up = LR.upper.u,



                    mu_adjusted_benefit = mle.b,
                    LR_mu_adjusted_low = LR.lower.b,
                    LR_mu_adjusted_up = LR.upper.b,


                    average_sigma_squared_unadjusted = sigma_squared_average_unadjusted,
                    sigma_squared_average_adjusted = sigma_squared_average_adjusted



        ))


      }


    } else {



      return(list(
        mu_unadjusted = mle.u,
        mu_adjusted_benefit = mle.b,

        average_sigma_squared_unadjusted = sigma_squared_average_unadjusted,
        sigma_squared_average_adjusted = sigma_squared_average_adjusted
      ))


    }


    #Adjustment for ORB in harmful outcome
  } else if (outcome == "harm"){

    #Adjusted log-likelihood function for harmful outcome to be maximized
    f.adj.h <- function(mu, logRR, sigma_squared, sigma_squared_imputed, tau_squared) {


      -(1/2)*sum(log(sigma_squared + tau_squared) + ((logRR - mu)^2)/(sigma_squared + tau_squared)) +

        if (length(sigma_squared_imputed)>0){

          sum(log(pnorm(mu/sqrt(sigma_squared_imputed + tau_squared))))

        } else {

          0
        }
    }

    #Maximize log-likelihoood

    fit.adj.h <- optim(init_params, f.adj.h, logRR = logRR, tau_squared=tau_squared, sigma_squared = sigma_squared, sigma_squared_imputed = sigma_squared_imputed,
                       method = method,
                       lower = lower,
                       upper = upper,
                       control = list(fnscale = -1),

                       hessian = my.hessian)


    #Adjusted mu and tau_squared estimates
    mle.h <- fit.adj.h$par




    #LIKELIHOOD RATIO CONFIDENCE INTERVALS

    if (LR.CI){


      #Unadjusted
      z <- qchisq(1-alpha_harm, df=1) #3.841

      #Unadjsuted
      LRstat.u <- function(mu){

        2*(f.u(mle.u, logRR, sigma_squared, tau_squared) - f.u(mu, logRR, sigma_squared, tau_squared))

      }

      #From prof Held's book likelihood inference
      LR.lower.u <- uniroot(Vectorize(function(mu){LRstat.u(mu) - z}),
                            interval=c(-10,mle.u))$root

      LR.upper.u <- uniroot(Vectorize(function(mu){LRstat.u(mu) - z}),
                            interval=c(mle.u, 10))$root

      #harm adjusted
      LRstat.h <- function(mu){

        2*(f.adj.h(mle.h, logRR, sigma_squared, sigma_squared_imputed, tau_squared) - f.adj.h(mu, logRR, sigma_squared, sigma_squared_imputed, tau_squared))

      }

      LR.lower.h <- uniroot(Vectorize(function(mu){LRstat.h(mu) - z}),
                            interval=c(-10,mle.u))$root

      LR.upper.h <- uniroot(Vectorize(function(mu){LRstat.h(mu) - z}),
                            interval=c(mle.u, 10))$root


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
                    LR_mu_unadjusted_low = LR.lower.u,
                    LR_mu_unadjusted_up = LR.upper.u,


                    CI_unadjusted_low_WALD = ci.u[1],
                    CI_unadjusted_up_WALD = ci.u[2],

                    mu_adjusted_harm = mle.h,
                    LR_mu_adjusted_low = LR.lower.h,
                    LR_mu_adjusted_up = LR.upper.h,

                    average_sigma_squared_unadjusted = sigma_squared_average_unadjusted,
                    sigma_squared_average_adjusted = sigma_squared_average_adjusted,


                    CI_adjusted_harm_low_WALD = ci.u.adj.h[1],
                    CI_adjusted_harm_up_WALD = ci.u.adj.h[2]

        ))

      } else {


        return(list(mu_unadjusted = mle.u,
                    LR_mu_unadjusted_low = LR.lower.u,
                    LR_mu_unadjusted_up = LR.upper.u,

                    mu_adjusted_harm = mle.h,
                    LR_mu_adjusted_low = LR.lower.h,
                    LR_mu_adjusted_up = LR.upper.h,


                    average_sigma_squared_unadjusted = sigma_squared_average_unadjusted,
                    sigma_squared_average_adjusted = sigma_squared_average_adjusted



        ))

      }

    }else {

      return(list(
        mu_unadjusted = mle.u,
        mu_adjusted_harm = mle.h,

        average_sigma_squared_unadjusted = sigma_squared_average_unadjusted,
        sigma_squared_average_adjusted = sigma_squared_average_adjusted
      ))
    }


  } else {

    return("invalid outcome input")
  }


}
