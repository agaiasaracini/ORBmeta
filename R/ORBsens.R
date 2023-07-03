

ORBsens <- function(p_vals=NULL,
                      eta_vals=NULL,
                      a=NULL, c=NULL,
                      mu1=NULL, mu2=NULL, sd1=NULL, sd2=NULL,
                      n1, n2,
                      ntot,
                      sign_level,
                      init_param,
                      outcome){

  p1 <- p_vals[1]
  p2 <- p_vals[2]
  eta1 <- eta_vals[1]
  eta2 <- eta_vals[2]


  #Significance
  z_alpha <- qnorm(1-sign_level/2)

  if (!is.null(a) & !is.null(c)){
  #Reported outcomes
  #Indecies where we don't have low/high, i.e., reported outcomes
  #If we turn C into numeric, the low/high (unreported) become NA
  #Where do we have low and were do we have high
  Rep_index <- which(!is.na(as.numeric(a)))
  HR_index <- which(a == "high")
  LR_index <- which(a == "low")

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

    a_0_both <- which(a_rep == 0.5 & c_rep == 0.5)

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
  n_HR <- as.numeric(ntot[HR_index])
  n_LR <- as.numeric(ntot[LR_index])

  #logRR and standard error needed to evaluate if the study is significant or not
  #Outcome is significant if |logRR| > z_alpha*sigma
  y <- log((a_rep*n2_rep)/(c_rep*n1_rep))
  s <- sqrt(((n1_rep - a_rep)/(n1_rep*a_rep)) + ((n2_rep - c_rep)/(n2_rep*c_rep)))
  s2 <- s^2

  #Check which are significant and which are not
  #we divide into the reported y (logRR) which are significant and those which are not
  y_S <- y[which(abs(y) > z_alpha*s)]
  y_s <- y[which(abs(y) <= z_alpha*s)]

  #Which are positive and which are negative? Needed for my adjustment
  y_P <- y[which(y>0)]
  y_p <- y[which(y<=0)]

  #How many are pos and neg?
  N_rep_P <- length(y_P)
  N_rep_p <- length(y_p)

  #How many studies are reported and non significant
  N_rep_s <- length(y_s)



  #Imputed values of sigma squared for the unreported studies
  #K value based on the reported studies
  k <- sum(1/(s2))/sum((n1_rep + n2_rep))

  } else if (!is.null(mu1) & !is.null(mu2) & !is.null(sd1) & !is.null(sd2)){


    Rep_index <- which(!is.na(as.numeric(mu1)))
    HR_index <- which(mu1 == "high")
    LR_index <- which(mu1 == "low")

    #Unreported study sizes, we might have info from n1,n2 or just the total
    n_HR <- as.numeric(ntot[HR_index])
    n_LR <- as.numeric(ntot[LR_index])

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
    y <- mu1_rep - mu2_rep
    s <- sqrt((as.numeric(sd1_rep)^2)/(as.numeric(n1_rep)) + (as.numeric(sd2_rep)^2)/(as.numeric(n2_rep)))
    s2 <- s^2

    #check which are significant and which are not
    #we divide into the reported y (logRR) which are significant and those which are not
    y_S <- y[which(abs(y) > z_alpha*s)]
    y_s <- y[which(abs(y) <= z_alpha*s)]

    N_rep_s <- length(y_s)

    #Which are positive and which are negative? Needed for my adjustment
    y_P <- y[which(y>0)]
    y_p <- y[which(y<=0)]

    #How many are pos and neg?
    N_rep_P <- length(y_P)
    N_rep_p <- length(y_p)

    k <- sum(1/(s2))/sum((n1_rep + n2_rep))

  } else {

    return("Error: invalid inputs. Input either a and c values or y1,sd1 and y2,sd2 values.")
  }

  #UNADJUSTED log-likelihood
  loglik.unadjusted <- function(param, y, s2, N_rep, N_rep_s){

    theta <- param[1]
    t2 <- param[2]
    alpha <- param[3]
    beta <- param[4]

    f <- (1/sqrt(2*pi*(s2 + t2)))*exp((-1/2)*(((y - theta)^2)/(s2 + t2))) #reported studies

    sum(log(f)) +  N_rep*log(alpha) + #REPORTED
      N_rep_s*log(beta) #REPORTED AND NOT SIGNIFICANT

  }

  init_param <- init_param

  fit.unadjusted <- optim(init_param, loglik.unadjusted,  y=y, s2=s2,
                          N_rep=N_rep, N_rep_s=N_rep_s,
                          #method = "BFGS",
                          method =  "L-BFGS-B",
                          lower = c(-10, 0.01, 0.01, 0.01),
                          upper = c(10, 1,1,1),
                          #method="Nelder-Mead",
                          control = list(fnscale = -1),
                          hessian=TRUE)


  mle.unadjusted <- fit.unadjusted$par[1]

  a <- sign_level

  #Unsdjusted
  fisher_info.unadjusted <- solve(-fit.unadjusted$hessian)
  s.unadjusted <- sqrt(diag(fisher_info.unadjusted)[1])
  ci.unadjusted <- fit.unadjusted$par[1] + qnorm(c(a/2, 1-a/2)) * s.unadjusted
  width.unadjusted <- abs(ci.unadjusted[1]-ci.unadjusted[2])

  #ADJUSTMENT

  if (outcome == "benefit"){

  if (length(HR_index)>0 & length(LR_index)>0){ #We have both HR and LR

    #Imputed variances for the HR studies
    s2_imp_HR <- 1/(k*n_HR)
    #Imputed variances for th LR studies
    s2_imp_LR <- 1/(k*n_LR)

  loglik <- function(param, p1, p2, y, s2, s2_imp_HR, s2_imp_LR, N_rep, N_rep_s){

    #p1, p2 parameters which are assumed to be known (we can do a sensitivity analysis)
    #y the logRR of the reported studies
    #s2 the sigma squared of the reported studies
    #s2_imp_HR, s2_imp_LR the imputed variances of the HR and LR studies
    #n_rep number of reportes studies
    #n_rep_s number of reported and non significant studies


    #parameters to be estimated with maximum likelihood
    theta <- param[1]
    t2 <- param[2]
    alpha <- param[3]
    beta <- param[4]

    f <- (1/sqrt(2*pi*(s2 + t2)))*exp((-1/2)*(((y - theta)^2)/(s2 + t2))) #reported studies

    Q_HR <- pnorm((z_alpha*sqrt(s2_imp_HR) - theta)/(sqrt(s2_imp_HR + t2))) -
      pnorm((-z_alpha*sqrt(s2_imp_HR) - theta)/(sqrt(s2_imp_HR + t2)))

    Q_LR <- pnorm((z_alpha*sqrt(s2_imp_LR) - theta)/(sqrt(s2_imp_LR + t2))) -
      pnorm((-z_alpha*sqrt(s2_imp_LR) -  theta)/(sqrt(s2_imp_LR + t2)))

    sum(log(f)) +
    N_rep*log(alpha) + #REPORTED
    N_rep_s*log(beta) +  #REPORTED AND NOT SIGNIFICANT
    sum(log(p1*alpha*(Q_HR)*(1-beta) + (1-p2)*(1-alpha))) + #HIGH RISK
    sum(log((1-p1)*alpha*(Q_LR)*(1-beta) + p2*(1-alpha))) #LOW RISK

  }

  init_param <- init_param

  fit <- optim(init_param, loglik, p1=p1, p2=p2, y=y, s2=s2, s2_imp_HR=s2_imp_HR, s2_imp_LR=s2_imp_LR, N_rep=N_rep, N_rep_s=N_rep_s,
               method =  "L-BFGS-B",
               lower = c(-10, 0.001, 0.001, 0.001),
               upper = c(10, 20,1,1),
               control = list(fnscale = -1),
               hessian=TRUE)

  mle <- fit$par[1]
  val <- fit$value


  } else if (length(HR_index)>0 & length(LR_index)==0){ #We have HR but not LR



  #Imputed variances for the HR studies
  s2_imp_HR <- 1/(k*n_HR)

  loglik <- function(param, p1, p2, y, s2, s2_imp_HR, N_rep, N_rep_s){

    #p1, p2 parameters which are assumed to be known (we can do a sensitivity analysis)
    #y the logRR of the reported studies
    #s2 the sigma squared of the reported studies
    #s2_imp_HR, s2_imp_LR the imputed variances of the HR and LR studies
    #n_rep number of reportes studies
    #n_rep_s number of reported and non significant studies


    #parameters to be estimated with maximum likelihood
    theta <- param[1]
    t2 <- param[2]
    alpha <- param[3]
    beta <- param[4]

    f <- (1/sqrt(2*pi*(s2 + t2)))*exp((-1/2)*(((y - theta)^2)/(s2 + t2))) #reported studies

    Q_HR <- pnorm((z_alpha*sqrt(s2_imp_HR) - theta)/(sqrt(s2_imp_HR + t2))) -
      pnorm((-z_alpha*sqrt(s2_imp_HR) - theta)/(sqrt(s2_imp_HR + t2)))


    sum(log(f)) +
      N_rep*log(alpha) + #REPORTED
      N_rep_s*log(beta) +  #REPORTED AND NOT SIGNIFICANT
      sum(log(p1*alpha*(Q_HR)*(1-beta) + (1-p2)*(1-alpha))) #HIGH RISK

  }

  init_param <- init_param

  fit <- optim(init_param, loglik, p1=p1, p2=p2, y=y, s2=s2, s2_imp_HR=s2_imp_HR, N_rep=N_rep, N_rep_s=N_rep_s,
               method =  "L-BFGS-B",
               control = list(fnscale = -1),
               lower = c(-5, 0.05, 0.05, 0.05),
               upper = c(10, 10, 1, 1),
               hessian=TRUE)

  mle <- fit$par[1]
  val <- fit$value

  } else if (length(HR_index)==0 & length(LR_index)>0){ #We have LR bur not HR



    #Imputed variances for th LR studies
    s2_imp_LR <- 1/(k*n_LR)

    loglik <- function(param, p1, p2, y, s2, s2_imp_LR, N_rep, N_rep_s){

      #p1, p2 parameters which are assumed to be known (we can do a sensitivity analysis)
      #y the logRR of the reported studies
      #s2 the sigma squared of the reported studies
      #s2_imp_HR, s2_imp_LR the imputed variances of the HR and LR studies
      #n_rep number of reportes studies
      #n_rep_s number of reported and non significant studies


      #parameters to be estimated with maximum likelihood
      theta <- param[1]
      t2 <- param[2]
      alpha <- param[3]
      beta <- param[4]

      f <- (1/sqrt(2*pi*(s2 + t2)))*exp((-1/2)*(((y - theta)^2)/(s2 + t2))) #reported studies


      Q_LR <- pnorm((z_alpha*sqrt(s2_imp_LR) - theta)/(sqrt(s2_imp_LR + t2))) -
        pnorm((-z_alpha*sqrt(s2_imp_LR) -  theta)/(sqrt(s2_imp_LR + t2)))

      sum(log(f)) +
        N_rep*log(alpha) + #REPORTED
        N_rep_s*log(beta) +  #REPORTED AND NOT SIGNIFICANT
        sum(log((1-p1)*alpha*(Q_LR)*(1-beta) + p2*(1-alpha))) #LOW RISK

    }

    init_param <- init_param

    fit <- optim(init_param, loglik, p1=p1, p2=p2, y=y, s2=s2, s2_imp_LR=s2_imp_LR, N_rep=N_rep, N_rep_s=N_rep_s,
                 method =  "L-BFGS-B",
                 lower = c(-10, 0.01, 0.01, 0.01),
                 upper = c(10, 1,1,1),
                 control = list(fnscale = -1),
                 hessian=TRUE)

    mle <- fit$par[1]
    val <- fit$value


      } else { #No ORB

        loglik <- function(param, p1, p2, y, s2, N_rep, N_rep_s){

          #p1, p2 parameters which are assumed to be known (we can do a sensitivity analysis)
          #y the logRR of the reported studies
          #s2 the sigma squared of the reported studies

          #n_rep number of reportes studies
          #n_rep_s number of reported and non significant studies


          #parameters to be estimated with maximum likelihood
          theta <- param[1]
          t2 <- param[2]
          alpha <- param[3]
          beta <- param[4]

          f <- (1/sqrt(2*pi*(s2 + t2)))*exp((-1/2)*(((y - theta)^2)/(s2 + t2))) #reported studies

          sum(log(f)) +
            N_rep*log(alpha) + #REPORTED
            N_rep_s*log(beta)  #REPORTED AND NOT SIGNIFICANT


  }


  init_param <- init_param

  fit <- optim(init_param, loglik, p1=p1, p2=p2, y=y, s2=s2, N_rep=N_rep, N_rep_s=N_rep_s,
               #method = "BFGS",
               method =  "L-BFGS-B",
               lower = c(-5, 0.01, 0.01, 0.01),
               upper = c(10, 10,1,1),
               control = list(fnscale = -1),
               hessian=TRUE)


  mle <- fit$par[1]
  val <- fit$value

      }

  #WALD CONFIDENCE INTERVALS
  a <- sign_level #for harm Copas et al use 99% conf level
  #Unadjusted
  #fisher_info.unadjusted <- solve(-fit.unadjusted$hessian)
  #s.unadjusted <- sqrt(diag(fisher_info.unadjusted)[1])
  #ci.unadjusted <- fit.unadjusted$par[1] + qnorm(c(a/2, 1-a/2)) * s.unadjusted
  #width.unadjusted <- abs(ci.unadjusted[1]-ci.unadjusted[2])

  #Adjusted
  fisher_info <- solve(-fit$hessian)
  s <- sqrt(diag(fisher_info)[1])
  ci <- fit$par[1] + qnorm(c(a/2, 1-a/2)) * s
  width <- abs(ci[1]-ci[2])


  return(list(
              mle.unadjusted = mle.unadjusted,
              se.unadjusted = s.unadjusted,
              width.unadjusted = width.unadjusted,
              ci.unadjusted = ci.unadjusted,

              mle.adjusted = mle,
              se.adjusted = s,
              width.adjusted = width,
              ci.adjusted = ci))

  } else if (outcome == "harm"){




    if (length(HR_index)>0 & length(LR_index)>0){ #We have both HR and LR

      #Imputed variances for the HR studies
      s2_imp_HR <- 1/(k*n_HR)
      #Imputed variances for th LR studies
      s2_imp_LR <- 1/(k*n_LR)

      loglik <- function(param, eta1, eta2, y, s2, s2_imp_HR, s2_imp_LR, N_rep, N_rep_P){

        #p1, p2 parameters which are assumed to be known (we can do a sensitivity analysis)
        #y the logRR of the reported studies
        #s2 the sigma squared of the reported studies
        #s2_imp_HR, s2_imp_LR the imputed variances of the HR and LR studies
        #n_rep number of reportes studies
        #n_rep_s number of reported and non significant studies


        #parameters to be estimated with maximum likelihood
        theta <- param[1]
        t2 <- param[2]
        alpha <- param[3]
        gamm <- param[4]

        f <- (1/sqrt(2*pi*(s2 + t2)))*exp((-1/2)*(((y - theta)^2)/(s2 + t2))) #reported studies

        Q_HR <- pnorm(theta/(sqrt(s2_imp_HR + t2)))

        Q_LR <- pnorm(theta/(sqrt(s2_imp_LR + t2)))

        sum(log(f)) +
          N_rep*log(alpha) + #REPORTED
          N_rep_P*log(gamm) +  #REPORTED AND NOT SIGNIFICANT
          sum(log(eta1*alpha*(Q_HR)*(1-gamm) + (1-eta2)*(1-alpha))) + #HIGH RISK
          sum(log((1-eta1)*alpha*(Q_LR)*(1-gamm) + eta2*(1-alpha))) #LOW RISK

      }

      init_param <- init_param

      fit <- optim(init_param, loglik, eta1=eta1, eta2=eta2, y=y, s2=s2, s2_imp_HR=s2_imp_HR, s2_imp_LR=s2_imp_LR, N_rep=N_rep, N_rep_P=N_rep_P,
                   method =  "L-BFGS-B",
                  #method="Nelder-Mead",
                   lower = c(-10, 0.001, 0.001, 0.001),
                   upper = c(10, 20,1,1),
                   control = list(fnscale = -1),
                   hessian=TRUE)

      mle <- fit$par[1]
      val <- fit$value


    } else if (length(HR_index)>0 & length(LR_index)==0){ #We have HR but not LR



      #Imputed variances for the HR studies
      s2_imp_HR <- 1/(k*n_HR)

      loglik <- function(param, eta1, eta2, y, s2, s2_imp_HR, N_rep, N_rep_P){

        #p1, p2 parameters which are assumed to be known (we can do a sensitivity analysis)
        #y the logRR of the reported studies
        #s2 the sigma squared of the reported studies
        #s2_imp_HR, s2_imp_LR the imputed variances of the HR and LR studies
        #n_rep number of reportes studies
        #n_rep_s number of reported and non significant studies


        #parameters to be estimated with maximum likelihood
        theta <- param[1]
        t2 <- param[2]
        alpha <- param[3]
        gamm <- param[4]

        f <- (1/sqrt(2*pi*(s2 + t2)))*exp((-1/2)*(((y - theta)^2)/(s2 + t2))) #reported studies

        Q_HR <- pnorm(theta/(sqrt(s2_imp_HR + t2)))


        sum(log(f)) +
          N_rep*log(alpha) + #REPORTED
          N_rep_P*log(gamm) +  #REPORTED AND NOT SIGNIFICANT
          sum(log(eta1*alpha*(Q_HR)*(1-gamm) + (1-eta2)*(1-alpha)))  #HIGH RISK


      }

      init_param <- init_param

      fit <- optim(init_param, loglik, eta1=eta1, eta2=eta2, y=y, s2=s2, s2_imp_HR=s2_imp_HR, N_rep=N_rep, N_rep_P=N_rep_P,
                   method =  "L-BFGS-B",
                   control = list(fnscale = -1),
                   lower = c(-5, 0.05, 0.05, 0.05),
                   upper = c(10, 10, 1, 1),
                   hessian=TRUE)

      mle <- fit$par[1]
      val <- fit$value

    } else if (length(HR_index)==0 & length(LR_index)>0){ #We have LR bur not HR



      #Imputed variances for th LR studies
      s2_imp_LR <- 1/(k*n_LR)

      loglik <- function(param, eta1, eta2, y, s2, s2_imp_LR, N_rep, N_rep_P){

        #p1, p2 parameters which are assumed to be known (we can do a sensitivity analysis)
        #y the logRR of the reported studies
        #s2 the sigma squared of the reported studies
        #s2_imp_HR, s2_imp_LR the imputed variances of the HR and LR studies
        #n_rep number of reportes studies
        #n_rep_s number of reported and non significant studies


        #parameters to be estimated with maximum likelihood
        theta <- param[1]
        t2 <- param[2]
        alpha <- param[3]
        gamm <- param[4]

        f <- (1/sqrt(2*pi*(s2 + t2)))*exp((-1/2)*(((y - theta)^2)/(s2 + t2))) #reported studies


        Q_LR <- pnorm(theta/(sqrt(s2_imp_LR + t2)))

        sum(log(f)) +
          N_rep*log(alpha) + #REPORTED
          N_rep_P*log(gamm) +  #REPORTED AND NOT SIGNIFICANT
          sum(log((1-p1)*alpha*(Q_LR)*(1-gamm) + eta2*(1-alpha))) #LOW RISK

      }

      init_param <- init_param

      fit <- optim(init_param, loglik, eta1=eta1, eta2=eta2, y=y, s2=s2, s2_imp_LR=s2_imp_LR, N_rep=N_rep, N_rep_P=N_rep_P,
                   method =  "L-BFGS-B",
                   lower = c(-10, 0.01, 0.01, 0.01),
                   upper = c(10, 1,1,1),
                   control = list(fnscale = -1),
                   hessian=TRUE)

      mle <- fit$par[1]
      val <- fit$value


    } else { #No ORB

      loglik <- function(param, eta1, eta2, y, s2, N_rep, N_rep_P){

        #p1, p2 parameters which are assumed to be known (we can do a sensitivity analysis)
        #y the logRR of the reported studies
        #s2 the sigma squared of the reported studies

        #n_rep number of reportes studies
        #n_rep_s number of reported and non significant studies


        #parameters to be estimated with maximum likelihood
        theta <- param[1]
        t2 <- param[2]
        alpha <- param[3]
        gamm <- param[4]

        f <- (1/sqrt(2*pi*(s2 + t2)))*exp((-1/2)*(((y - theta)^2)/(s2 + t2))) #reported studies

        sum(log(f)) +
          N_rep*log(alpha) + #REPORTED
          N_rep_P*log(gamm)  #REPORTED AND NOT SIGNIFICANT


      }


      init_param <- init_param

      fit <- optim(init_param, loglik, eta1=eta1, eat2=eta2, y=y, s2=s2, N_rep=N_rep, N_rep_P=N_rep_P,
                   #method = "BFGS",
                   method =  "L-BFGS-B",
                   lower = c(-5, 0.01, 0.01, 0.01),
                   upper = c(10, 10,1,1),
                   control = list(fnscale = -1),
                   hessian=TRUE)


      mle <- fit$par[1]
      val <- fit$value

    }

    #WALD CONFIDENCE INTERVALS
    a <- sign_level #for harm Copas et al use 99% conf level
    #Unadjusted
    #fisher_info.unadjusted <- solve(-fit.unadjusted$hessian)
    #s.unadjusted <- sqrt(diag(fisher_info.unadjusted)[1])
    #ci.unadjusted <- fit.unadjusted$par[1] + qnorm(c(a/2, 1-a/2)) * s.unadjusted
    #width.unadjusted <- abs(ci.unadjusted[1]-ci.unadjusted[2])

    #Adjusted
    fisher_info <- solve(-fit$hessian)
    s <- sqrt(diag(fisher_info)[1])
    ci <- fit$par[1] + qnorm(c(a/2, 1-a/2)) * s
    width <- abs(ci[1]-ci[2])


    return(list(
      mle.unadjusted = mle.unadjusted,
      se.unadjusted = s.unadjusted,
      width.unadjusted = width.unadjusted,
      ci.unadjusted = ci.unadjusted,

      mle.adjusted = mle,
      se.adjusted = s,
      width.adjusted = width,
      ci.adjusted = ci))



  } else {

    return("invalid outcome input")
}

}



