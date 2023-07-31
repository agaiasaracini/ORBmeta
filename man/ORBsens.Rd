\name{ORBsens}
\alias{ORBsens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ORB sens
}
\description{
Sensitivity Analaysis of Copas et al. 2014 methodlogy plus extension
}
\usage{
ORBsens(p_vals = NULL, eta_vals = NULL, a = NULL, c = NULL, mu1 = NULL, mu2 = NULL, sd1 = NULL, sd2 = NULL, n1, n2, ntot, sign_level, init_param, outcome, opt_method, lower = c(-10, 0.05, 0.05, 0.05), upper = c(10, 10, 1, 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p_vals}{
%%     ~~Describe \code{p_vals} here~~
}
  \item{eta_vals}{
%%     ~~Describe \code{eta_vals} here~~
}
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{c}{
%%     ~~Describe \code{c} here~~
}
  \item{mu1}{
%%     ~~Describe \code{mu1} here~~
}
  \item{mu2}{
%%     ~~Describe \code{mu2} here~~
}
  \item{sd1}{
%%     ~~Describe \code{sd1} here~~
}
  \item{sd2}{
%%     ~~Describe \code{sd2} here~~
}
  \item{n1}{
%%     ~~Describe \code{n1} here~~
}
  \item{n2}{
%%     ~~Describe \code{n2} here~~
}
  \item{ntot}{
%%     ~~Describe \code{ntot} here~~
}
  \item{sign_level}{
%%     ~~Describe \code{sign_level} here~~
}
  \item{init_param}{
%%     ~~Describe \code{init_param} here~~
}
  \item{outcome}{
%%     ~~Describe \code{outcome} here~~
}
  \item{opt_method}{
%%     ~~Describe \code{opt_method} here~~
}
  \item{lower}{
%%     ~~Describe \code{lower} here~~
}
  \item{upper}{
%%     ~~Describe \code{upper} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (p_vals = NULL, eta_vals = NULL, a = NULL, c = NULL,
    mu1 = NULL, mu2 = NULL, sd1 = NULL, sd2 = NULL, n1, n2, ntot,
    sign_level, init_param, outcome, opt_method, lower = c(-10,
        0.05, 0.05, 0.05), upper = c(10, 10, 1, 1))
{
    p1 <- p_vals[1]
    p2 <- p_vals[2]
    eta1 <- eta_vals[1]
    eta2 <- eta_vals[2]
    method <- opt_method
    lower <- lower
    upper <- upper
    z_alpha <- qnorm(1 - sign_level/2)
    if (!is.null(a) & !is.null(c)) {
        Rep_index <- which(!is.na(as.numeric(a)))
        HR_index <- which(a == "high")
        LR_index <- which(a == "low")
        a_rep <- as.numeric(a[Rep_index])
        c_rep <- as.numeric(c[Rep_index])
        n1_rep <- as.numeric(n1[Rep_index])
        n2_rep <- as.numeric(n2[Rep_index])
        if (length(a_rep == 0) | length(c_rep == 0) > 0) {
            a_0_index <- c(which(a_rep == 0), which(c_rep ==
                0))
            a_rep[a_0_index] <- a_rep[a_0_index] + 0.5
            c_rep[a_0_index] <- c_rep[a_0_index] + 0.5
            n1_rep[a_0_index] <- n1_rep[a_0_index] + 0.5
            n2_rep[a_0_index] <- n2_rep[a_0_index] + 0.5
            a_0_both <- which(a_rep == 0.5 & c_rep == 0.5)
            if (length(a_0_both) > 0) {
                a_rep <- a_rep[-a_0_both]
                c_rep <- c_rep[-a_0_both]
                n1_rep <- n1_rep[-a_0_both]
                n2_rep <- n2_rep[-a_0_both]
            }
            else {
                a_rep <- a_rep
                c_rep <- c_rep
                n1_rep <- n1_rep
                n2_rep <- n2_rep
            }
        }
        else {
            a_rep <- a_rep
            c_rep <- c_rep
            n1_rep <- n1_rep
            n2_rep <- n2_rep
        }
        N_rep <- length(Rep_index)
        n_HR <- as.numeric(ntot[HR_index])
        n_LR <- as.numeric(ntot[LR_index])
        y <- log((a_rep * n2_rep)/(c_rep * n1_rep))
        s <- sqrt(((n1_rep - a_rep)/(n1_rep * a_rep)) + ((n2_rep -
            c_rep)/(n2_rep * c_rep)))
        s2 <- s^2
        y_S <- y[which(abs(y) > z_alpha * s)]
        y_s <- y[which(abs(y) <= z_alpha * s)]
        y_P <- y[which(y > 0)]
        y_p <- y[which(y <= 0)]
        N_rep_P <- length(y_P)
        N_rep_p <- length(y_p)
        N_rep_s <- length(y_s)
        k <- sum(1/(s2))/sum((n1_rep + n2_rep))
    }
    else if (!is.null(mu1) & !is.null(mu2) & !is.null(sd1) &
        !is.null(sd2)) {
        Rep_index <- which(!is.na(as.numeric(mu1)))
        HR_index <- which(mu1 == "high")
        LR_index <- which(mu1 == "low")
        n_HR <- as.numeric(ntot[HR_index])
        n_LR <- as.numeric(ntot[LR_index])
        mu1_rep <- as.numeric(mu1[Rep_index])
        mu2_rep <- as.numeric(mu2[Rep_index])
        sd1_rep <- as.numeric(sd1[Rep_index])
        sd2_rep <- as.numeric(sd2[Rep_index])
        n1_rep <- as.numeric(n1[Rep_index])
        n2_rep <- as.numeric(n2[Rep_index])
        N_rep <- length(Rep_index)
        y <- mu1_rep - mu2_rep
        s <- sqrt((as.numeric(sd1_rep)^2)/(as.numeric(n1_rep)) +
            (as.numeric(sd2_rep)^2)/(as.numeric(n2_rep)))
        s2 <- s^2
        y_S <- y[which(abs(y) > z_alpha * s)]
        y_s <- y[which(abs(y) <= z_alpha * s)]
        N_rep_s <- length(y_s)
        y_P <- y[which(y > 0)]
        y_p <- y[which(y <= 0)]
        N_rep_P <- length(y_P)
        N_rep_p <- length(y_p)
        k <- sum(1/(s2))/sum((n1_rep + n2_rep))
    }
    else {
        return("Error: invalid inputs. Input either a and c values or y1,sd1 and y2,sd2 values.")
    }
    loglik.unadjusted <- function(param, y, s2, N_rep, N_rep_s) {
        theta <- param[1]
        t2 <- param[2]
        alpha <- param[3]
        beta <- param[4]
        f <- (1/sqrt(2 * pi * (s2 + t2))) * exp((-1/2) * (((y -
            theta)^2)/(s2 + t2)))
        sum(log(f)) + N_rep * log(alpha) + N_rep_s * log(beta)
    }
    init_param <- init_param
    if (method == "L-BFGS-B") {
        fit.unadjusted <- optim(init_param, loglik.unadjusted,
            y = y, s2 = s2, N_rep = N_rep, N_rep_s = N_rep_s,
            method = method, lower = lower, upper = upper, control = list(fnscale = -1),
            hessian = TRUE)
    }
    else {
        fit.unadjusted <- optim(init_param, loglik.unadjusted,
            y = y, s2 = s2, N_rep = N_rep, N_rep_s = N_rep_s,
            method = method, control = list(fnscale = -1), hessian = TRUE)
    }
    mle.unadjusted <- fit.unadjusted$par[1]
    a <- sign_level
    fisher_info.unadjusted <- solve(-fit.unadjusted$hessian)
    s.unadjusted <- sqrt(diag(fisher_info.unadjusted)[1])
    ci.unadjusted <- fit.unadjusted$par[1] + qnorm(c(a/2, 1 -
        a/2)) * s.unadjusted
    width.unadjusted <- abs(ci.unadjusted[1] - ci.unadjusted[2])
    if (outcome == "benefit") {
        if (length(HR_index) > 0 & length(LR_index) > 0) {
            s2_imp_HR <- 1/(k * n_HR)
            s2_imp_LR <- 1/(k * n_LR)
            loglik <- function(param, p1, p2, y, s2, s2_imp_HR,
                s2_imp_LR, N_rep, N_rep_s) {
                theta <- param[1]
                t2 <- param[2]
                alpha <- param[3]
                beta <- param[4]
                f <- (1/sqrt(2 * pi * (s2 + t2))) * exp((-1/2) *
                  (((y - theta)^2)/(s2 + t2)))
                Q_HR <- pnorm((z_alpha * sqrt(s2_imp_HR) - theta)/(sqrt(s2_imp_HR +
                  t2))) - pnorm((-z_alpha * sqrt(s2_imp_HR) -
                  theta)/(sqrt(s2_imp_HR + t2)))
                Q_LR <- pnorm((z_alpha * sqrt(s2_imp_LR) - theta)/(sqrt(s2_imp_LR +
                  t2))) - pnorm((-z_alpha * sqrt(s2_imp_LR) -
                  theta)/(sqrt(s2_imp_LR + t2)))
                sum(log(f)) + N_rep * log(alpha) + N_rep_s *
                  log(beta) + sum(log(p1 * alpha * (Q_HR) * (1 -
                  beta) + (1 - p2) * (1 - alpha))) + sum(log((1 -
                  p1) * alpha * (Q_LR) * (1 - beta) + p2 * (1 -
                  alpha)))
            }
            init_param <- init_param
            if (method == "L-BFGS-B") {
                fit <- optim(init_param, loglik, p1 = p1, p2 = p2,
                  y = y, s2 = s2, s2_imp_HR = s2_imp_HR, s2_imp_LR = s2_imp_LR,
                  N_rep = N_rep, N_rep_s = N_rep_s, method = method,
                  lower = lower, upper = upper, control = list(fnscale = -1),
                  hessian = TRUE)
            }
            else {
                fit <- optim(init_param, loglik, p1 = p1, p2 = p2,
                  y = y, s2 = s2, s2_imp_HR = s2_imp_HR, s2_imp_LR = s2_imp_LR,
                  N_rep = N_rep, N_rep_s = N_rep_s, method = method,
                  control = list(fnscale = -1), hessian = TRUE)
            }
            mle <- fit$par[1]
            val <- fit$value
        }
        else if (length(HR_index) > 0 & length(LR_index) == 0) {
            s2_imp_HR <- 1/(k * n_HR)
            loglik <- function(param, p1, p2, y, s2, s2_imp_HR,
                N_rep, N_rep_s) {
                theta <- param[1]
                t2 <- param[2]
                alpha <- param[3]
                beta <- param[4]
                f <- (1/sqrt(2 * pi * (s2 + t2))) * exp((-1/2) *
                  (((y - theta)^2)/(s2 + t2)))
                Q_HR <- pnorm((z_alpha * sqrt(s2_imp_HR) - theta)/(sqrt(s2_imp_HR +
                  t2)))
                -pnorm((-z_alpha * sqrt(s2_imp_HR) - theta)/(sqrt(s2_imp_HR +
                  t2)))
                sum(log(f)) + N_rep * log(alpha) + N_rep_s *
                  log(beta) + sum(log(p1 * alpha * (Q_HR) * (1 -
                  beta) + (1 - p2) * (1 - alpha)))
            }
            init_param <- init_param
            if (method == "L-BFGS-B") {
                fit <- optim(init_param, loglik, p1 = p1, p2 = p2,
                  y = y, s2 = s2, s2_imp_HR = s2_imp_HR, N_rep = N_rep,
                  N_rep_s = N_rep_s, method = method, control = list(fnscale = -1),
                  lower = lower, upper = upper, hessian = TRUE)
            }
            else {
                fit <- optim(init_param, loglik, p1 = p1, p2 = p2,
                  y = y, s2 = s2, s2_imp_HR = s2_imp_HR, N_rep = N_rep,
                  N_rep_s = N_rep_s, method = method, control = list(fnscale = -1),
                  hessian = TRUE)
            }
            mle <- fit$par[1]
            val <- fit$value
        }
        else if (length(HR_index) == 0 & length(LR_index) > 0) {
            s2_imp_LR <- 1/(k * n_LR)
            loglik <- function(param, p1, p2, y, s2, s2_imp_LR,
                N_rep, N_rep_s) {
                theta <- param[1]
                t2 <- param[2]
                alpha <- param[3]
                beta <- param[4]
                f <- (1/sqrt(2 * pi * (s2 + t2))) * exp((-1/2) *
                  (((y - theta)^2)/(s2 + t2)))
                Q_LR <- pnorm((z_alpha * sqrt(s2_imp_LR) - theta)/(sqrt(s2_imp_LR +
                  t2))) - pnorm((-z_alpha * sqrt(s2_imp_LR) -
                  theta)/(sqrt(s2_imp_LR + t2)))
                sum(log(f)) + N_rep * log(alpha) + N_rep_s *
                  log(beta) + sum(log((1 - p1) * alpha * (Q_LR) *
                  (1 - beta) + p2 * (1 - alpha)))
            }
            init_param <- init_param
            if (method == "L-BFGS-B") {
                fit <- optim(init_param, loglik, p1 = p1, p2 = p2,
                  y = y, s2 = s2, s2_imp_LR = s2_imp_LR, N_rep = N_rep,
                  N_rep_s = N_rep_s, method = method, lower = lower,
                  upper = upper, control = list(fnscale = -1),
                  hessian = TRUE)
            }
            else {
                fit <- optim(init_param, loglik, p1 = p1, p2 = p2,
                  y = y, s2 = s2, s2_imp_LR = s2_imp_LR, N_rep = N_rep,
                  N_rep_s = N_rep_s, method = method, control = list(fnscale = -1),
                  hessian = TRUE)
            }
            mle <- fit$par[1]
            val <- fit$value
        }
        else {
            loglik <- function(param, p1, p2, y, s2, N_rep, N_rep_s) {
                theta <- param[1]
                t2 <- param[2]
                alpha <- param[3]
                beta <- param[4]
                f <- (1/sqrt(2 * pi * (s2 + t2))) * exp((-1/2) *
                  (((y - theta)^2)/(s2 + t2)))
                sum(log(f)) + N_rep * log(alpha) + N_rep_s *
                  log(beta)
            }
            init_param <- init_param
            if (method == "L-BFGS-B") {
                fit <- optim(init_param, loglik, p1 = p1, p2 = p2,
                  y = y, s2 = s2, N_rep = N_rep, N_rep_s = N_rep_s,
                  method = method, lower = lower, upper = upper,
                  control = list(fnscale = -1), hessian = TRUE)
            }
            else {
                fit <- optim(init_param, loglik, p1 = p1, p2 = p2,
                  y = y, s2 = s2, N_rep = N_rep, N_rep_s = N_rep_s,
                  method = method, control = list(fnscale = -1),
                  hessian = TRUE)
            }
            mle <- fit$par[1]
            val <- fit$value
        }
        a <- sign_level
        fisher_info <- solve(-fit$hessian)
        s <- sqrt(diag(fisher_info)[1])
        ci <- fit$par[1] + qnorm(c(a/2, 1 - a/2)) * s
        width <- abs(ci[1] - ci[2])
        return(list(mle.unadjusted = mle.unadjusted, se.unadjusted = s.unadjusted,
            width.unadjusted = width.unadjusted, ci.unadjusted = ci.unadjusted,
            mle.adjusted = mle, se.adjusted = s, width.adjusted = width,
            ci.adjusted = ci))
    }
    else if (outcome == "harm") {
        if (length(HR_index) > 0 & length(LR_index) > 0) {
            s2_imp_HR <- 1/(k * n_HR)
            s2_imp_LR <- 1/(k * n_LR)
            loglik <- function(param, eta1, eta2, y, s2, s2_imp_HR,
                s2_imp_LR, N_rep, N_rep_P) {
                theta <- param[1]
                t2 <- param[2]
                alpha <- param[3]
                gamm <- param[4]
                f <- (1/sqrt(2 * pi * (s2 + t2))) * exp((-1/2) *
                  (((y - theta)^2)/(s2 + t2)))
                Q_HR <- pnorm(theta/(sqrt(s2_imp_HR + t2)))
                Q_LR <- pnorm(theta/(sqrt(s2_imp_LR + t2)))
                sum(log(f)) + N_rep * log(alpha) + N_rep_P *
                  log(gamm) + sum(log(eta1 * alpha * (Q_HR) *
                  (1 - gamm) + (1 - eta2) * (1 - alpha))) + sum(log((1 -
                  eta1) * alpha * (Q_LR) * (1 - gamm) + eta2 *
                  (1 - alpha)))
            }
            init_param <- init_param
            if (method == "L-BFGS-B") {
                fit <- optim(init_param, loglik, eta1 = eta1,
                  eta2 = eta2, y = y, s2 = s2, s2_imp_HR = s2_imp_HR,
                  s2_imp_LR = s2_imp_LR, N_rep = N_rep, N_rep_P = N_rep_P,
                  method = method, lower = lower, upper = upper,
                  control = list(fnscale = -1), hessian = TRUE)
            }
            else {
                fit <- optim(init_param, loglik, eta1 = eta1,
                  eta2 = eta2, y = y, s2 = s2, s2_imp_HR = s2_imp_HR,
                  s2_imp_LR = s2_imp_LR, N_rep = N_rep, N_rep_P = N_rep_P,
                  method = method, control = list(fnscale = -1),
                  hessian = TRUE)
            }
            mle <- fit$par[1]
            val <- fit$value
        }
        else if (length(HR_index) > 0 & length(LR_index) == 0) {
            s2_imp_HR <- 1/(k * n_HR)
            loglik <- function(param, eta1, eta2, y, s2, s2_imp_HR,
                N_rep, N_rep_P) {
                theta <- param[1]
                t2 <- param[2]
                alpha <- param[3]
                gamm <- param[4]
                f <- (1/sqrt(2 * pi * (s2 + t2))) * exp((-1/2) *
                  (((y - theta)^2)/(s2 + t2)))
                Q_HR <- pnorm(theta/(sqrt(s2_imp_HR + t2)))
                sum(log(f)) + N_rep * log(alpha) + N_rep_P *
                  log(gamm) + sum(log(eta1 * alpha * (Q_HR) *
                  (1 - gamm) + (1 - eta2) * (1 - alpha)))
            }
            init_param <- init_param
            if (method == "L-BFGS-B") {
                fit <- optim(init_param, loglik, eta1 = eta1,
                  eta2 = eta2, y = y, s2 = s2, s2_imp_HR = s2_imp_HR,
                  N_rep = N_rep, N_rep_P = N_rep_P, method = method,
                  control = list(fnscale = -1), lower = lower,
                  upper = upper, hessian = TRUE)
            }
            else {
                fit <- optim(init_param, loglik, eta1 = eta1,
                  eta2 = eta2, y = y, s2 = s2, s2_imp_HR = s2_imp_HR,
                  N_rep = N_rep, N_rep_P = N_rep_P, method = method,
                  control = list(fnscale = -1), hessian = TRUE)
            }
            mle <- fit$par[1]
            val <- fit$value
        }
        else if (length(HR_index) == 0 & length(LR_index) > 0) {
            s2_imp_LR <- 1/(k * n_LR)
            loglik <- function(param, eta1, eta2, y, s2, s2_imp_LR,
                N_rep, N_rep_P) {
                theta <- param[1]
                t2 <- param[2]
                alpha <- param[3]
                gamm <- param[4]
                f <- (1/sqrt(2 * pi * (s2 + t2))) * exp((-1/2) *
                  (((y - theta)^2)/(s2 + t2)))
                Q_LR <- pnorm(theta/(sqrt(s2_imp_LR + t2)))
                sum(log(f)) + N_rep * log(alpha) + N_rep_P *
                  log(gamm) + sum(log((1 - p1) * alpha * (Q_LR) *
                  (1 - gamm) + eta2 * (1 - alpha)))
            }
            init_param <- init_param
            if (method == "L-BFGS-B") {
                fit <- optim(init_param, loglik, eta1 = eta1,
                  eta2 = eta2, y = y, s2 = s2, s2_imp_LR = s2_imp_LR,
                  N_rep = N_rep, N_rep_P = N_rep_P, method = method,
                  lower = lower, upper = upper, control = list(fnscale = -1),
                  hessian = TRUE)
            }
            else {
                fit <- optim(init_param, loglik, eta1 = eta1,
                  eta2 = eta2, y = y, s2 = s2, s2_imp_LR = s2_imp_LR,
                  N_rep = N_rep, N_rep_P = N_rep_P, method = method,
                  control = list(fnscale = -1), hessian = TRUE)
            }
            mle <- fit$par[1]
            val <- fit$value
        }
        else {
            loglik <- function(param, eta1, eta2, y, s2, N_rep,
                N_rep_P) {
                theta <- param[1]
                t2 <- param[2]
                alpha <- param[3]
                gamm <- param[4]
                f <- (1/sqrt(2 * pi * (s2 + t2))) * exp((-1/2) *
                  (((y - theta)^2)/(s2 + t2)))
                sum(log(f)) + N_rep * log(alpha) + N_rep_P *
                  log(gamm)
            }
            init_param <- init_param
            if (method == "L-BFGS-B") {
                fit <- optim(init_param, loglik, eta1 = eta1,
                  eat2 = eta2, y = y, s2 = s2, N_rep = N_rep,
                  N_rep_P = N_rep_P, method = method, lower = lower,
                  upper = upper, control = list(fnscale = -1),
                  hessian = TRUE)
            }
            else {
                fit <- optim(init_param, loglik, eta1 = eta1,
                  eat2 = eta2, y = y, s2 = s2, N_rep = N_rep,
                  N_rep_P = N_rep_P, method = method, control = list(fnscale = -1),
                  hessian = TRUE)
            }
            mle <- fit$par[1]
            val <- fit$value
        }
        a <- sign_level
        fisher_info <- solve(-fit$hessian)
        s <- sqrt(diag(fisher_info)[1])
        ci <- fit$par[1] + qnorm(c(a/2, 1 - a/2)) * s
        width <- abs(ci[1] - ci[2])
        return(list(mle.unadjusted = mle.unadjusted, se.unadjusted = s.unadjusted,
            width.unadjusted = width.unadjusted, ci.unadjusted = ci.unadjusted,
            mle.adjusted = mle, se.adjusted = s, width.adjusted = width,
            ci.adjusted = ci))
    }
    else {
        return("invalid outcome input")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
