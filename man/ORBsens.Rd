\name{ORBsens}
\alias{ORBsens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to adjust for outcome reporting bias (ORB) in a random-effect meta-analysis and conduct sensitivity analyses, according to the methodology developed by Copas et al. (2014) with extensions.
}
\description{
An implementation of the Copas et al. (2014) method to adjust for outcome reporting bias (ORB) in meta-anlyses of clinical trials for the random-effect model for beneficial outcomes, and extention to harmful outcomes (see Details). The \code{ORBsens} function can take in binary outcome data (i.e., the user should provide \code{a}, \code{c}, \code{n1}, \code{n2}), as well as means and standard errors (i.e., the user should provide \code{mu1}, \code{mu2}, \code{sd1}, \code{sd2}, \code{n1}, \code{n2}). The outcome variable is reported for some studies, and unreported for others. The unreported studies should have entries "high", i.e., high risk of bias, or "low", i.e., low risk of bias. The risk of bias classification should be done according to the ORBIT methodology. For the method to work, the sample sizes are required for all studies (reported and unreported). The function includes uncertainty in the risk of bias classification through the \code{p_vals} or \code{eta_vals} for beneficial and harmful outcomes respectively (see Details).
}
\usage{
ORBsens(
p_vals = NULL,
eta_vals = NULL,
a = NULL,
c = NULL,
mu1 = NULL,
mu2 = NULL,
sd1 = NULL,
sd2 = NULL,
n1,
n2,
ntot,
sign_level,
init_param,
outcome,
opt_method,
lower = c(-10, 0.05, 0.05, 0.05),
upper = c(10, 10, 1, 1)
)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p_vals}{A vector of size 2 specifiying the probabilities \eqn{\rho_1} and \eqn{\rho_2} (see Details).}

  \item{eta_vals}{A vector of size 2 specifiying the probabilities \eqn{\eta_1} and \eqn{\eta_2} (see Details).}

\item{a}{A vector of observed counts for the outcome of interest in the treatment group for each study in the meta-analysis. For the studies which do not report the outcome, the entries are either "high" or "low".}

\item{c}{A vector of observed counts for the outcome of interest in the control group for each study in the meta-analysis. For the studies which do not report the outcome, the entries are either "high" or "low".}

\item{mu1}{A vector of observed means for the outcome of interest in the treatment group for each study in the meta-analysis. For the studies which do not report the outcome, the entries are either "high" or "low".}

\item{mu2}{A vector of observed means for the outcome of interest in the control group for each study in the meta-analysis. For the studies which do not report the outcome, the entries are either "high" or "low".}

\item{sd1}{A vector of observed standard errors for the mean in the treatment group for each study in the meta-analysis. For the studies which do not report outcome, the entries are either "high" or "low".}

\item{sd2}{A vector of observed standard errors for the mean in the control group for each study in the meta-analysis. For the studies which do not report the outcome, the entries are either "high" or "low".}

\item{n1}{A vector of sample sizes in the treatment group for each study in the meta-analysis. The sample sizes need to be reported for all studies in the meta-analysis.}

\item{n2}{A vector of sample sizes in the control group for each study in the meta-analysis. The sample sizes need to be reported for all studies in the meta-analysis.}

  \item{ntot}{
A vector of total sample sizes for each study in the meta-analysis. The sample sizes need to be reported for all studies in the meta-analysis.
}
  \item{sign_level}{
The value for the desired confidence level if the outcome of interest is beneficial, e.g., \eqn{\alpha=0.05} for beneficial outcomes, and \eqn{\alpha=0.01} for harmful outcomes according to Copas et al. (2019).
}
  \item{init_param}{
A vector of length 4, containing the values to be used for initialization of the joint optimization. For beneficial outcomes these are:  \code{init_param = c ( \eqn{\mu}, \eqn{\tau^2}, \eqn{\alpha}, \eqn{\beta} )}.For harmful outcomes these are: \code{init_param = c( \eqn{\mu}, \eqn{\tau^2}, \eqn{\alpha}, \eqn{\gamma} )} See Details for definitions of \eqn{\alpha}, \eqn{\beta} and \eqn{\gamma}.
}
  \item{outcome}{
A character indicating whether the outcome of interest is beneficial ("benefit") or harmful ("harm").
}
  \item{opt_method}{
A character specifying the method used in the \eqn{\texttt{optim()}} function for the joint optimization of the above parameters. Reccomended is "L-BFGS-B".
}
  \item{lower}{
A vector of length 4 containing the lower bounds for the optimization problem for \eqn{\mu, \tau^2, \alpha, \beta} or \eqn{\mu, \tau^2, \alpha, \gamma}. Default is \code{c(-10, 0.05, 0.05, 0.05)}. The lower bounds are specified only if the "L-BFGS-B" optimization method is used.
}
  \item{upper}{
A vector of length 4 containing the upper bounds for the optimization problem for \eqn{\mu, \tau^2, \alpha, \beta} or \eqn{\mu, \tau^2, \alpha, \gamma}. Default is \code{c(10, 10, 1, 1)}. The upper bounds are specified only if the "L-BFGS-B" optimization method is used.
}

}
\details{

For a beneficial outcome, the log-likelihood to be maximized to obtain an ORB-adjusted estimate of the treatment effect \eqn{\mu} is

\deqn{

l_{\rho_{1}, \rho_{2}} (\mu, \tau^2, \alpha, \beta) = \\
\sum_{ i \in \text{Rep}} \log  [f_{i}(y_{i})] +
N_{R} \log (\alpha) +
N_{R \cap s} \log (\beta) + \\
\sum_{i \in \text{HR}} \log [\rho_{1} \alpha Q_{i} (1 - \beta) + (1 - \rho_{2})(1 - \alpha) ] + \\
\sum_{i \in \text{LR}} \log [(1-\rho_{1})\alpha Q_{i}(1-\beta) + \rho_{2}(1-\alpha) ]


}

where \eqn{\mu} the global treatment effect, \eqn{\tau^2} is the the between-study heterogeneity, \eqn{\alpha} is the probability of an outcome being measured, and \eqn{\beta} is the probability an outcome being reported given that it was measured and non-significant. The \eqn{\rho_{1}, \rho_{2}} reliability parameters are treated as fixed. \eqn{\rho_{1}} is defined as the probability that an outcome is classified as high risk of bias given that it was measured, non-significant, and unreported. \eqn{\rho_{2}} is the probability that an outcome is classified as low risk of bias given that it was not measured and was not reported. \eqn{Q_i} is the probability of the observed treatment effect \eqn{y_i} being significant, while \eqn{f_i(y_i)} is the pdf of a normal distribution. See Copas et al. (2014) for more details and for the derivation of the log-likelihood function.

For a harmful outcome, the log-likelihood to be maximized is

\deqn{

L_{\eta_{1}, \eta_{2}} (\mu, \tau^2, \alpha, \gamma)  = \\
\sum_{i \in \text{Rep}} \log [f_{i}(y_{i}) ] + N_{\text{R}} \log (\alpha) + N_{\text{R} \cap P } \log (\gamma) +\\
\sum_{i \in \text{LR}} \log [\alpha Q_{h_i} (1-\gamma) (1-\eta_1) + (1-\alpha)\eta_2 ] + \\
\sum_{i \in \text{HR}} \log [\alpha Q_{h_i} (1-\gamma) \eta_1 + (1-\alpha)(1-\eta_2) ]

}

where \eqn{\mu} the global treatment effect,\eqn{\tau^2} is the the between-study heterogeneity, \eqn{\alpha} is the probability of an outcome being measured, and \eqn{\gamma} is the probability an outcome being reported given that it was measured and positive. The \eqn{\eta_{1}, \eta_{2}} reliability parameters are treated as fixed. \eqn{\eta_{1}} is defined as the probability that an outcome is classified as high risk of bias given that it was measured, positive, and unreported. \eqn{\eta_{2}} is the probability that an outcome is classified as low risk of bias given that it was not measured and was not reported.


}
\value{

Returns a list with the following elements

\item{mle.unadjusted}{Unadjusted MLE of the global treatment effect \eqn{\mu}}

\item{se.unadjusted}{Unadjusted standard error of the MLE for the global treatment effect \eqn{\mu}}

\item{width.unadjusted}{Unadjusted width of the Wald CI for the global treatment effect \eqn{\mu}}

\item{ci.unadjusted}{Unadjusted Wald CI for the global treatment effect \eqn{\mu}}

\item{mle.adjusted}{ORB-adjusted MLE of the global treatment effect \eqn{\mu}}

\item{se.adjusted}{ORB-adjusted standard error of the MLE for the global treatment effect \eqn{\mu}}

\item{width.adjusted}{ORB-adjusted width of the Wald CI for the global treatment effect \eqn{\mu}}

\item{ci.unadjusted}{ORB-adjusted Wald CI for the global treatment effect \eqn{\mu}}


}
\references{

Copas, J., Dwan, K., Kirkham, J., and Williamson, P. (2014). A model-based correction for outcome reporting bias in meta-analysis. Biostatistics, 15, 370–383.

}
\author{
Alessandra Gaia Saracini under the supervision of Prof. Dr. Leonhard Held.
}
\note{
R package created as part of Master's Thesis in Statistics at ETH Zurich.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

#Sensitivity analysis contour plot for the ORB-adjusted estiamte of the treatment effect of a beneficial outcome from PPHData.

p1.grid <- seq(0.7, 1, by=0.1)
p2.grid <- seq(0.7, 1, by=0.1)
grid <- expand.grid(p1 = p1.grid, p2 = p2.grid)
values.ex1 <- matrix(data = apply(grid,
                                  1,
                                  function(p) ORBmeta:::ORBsens(p,
                                                      a=PPHData$`PPH (C)`,
                                                      c=PPHData$`PPH (T)`,
                                                      n1=PPHData$`Sample (C)`,
                                                      n2=PPHData$`Sample (T)`,
                                                      ntot=PPHData$`Sample Total`,
                                                      sign_level=0.05,
                                                      init_param=c(0.7, 0.6, 0.9, 0.9),
                                  										outcome="benefit",
                                  										opt_method="L-BFGS-B")$mle.adjusted),
                                  nrow = length(p1.grid))




filled.contour(x=p1.grid,
               y=p2.grid,
               z=values.ex1,
               levels = seq(-0.7, -0.6, by=0.005),
               xlab = expression(paste(rho[1], sep = "")),
               ylab = expression(paste(rho[2], sep = "")),
               plot.axes = {
                 axis(1)
                 axis(2)
                 contour(x=p1.grid,
                         y=p2.grid,
                         z=values.ex1,
                         levels = seq(-0.7, -0.6, by=0.005),
                         add = TRUE, lwd = 1)
               },
							 main="ORB-adjusted Treatment Effect")


#Sensitivity analysis contour plot for the ORB-adjusted estiamte of the treatment effect of a harmful outcome from EpilepsyData.

eta1.grid <- seq(0.7, 1, by=0.1)
eta2.grid <- seq(0.7, 1, by=0.1)
grid <- expand.grid(eta1 = eta1.grid, eta2 = eta2.grid)
values <- matrix(data = apply(grid,
                              1,
                              function(eta) ORBmeta:::ORBsens(p_vals=NULL,
                                                    eta,
                                                    a=EpilepsyData$`Somnolence (T)`,
                                                    c=EpilepsyData$`Somnolence (C)`,
                                                    n1=EpilepsyData$`Sample (T)`,
                                                    n2=EpilepsyData$`Sample (C)`,
                                                    ntot=as.numeric(EpilepsyData$`Sample (T)`) + as.numeric(EpilepsyData$`Sample (C)`),
                              											outcome="harm",
                                                    sign_level=0.05,
                                                    init_param=c(0.7, 0.6, 0.9, 0.9),
                              											opt_method="L-BFGS-B",
                              											lower = c(-10, 0.001, 0.01, 0.01),
                                                    upper = c(10, 10,0.99,0.99))$mle.adjusted),
                               nrow = length(eta1.grid))


filled.contour(x=eta1.grid,
               y=eta2.grid,
               z=values,
               levels = seq(0.79, 0.802, by=0.001),
               xlab = expression(paste(eta[1], sep = "")),
               ylab = expression(paste(eta[2], sep = "")),
               plot.axes = {
                 axis(1)
                 axis(2)
                 contour(x=eta1.grid,
                         y=eta2.grid,
                         z=values,
                         levels = seq(0.68, 0.81, by=0.001),
                         add = TRUE, lwd = 1)
               },
							 main="ORB-adjusted Treatment Effect")





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
