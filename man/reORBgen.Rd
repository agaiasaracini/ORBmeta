\name{reROBgen}
\alias{reROBgen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to adjust for outcome reporting bias (ORB) in an outcome of a meta-analysis of clinical trials
}
\description{
An implementation of the Copas et al. (2019) method to adjust for outcome reporting bias (ORB) in meta-anlyses of clinical trials with extention to different selection functions to model the ORB mechanism (see details). The reORBgen function can take in binary outcome data (i.e., the user shoul provide a, c, n1, n2), means and standard errors (i.e., the user should provide mu1, mu2, sd1, sd2, n1, n2), or a treatment effect assumed to be normally distributed (i.e., the user should provide y, s, n1, n2). The outcome variable is reported for some studies, and unreported for others. The studies for which the observed outcome and standard error are missing should have entries "high", i.e., high risk of bias, or "low", i.e., low risk of bias. For the method to work, the sample sizes are required for the studies which do not report the observed treatment effect. The risk of bias classification should be done according to the ORBIT methodology (cite website). The function can adjust both beneficial and harmful outcomes.
}
\usage{
(a=NULL,
c=NULL,
mu1=NULL,
mu2=NULL,
sd1=NULL,
sd2=NULL,
y=NULL,
s=NULL,
n1,
n2,
outcome,
init_param,
alpha_ben=NULL,
alpha_ben_one.sided = TRUE,
alpha_harm=NULL,
true.SE=NULL,
LR.CI = FALSE,
Wald.CI = FALSE,
selection.benefit = "Copas.oneside",
weight_param=15,
opt_method="L-BFGS-B",
lower = c(-5, 0.00001),
upper = c(5,5),
low.risk = FALSE
}
%- maybe also 'usage' for other objects documented here.
\arguments{

\item{a}{A vector of observed counts for the outcome of interest in the treatment group for each study in the meta-analysis. For the studies which do not report the counts in the outcome of interest, the entries are either "high" or "low".}

\item{c}{A vector of observed counts for the outcome of interest in the control group for each study in the meta-analysis. For the studies which do not report the counts in the outcome of interest, the entries are either "high" or "low".}

\item{mu1}{A vector of observed means for the outcome of interest in the treatment group for each study in the meta-analysis. For the studies which do not report the means in the outcome of interest, the entries are either "high" or "low".}

\item{mu2}{A vector of observed means for the outcome of interest in the control group for each study in the meta-analysis. For the studies which do not report the means in the outcome of interest, the entries are either "high" or "low".}

\item{sd1}{A vector of observed standard errors for the mean in the treatment group for each study in the meta-analysis. For the studies which do not report the mean/standard error in the outcome of interest, the entries are either "high" or "low".}

\item{sd2}{A vector of observed standard errors for the mean in the control group for each study in the meta-analysis. For the studies which do not report the mean/standard error in the outcome of interest, the entries are either "high" or "low".}

\item{y}{A vector of observed treatment effects for the outcome of interest. For the studies which do not report the observed treatment effect in the outcome of interest, the entries are either "high" or "low".}

\item{s}{A vector of observed standard errors for treatment effect. For  the studies which do not report the treatment effect/standard error, the entries are either "high" or "low". }

\item{n1}{A vector of sample sizes in the treatment group for each study in the meta-analysis. The sample sizes need to be reported for all studies in the meta-analysis.}

\item{n2}{A vector of sample sizes in the control group for each study in the meta-analysis. The sample sizes need to be reported for all studies in the meta-analysis.}

\item{outcome}{A character indicating whether the outcome of interest is beneficial ("benefit") or harmful ("harm").}

\item{init_param}{A vector of length two, containing the values to be used for initialization of the joint optimization of the global treatment effect $mu$ and the between-study heterogeneity $tau^2$ of the random effect model, e.g., initi_param = c(0.1, 0.05).}

\item{alpha_ben}{The value for the desired confidence level if the outcome of interest is beneficial, e.g., $alpha=0.05$ in Copas et al. (2019).}

\item{alpha_ben_one.sided}{If TRUE, the one-sided significance threshold qnorm(1-alpha) is used for the selection functions (see details); if FALSE, the two-sided one is used qnorm(1-alpha/2).}

\item{alpha_harm}{The value for the desired confidence level if the outcome of interest is beneficial, e.g., $alpha=0.01$ in Copas et al. (2019).}

\item{true.SE}{A vector of the observed standard errors of the treatment effects for all studies in the meta-analysis, including the ones for which we have "high" or "low", if known. Of note, this works if the outcome data is provided as an observed treatment effect y.}

\item{LR.CI}{If TRUE calculates profile likelihood confidence intervals based on the likelihood ratio statistic.}

\item{Wald.CI}{If TRUE calculates the Wald confidence intervals from the Hessian matrix of the optimization problem.}

\item{selection.benefit}{A character specifying the selection function to be used for the ORB-adjustment of a beneficial outcome. The options are "Copas.twoside", "Copas.oneside" (default), "Neg.exp.piecewise", and "Sigmoid.cont". See details.}

\item{weight_param}{A positive number to be included if the selection function from selection.benefit is "Neg.exp.piecewise" or "Sigmoid.cont". See details.}

\item{opt_method}{A character specifying the method used in the optim() function for the joint optimization of the global treatment effect $mu$ and the heterogeneity variance $tau^2$. For numerical reasons, the reccomended is "L-BFGS-B" when using "Neg.exp.piecewise" or "Sigmoid.cont" as selection.benefit, otherwise "Nelder-Mead".}

\item{lower}{A vector of length two containing the lower bounds for the optimization problem, if "L-BFGS-B" is used. The first is the lower bound for the global treament effect $mu$ and the second for the heterogeneity variance $tau^2$. }

\item{upper}{A vector of length two containing the upper bounds for the optimization problem, if "L-BFGS-B" is used. The first is the upper bound for the global treament effect $mu$ and the second for the heterogeneity variance $tau^2$.}

\item{low.risk}{If TRUE, the studies for which the missing outcome is classified as low risk of bias, i.e., "low", are also included in the adjustment. If FALSE, only the "high" are used, and the "low" are disregarded. It is reccomended to use both "high" and "low" when selection.function is "Neg.exp.piecewise" or "Sigmoid.cont", and only "high" when selection.function is "Copas.twoside" or "Copas.oneside", as originally proposed by Copas et al. (2019). See details.}

}

\details{

For a beneficial outcome, difference underlying missing data mechanisms can be assumed, which result in a different weight attributed to the ORB-adjusted log-likelihood.

The "Copas.twoside" assumes, among the unreported studies classified as high risk of bias ("high"), the following probability of unreporting \eqn{\Omega_0(y_i, z_{\alpha}\sigma_i)} , where \eqn{z_{\alpha} = \text{qnorm}(1-\alpha/2)}. The underlying ORB mechanism proposed by Copas et al. (2019) assumes that, for a beneficial outcome, an unrepored study classified as "high", is missing because it was not significant, assuming two sided siginficance.

\deqn{\Omega_0(y_i, z_{\alpha}\sigma_i) = \begin{cases}0 &  y_i \in \{ -\infty, -z_{\alpha}\sigma_i \} \cup \{ z_{\alpha}\sigma_i, + \infty \}  \\
1 & y_i \in  \{ -z_{\alpha}\sigma_i, +z_{\alpha}\sigma_i \} \end{cases}}

The "Copas.oneside" is an adaptation of the above, where the one-sided selection function is used, under the assumption that for a positive direction of treatment, a significant effect but in the "wrong" direction is very unlikely to be repotred. Hence the probability of unreporting \eqn{\Omega_1(y_i, z_{\alpha}\sigma_i)} is:

\deqn{  \Omega_1(y_i, z_{\alpha}\sigma_i) =  \begin{cases}1 & y_i \in \{ - \infty, z_{\alpha}\sigma_i \} \\
0 & y_i \in \{ z_{\alpha}\sigma_i, + \infty \} \end{cases}  }


}
\value{

Returns a list containing the following results:

\item{mu_unadjusted}{Maximum likelihood (ML) estimate of the treatment effect, unadjusted for ORB}

\item{LR_mu_unadjusted_low}{Lower bound of likelihood ratio (LR) confidence interval (CI) for unadjusted estimate}

\item{LR_mu_unadjusted_up}{Upper bound of likelihood ratio (LR) confidence interval (CI) for unadjusted estimate}

\item{CI_unadjusted_low_WALD}{}

\item{CI_unadjusted_up_WALD}{}

\item{mu_adjusted_benefit}{ML estimate of the treatment effect, adjusted for ORB, assuming beneficial outcome}

\item{mu_adjusted_harm}{}

\item{LR_mu_adjusted_low}{Lower bound of likelihood ratio (LR) confidence interval (CI) for adjusted estimate}

\item{LR_mu_adjusted_up}{Upper bound of likelihood ratio (LR) confidence interval (CI) for adjusted estimate}

\item{tau_squared_unadjusted}{ML estimate of heterogeneity parameter \eqn{\tau^2}, from unadjusted log-likelihood function}

\item{tau_squared_adjusted}{ML estimate of heterogeneity parameter \eqn{\tau^2}, from adjusted log-likelihood function}

\item{tau_squared_unadjusted_REML}{}

\item{tau_squared_adjusted_REML}{}

\item{average_sigma_squared_unadjusted}{}

\item{sigma_squared_average_adjusted}{}

\item{CI_adjusted_benefit_low_WALD}{}

\item{CI_adjusted_benefit_up_WALD}{}

%% ...
}
\references{




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

Generate a meta-analysis dataset with \code{simualte_meta_data} with ORB.
dat <- simulate_meta_data(n_studies=10, mu=0.1, tau_squared = 0.002, seed=123456789, gamma=1.5)$meta_data_missing

Apply the reORBgen function.
reORBgen(y=dat$y, s=dat$sd, n1=dat$`N treatment`, n2=dat$`N control`, outcome="benefit", init_param = c(0.1, 0.005), alpha=0.05, true.SE=NULL, LR.CI=TRUE)

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
