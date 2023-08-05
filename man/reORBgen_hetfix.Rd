\name{reORBgen_hetfix}
\alias{reORBgen_hetfix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Version of the reORBgen() function when the heterogeneity variance \eqn{\tau^2} is known.
}
\description{

An implementation of the Copas et al. (2019) method to adjust for outcome reporting bias (ORB) in meta-anlyses of clinical trials for the random-effectt model with extention to different selection functions to model the ORB mechanism, when treating the \eqn{\tau^2} as a known and fixed parameter. Hence, all optimization methods are one-dimesional only for \eqn{\mu}.

}
\usage{

reORBgen_hetfixed(
a = NULL,
c = NULL,
mu1 = NULL,
mu2 = NULL,
sd1 = NULL,
sd2 = NULL,
y = NULL,
s = NULL,
n1,
n2,
outcome,
init_param,
alpha_ben = NULL,
alpha_ben_one.sided = TRUE,
alpha_harm = NULL,
true.SE = NULL,
LR.CI = FALSE,
Wald.CI = FALSE,
selection.benefit = "Copas.oneside",
weight_param = 15,
opt_method = "Brent",
lower = -10,
upper = 10,
low.risk = FALSE,
tau_squared_fixed
)


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

\item{init_param}{A numer value, to be used for initialization of the one-dimensional optimization of the global treatment effect \eqn{\mu}, e.g., initi_param = 0.5.}

\item{alpha_ben}{The value for the desired confidence level if the outcome of interest is beneficial, e.g., \eqn{\alpha=0.05} in Copas et al. (2019).}

\item{alpha_ben_one.sided}{If TRUE, the one-sided significance threshold \eqn{\text{qnorm}(1-\alpha)} is used for the selection functions (see details); if FALSE, the two-sided one is used \eqn{\text{qnorm}(1-\alpha/2)}. Only needed for beneficial outcomes. Reccomended to use alpha_bet_one.sided = FALSE for all selection functions.}

\item{alpha_harm}{The value for the desired confidence level if the outcome of interest is beneficial, e.g., \eqn{\alpha=0.01} in Copas et al. (2019).}

\item{true.SE}{A vector of the observed standard errors of the treatment effects for all studies in the meta-analysis, including the ones for which we have "high" or "low", if known. Of note, this works if the outcome data is provided as an observed treatment effect y.}

\item{LR.CI}{If TRUE calculates likelihood ratio confidence intervals based on the likelihood ratio statistic.}

\item{Wald.CI}{If TRUE calculates the one-dimensional Wald confidence intervals.}

\item{selection.benefit}{A character specifying the selection function to be used for the ORB-adjustment of a beneficial outcome. The options are "Copas.twoside", "Copas.oneside" (default), "Neg.exp.piecewise", and "Sigmoid.cont". See details. Only needed for beneficial outcomes.}

\item{weight_param}{A positive number to be included if the selection function from selection.benefit is "Neg.exp.piecewise" or "Sigmoid.cont". See details. Only needed for beneficial outcomes.}

\item{opt_method}{A character specifying the method used in the \eqn{\texttt{optim()}} function for the one-dimensional optimization of the global treatment effect \eqn{\mu}. The default is "Brent".}

\item{lower}{A numerber value for the lower bound for the optimization problem for \eqn{\mu}.}

\item{upper}{A numerber value for the upper bound for the optimization problem for \eqn{\mu}.}

\item{low.risk}{If TRUE, the studies for which the missing outcome is classified as low risk of bias, i.e., "low", are also included in the adjustment. If FALSE, only the "high" are used, and the "low" are disregarded. It is reccomended to use both "high" and "low" when selection.function is "Neg.exp.piecewise" or "Sigmoid.cont", and only "high" when selection.function is "Copas.twoside" or "Copas.oneside", as originally proposed by Copas et al. (2019). See details. Only needed for beneficial outcomes.}

\item{tau_squared_fixed}{A number to be used for the true underlying heterogeneity variance \eqn{\tau^2}}

}
\details{

See details of reORBgen()

}
\value{

Returns a list containing the following results:

\item{mu_unadjusted}{Maximum likelihood (ML) estimate of the treatment effect, unadjusted for ORB}

\item{LR_mu_unadjusted_low}{Lower bound of likelihood ratio (LR) confidence interval (CI) for unadjusted estimate}

\item{LR_mu_unadjusted_up}{Upper bound of likelihood ratio (LR) confidence interval (CI) for unadjusted estimate}

\item{CI_unadjusted_low_WALD}{Lower bound of Wald CI for unadjusted estimate}

\item{CI_unadjusted_up_WALD}{Upper bound of Wald CI for unadjusted estimate}

\item{mu_adjusted_benefit}{ML estimate of the treatment effect, adjusted for ORB, assuming beneficial outcome}

\item{mu_adjusted_harm}{ML estimate of the treatment effect, adjusted for ORB, assuming harmful outcome}

\item{LR_mu_adjusted_low}{Lower bound of likelihood ratio (LR) confidence interval (CI) for adjusted estimate}

\item{LR_mu_adjusted_up}{Upper bound of likelihood ratio (LR) confidence interval (CI) for adjusted estimate}

\item{average_sigma_squared_unadjusted}{Average \eqn{\sigma_{i}^2} across the reported studies}

\item{sigma_squared_average_adjusted}{Average \eqn{\sigma_{i}^2} across the reported and unreported studies, where the \eqn{\sigma_{i}^2} for unreported studies is imputed according to the method of Copas et al. (2019)}

\item{CI_adjusted_benefit_low_WALD}{Lower bound of Wald CI for beneficial outcome, adjusted for ORB}

\item{CI_adjusted_benefit_up_WALD}{Upper bound of Wald CI for beneficial outcome, adjusted for ORB}

\item{CI_adjusted_harm_low_WALD}{Lower bound of Wald CI for harmful outcome, adjusted for ORB}

\item{CI_adjusted_harm_up_WALD}{Upper bound of Wald CI for harmful outcome, adjusted for ORB}

%% ...


}
\references{

Copas, J., Marson, A., Williamson, P., and Kirkham, J. (2019). Model-based sensitivity analysis for outcome reporting bias in the meta analysis of benefit and harm outcomes. Statistical Methods in Medical Research, 28, 889–903.


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


#Generate a meta-analysis dataset for beneficial outcome with ORB using the ORBmeta:::simualte_meta_data() function in this package .

set.seed(123)

dat <- ORBmeta:::simulate_meta_data(
n_studies = 15,
mu = 0.1,
n_treatment = 50,
n_control = 50,
tau_squared = 0.002,
gamma = 1.5)

#Apply the reORBgen function to the simulated beneficial outcome, using the true underlying heterogeneity which is known.

ORBmeta:::reORBgen_hetfixed(
a = NULL,
c = NULL,
mu1 = NULL,
mu2 = NULL,
sd1 = NULL,
sd2 = NULL,
y = dat$y,
s = dat$s,
n1 = dat$n_treatment,
n2 = dat$n_control,
outcome = "benefit",
alpha_ben = 0.05,
alpha_ben_one.sided = FALSE,
init_param = 0.1,
true.SE = NULL,
LR.CI = TRUE,
Wald.CI = TRUE,
selection.benefit = "Neg.exp.piecewise",
weight_param = 7,
opt_method = "Brent",
lower = -10,
upper = 10,
low.risk = TRUE,
tau_squared_fixed = 0.002)



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
