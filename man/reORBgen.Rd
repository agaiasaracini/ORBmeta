\name{reROBgen}
\alias{reROBgen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to adjust for outcome reporting bias (ORB) in a random effects meta-analysis, according to the methodology developed by Copas et al. (2019) with extensions.
}
\description{
An implementation of the Copas et al. (2019) method to adjust for ORB in meta-anlyses of clinical trials for the random effects model with extention to different selection functions to model the ORB mechanism (see Details). The \code{reORBgen} function can take in binary outcome data (i.e., the user should provide \code{a}, \code{c}, \code{n1}, \code{n2}), means and standard errors (i.e., the user should provide \code{mu1}, \code{mu2}, \code{sd1}, \code{sd2}, \code{n1}, \code{n2}), or a treatment effect assumed to be normally distributed (i.e., the user should provide \code{y}, \code{s}, \code{n1}, \code{n2}). The outcome variable is reported for some studies, and unreported for others. The unreported studies should have entries "high", i.e., high risk of bias, or "low", i.e., low risk of bias. The risk of bias classification should be done according to the Outcome Reporting Bias in Trials (ORBIT) methodology. For the method to work, the sample sizes are required for all studies (reported and unreported). The function can adjust for ORB in both beneficial and harmful outcomes.
}
\usage{
reORBgen(
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
opt_method = "L-BFGS-B",
lower = c(-5, 0.00001),
upper = c(5,5),
low.risk = FALSE
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

\item{init_param}{A vector of length two, containing the values to be used for initialization of the joint optimization of the global treatment effect \eqn{\mu} and the between-study heterogeneity \eqn{\tau^2} of the random effect model, e.g., initi_param = c(0.1, 0.05).}

\item{alpha_ben}{The value for the desired confidence level if the outcome of interest is beneficial, e.g., \eqn{\alpha=0.05} in Copas et al. (2019).}

\item{alpha_ben_one.sided}{If TRUE, the one-sided significance threshold \eqn{\text{qnorm}(1-\alpha)} is used for the selection functions (see Details); if FALSE, the two-sided one is used \eqn{\text{qnorm}(1-\alpha/2)}. Only needed for beneficial outcomes. Reccomended to use alpha_bet_one.sided = FALSE for all selection functions.}

\item{alpha_harm}{The value for the desired confidence level if the outcome of interest is harmful, e.g., \eqn{\alpha=0.01} in Copas et al. (2019).}

\item{true.SE}{A vector of the observed standard errors of the treatment effects for all studies in the meta-analysis, including the ones for which we have "high" or "low", if known. Of note, this works if the outcome data is provided as an observed treatment effect y.}

\item{LR.CI}{If TRUE calculates profile likelihood (PL) confidence intervals (CIs) based on the likelihood ratio statistic.}

\item{Wald.CI}{If TRUE calculates the Wald confidence intervals from the Hessian matrix of the optimization problem.}

\item{selection.benefit}{A character specifying the selection function to be used for the ORB-adjustment of a beneficial outcome. The options are "Copas.twoside", "Copas.oneside" (default), "Neg.exp.piecewise", and "Sigmoid.cont". See Details. Only needed for beneficial outcomes.}

\item{weight_param}{A positive number to be included if the selection function from selection.benefit is "Neg.exp.piecewise" or "Sigmoid.cont". See details. Only needed for beneficial outcomes.}

\item{opt_method}{A character specifying the method used in the \eqn{\texttt{optim()}} function for the joint optimization of the global treatment effect \eqn{\mu} and the heterogeneity variance \eqn{\tau^2}. For numerical reasons, the reccomended is "L-BFGS-B" when using "Neg.exp.piecewise" or "Sigmoid.cont" as selection.benefit, otherwise "Nelder-Mead".}

\item{lower}{A vector of length two containing the lower bounds for the optimization problem, if "L-BFGS-B" is used. The first is the lower bound for the global treament effect \eqn{\mu} and the second for the heterogeneity variance \eqn{\tau^2}. }

\item{upper}{A vector of length two containing the upper bounds for the optimization problem, if "L-BFGS-B" is used. The first is the upper bound for the global treament effect \eqn{\mu} and the second for the heterogeneity variance \eqn{\tau^2}.}

\item{low.risk}{If TRUE, the studies for which the missing outcome is classified as low risk of bias, i.e., "low", are also included in the adjustment. If FALSE, only the "high" are used, and the "low" are disregarded. It is reccomended to use both "high" and "low" when selection.function is "Neg.exp.piecewise" or "Sigmoid.cont", and only "high" when selection.function is "Copas.twoside" or "Copas.oneside", as originally proposed by Copas et al. (2019). See details. Only needed for beneficial outcomes.}

}

\details{

For a beneficial outcome, different underlying missing data mechanisms can be assumed, which result in a different weight attributed to the ORB-adjusted log-likelihood.

The "Copas.twoside" from Copas et al. (2019) assumes, among the unreported studies classified as high risk of bias ("high"), the probability of unreporting \eqn{\Omega_0(y_i, z_{\alpha}\sigma_i)} , where \eqn{z_{\alpha} = \text{qnorm}(1-\alpha/2)}. The underlying ORB mechanism proposed by Copas et al. (2019) assumes that, for a beneficial outcome, an unrepored study classified as "high" is missing because it was not significant, assuming two-sided siginficance. Using "Copas.twoside" implies that only the "high" unreported studies are used, thus low.risk = FALSE.

\deqn{\Omega_0(y_i, z_{\alpha}\sigma_i) = \begin{cases}0 &  y_i \in \{ -\infty, -z_{\alpha}\sigma_i \} \cup \{ z_{\alpha}\sigma_i, + \infty \}  \\
1 & y_i \in  \{ -z_{\alpha}\sigma_i, +z_{\alpha}\sigma_i \} \end{cases}}

The "Copas.oneside" is an adaptation of the above, where the one-sided probability of unreporting \eqn{\Omega_1(y_i, z_{\alpha}\sigma_i)} is used, under the assumption that, for a positive direction of treatment, a significant effect, but in the "wrong"/harmful direction, is very unlikely to be repotred. Using "Copas.oneisde" implies that only the "high" unreported studies are used, thus low.risk = FALSE.

\deqn{  \Omega_1(y_i, z_{\alpha}\sigma_i) =  \begin{cases}1 & y_i \in \{ - \infty, z_{\alpha}\sigma_i \} \\
0 & y_i \in \{ z_{\alpha}\sigma_i, + \infty \} \end{cases}  }

The "Neg.exp.piecewise" is a generalization of the above piece-wise probabilities of unreporting, which assumes that, among all unerported studies ("high" and "low"), a non-significant one has probability of unreporing = 1, and a significant one (one-sided) is a negative exponential function with rate parameter \eqn{\rho > 0}, instead of being 0 as the above. The function, \eqn{\Omega_2(y_i, z_{\alpha}\sigma_i, \lambda)} is still piece-wise, but smooth after the \eqn{z_{\alpha}\sigma_i} threshold for significance. Using "Neg.exp.piecewise" implies that both "high" and "low" unreported studies are used, thus low.risk = TRUE.

\deqn{\Omega_2(y_i, z_{\alpha}\sigma_i, \lambda) =  \begin{cases}1 & y_i \in \{ - \infty, z_{\alpha}\sigma_i \} \\
e^{-\rho \cdot (y_i - z_{\alpha}\sigma_i)} & y_i \in \{ z_{\alpha}\sigma_i, + \infty \} \end{cases} }


The "Sigmoid.cont" is a further generalizaiton, in which the probability of an unreported study not reporting the outcome  is a continous decreasing function of the observed treatment effect with rate parameter \eqn{\lambda > 0}. At the threshold for one-sided significance \eqn{z_{\alpha}\sigma_i} the probability of unreporting is 0.5. Using "Sigmoid.cont" implies that both "high" and "low" unreported studies are used, thus low.risk = TRUE.

\deqn{\Omega_3(y_i, z_{\alpha}\sigma_i, \lambda) =  \frac{1}{1+ e^{\lambda(y_i- z_{\alpha}\sigma_i)}}}


For a harmful outcome, the following missing data mechanism is assumed, as proposed originally in Copas et al. (2019): among the unreported studies at high risk of bias ("high"), an unreported study is missing if its treatment effect was positive, assuming that a positive value means that the treatment has a negative/harmful effect on that outcome. The probability of unreporting is thus

\deqn{  \Omega_4(y_i) =  \begin{cases}0 & y_i > 0 \\
1 & y_i < 0 \end{cases}  }


}
\value{

Returns a list containing the following results:

\item{mu_unadjusted}{Maximum likelihood (ML) estimate of the treatment effect, unadjusted for ORB}

\item{LR_mu_unadjusted_low}{Lower bound of likelihood ratio (LR) confidence interval (CI) for unadjusted estimate}

\item{LR_mu_unadjusted_up}{Upper bound of the LR CI for unadjusted estimate}

\item{CI_unadjusted_low_WALD}{Lower bound of Wald CI for unadjusted estimate}

\item{CI_unadjusted_up_WALD}{Upper bound of Wald CI for unadjusted estimate}

\item{mu_adjusted_benefit}{ML estimate of the treatment effect, adjusted for ORB, assuming beneficial outcome}

\item{mu_adjusted_harm}{ML estimate of the treatment effect, adjusted for ORB, assuming harmful outcome}

\item{LR_mu_adjusted_low}{Lower bound of the LR CI for adjusted estimate}

\item{LR_mu_adjusted_up}{Upper bound of the LR CI for adjusted estimate}

\item{tau_squared_unadjusted}{ML estimate of heterogeneity parameter \eqn{\tau^2}, from unadjusted log-likelihood function}

\item{tau_squared_adjusted}{ML estimate of heterogeneity parameter \eqn{\tau^2}, from adjusted log-likelihood function}

\item{tau_squared_unadjusted_REML}{REML estimate of heterogeneity parameter \eqn{\tau^2}, from unadjusted log-likelihood function}

\item{tau_squared_adjusted_REML}{REML estimate of heterogeneity parameter \eqn{\tau^2}, from adjusted log-likelihood function}

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

#Apply the reORBgen function to the simulated beneficial outcome.

ORBmeta:::reORBgen(
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
init_param = c(0.1, 0.005),
true.SE = NULL,
LR.CI = TRUE,
Wald.CI = TRUE,
selection.benefit = "Neg.exp.piecewise",
weight_param = 7,
opt_method = "L-BFGS-B",
lower = c(-5, 0.00001),
upper = c(5,5),
low.risk = TRUE)


#Use a beneficial outcome from the EpilepsyData.


ORBmeta:::reORBgen(
a = EpilepsyData$`Seizure Free (T)`,
c = EpilepsyData$`Seizure Free (C)`,
mu1 = NULL,
mu2 = NULL,
sd1 = NULL,
sd2 = NULL,
y = NULL,
s = NULL,
n1 = EpilepsyData$`Sample (T)`,
n2 = EpilepsyData$`Sample (C)`,
outcome = "benefit",
init_param = c(0.5, 0.1),
alpha_ben = 0.05,
alpha_ben_one.sided = FALSE,
alpha_harm = NULL,
true.SE = NULL,
LR.CI = TRUE,
Wald.CI = FALSE,
selection.benefit = "Neg.exp.piecewise",
weight_param = 7,
opt_method = "L-BFGS-B",
lower = c(-5, 0.00001),
upper = c(5,5),
low.risk = TRUE)


#Use a harmful outcome from the EpilepsyData.


ORBmeta:::reORBgen(
a = EpilepsyData$`Dizziness (T)`,
c = EpilepsyData$`Dizziness (C)`,
mu1 = NULL,
mu2 = NULL,
sd1 = NULL,
sd2 = NULL,
y = NULL,
s = NULL,
n1 = EpilepsyData$`Sample (T)`,
n2 = EpilepsyData$`Sample (C)`,
outcome = "harm",
init_param = c(0.5, 0.1),
alpha_ben = NULL,
alpha_ben_one.sided = FALSE,
alpha_harm = 0.01,
true.SE = NULL,
LR.CI = TRUE,
Wald.CI = FALSE,
opt_method = "L-BFGS-B",
lower = c(-5, 0.00001),
upper = c(5,5),
low.risk = FALSE)


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
