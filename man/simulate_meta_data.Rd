\name{simulate_meta_data}
\alias{simulate_meta_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate Meta-Analysis Data with ORB
}
\description{
The \code{simulate_meta_data} function generates simulated data for a meta-analysis for a beneficial outcome of a randomized controlled trials (RCTs) so as to include ORB. It creates a dataset of studies with the following variables: sample sizes in treatment and control, observed treatment effects, standard errors, and one-sided p-values.

First, it generates a meta-analysis according to the random-effect model, where all studies have the same sample sizes for treatment and control. Then, it removes studies to simulate ORB and replaces some treatment effects and standard errors with high risk of bias ("high") or low risk of bias ("low").

Unreported studies at high risk of bias ("high") are removed according to the following selection function: \deqn{exp(-4 \cdot p_i^\gamma)} used in Copas and Henmi (2010), representing the probability of reporting the outcome as a function of the one-sided p-value: \deqn{p_i = \Phi(- \frac{y_i}{\sigma_i})}.

Among the studies which are left reported, a number between 1 and the number of "high" ones is chosen to classified as low risk of bias ("low").

If there are less than 3 reported studies, the function chooses at random among the "high" and "low" ones some to be turned back into reported, so as to guarantee that at least 3 are reported.


}
\usage{
simulate_meta_data(
n_studies,
mu,
tau_squared,
n_treatment=50,
n_control=50,
gamma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{n_studies}{
The number of studies in the meta-analysis.
}
\item{mu}{
The true treatment effect \eqn{\mu}.
}
\item{tau_squared}{
The true heterogeneity parameter \eqn{\tau^2}.
}
\item{n_treatment=50}{
The sample size in the control. This will be the same across all studies.
}
\item{n_control=50}{
The sample size in the treatment. This will be the same across all studies.
}

\item{gamma}{
The parameter for ORB modeling, controlling the probability of studies reporting "high". Higher values of gamma result in a higher probability of reporting "high". Based on Copas and Henmi (2010).
}


}
\details{

The y treatment effect is assumed to for a beneficial outcomes, where a positive/significant value indicates a beneficial effect of treatment on that outcome.

}
\value{

The function returns a dataset with the following columns
\item{study}{
Study IDs, i.e., a number from 1 to the number of studies chosen for the simualted meta-analysis.
}
\item{n_control}{
Sizes of control groups, equal across all studies.
}

\item{n_treatment}{
Sizes of treatment groups, equal across all studies.
}

\item{y}{Observed treatment effect. It will be reported for some studies, and "high" or "low" for others.}

\item{s}{Observed standard error. It will be reported for some studies, and "high" or "low" for others.}

\item{p_value}{The one-sided p-values used to generate "high" studies.}

\item{s_true}{The true standard errors before ORB is simualted and some are removed.}

}
\references{


Copas, J. and Henmi, M (2010).Confidence intervals for random effects meta-analysis and robustness to publication bias. Statistics in Medicine, 29 (29), 2969-3068.

Copas, J., Marson, A., Williamson, P., and Kirkham, J. (2019). Model-based sensitivity analysis for outcome reporting bias in the meta analysis of benefit and harm outcomes. Statistical Methods in Medical Research, 28, 889–903



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


#Generate a meta-analysis dataset for beneficial outcome with ORB.

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
