\name{simulate_meta_data}
\alias{simulate_meta_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate Meta-Analysis Data with ORB
}
\description{
The \code{simulate_meta_data} function generates simulated data for a meta-analysis of randomized controlled trials (RCTs) with the option to include Outcome Reporting Bias (ORB). It creates a dataset of studies reporing the following characteristics: sample sizes, observed treatment effects, standard errors, and p-values. ORB is simulated according to the following selection function: \deqn{exp(-4 \cdot p_i^\gamma)} used in Copas and Henmi (2010), representing the probability of reporting the outcome as a function of the p-value. Unreported studies, removed acccording to this missing data mechanism, are marked as "high", i.e., "high risk of bias", so as to apply the ORBIT methdology for ORB adjustment used in Copas et al. (2019). The function returns two datasets: one with all the generated data before the ORB simulation, i.e., meta_data, and another with missing data due to ORB,  i.e, meta_data_missing.
}
\usage{
simulate_meta_data(
n_studies,
mu,
tau_squared,
min_n = 50,
max_n = 500,
seed,
gamma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{n_studies}{
The number of studies in the meta-analysis.
}
\item{mu}{
The true treatment effect mean.
}
\item{tau_squared}{
The true heterogeneity parameter.
}
\item{min_n}{
The minimum sample size per treatment arm. Default is 50.
}
\item{max_n}{
The maximum sample size per treatment arm. Default is 500.
}
\item{seed}{
The seed for random number generation. If not provided, the seed will be random.
}
\item{gamma}{
The parameter for ORB modeling, controlling the probability of studies reporting "high". Higher values of gamma result in a higher probability of reporting "high". Based on
}


}
\details{
The function uses the \code{biostatUZH} package for some of the calculations. If the package is not installed, an error will be displayed.

The function ensures that the generated dataset has at least three studies with reported values, as having fewer studies can result in unstable calculations and unrealistic meta-analyses.

The returned datasets include formatted columns for better readability.
}
\value{

The function returns a list with two elements:
\item{meta_data}{
The complete dataset of generated meta-analysis data, including sample sizes, observed treatment effects, standard errors, p-values, etc.
}
\item{meta_data_missing}{
The dataset with missing data due to ORB, where some studies report "high" instead of actual values based on their individual replacement probabilities.
}


}
\references{


Copas, J. and Henmi, M (2010).Confidence intervals for random effects meta-analysis and robustness to publication bias. Statistics in Medicine, 29 (29), 2969-3068.

Copas, J., Marson, A., Williamson, P., and Kirkham, J. (2019). Model-based sensitivity analysis for outcome reporting bias in the meta analysis of benefit and harm outcomes. Statistical Methods in Medical Research, 28, 889–903



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


Simulate a meta-analysis dataset without ORB
meta_data <- simulate_meta_data(n_studies = 50, mu = 0.5, tau_squared = 0.1)

Simulate a meta-analysis dataset with ORB
meta_data_with_orb <- simulate_meta_data(n_studies = 50, mu = 0.5, tau_squared = 0.1, gamma = 2)

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
