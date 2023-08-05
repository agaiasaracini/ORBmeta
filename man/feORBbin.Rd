\name{feORBbin}
\alias{feORBbin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to adjust for outcome reporting bias (ORB) in a fixed-effect meta-analysis, according to the methodology developed by Copas et al. (2019).
}
\description{
An implementation of the Copas et al. (2019) method to adjust for outcome reporting bias (ORB) in meta-anlyses of clinical trials for the fixed-effect model. The function takes in binary outcome data and calculates the log RR. The studies for which the observed outcome and standard error are missing should have entries "high", i.e., high risk of bias, or "low", i.e., low risk of bias. For the method to work, the sample sizes are required for the studies which do not report the observed treatment effect. The risk of bias classification should be done according to the ORBIT methodology. The function can adjust both beneficial and harmful outcomes.
}
\usage{
feORBbin(
a,
c,
n1,
n2,
outcome,
init_param,
alpha_ben,
alpha_harm,
true.SE=NULL
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

 \item{a}{A vector of observed counts for the outcome of interest in the treatment group for each study in the meta-analysis. For the studies which do not report the counts in the outcome of interest, the entries are either "high" or "low".}

\item{c}{A vector of observed counts for the outcome of interest in the control group for each study in the meta-analysis. For the studies which do not report the counts in the outcome of interest, the entries are either "high" or "low".}

\item{n1}{A vector of sample sizes in the treatment group for each study in the meta-analysis. The sample sizes need to be reported for all studies in the meta-analysis.}

\item{n2}{A vector of sample sizes in the control group for each study in the meta-analysis. The sample sizes need to be reported for all studies in the meta-analysis.}

\item{outcome}{A character indicating whether the outcome of interest is beneficial ("benefit") or harmful ("harm").}

\item{init_param}{A number to be used for the initialization of the optimization of \eqn{\mu}, i.e., the true treatment effect in the fixed effect model}

\item{alpha_ben}{The value for the desired confidence level if the outcome of interest is beneficial, e.g., \eqn{\alpha=0.05} in Copas et al. (2019).}

\item{alpha_harm}{The value for the desired confidence level if the outcome of interest is beneficial, e.g., \eqn{\alpha=0.01} in Copas et al. (2019).}

\item{true.SE}{A vector of the observed standard errors of the treatment effects (log RR) for all studies in the meta-analysis, including the ones for which we have "high" or "low", if known. Of note, this means that the standard errors of the log RR need to be provided.}

}
\details{

The function implements the Copas et al. (2019) methodology for the fixed-effect model, with the only difference that, for a beneficial outcome classified as high risk of bias, the following function for the probability of unreporting is used:

\deqn{  \Omega_1(y_i, z_{\alpha}\sigma_i) =  \begin{cases}1 & y_i \in \{ - \infty, z_{\alpha}\sigma_i \} \\
0 & y_i \in \{ z_{\alpha}\sigma_i, + \infty \} \end{cases}  }

whereas in Copas et al. (2019), the following was used:

\deqn{\Omega_0(y_i, z_{\alpha}\sigma_i) = \begin{cases}0 &  y_i \in \{ -\infty, -z_{\alpha}\sigma_i \} \cup \{ z_{\alpha}\sigma_i, + \infty \}  \\
1 & y_i \in  \{ -z_{\alpha}\sigma_i, +z_{\alpha}\sigma_i \} \end{cases}}

Here, \eqn{z_{\alpha} = \texttt{qnorm}(1-\alpha/2)}.

For a harmful outcome, that is unreported and classified as high risk of bias it is assumed that the probability of an outcome to be unreported is:

\deqn{  \Omega_4(y_i) =  \begin{cases}0 & y_i > 0 \\
1 & y_i < 0 \end{cases}  }



}
\value{

\item{logRR_unadjusted}{MLE of the logRR common effect, unadjusted for ORB}
\item{CI_unadjusted_low}{Lower bound of the likelihood ratio (LR) confidence interaval (CI) for the unadjusted estimate}
\item{CI_unadjusted_up}{Upper bound of the LR CI for the unadjusted estimate}
\item{CI_unadjusted_low_WALD}{Lower bound of the Wald CI for the unadjusted estimate}
\item{CI_unadjusted_up_WALD}{Upper bound of the Wald CI for the unadjusted estimate}
\item{logRR_adjusted_benefit}{MLE of the logRR common effect, adjusted for ORB when the outcome is beneficial}
\item{logRR_adjusted_harm}{MLE of the logRR common effect, adjusted for ORB when the outcome is harmful}
\item{CI_adjusted_benefit_low}{Lower bound of the likelihood ratio (LR) confidence interaval (CI) for the adjusted estimate of a beneficial outcome}
\item{CI_adjusted_benefit_up}{Upper bound of the likelihood ratio (LR) confidence interaval (CI) for the adjusted estimate of a beneficial outcome}
\item{CI_adjusted_harm_low}{Lower bound of the likelihood ratio (LR) confidence interaval (CI) for the adjusted estimate of a harmful outcome}
\item{CI_adjusted_harm_up}{Upper bound of the likelihood ratio (LR) confidence interaval (CI) for the adjusted estimate of a harmful outcome}
\item{CI_adjusted_benefit_low_WALD}{Lower bound of the Wald CI for the adjusted estimate of a beneficial outcome}
\item{CI_adjusted_benefit_up_WALD}{Upper bound of the Wald CI for the adjusted estimate of a beneficial outcome}
\item{CI_adjusted_harm_low_WALD}{Lower bound of the Wald CI for the adjusted estimate of a harmful outcome}
\item{CI_adjusted_harm_up_WALD}{Upper bound of the Wald CI for the adjusted estimate of a harmful outcome}
\item{logRR}{Individual logRR of the observed treatment effects for each study}
\item{p_value}{Individual one-sided p-values for each study}

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

data(EpilepsyData)

#Apply the reORBgen function to a beneficial outcome.

ORBmeta:::feORBbin(
a=EpilepsyData$`Seizure Free (T)`,
c=EpilepsyData$`Seizure Free (C)`,
n1=EpilepsyData$`Sample (T)`,
n2=EpilepsyData$`Sample (C)`,
outcome="benefit",
alpha_ben=0.05,
init_param = 0.1,
true.SE=NULL
)

#Apply the reORBgen function to a harmful outcome.

ORBmeta:::reORBgen(
a=EpilepsyData$`Nausea/vomiting (T)`,
c=EpilepsyData$`Nausea/vomiting (C)`,
n1=EpilepsyData$`Sample (T)`,
n2=EpilepsyData$`Sample (C)`,
outcome="harm",
alpha_harm=0.01,
init_param = 0.1,
true.SE=NULL
)

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
