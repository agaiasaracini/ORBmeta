\name{reROBgen}
\alias{reROBgen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to adjust for outcome reporting bias (ORB) in an outcome of a meta-analysis of clinical trials
}
\description{
An implementation of the Copas et al. (2019) method to adjust for outcome reporting bias (ORB) in meta-anlyses of clinical trials. The function takes in outcome data from a variable that is assumed to be normally distributed (i.e., observed treatment effect, standard error, and sample size of the treatment arms). The outcome variable is reported for some studies, and unreported for others. The studies for which the observed outcome and standard error are missing should have entries "high", i.e., high risk of bias. For the method to work, the sample sizes are required for the studies which do not report the observed treatment effect. The risk of bias classification should be done according to the ORBIT methodology (cite website).
}
\usage{
reROBgen(
y,
s,
n1,
n2,
outcome,
init_param,
alpha,
true.SE=NULL,
LR.CI = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{A vector of observed treatment effects. For the studies which do not report the outcome of interest, the entry is "high".}
  \item{s}{A vector of corresponding standard errors. For  the studies which do not report the outcome, the entry is "high". }
  \item{n1}{A vector of the corresponding sample size of the treatment group. Needs to be provided even for the studies which do not report the outcome.}
  \item{n2}{A vector of the corresponfing sample size of the control group. Needs to be provided even for the studies which do not report the outcome.}
  \item{outcome}{A character indicating whether the outcome of interest is beneficial ("benefit") or harmful ("harm").}
  \item{init_param}{A vector of length two, containing the values to be used for initialization of the joint optimization of the gloval treatment effect $mu$ and the between-study heterogeneity $tau^2$ of the random effect model, e.g., initi_param = c(0.1, 0.05).}
  \item{alpha}{The value for the desired confidence level. Default $alpha=0.05$.}
  \item{true.SE=NULL}{Takes in a vector with the true standard errors of all studies, including the missing ones, if known. }
  \item{LR.CI = TRUE}{Calcultes likelihood ratio (LR) confidence intervals (CI)}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

%%\item{comp1 }{Description of 'comp1'}
%%\item{comp2 }{Description of 'comp2'}
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
function (x)
{
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
