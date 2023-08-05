\name{WeightData}
\alias{WeightData}
\docType{data}
\title{
Meta-analysis of the effect of school feeding on weight gain, in the presence of outcome reporting bias (ORB) and risk of bias classification done according to the ORBIT methodology.
}
\description{
The data is originally from the meta-analysis study by  Kristjansson et al. (2007), and is used in Copas et al. (2014) to shwcase their ORB adjustment method. The study investigates the effect of school feeding for improving the phsysical and psychosocial health of disadvantaged students. The meta-analysis includes 3 studies which report the mean weight gain (seen as a beneficial outcome) and observed standard deviation in the treatment and control groups. All 3 studies which do not report the outcome are classified as high risk of bias ("high") by the ORBIT methodology.
}
\usage{data("WeightData")}
\format{
  A data frame with 6 observations on the following 7 variables.
  \describe{
    \item{\code{n1}}{Sample size in the treatment group}
    \item{\code{n2}}{Sample size in the control group}
    \item{\code{totn}}{Total sample size}
    \item{\code{mean1}}{Average weight gain in the treatment group}
    \item{\code{sd1}}{Observed standard deviation in the treatment group}
    \item{\code{mean2}}{Average weight gain in the control group}
    \item{\code{sd2}}{Observed standard deviation in the control group}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{

Copas, J., Dwan, K., Kirkham, J., and Williamson, P. (2014). A model-based correction for outcome reporting bias in meta-analysis. Biostatistics, 15, 370–383.

Kristjansson, B., Petticrew, M., MacDonald, B., Krasevec, J., Janzen, L., Greenhalgh, T., Wells, G., MacGowan, J., Farmer, A., Shea, B., Tugwell, P., and Welch, V. (2007). School feeding for improving the physical and psychosocial health of disadvantaged students. Cochrane Database of Systematic Reviews

}
\examples{

data(WeightData)

#Apply the ORBsens function to the data

ORBmeta:::ORBsens(
p_vals = c(0.7, 0.8),
eta_vals = NULL,
a = NULL,
c = NULL,
mu1 = WeightData$mean1,
mu2 = WeightData$mean2,
sd1 = WeightData$sd1,
sd2 = WeightData$sd2,
n1 = WeightData$n1,
n2 = WeightData$n2,
ntot = WeightData$totn,
sign_level = 0.05,
init_param = c(0.5, 0.6, 0.7, 0.7),
outcome = "benefit",
opt_method = "L-BFGS-B",
lower = c(-1, 0.001, 0.01, 0.01),
upper= c(1, 5, 0.99, 0.99))



}
\keyword{datasets}
