\name{PPHData}
\alias{PPHData}
\docType{data}
\title{
Meta-analysis of prophylactic use of ergot alkaloids for prevention of post partum hemmorage (PPH), in the presence of outcome reporting bias (ORB) and risk of bias classification done according to the ORBIT methodology.
}
\description{
The data is originally from the meta-analysis study by Liabsuetrakul et al. (2007), and is used in Copas et al. (2014) to shwcase their ORB adjustment method. The study investigates the effectiveness and safety of prophylactic use of ergot alkaloids in the third stage of labor compared with non-uterotonic agents for the prevention of PPH. The study compares control to treatment and treats PPH as a beneficial outcome. The meta-analysis includes 4 studies which report the number of cases of PPH in the control and treatment groups. Out of the 5 studies which do not report the outcome, 4 are classified as low risk of bias ("low"), and 1 as high risk of bias ("high") by the ORBIT methodology.
}
\usage{data("PPHData")}
\format{
  A data frame with 9 observations on the following 5 variables.
  \describe{
    \item{\samp{Sample (C)}}{Sample size in the control group}
    \item{\samp{Sample (T)}}{Sample size in the treatment group}
    \item{\samp{Sample Total}}{Total sample size in the studies}
    \item{\samp{PPH (C)}}{Cases of PPH in the control group}
    \item{\samp{PPH (T)}}{Cases of PPH in the treatment group}
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

Liabsuetrakul, T., Choobun, T., Peeyananjarassri, K., and Islam, M. (2007). Prophylactic use of ergot alkaloids in the third stage of labour. Cochrane Database of Systematic Reviews.

}
\examples{

data(PPHData)

Apply the ORBsens function to the data

ORBsens(
p_vals = c(0.7, 0.7),
eta_vals = NULL,
a = PPHData$`Sample (C)`,
c = PPHData$`Sample (T)`,
mu1 = NULL,
mu2 = NULL,
sd1 = NULL,
sd2 = NULL,
n1 = PPHData$`Sample (C)`,
n2 = PPHData$`Sample (T)`,
ntot = PPHData$`Sample Total`,
sign_level = 0.05,
init_param = c(0.7, 0.6, 0.9, 0.9),
outcome = "benefit",
opt_method = "L-BFGS-B",
lower = c(-10, 0.05, 0.05, 0.05),
upper = c(10, 10,1,1))



}
\keyword{datasets}
