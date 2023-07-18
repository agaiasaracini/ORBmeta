\name{EpilepsyData}
\alias{EpilepsyData}
\docType{data}
\title{
Epilepsy Data with ORB
}
\description{
The data is originally from the meta-analysis study by Bresnah et al. (2019), and is used in Copas et al. (2019) to shwcase their ORB adjustment method. The study investigates the effect of Topiramate, an antiepileptic drug first marketed in 1996, when used as an add-on treatment for drug-resistant focal epilepsy. The meta-analysis consists of 12 studies and 14 outcomes, 2 beneficial ones and 12 harmful ones. Each unreported outcome for each study is classified as high risk of bias ("high") or low risk of bias ("low") by the ORBIT methdology.
}
\usage{data("EpilepsyData")}
\format{

  \describe{
    \item{\samp{Sample (T)}}{Sample size in the treatment group}
    \item{\samp{Sample (C)}}{Sample size in the control group}
    \item{\samp{50\% Seizure Reduction (T)}}{Benficial outcome: cases of 50\% seizure reduction in the treatment group}
    \item{\samp{50\% Seizure Reduction (C)}}{Benficial outcome: cases of 50\% seizure reduction in the control group}
    \item{\samp{Seizure Free (T)}}{Beneficial outcome: cases of seizure freedom in the treatment group}
    \item{\samp{Seizure Free (C)}}{Beneficial outcome: cases of seizure freedom in the control group}
    \item{\samp{Trt withdrawl (T)}}{Harmful outcome: cases of treatment withdrawl in the treatment group}
    \item{\samp{Trt withdrawl (C)}}{Harmful outcome: cases of treatment withdrawl in the control group}
    \item{\samp{Dizziness (T)}}{Harmful outcome: cases of diziness in the treatment group}
    \item{\samp{Dizziness (C)}}{Harmful outcome: cases of diziness in the control group}
    \item{\samp{Headache (T)}}{Harmful outcome: cases of headache in the treatment group}
    \item{\samp{Headache (C)}}{Harmful outcome: cases of headache in the control group}
    \item{\samp{Nausea/vomiting (T)}}{Harmful outcome: cases of nausea/vomiting in the treatment group}
    \item{\samp{Nausea/vomiting (C)}}{cases of nausea/vomiting in the control group}
    \item{\samp{Paraesthesias (T)}}{Harmful outcome: cases of paraesthesias in the treatment group}
    \item{\samp{Paraesthesias (C)}}{Harmful outcome: cases of paraesthesias in the control group}
    \item{\samp{Weight loss (T)}}{Harmful outcome: cases of weight loss in the treatment group}
    \item{\samp{Weight loss (C)}}{Harmful outcome: cases of weight loss in the control group}
    \item{\samp{Fatigue (T)}}{Harmful outcome: cases of fatigue in the treatment group}
    \item{\samp{Fatigue (C)}}{Harmful outcome: cases of fatigue in the control group}
    \item{\samp{Somnolence (T)}}{Harmful outcome: cases of somnolence in the treatment group}
    \item{\samp{Somnolence (C)}}{Harmful outcome: cases of somnolence in the control group}
    \item{\samp{Concentration impairment (T)}}{Harmful outcome: cases of concentration impairment in the treatment group}
    \item{\samp{Concentration impairment (C)}}{Harmful outcome: cases of concentration impairment in the control group}
    \item{\samp{Speech difficulty (T)}}{Harmful outcome: cases of speech difficulty in the treatment group}
    \item{\samp{Speech difficulty (C)}}{Harmful outcome: cases of speech difficulty in the control group}
    \item{\samp{Thinking abnormality (T)}}{Harmful outcome: cases of thinking abnormality in the treatment group}
    \item{\samp{Thinking abnormality (C)}}{Harmful outcome: cases of thinking abnormality in the control group}
    \item{\samp{Ataxia (T)}}{Harmful outcome: cases of ataxia in the treatment group}
    \item{\samp{Ataxia (C)}}{Harmful outcome: cases of ataxia in the control group}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{

Bresnah, R., Hounsome, J., Jette, N., Hutton, J., and Marson, A. (2019). Topiramate addâon therapy for drugâresistant focal epilepsy. Cochrane Library.

Copas, J., Marson, A., Williamson, P., and Kirkham, J. (2019). Model-based sensitivity analysis for outcome reporting bias in the meta analysis of benefit and harm outcomes. Statistical Methods in Medical Research, 28, 889–903.


}
\examples{

data(EpilepsyData)

Apply the reORBgen function to a beneficial outcome.

reORBgen(
a=EpilepsyData$`50% Seizure Reduction (T)`,
c=EpilepsyData$`50% Seizure Reduction (C)`,
n1=EpilepsyData$`Sample (T)`,
n1=EpilepsyData$`Sample (C)`,
outcome="benefit",
alpha_ben=0.05,
alpha_ben_one.sided = FALSE,
init_param = c(0.1, 0.005),
true.SE=NULL,
LR.CI=TRUE,
Wald.CI = TRUE,
selection.benefit = "Neg.exp.piecewise",
weight_param=7,
opt_method="L-BFGS-B",
lower = c(-5, 0.00001),
upper = c(5,5),
low.risk = TRUE)

Apply the reORBgen function to a harmful outcome.

reORBgen(
a=EpilepsyData$`Nausea/vomiting (T)`,
c=EpilepsyData$`Nausea/vomiting (C)`,
n1=EpilepsyData$`Sample (T)`,
n1=EpilepsyData$`Sample (C)`,
outcome="harm",
alpha_harm=0.01,
alpha_ben_one.sided = FALSE,
init_param = c(0.1, 0.005),
true.SE=NULL,
LR.CI=TRUE,
Wald.CI = TRUE,
opt_method="L-BFGS-B",
lower = c(-5, 0.00001),
upper = c(5,5),
low.risk = FALSE)


}
\keyword{datasets}
