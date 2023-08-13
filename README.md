# ORBmeta


The package implements the methodology developed by Copas et al. (2019) and Copas et al. (2014) as well as extensions, to adjust for outcome reportin bias (ORB) in meta-analyses. 

# Installing ORBmeta

In the console follow these steps:

install.packages("devtools")

library(devtools)

devtools::install_github("agaiasaracini/ORBmeta")

library(ORBmeta)

# Functions
The package contains the function reORBgen(), which uses the random effects model to implement the methodology of Copas et al. (2019) and extensions based on generalized functions that describe the missing data mechanism. The function can be used of beneficial and harmful outcomes, and can take in: i) binary data, from which log RR are calculated as the treatment effect, ii) means, from which differences in means are computed, iii) a normally distributed treatment effect. The function reORBgen_hetfixed() is the same but assumes that the heterogeneity variance is a known and fixed parameter, and hence only optimizes with respect to the global treatment effect. The function feORBbin() implements the Copas et al. (2019) methodology for the fixed effect model, without extensions, and for binary outcome data only. The ORBsens() function performs a sensitivity analysis for ORB adjustment based on the probability that the Outcome Reporting Bias in Trial (ORBIT) assessment into high/low risk of bias is correct. This is done by implementing the method of Copas et al. (2014), which was constructed for beneficial outcomes, with its extension to harmful outcomes. The function works for both binary data and means. Finally, the function simulate_meta_data() simulates a meta-analysis in the presence of both low and high risk of bias study outcomes, 
in the random effects model context. This function can be used to generate a meta-analysis in the presence of ORB and subsequently apply the functions reORBgen(), reORBgen_hetfixed(), and ORBsens().

# Datasets
The package also contains some example datasets where ORB is present and the ORBIT methdology is used to classify the unreported data in high/low risk of bias. EpilepsyData contains multiple outcomes, both harmful and beneficial, and was used in Copas et al. (2019), 
PPHData and WeightData contain one beneficial outcome and were originally used in Copas et al. (2014).

# References

Copas, J., Marson, A., Williamson, P., and Kirkham, J. (2019). Model-based sensitivity analysis for outcome reporting bias in the meta analysis of benefit and harm outcomes. Statistical Methods in Medical Research, 28, 889–903. 

Copas, J., Dwan, K., Kirkham, J., and Williamson, P. (2014). A model-based correction for outcome reporting bias in meta-analysis. Biostatistics, 15, 370–383. 
