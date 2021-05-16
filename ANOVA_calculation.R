## R code template for analyses reported in Li, Koester, and Lachance et al iScience 2021
## DOI:https://doi.org/10.1016/j.isci.2021.102508


# load packages -----------------------------------------------------------
library(ez)


# calculate ANOVA statistics ----------------------------------------------

# two different methods for calculating an ANOVA presented below #

# option 1: calculating ANOVA using the basic R function aov() #
#   A = Dataset e.g. "Supplemental Table 3"
#   B = Response variable e.g. "whole gut transit time"
#   X, Y, Z, etc. = Independent variables e.g. "Microbiota", "Diet", Genotype", etc.
res.aov = aov(B ~ X + Y + Z, data = A)
summary(res.aov)

# option 2: calculating ANOVA using the function ezANOVA from the ez package #
#   A = Dataset e.g. "Supplemental Table 3"
#   B = Response variable e.g. "whole gut transit time"
#   C = Variable describing individual subject e.g. MouseID
#   X, Y, Z, etc. = Independent variables e.g. "Microbiota", "Diet", Genotype", etc.
ezANOVA(data = A, dv = B, wid = C, between = c(X, Y, Z))$ANOVA

