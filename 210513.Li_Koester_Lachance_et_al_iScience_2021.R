## R code template for analyses reported in Li, Koester, and Lachance et al iScience 2021
## DOI:https://doi.org/10.1016/j.isci.2021.102508


# load packages -----------------------------------------------------------
library(ggplot2)
library(pheatmap)
library(reshape2)
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



# create Loess curves from FITC-dextran motility assay data ---------------

# Do not read in table S3 directly from Excel. Instead use the .txt file posted on GitHub. This file contains the same information, but is in a format that is R readable.
# Make sure the .txt file is saved in the same folder as this code; otherwise, you may get an error when importing the dataset. Alternatively, you could set working directory with setwd().

# read in Supplemental Table 3 data #
dat.FITC = read.table(file = "Table_S3_For_GitHub_v210415.txt", header = T, sep = "\t")

# edit the code that follows this comment block in order to define the specific cohorts you would like to analyze #
#   microbiota options = "GF", "BSH-low-2", "BSH-low-1", "BSH-high-1", "BSH-high-2", "SPF", "BSH-low-1 + BSH-high-1"
#   genotype options = "Swiss-Webster wild-type", "C57BL/6 wild-type", "C57BL/6 Ret+/-"
#   diet options = "Turmeric", "Bland", "Bland-->Turmeric", "Turmeric-->Bland"
#   FITC_time options = 2, 4
v.Microbiota = "GF"
v.Genotype = "Swiss-Webster wild-type"
v.Diet = "Turmeric"
v.FITC_Time = 2

# subset FITC data for above microbiota, genotype, and diet variables defined by the user #
dat.FITC.sub = subset(dat.FITC, Microbiota == v.Microbiota & Genotype == v.Genotype & Diet == v.Diet & Interval.between.FITC.gavage.and.euthanasia == v.FITC_Time)

# calculate mean FITC values for each gut segment, which we will index 1-12 #
# FITC data in columns 10 through 21 in dat.FITC.sub, which is why we calculate colMeans() on c(10:21)
dat.FITC.sub.colMeans = colMeans(dat.FITC.sub[,c(10:21)])
dat.FITC.sub.colMeans = as.data.frame(cbind("Index" = c(1:12), "FITC_Mean" = dat.FITC.sub.colMeans))

# calculate smoothed Loess curves using user defined span variable #
# Our Loess curves were generated with a span of 0.5 which means half of the data points were used in the calculation. A lower span will make the curve more jagged and a larger span will make a smoother curve.
v.span = 0.50
loess_FITC = loess(FITC_Mean ~ Index, data = dat.FITC.sub.colMeans, span = v.span)
smoothed_FITC = predict(loess_FITC, se=T)
dat.graph_loess = cbind("Index" = c(1:12),
                        "ID" = "Loess",
                        melt(smoothed_FITC$fit, value.name = "Value"),
                        melt(smoothed_FITC$se.fit, value.name = "SEM"))

# plot Loess curve calculated above with individual mice #
dat.FITC.sub_transpose = t(dat.FITC.sub[,c(10:21)])  # transpose FITC columns in dat.FITC.sub so we can plot the individual mice
colnames(dat.FITC.sub_transpose) = dat.FITC.sub$Mouse.ID  # define column names as mouse.ID
rownames(dat.FITC.sub_transpose) = c(1:12)  # define row names as c(1:12) which corresponds to segment index

# use melt() to create flat file of individual mice; columns correspond to Mouse.ID and FITC value #
dat.mouse.FITC = melt(dat.FITC.sub_transpose, varnames = c("Index", "ID"), value.name = "Value")

# create column for standard errors in order to merge with Loess data for the cohort
dat.mouse.FITC$SEM = NA

# combine Loess FITC dataframe with mouse FITC dataframe
dat.graph_FITC = rbind(dat.graph_loess, dat.mouse.FITC)

# generate plot of the Loess curve (with error bars) for the entire cohort along with individual mouse readouts #
#   (subset dataframe given to the error bar function for just the loess data)
ggplot(dat.graph_FITC) +
  geom_line(aes(x = Index, y = Value, color = ID), size = 1) +
  geom_errorbar(data = subset(dat.graph_FITC, ID == "Loess"), aes(x = Index, ymin=Value-SEM, ymax=Value+SEM, color = ID), width=.5, position=position_dodge(0.9), size = 1) + 
  xlab("Segment") + ylab("% FITC") + scale_x_continuous(breaks = dat.graph_FITC$Index) + ggtitle("Proportion FITC Measured in Intestinal Segments")  +
  theme(panel.background = element_rect(fill = "white", colour = "black"), legend.justification=c(1,0.01), legend.title = element_blank(), plot.title = element_text(hjust = 0.5), legend.key=element_blank())
# use ggsave to save the file #



# calculate geometric means -----------------------------------------------

geo_mean <- function(v.FITC.indiv)  {
  return (sum(v.FITC.indiv * c(1:12)))
}

# examples using the data structure dat.FITC.sub_transpose created above; pay attention to the orientation of your table #
v.FITC.indiv = dat.FITC.sub_transpose[,1] 
geo_mean(v.FITC.indiv = v.FITC.indiv)

v_geoMeans = rep(NA, ncol(dat.FITC.sub_transpose))
for (i_col in 1:ncol(dat.FITC.sub_transpose))  {
  v_geoMeans[i_col] = geo_mean(v.FITC.indiv = dat.FITC.sub_transpose[,i_col])
}

