## R code template for analyses reported in Li, Koester, and Lachance et al iScience 2021
## DOI:https://doi.org/10.1016/j.isci.2021.102508


# calculate geometric means -----------------------------------------------

geo_mean <- function(v.FITC.indiv)  {
  return (sum(v.FITC.indiv * c(1:12)))
}

# examples using the data structure dat.FITC.sub_transpose created using the Loess curve generation script; pay attention to the orientation of your table #
v.FITC.indiv = dat.FITC.sub_transpose[,1] 
geo_mean(v.FITC.indiv = v.FITC.indiv)

v_geoMeans = rep(NA, ncol(dat.FITC.sub_transpose))
for (i_col in 1:ncol(dat.FITC.sub_transpose))  {
  v_geoMeans[i_col] = geo_mean(v.FITC.indiv = dat.FITC.sub_transpose[,i_col])
}

