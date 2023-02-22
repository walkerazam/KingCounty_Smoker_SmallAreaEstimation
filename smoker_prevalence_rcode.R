## R-Code Appendix

# Imports
library(SUMMER)
if (!isTRUE(requireNamespace("INLA", quietly = TRUE))) {
  install.packages("INLA",
                   repos = c(getOption("repos"),
                             INLA="https://inla.r-inla-download.org/R/stable"),
                   dep=TRUE)
}
library(sf) # Load sf for spatial analysis
library(prioritizr) # Allows us to create an adjacency matrix
library(survey)
library(ggplot2)


## Loading data and pre-processing
data(BRFSS)
# Dropping missing smoker1 rows or missing HRA codes
BRFSS <- subset(BRFSS, !is.na(BRFSS$smoker1))
BRFSS <- subset(BRFSS, !is.na(BRFSS$hracode))
data(KingCounty)
# cast as spatial dataframe
KingCounty <- st_as_sf(KingCounty)
# compute adjacency matrix
mat <- adjacency_matrix(KingCounty)
# Setting row and col name to HRA names in our neighbor
# matrix
colnames(mat) <- rownames(mat) <- KingCounty$HRA2010v2_
mat <- as.matrix(mat[1:dim(mat)[1], 1:dim(mat)[1]])

## Q1: Naive Direct Estimates

smoothed <- smoothSurvey(
  data = BRFSS, geo = KingCounty, Amat = NULL, responseType = "binary",
  responseVar = "smoker1", strataVar = NULL, weightVar = NULL, 
  regionVar = "hracode",
  clusterVar = NULL, CI = 0.95
)
# assigning a dataframe for mapping purposes:
naive_df <- smoothed$HT
# find standard error by taking the square root of the variance
naive_df$naive_se <- sqrt(smoothed$HT$HT.var)
# mapping mean posterior estimates for the Naive Direct Estimates
mapPlot(
  data = naive_df, geo = KingCounty,
  variables = c("HT.est"),
  labels = c("Naive Direct Estimates"), by.data = "region", by.geo = "HRA2010v2_", legend.label="Estimate Value"
)
# mapping standard errors
mapPlot(
  data = naive_df, geo = KingCounty,
  variables = c("naive_se"),
  labels = c("Naive Standard Errors"), by.data = "region", by.geo = "HRA2010v2_", legend.label="Standard Error"
)

## Question 2: Weighted Estimation

# Design object using strata and weights from BRFSS (rwt_llcp is final weights)
design <- svydesign(
  ids = ~1, weights = ~rwt_llcp,
  strata = ~strata, data = BRFSS
)
# getting direct estimates
direct <- svyby(~smoker1, ~hracode, design, svymean)
data(KingCounty)
# mapping weighted estimates
mapPlot(
  data = direct, geo = KingCounty,
  variables = c("smoker1"),
  labels = c("Weighted Estimates"), by.data = "hracode", by.geo = "HRA2010v2_", legend.label = "Weighted Estimates"
)
# mapping the standard errors
mapPlot(
  data = direct, geo = KingCounty,
  variables = c("se"),
  labels = c("Weighted Standard Errors"), by.data = "hracode", by.geo = "HRA2010v2_",  legend.label = "Standard Errors"
)

## Question 3: Comparing Naive and Weighted Estimates
# merging the naive and weighted dataframes for plotting
merged <-  merge(naive_df, direct, by.x=c("region"), by.y=c('hracode'))
# plotting our results
ggplot(merged, aes(x = HT.est, y = smoker1)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  ggtitle("Naive vs. Weighted Estimates") +
  xlab("Naive Estimates") +
  ylab("Weighted Estimates") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlim(0, 0.3) + ylim(0, 0.3)
ggplot(merged, aes(x = naive_se, y = se)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  ggtitle("Naive vs. Weighted Standard Errors") +
  xlab("Naive SE") +
  ylab("Weighted SE") +
  theme(plot.title = element_text(hjust = 0.5)) + xlim(0, 0.1) + ylim(0, 0.1)

## Question 4: Smoothed Naive Estimation

# Specifying Amat = mat to allow for BYM2 random effects
smoothed_bym2 <- smoothSurvey(
  data = BRFSS, geo = KingCounty, Amat = mat, responseType = "binary", 
  responseVar = "smoker1", strataVar = NULL, weightVar = NULL,
  regionVar = "hracode", clusterVar = NULL,
  CI = 0.95
)
# extracting posterior medians from smoothed_bym2$smooth$median
mapPlot(
  data = smoothed_bym2$smooth, geo = KingCounty,
  variables = c("median"),
  labels = c("Smoothed Posterior Medians"), by.data = "region", by.geo = "HRA2010v2_", legend.label = 'Posterior Medians'
)
# taking the square root of the variance
smoothed_bym2$smooth$std <- sqrt(smoothed_bym2$smooth$var)
# plotting
mapPlot(
  data = smoothed_bym2$smooth, geo = KingCounty,
  variables = c("std"),
  labels = c("Smoothed Posterior Standard Deviations"), by.data = "region", by.geo = "HRA2010v2_", legend.label = 'Standard Deviation'
)


## Question 5: Comparing Naive and Smoothed Binomial Estimates

# merging naive and smoothed results
merged2 <-  merge(naive_df, smoothed_bym2$smooth, by = 'region')
# plotting...
ggplot(merged2, aes(x = HT.est, y = median)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  ggtitle("Naive vs. Smoothed Estimates") +
  xlab("Naive Estimates") +
  ylab("Smoothed Estimates") +
  theme(plot.title = element_text(hjust = 0.5)) + xlim(0, 0.25) + ylim(0, 0.25)
# using square root of variance
merged2$smoothed_se <- sqrt(merged2$var)
# plotting
ggplot(merged2, aes(x = naive_se, y = smoothed_se)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  ggtitle("Naive vs. Smoothed Standard Errors") +
  xlab("Naive SE") +
  ylab("Smoothed SE") +
  theme(plot.title = element_text(hjust = 0.5)) + xlim(0, 0.04) + ylim(0, 0.04)

## Question 6: Smoothed Weighted Estimations

# Specifying Amat = mat to allow for BYM2 random effects
FHmodel <- smoothSurvey(
  data = BRFSS, geo = KingCounty, Amat = mat, responseType = "binary", 
  responseVar = "smoker1", strataVar = "strata", weightVar = "rwt_llcp",
  regionVar = "hracode", clusterVar = "~1", CI = 0.95
)
# extracting posterior medians: FHmodel$smooth$median
mapPlot(
  data = FHmodel$smooth, geo = KingCounty,
  variables = c("median"),
  labels = c("Smoothed Weighted Posterior Medians"), by.data = "region", by.geo = "HRA2010v2_", legend.label = 'Posterior Medians'
)
# getting standard deviation by taking the sqrt of variance
FHmodel$smooth$std <- sqrt(FHmodel$smooth$var)
# plotting...
mapPlot(
  data = FHmodel$smooth, geo = KingCounty,
  variables = c("std"),
  labels = c("Smoothed Weighted Standard Deviations"), by.data = "region", by.geo = "HRA2010v2_", legend.label = 'Standard Deviation'
)

## Question 7: Comparing Weighted and Smoothed Weighted Estimates

merged3 <-  merge(direct, FHmodel$smooth, by.x = c('hracode'), by.y = c('region') )
ggplot(merged3, aes(x = smoker1, y = median)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  ggtitle("Weighted vs. Smoothed Weighted Estimates") +
  xlab("Weighted Estimates") +
  ylab("Smoothed Weighted Estimates") +
  theme(plot.title = element_text(hjust = 0.5)) + xlim(0, 0.35) + ylim(0, 0.35)

merged3$std <- sqrt(FHmodel$smooth$var)
ggplot(merged3, aes(x = se, y = std)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  ggtitle("Weighted vs. Smoothed Weighted Standard Errors") +
  xlab("Weighted SE") +
  ylab("Smoothed Weighted SE") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlim(0, 0.09) + ylim(0, 0.09)

## Question 8: Summarize the HRA variation in smoking prevalence across King County

library(gridExtra)
## plot of all maps
p1 <-  mapPlot(
  data = naive_df, geo = KingCounty,
  variables = c("HT.est"),
  labels = c("Naive Direct Estimates"), by.data = "region", by.geo = "HRA2010v2_", legend.label="Estimate Value"
)
p2 <- mapPlot(
  data = direct, geo = KingCounty,
  variables = c("smoker1"),
  labels = c("Weighted Estimates"), by.data = "hracode", by.geo = "HRA2010v2_", legend.label = "Weighted Estimates"
)
p3 <- mapPlot(
  data = smoothed_bym2$smooth, geo = KingCounty,
  variables = c("median"),
  labels = c("Smoothed Posterior Medians"), by.data = "region", by.geo = "HRA2010v2_", legend.label = 'Posterior Medians'
)
p4 <- mapPlot(
  data = FHmodel$smooth, geo = KingCounty,
  variables = c("median"),
  labels = c("Smoothed Weighted Posterior Medians"), by.data = "region", by.geo = "HRA2010v2_", legend.label = 'Posterior Medians'
)
grid.arrange(grobs = list(p1, p2, p3, p4), ncol = 2)

# plot of all errors
e1 <-  mapPlot(
  data = naive_df, geo = KingCounty,
  variables = c("naive_se"),
  labels = c("Naive Standard Errors"), by.data = "region", by.geo = "HRA2010v2_", legend.label="Error Value"
)
e2 <- mapPlot(
  data = direct, geo = KingCounty,
  variables = c("se"),
  labels = c("Weighted Standard Errors"), by.data = "hracode", by.geo = "HRA2010v2_", legend.label = "Error Value"
)
e3 <- mapPlot(
  data = smoothed_bym2$smooth, geo = KingCounty,
  variables = c("std"),
  labels = c("Smoothed Posterior Standard Deviation"), by.data = "region", by.geo = "HRA2010v2_", legend.label = 'Error Value'
)
e4 <- mapPlot(
  data = FHmodel$smooth, geo = KingCounty,
  variables = c("std"),
  labels = c("Smoothed Weighted Posterior Standard Deviations"), by.data = "region", by.geo = "HRA2010v2_", legend.label = 'Error Value'
)
grid.arrange(grobs = list(e1, e2, e3, e4), ncol = 2)
