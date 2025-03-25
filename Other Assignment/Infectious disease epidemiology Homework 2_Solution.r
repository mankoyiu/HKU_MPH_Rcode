#Fraser et al NEJM 1977
# Figure 1
onsetDays = 3:29
onsetCounts = c(1,5,10,24,21,22,20,14,9,7,7,4,3);
incuPeriods = c(2,rep(3,5),rep(4,10),rep(5,24),rep(6,21),rep(7,22),
                rep(8,20),rep(9,14),rep(10,9),rep(11,7),rep(12,7),rep(13,4),rep(14,3));
# 1) "fitdistr"
library(MASS);
fitGamma = fitdistr(incuPeriods,"gamma");
fitGamma
fitMean = as.numeric(fitGamma$estimate[1]/fitGamma$estimate[2]);
fitMean
fitVar = as.numeric(fitGamma$estimate[1]/fitGamma$estimate[2]/fitGamma$estimate[2]);
fitVar

# 95th percentile
qgamma(.95, fitGamma$estimate[1], fitGamma$estimate[2])

# Histogram
hist(incuPeriods, prob=TRUE, breaks=0:15+0.5, main="Incubation period distribution", xlab="Incubation period",
ylab="Probability")

# Fitted curve
curve(dgamma(x,fitGamma$estimate[1],fitGamma$estimate[2]),add=TRUE,col="blue")

# 2) "optim"
funNegLogLikelihood = function(x)
{
  vshape = x[1];
  vscale = x[2];
  logLikelihood = 0;
  for(ii in 1:length(incuPeriods))
  {
    logLikelihood = logLikelihood + log(dgamma(incuPeriods[ii],shape=vshape,scale=vscale));
  }
  return(-logLikelihood);
}
optimGamma = optim(c(9,1),funNegLogLikelihood);
optimGamma
optimMean = as.numeric(optimGamma$par[1]*optimGamma$par[2]);
optimMean
optimVar = as.numeric(optimGamma$par[1]*optimGamma$par[2]*optimGamma$par[2]);
optimVar
