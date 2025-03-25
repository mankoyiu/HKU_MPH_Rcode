#==============================================================================
# CMED 6211: Infectious Disease Epidemiology
# Solution for Homework 5 
# Mar 5, 2018
#==============================================================================
rm(list = ls());

library(MASS)
library(binom)

# Question 1.

# t1 = Day 30
# t2 = Day 60

numSubject_t1 <- 195
numPositive_t1 <- 30
numNegative_t1 <- numSubject_t1-numPositive_t1

numSubject_t2 <- 175
numPositive_t2 <- 86
numNegative_t2 <- numSubject_t2-numPositive_t2

serosurvey <- data.frame(subject.no=1:numSubject_t1, positive_t1=c(rep(1,numPositive_t1),rep(0,numNegative_t1)),
                         positive_t2=c(rep(1,numPositive_t2),rep(0,numNegative_t2), rep(NA,numSubject_t1-numSubject_t2)))

# 1.1 Use the binomial distribution to provide MLE estimate and the CIs

# (1) Compute the negative log likelihood of a data vector y given a guess at the parameter theta
NLL <- function(theta,y) {
  NLL <- 0
  N <- length(y)
  for (i in 1:N) {
   if (y[i]==1) {p <- theta} # seroposistive
    if (y[i]==0) {p <- 1-theta} # seronegative
    NLL <- NLL + log(p)
    }
 -NLL
  }

# (2) Use a built-in optimizer to find the parameter value theta that minimizes the negative log likelihood
out <- optim(par=0.5, fn=NLL, method="Brent", lower=0.00001, upper=0.99999, y=serosurvey$positive_t1, hessian = TRUE)
p1 <- theta <- out$par # p1 = 0.15
out <- optim(par=0.5, fn=NLL, method="Brent", lower=0.00001, upper=0.99999, y=serosurvey$positive_t2[1:numSubject_t2], hessian = TRUE)
p2 <- theta <- out$par # p2 = 0.49

# (3) 95% confidence intervals of the binomial distribution - using all 8 different methods
p1_conf <- binom.confint(30, 195, conf.level = 0.95, methods = "all") # Day 30
p2_conf <- binom.confint(86, 175, conf.level = 0.95, methods = "all") # Day 60

#########################################################################################
# 1.1 Method 2 - Providing MLE

#a) - Day 30
x=30
n=195

#Building the likelihood function
likelhd <- function(p) dbinom(x,size=n,prob=p)
plot(likelhd,0,1,xlab="p",ylab="l(p)",main="Binomial likelihood - Day 30")

#maximum likelihood estimate
opt<-optimize(likelhd,c(0,1),maximum=TRUE)
opt$maximum

#b) The same for Day 60

#==============================================================================
# 1.2 Estimate cumulative infection attack rate on days 30 and 60

# 90% of infected individuals become seropositive

cIAR_t1 <- p1/0.9 # 0.17 
cIAR_t2 <- p2/0.9 # 0.55 

act_numPositive_t1 <- numPositive_t1/0.9
act_numPositive_t2 <- numPositive_t2/0.9

# 95% confidence intervals of the binomial distribution - using all 8 different methods
cIAR_t1_conf <- binom.confint(act_numPositive_t1, 195, conf.level = 0.95, methods = "all") # Day 30
cIAR_t2_conf <- binom.confint(act_numPositive_t2, 175, conf.level = 0.95, methods = "all") # Day 60

#==============================================================================
# 1.3 Estimate the infection attack rate between days 30 and 60
cIAR_t2 - cIAR_t1 # 0.38

#===============================================================================

# Question 2.

# 1. Data from the HK Center of Health Protection
weeklyILI <- read.csv("Data for Q2.csv",header=TRUE)

# 2. Observed ILI cases per 1,000 consultations starting from week 5
numWeek <- length(weeklyILI$Week)
#observed <- weeklyILI[1:numWeek,]

# 3. Expected ILI cases per 1,000 consultations starting from week 5

# Method 3.1: Use 'SMA' (simple moving average) built-in function
weeklyILI$expected <- c(NA, SMA(weeklyILI$ILI[1:(numWeek-1)], n=4))

# Method 3.2: Create a simple for-loop function to calculate the mean ILI over the past 4 weeks
weeklyILI$expected <- NA

for (i in 5:numWeek)
{
  weeklyILI$expected[i] <- mean(weeklyILI$ILI[(i-4):(i-1)])
}


# 4. 95% CI of the observed incidence
weeklyILI$LL <- qpois(0.025,weeklyILI$ILI) #Lower Limit
weeklyILI$UL <- qpois(0.975,weeklyILI$ILI) #Upper Limit

# 5. Alert generated when expected incidence < lower limit of the observed incidence
alert <- weeklyILI$Week[weeklyILI$expected<weeklyILI$LL]
alert # Alert generated in Week 9, 10, and 40, and alarm issued at week 10.

# Plot 
plot(weeklyILI$Week,weeklyILI$expected,type='l',col='blue',lwd=2,
     xlim=c(0,numWeek),ylim=c(10,75),
     xlab='Week',ylab='ILI cases per 1,000 consultations',
     main='Weekly syndromic surveillance') #Expected ILI
lines(weeklyILI$ILI,lwd=2,col='red') #Observed ILI
lines(weeklyILI$LL,lwd=2,lty=2) #95% CI lower limit
lines(weeklyILI$UL,lwd=2,lty=2) #95% CI upper limit

legend("topright",c('Expected incidence','Observed incidence','95% CI of observed incidence'),
       lty=c(1,1,2),lwd=2,col=c('blue','red','black'),bty='n',cex=0.8)
	   
arrows(9,40,9,36,length=0.2,angle=20,lwd=2,col='orange') #An arrow pointing at week 9
arrows(10,39,10,35,length=0.2,angle=20,lwd=2,col='orange') #An arrow pointing at week 10
arrows(40,41,40,37,length=0.2,angle=20,lwd=2,col='orange') #An arrow pointing at week 40

colnames(weeklyILI) <- c("Week", "Observed", "Expected", "95% Lower Limit", "95% Upper Limit")
write.csv(weeklyILI, file="Weekly ILI.csv", row.names=F)

###
#END of the script

