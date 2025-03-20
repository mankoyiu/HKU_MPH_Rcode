
# clear old files
rm(list=ls(all=TRUE))

# install statistical packages for meta-analysis
install.packages("meta", dependencies=TRUE)
library(meta)


# read summary estimates
getwd()
setwd("C:/")
colon <- read.csv("colon.csv", header=TRUE)
colon2 <- read.csv("colon2.csv", header=TRUE)
colon2revised <- read.csv("colon2revised.csv", header=TRUE)

colon <- colon[order(-colon$year,colon$studlab),] 
head(colon)
colon2 <- colon2[order(-colon$year,colon$studlab),] 
head(colon2)
colon2revised <- colon2revised[order(-colon$year,colon$studlab),] 
head(colon2revised)


# run meta-analysis for binary outcomes
# Fixed effect OR (inverse variance weighting) for all studies included, OR=1.54 using metabin

allsound<-NA
allsound <- metacont(vn, vscore, vsd, cn, cscore, csd, studlab,
                 data=colon, median.e=vmedian, q1.e=vq1, q3.e=vq3, min.e=vmin, max.e=vmax, 
                 median.c=cmedian, q1.c=cq1, q3.c=cq3, min.c=cmin, max.c=cmax, 
                 method.mean = "Luo", method.sd = "Shi",  label.e="With Visual Distraction", label.c="Usual")
allsound

# draw funnel plot
windows(height=10, width=10)
funnel(allsound)

# draw forest plot
windows(height=10, width=10)
forest(allsound,overall=T)
forest(allsound,comb.r=FALSE)
forest(allsound,comb.f=FALSE)

# comb.r - indicate whether random effect should be plotted; if comb.r=FALSE - only fixed effect is plotted
# comb.f - indicate whether fixed effect should be plotted; if comb.f=FALSE - only random effect is plotted
# pooled.events  indicating whether total number of events should be given in the figure.
 

# stratified analysis (by funding)
# Fixed effect OR (inverse variance weighting) for industry funded ones (pharm==1) of all studies included, OR=0.89 using metabin

allnosound<-NA
allnosound <- metacont(vn, vscore, vsd, cn, cscore, csd, studlab,
                     data=colon2, median.e=vmedian, q1.e=vq1, q3.e=vq3, min.e=vmin, max.e=vmax, 
                     median.c=cmedian, q1.c=cq1, q3.c=cq3, min.c=cmin, max.c=cmax, 
                     method.mean = "Luo", method.sd = "Shi", sm="MD", label.e="With Visual Distraction", label.c="Usual")
allnosound

windows(height=10, width=10)
funnel(allnosound)

# draw forest plot
windows(height=10, width=10)
forest(allnosound,overall=T)
forest(allnosound,comb.r=FALSE)
forest(allnosound,comb.f=FALSE)

allnosound2<-NA
allnosound2 <- metacont(vn, vscore, vsd, cn, cscore, csd, studlab,
                       data=colon2revised, median.e=vmedian, q1.e=vq1, q3.e=vq3, min.e=vmin, max.e=vmax, 
                       median.c=cmedian, q1.c=cq1, q3.c=cq3, min.c=cmin, max.c=cmax, 
                       method.mean = "Luo", method.sd = "Shi", sm="MD", label.e="With Visual Distraction", label.c="Usual")
allnosound2

windows(height=10, width=10)
funnel(allnosound2)

# draw forest plot
windows(height=10, width=10)
forest(allnosound2,overall=T)
forest(allnosound2,comb.r=FALSE)
forest(allnosound2,comb.f=FALSE)
