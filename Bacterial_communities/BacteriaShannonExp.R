# Shannon exp bacteria
# Load library
library(psych)
library(dunn.test)

# Load the data
BacteriaExpH<-read.table("ShannonExpBacteria.csv",h=T, sep=",")
# Visualize
BacteriaExpH
# Set as factor the landcover
BacteriaExpH$ProjectFocus<-as.factor(BacteriaExpH$ProjectFocus)
# Verify
summary(BacteriaExpH)
# Calculate average
tapply(BacteriaExpH$Shannon,BacteriaExpH$ProjectFocus, mean)
# Calculate STD
tapply(BacteriaExpH$Shannon, BacteriaExpH$ProjectFocus, sd)
# Or smply
describeBy(BacteriaExpH$Shannon, group=BacteriaExpH$ProjectFocus)

# KW test and Dunn's test
dunn.test(BacteriaExpH$Shannon, BacteriaExpH$ProjectFocus, method="bonferroni")


