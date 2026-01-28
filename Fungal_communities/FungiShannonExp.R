# Shannon exp funghi
# Load library
library(psych)
library(dunn.test)

# Load the data
FunghiExpH<-read.table("ShannonExpFungi.csv",h=T, sep=",")
# Visualize
FunghiExpH
# Set as factor the landcover
FunghiExpH$SampleGroup<-as.factor(FunghiExpH$SampleGroup)
# Verify
summary(FunghiExpH)
# Calculate average
tapply(FunghiExpH$ShannonExp,FunghiExpH$SampleGroup, mean)
# Calculate STD
tapply(FunghiExpH$ShannonExp, FunghiExpH$SampleGroup, sd)
# Or smply
describeBy(FunghiExpH$ShannonExp, group=FunghiExpH$SampleGroup)

# KW test and Dunn's test
dunn.test(FunghiExpH$ShannonExp, FunghiExpH$SampleGroup, method="bonferroni")


