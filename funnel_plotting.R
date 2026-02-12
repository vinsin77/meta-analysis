#####FUNEL PLOT##########
#https://www.metafor-project.org/doku.php/plots:forest_plot
##Libraries
library(Matrix)
library(metadat)
library(numDeriv)
library(metafor)
library(ggplot2)
##input data

dat <- read.csv("test_data.csv", header = T, sep=",")
#dat
dat1 <- escalc(measure="OR", ai=case_MT, bi=case_WT, ci=control_MT, di=control_WT, data=dat)
#dat1
res1 <- rma(yi, vi, data=dat1, slab=paste(author, year, sep=", ")) #random-effect
res1 <- rma(yi, vi, data=dat1, method="FE", slab=paste(author, year, sep=", "))#fixed-effect
#res1 
### draw a standard funnel plot
funnel(res1)
### show risk ratio values on x-axis (log scale)
funnel(res1, refline=0, ylim=c(0, 1), pch=1, cex=2)
#funnel(res1, refline=0, ylim=c(0, 0.5), pch=1, cex=1.5, label="out")#pch=1 makes empty circle shape
title("Funnel Plot with Pseudo 95% Confidence Limits")
#layout(matrix(c(1, 2, 0, 3), nrow = 2, byrow = TRUE))
###EGGER'S TEST####
#regtest(res1) ## use of a standard mixed-effects meta-regression model is described for the purposes of the regression test for funnel plot asymmetry. 
regtest(res1, model="lm") ## classical test
regtest(res1, model="rma", predictor="sei") ##mixed effects

##Begg's test:
regtest(res1, model="rma")
