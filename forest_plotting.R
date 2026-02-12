####################      FOREST PLOT       ##########
#https://www.metafor-project.org/doku.php/plots:forest_plot
##Libraries
library(Matrix)
library(metadat)
library(numDeriv)
library(metafor)
library(ggplot2)
library(meta)
##input data

dat <- read.csv("test_data.csv", header = T, sep=",")
colnames(dat)

#no: sorting number
#author: the first author of the study
#year: publication year
#Ca_Co_total: the sum of total cases and control numbers
#Ca: total case number
#Co: total control number
#case_WT: wild-type homozygote(MM) in cases
#case_het: heterozygote (Mm) in cases
#case_hom: mutant-type homozygote (mm) in cases
#control_WT: wild-type homozygote(MM) in controls
#control_het: heterozygote (Mm) in controls
#control_hom: mutant-type homozygote (mm) in controls


#I use 3 genetic models: allelic, dominant and recessive:(https://pmc.ncbi.nlm.nih.gov/articles/PMC4390695/)
#a major allele (M) and a minor allele (m).Thus, the genotype can be a major allele homozygote (MM), a heterozygote (Mm) or a minor allele homozygote (mm)
#Allelic model: Mutated allele  vs Wild-type alleles 
#Dominant model: compares  Mm + mm versus MM (in cases: case_hom + case_het vs. case_WT) (in controls: control_hom + control_het vs. control_WT)
#Recessive model: compares mm versus MM + Mm (in cases: case_hom vs. case_WT + case_het) (in controls: control_hom vs. control_WT + control_het)

library(dplyr)

dat <- dat %>%
  # Create new columns for AM genetic model
  mutate(case_AM_WT = (case_WT * 2) + case_het,
         case_AM_mut = (case_hom * 2) + case_het,
         control_AM_WT = (control_WT * 2) + control_het,
         control_AM_mut = (control_hom * 2) + control_het) %>%
  # Create new columns for DM genetic model
  mutate(case_DM_WT = case_WT,
         case_DM_mut = case_het + case_hom,
         control_DM_WT = control_WT,
         control_DM_mut = control_het + control_hom) %>%
  # Create new columns for RM genetic model
  mutate(case_RM_WT = case_WT + case_het,
         case_RM_mut = case_hom,
         control_RM_WT = control_WT + control_het,
         control_RM_mut = control_hom)

# View the updated dataframe
head(dat)


###################################[[1]]forest flot for allelic model ######################################
#
#
dat1 <- escalc(measure="OR", ai=case_AM_mut, bi=case_AM_WT, ci=control_AM_mut, di=control_AM_WT, data=dat)
#res1 
res1 <- rma(yi, vi, data=dat1, slab=paste(author, year, sep=", ")) #random-effect OR
res1 <- rma(yi, vi, data=dat1, method="FE", slab=paste(author, year, sep=", "))#fixed-effect

### total number of studies
k <- nrow(dat1)

### generate point sizes
psize <- weights(res1)
psize <- 1.2 + (psize - min(psize)) / (max(psize) - min(psize))

### get the weights and format them as will be used in the forest plot
w_num <- weights(res1)
ord <- order(w_num, decreasing = TRUE)
dat1_ord <- dat1[ord, ]
res1_ord <- rma(
  yi, vi,
  data = dat1_ord,
  slab = paste(author, year, sep=", ")
)
dat_ord <- dat[ord, ]

#weights <- fmtx(weights(res1), digits=1)
sav <- forest(res1_ord, atransf=exp, showweights = TRUE, at=log(c(0.5, 1.0, 2.0)), xlim=c(-8,6),
              cex=1.1, mlab="", shade=TRUE, header = F, ylim = c(-4, res1$k + 3.5), ilab = cbind(
                paste0(dat_ord$case_AM_mut, "/", dat_ord$case_AM_WT),
                paste0(dat_ord$control_AM_mut, "/", dat_ord$control_AM_WT)
              ),)
par(xpd = NA)
op <- par(cex=1.1, font=2)

text(c(5.1), sav$ylim[2] - 1.5, c("Odds Ratio\n[95% CI]  "), cex=1)
#text(c(-9.5,-8,-5.9,-4.4), res$k+2.4, c("Events", "Total", "Events", "Total"), cex=1)
text(c(-3.4,-1.9), res1_ord$k+1.8, c("Case (m/M)", "Control (m/M)"), cex=1)

#text(c(-8.75,-5.25),     res$k+3.5, c("Case", "Control"), cex= 1)
text(c(3.9), sav$ylim[2] - 1.50, c("Weight\n(%)"), cex=1)## c(3) yatay lokasyon, k+1 olan dikey
text(c(-7.0), sav$ylim[2] - 1.70, "Author(s) and Year", cex=1)
par()
### add text with Q-value, dfs, p-value, and I^2 statistic
#"FE Model (Q = ",
text(-8, -0.9, pos=4, cex=1, bquote(paste(
  "RE Model (Q = ", .(fmtx(res1$QE, digits=2)),
  ", df = ", .(res1$k - res1$p), ", ",
  .(fmtp(res1$QEp, digits=3, pname="p", add0=TRUE, sep=TRUE, equal=TRUE)), "; ",
  I^2, " = ", .(fmtx(res1$I2, digits=1)), "%)")))

### add text for test of overall effect
text(-8, -1.8, pos=4, cex=1, bquote(paste("Test for overall effect: ",
                                           "Z=", .(fmtx(res1$zval, digits=2)), ", ",
                                           .(fmtp(res1$pval, digits=3, pname="P", add0=TRUE, equal=TRUE)))))




####################################[[2]]forest plot for dominant model#########################################
#
#
dat2 <- escalc(measure="OR", ai=case_DM_mut, bi=case_DM_WT, ci=control_DM_mut, di=control_DM_WT, data=dat)
#res2
res2 <- rma(yi, vi, data=dat2, slab=paste(author, year, sep=", ")) #random-effect OR
res2 <- rma(yi, vi, data=dat2, method="FE", slab=paste(author, year, sep=", "))#fixed-effect


### total number of studies
k <- nrow(dat2)

### generate point sizes
psize <- weights(res2)
psize <- 1.2 + (psize - min(psize)) / (max(psize) - min(psize))

### get the weights and format them as will be used in the forest plot
w_num2 <- weights(res2)
ord2   <- order(w_num2, decreasing = TRUE)

dat2_ord <- dat2[ord2, , drop = FALSE]

res2_ord <- rma(
  yi, vi,
  data = dat2_ord,
  slab = paste(author, year, sep=", ")
)

dat_ord2 <- dat[ord2, , drop = FALSE]

sav <- forest(
  res2_ord,
  atransf = exp,
  showweights = TRUE,
  at = log(c(0.5, 1.0, 2.0)),
  xlim = c(-8, 6),
  cex = 1.1,
  mlab = "",
  shade = TRUE,
  header = FALSE,
  ylim = c(-4, res2_ord$k + 3.5),
  ilab = cbind(
    paste0(dat_ord2$case_DM_mut, "/", dat_ord2$case_DM_WT),
    paste0(dat_ord2$control_DM_mut, "/", dat_ord2$control_DM_WT)
  )
)

par(xpd = NA)
op <- par(cex=1.1, font=2)

text(c(5.1), sav$ylim[2] - 1.5, c("Odds Ratio\n[95% CI]  "), cex=1)
#text(c(-9.5,-8,-5.9,-4.4), res$k+2.4, c("Events", "Total", "Events", "Total"), cex=1)
text(c(-3.0,-1.5), res2_ord$k+1.8, c("Case (m/M)", "Control (m/M)"), cex=1)

#text(c(-8.75,-5.25),     res$k+3.5, c("Case", "Control"), cex= 1)
text(c(3.9), sav$ylim[2] - 1.50, c("Weight\n(%)"), cex=1)## c(3) yatay lokasyon, k+1 olan dikey
text(c(-7.0), sav$ylim[2] - 1.70, "Author(s) and Year", cex=1)
par()
### add text with Q-value, dfs, p-value, and I^2 statistic
#"FE Model (Q = ",
text(-8, -0.9, pos=4, cex=1, bquote(paste(
  "FE Model (Q = ", .(fmtx(res2$QE, digits=2)),
  ", df = ", .(res2$k - res2$p), ", ",
  .(fmtp(res2$QEp, digits=3, pname="p", add0=TRUE, sep=TRUE, equal=TRUE)), "; ",
  I^2, " = ", .(fmtx(res2$I2, digits=1)), "%)")))

### add text for test of overall effect
text(-8, -1.8, pos=4, cex=1, bquote(paste("Test for overall effect: ",
                                          "Z=", .(fmtx(res2$zval, digits=2)), ", ",
                                          .(fmtp(res2$pval, digits=3, pname="P", add0=TRUE, equal=TRUE)))))



#################################[[3]]forest plot for recessive model#####################################
#
#
dat3 <- escalc(measure="OR", ai=case_RM_mut, bi=case_RM_WT, ci=control_RM_mut, di=control_RM_WT, data=dat)
#res3
res3 <- rma(yi, vi, data=dat3, slab=paste(author, year, sep=", ")) #random-effect OR
res3 <- rma(yi, vi, data=dat3, method="FE", slab=paste(author, year, sep=", "))#fixed-effect

### total number of studies
k <- nrow(dat3)

### generate point sizes
psize <- weights(res3)
psize <- 1.2 + (psize - min(psize)) / (max(psize) - min(psize))

### get the weights and format them as will be used in the forest plot
w_num3 <- weights(res3)
ord3   <- order(w_num3, decreasing = TRUE)

dat3_ord <- dat3[ord3, , drop = FALSE]

res3_ord <- rma(
  yi, vi,
  data = dat3_ord,
  slab = paste(author, year, sep=", ")
)

dat_ord3 <- dat[ord3, , drop = FALSE]

sav <- forest(
  res3_ord,
  atransf = exp,
  showweights = TRUE,
  at = log(c(0.5, 1.0, 2.0)),
  xlim = c(-8, 6),
  cex = 1.1,
  mlab = "",
  shade = TRUE,
  header = FALSE,
  ylim = c(-4, res3_ord$k + 3.5),
  ilab = cbind(
    paste0(dat_ord3$case_RM_mut, "/", dat_ord3$case_RM_WT),
    paste0(dat_ord3$control_RM_mut, "/", dat_ord3$control_RM_WT)
  )
)

par(xpd = NA)
op <- par(cex=1.1, font=2)

text(c(5.2), sav$ylim[2] - 1.2, c("Odds Ratio\n[95% CI]  "), cex=1)
#text(c(-9.5,-8,-5.9,-4.4), res$k+2.4, c("Events", "Total", "Events", "Total"), cex=1)
text(c(-5.5,-4.2), res3_ord$k+1.8, c("Case (m/M)", "Control (m/M)"), cex=1)

#text(c(-8.75,-5.25),     res$k+3.5, c("Case", "Control"), cex= 1)
text(c(3.8), sav$ylim[2] - 1.20, c("Weight\n(%)"), cex=1)## c(3) yatay lokasyon, k+1 olan dikey
text(c(-7.1), sav$ylim[2] - 1.70, "Author(s) and Year", cex=1)
par()
### add text with Q-value, dfs, p-value, and I^2 statistic
#"FE Model (Q = ",
text(-8, -0.9, pos=4, cex=1, bquote(paste(
  "FE Model (Q = ", .(fmtx(res3$QE, digits=2)),
  ", df = ", .(res3$k - res3$p), ", ",
  .(fmtp(res3$QEp, digits=3, pname="p", add0=TRUE, sep=TRUE, equal=TRUE)), "; ",
  I^2, " = ", .(fmtx(res3$I2, digits=1)), "%)")))

### add text for test of overall effect
text(-8, -2.3, pos=4, cex=1, bquote(paste("Test for overall effect: ",
                                          "Z=", .(fmtx(res3$zval, digits=2)), ", ",
                                          .(fmtp(res3$pval, digits=3, pname="P", add0=TRUE, equal=TRUE)))))
