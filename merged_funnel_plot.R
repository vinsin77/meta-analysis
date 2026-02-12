library(Matrix)
library(metadat)
library(numDeriv)
library(metafor)
library(ggplot2)
##input data
library(dplyr)
dat <- read.csv("PLA2G2A_rs11573156.csv", header = T, sep=",")


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

dat1 <- escalc(measure="OR", ai=case_AM_mut, bi=case_AM_WT, ci=control_AM_mut, di=control_AM_WT, data=dat)
dat2 <- escalc(measure="OR", ai=case_DM_mut, bi=case_DM_WT, ci=control_DM_mut, di=control_DM_WT, data=dat)
dat3 <- escalc(measure="OR", ai=case_RM_mut, bi=case_RM_WT, ci=control_RM_mut, di=control_RM_WT, data=dat)
#dat1
res1 <- rma(yi, vi, data=dat1, slab=paste(author, year, sep=", ")) #random-effect OR
res1 <- rma(yi, vi, data=dat1, method="FE", slab=paste(author, year, sep=", "))#fixed-effect
#res1 
res2 <- rma(yi, vi, data=dat2, slab=paste(author, year, sep=", ")) #random-effect OR
res2 <- rma(yi, vi, data=dat2, method="FE", slab=paste(author, year, sep=", "))#fixed-effect
#res2
res3 <- rma(yi, vi, data=dat3, slab=paste(author, year, sep=", ")) #random-effect OR
res3 <- rma(yi, vi, data=dat3, method="FE", slab=paste(author, year, sep=", "))#fixed-effect
#res3


library(grid)
library(gridGraphics)
library(gridExtra)
library(patchwork)
library(ggplotify)
#res1 <- rma(yi, vi, data=dat1, slab=paste(author, year, sep=", ")) 
funnel(res1)
funnel(res2)
funnel(res3)
# Create consistent aspect ratio (1:1 for equal width/height)
aspect_ratio <- 1

# Generate funnel plots and convert to ggplot with fixed aspect ratio
#funnel_plot_A:Allelic comparison; funnel_plot_B: dominant model; funnel_plot_C:recessive model
funnel_plot_A <- as.ggplot(~funnel(res1, refline=0, ylim=c(0, 2.1), pch=1, cex=2)) + 
  ggtitle("A") + theme(plot.title = element_text(face = "bold")) + coord_fixed(ratio = aspect_ratio)
funnel_plot_B <- as.ggplot(~funnel(res2, refline=0, ylim=c(0, 1.6), pch=1, cex=2)) + 
  ggtitle("B") + theme(plot.title = element_text(face = "bold")) + coord_fixed(ratio = aspect_ratio)
funnel_plot_C <- as.ggplot(~funnel(res3, refline=0, ylim=c(0, 0.85), xlim=c(-2.2, 2.2), pch=1, cex=2)) + 
  ggtitle("C") + theme(plot.title = element_text(face = "bold")) + coord_fixed(ratio = aspect_ratio)

# Combine the three funnel plots using patchwork, centering plot C
combined_plot <- (funnel_plot_A | funnel_plot_B) / (plot_spacer() | funnel_plot_C | plot_spacer())

# Apply custom layout with equal height distribution
combined_plot <- combined_plot + 
  plot_layout(heights = c(1, 1))  # Give equal height to both rows

# Display the combined plot
combined_plot


regtest(res1, model="lm") ## classical test
regtest(res2, model="lm")
regtest(res3, model="lm")

regtest(res1, model="rma")##Begg's test
regtest(res2, model="rma")
regtest(res3, model="rma")
