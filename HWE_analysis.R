######################    HWE  for control    #######
library(HardyWeinberg)
x <- c(MM=325, MN=233, NN=35)
#HW.test <- HWChisq(x,verbose=TRUE)
# Same test without continuity correction
HW.test <- HWChisq(x,cc=0,verbose=TRUE)
HW.exacttest <- HWExact(x, verbose = TRUE)# The exact test was preferred when the expected count under the Hardy-Weinberg law was less than five for any of the three genotypes
