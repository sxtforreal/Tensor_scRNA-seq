# Packages
library(tensr)
library(dplyr)
library(lattice)
library(ggplot2)
library(Rfast)

# Data
dat = readRDS('/hpc/group/jilab/xiaotan/data/filtered_tsr.RData')
#dat = readRDS('/Users/sunxiaotan/Desktop/Tensor/Covid data/filtered_tsr.RData')

# Combined matrix
combmat = do.call(cbind,dat)

# Retain genes that are expressed in at least 1% of cells:subjects
nonzero = rowSums(combmat != 0)
combmat = cbind(combmat, nonzero)
combmat = combmat[combmat[,ncol(combmat)] > 34,]
combmat = combmat[,-ncol(combmat)]

# Retain half of the genes with highest coefficient of variance (CV)
cv = rowcvs(combmat)
summary(cv)
combmat = cbind(combmat, cv)
combmat = combmat[combmat[,ncol(combmat)] > median(cv),]
combmat = combmat[,-ncol(combmat)]

# Standardization (gene-wise)
# Make sure genes are comparable --> control the effect from gene variation so we can isolate the effect from the sample variable.
combmat = t(combmat)
combmat = scale(combmat)
combmat = t(combmat)

# Transform list into array
ary <- array(unlist(combmat), dim = c(7421, 262, 13))
saveRDS(ary, file = "/hpc/group/jilab/xiaotan/data/covid_ary.RData")
#saveRDS(ary, file = "/Users/sunxiaotan/Desktop/scaled_filtered_ary.RData")

# HOOI
hooi_covid = hooi(ary, r = c(742, 26, 5), itermax = 5)
saveRDS(hooi_covid, file = "/hpc/group/jilab/xiaotan/data/hooi_covid.RData")
