library(rTensor)
library(tensr)
library(dplyr)
library(lattice)
library(ggplot2)
library(Rfast)
library(ggExtra)
library(findPC)
library(RColorBrewer)
library(scales)
library(cluster)
library(mclust)
dat = readRDS('/Users/sunxiaotan/Desktop/Tensor/Data/data.rds')
ary_scale1 = readRDS('/Users/sunxiaotan/Desktop/Tensor/Data/scaled_ary1.RData')

dat1 <- readRDS('/Users/sunxiaotan/Desktop/Tensor/Data/hooi/Scaling scheme 1/hooi1.RData')

### Functions
# Calculate the Frobenius Norm
Fnorm <- function(ary, G, U){
  approx = atrans(G, U)
  fnorm = fnorm(ary - approx)
  return(fnorm)
}

# Select optimal number of PCs for mode 2, scheme 1
selectPC_mod2 <- function(scaled_ary, tnsr, r2){
  ary = scaled_ary
  dim(ary) = c(7285*15, 209)
  score = ary %*% tnsr[["U"]][[2]]
  sd = apply(score, 2, sd)
  x = as.data.frame(sd)
  x$PC = rownames(x)
  x = x[order(x$sd, decreasing = T), ]
  num_PC = as.numeric(findPC::findPC(x$sd, r2))
  return(x$PC[1:num_PC])
}

# Select optimal number of PCs for mode 3, scheme 1
selectPC_mod3 <- function(scaled_ary, tnsr, r3){
  ary = scaled_ary
  dim(ary) = c(7285*209, 15)
  score = ary %*% tnsr[["U"]][[3]]
  sd = apply(score, 2, sd)
  x = as.data.frame(sd)
  x$PC = rownames(x)
  x = x[order(x$sd, decreasing = T), ]
  num_PC = as.numeric(findPC::findPC(x$sd, r3))
  return(x$PC[1:num_PC])
}

# Calculate entropy of accuracy
# 8 for kclust, 9 for hclust, 10 for mclust
H_acc <- function(x){
  p1x1 = sum(selected_PC$true_label == "Mi&HD" & selected_PC[,x] == 1)/sum(selected_PC[,x] == 1)
  p1x2 = sum(selected_PC$true_label == "Mi&HD" & selected_PC[,x] == 2)/sum(selected_PC[,x] == 2)
  p2x1 = sum(selected_PC$true_label == "Mod&Se" & selected_PC[,x] == 1)/sum(selected_PC[,x] == 1)
  p2x2 = sum(selected_PC$true_label == "Mod&Se" & selected_PC[,x] == 2)/sum(selected_PC[,x] == 2)
  return(-(p1x1*log(p1x1)+p1x2*log(p1x2)+p2x1*log(p2x1)+p2x2*log(p2x2))/2)
}

# Calculate entropy of purity
# 8 for kclust, 9 for hclust, 10 for mclust
H_pur <- function(x){
  p1x1 = sum(selected_PC$true_label == "Mi&HD" & selected_PC[,x] == 1)/sum(selected_PC$true_label == "Mi&HD")
  p1x2 = sum(selected_PC$true_label == "Mi&HD" & selected_PC[,x] == 2)/sum(selected_PC$true_label == "Mi&HD")
  p2x1 = sum(selected_PC$true_label == "Mod&Se" & selected_PC[,x] == 1)/sum(selected_PC$true_label == "Mod&Se")
  p2x2 = sum(selected_PC$true_label == "Mod&Se" & selected_PC[,x] == 2)/sum(selected_PC$true_label == "Mod&Se")
  return(-(p1x1*log(p1x1)+p1x2*log(p1x2)+p2x1*log(p2x1)+p2x2*log(p2x2))/2)
}

# Calculate median silhouette index
med_sil <- function(x) {
  ss = silhouette(selected_PC[,x], dist(selected_PC[,1:6]))
  return(median(ss[, 3]))
}

### Quantify selection of r1, r2, r3
fnorm_dat <- data.frame(matrix(ncol = 5, nrow = 12))
colnames(fnorm_dat) <- c('Scaling_scheme', 'r1', 'r2', 'r3', 'fnorm')
fnorm_dat[1:6,1] = 1
fnorm_dat[7:12,1] = 2
fnorm_dat[c(1:3,7:9),2] = 620
fnorm_dat[c(4:6,10:12),2] = 310
fnorm_dat[,4] = 5
fnorm_dat[c(1,4,7,10),3] = 20
fnorm_dat[c(2,5,8,11),3] = 10
fnorm_dat[c(3,6,9,12),3] = 5
fnorm_dat[1,5] = Fnorm(ary_scale1, dat1$G, dat1$U)
fnorm_dat[2,5] = Fnorm(ary_scale1, dat2$G, dat2$U)
fnorm_dat[3,5] = Fnorm(ary_scale1, dat3$G, dat3$U)
fnorm_dat[4,5] = Fnorm(ary_scale1, dat4$G, dat4$U)
fnorm_dat[5,5] = Fnorm(ary_scale1, dat5$G, dat5$U)
fnorm_dat[6,5] = Fnorm(ary_scale1, dat6$G, dat6$U)
fnorm_dat[7,5] = Fnorm(ary_scale2, dat7$G, dat7$U)
fnorm_dat[8,5] = Fnorm(ary_scale2, dat8$G, dat8$U)
fnorm_dat[9,5] = Fnorm(ary_scale2, dat9$G, dat9$U)
fnorm_dat[10,5] = Fnorm(ary_scale2, dat10$G, dat10$U)
fnorm_dat[11,5] = Fnorm(ary_scale2, dat11$G, dat11$U)
fnorm_dat[12,5] = Fnorm(ary_scale2, dat12$G, dat12$U)

fnorm_dat$kclust = NA
fnorm_dat$hclust = NA
fnorm_dat$mclust = NA

### Clustering
# Select optimal number of PCs
selectPC_mod2(ary_scale1, dat1, 20)

# K-Means (k = 2) -- HD&Mi vs Mod&Se
selected_PC = dat1[["U"]][[2]][,c(10, 8, 2, 3, 9, 13)]
colnames(selected_PC) = c("PC10", "PC8", "PC2", "PC3", "PC9", "PC13")
k <- kmeans(selected_PC, 2, nstart = 30, iter.max = 1000)
palette(alpha(brewer.pal(9, 'Set1'), 0.5))
pairs(selected_PC, col = k$clust, pch = 16)

# Agglomerative Hierarchical Clustering
# Choose ward method because it has the best performance on agglomerative coefficient
hc1 <- agnes(selected_PC, method = "ward")
pltree(hc1, cex = 0.6, hang = -1, main = "Dendrogram of agnes")
two_clusters <- cutree(hc1, k = 2)
table(two_clusters)
palette(alpha(brewer.pal(9, 'Set1'), 0.5))
pairs(selected_PC, col = two_clusters, pch = 16)

# Model-based clustering
#BIC = mclustBIC(selected_PC)
mod1 <- Mclust(selected_PC, G = 2)
summary(mod1, parameters = TRUE)
palette(alpha(brewer.pal(9, 'Set1'), 0.5))
pairs(selected_PC, col = mod1$classification, pch = 16)

### Metrics
true_label = as.data.frame(sub('.*:', '', colnames(dat[["Naive/Central Memory T cell_1:T/NK cell"]])))
colnames(true_label) = "true_label"
true_label$true_label[true_label$true_label == "Mod" | true_label$true_label == "Se"] <- "Mod&Se"
true_label$true_label[true_label$true_label == "Mi" | true_label$true_label == "HD"] <- "Mi&HD"

kclust = as.data.frame(k$clust)
colnames(kclust) = "kclust"

hclust = as.data.frame(two_clusters)
colnames(hclust) = "hclust"

mclust = as.data.frame(mod1$classification)
colnames(mclust) = "mclust"

selected_PC = cbind(selected_PC, true_label, kclust, hclust, mclust)

# Entropy of accuracy -- smaller better
H_acc(8) #k-means
H_acc(9) #hclust
H_acc(10) #mclust

# Entropy of purity -- smaller better
H_pur(8) #k-means
H_pur(9) #hclust
H_pur(10) #mclust

# Adjusted Rand Index
adjustedRandIndex(selected_PC$true_label, selected_PC$kclust)
adjustedRandIndex(selected_PC$true_label, selected_PC$hclust)
adjustedRandIndex(selected_PC$true_label, selected_PC$mclust)

# Median Silhouette index
med_sil(8)
med_sil(9)
med_sil(10)

### Plot results
d1 <- data.frame(PC5 = dat1[["U"]][[2]][,1], PC6 = dat1[["U"]][[2]][,2], severity = sub('.*:', '', colnames(dat[["Naive/Central Memory T cell_1:T/NK cell"]])))
d1$severity[d1$severity == "Mod" | d1$severity == "Se"] <- "Mod&Se"
d1$severity[d1$severity == "Mi" | d1$severity == "HD"] <- "Mi&HD"
p1 <- ggplot(d1, aes(x = PC5, y = PC6, col = severity)) + geom_point()
ggMarginal(p1, type="boxplot", groupColour = TRUE, groupFill = TRUE)
