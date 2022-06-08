library(dplyr)
library(ggplot2)
library(ggExtra)
library(RColorBrewer)
library(caret)
library(randomForest)
library(findPC)

dat = readRDS('/Users/sunxiaotan/Desktop/Tensor/Data/data.rds')
ary_scale1 = readRDS('/Users/sunxiaotan/Desktop/Tensor/Data/scaled_ary1.RData')
dat1 = readRDS('/Users/sunxiaotan/Desktop/Tensor/Data/hooi1.RData')
dat2 = readRDS('/Users/sunxiaotan/Desktop/Tensor/Data/hooi2.RData')
dat3 = readRDS('/Users/sunxiaotan/Desktop/Tensor/Data/hooi3.RData')
dat4 = readRDS('/Users/sunxiaotan/Desktop/Tensor/Data/hooi4.RData')
dat5 = readRDS('/Users/sunxiaotan/Desktop/Tensor/Data/hooi5.RData')
pca = readRDS('/Users/sunxiaotan/Desktop/Tensor/Data/pca.RData')

selectPC <- function(ary, tnsr, r2){
  ary = ary
  dim(ary) = c(7128*15, 209)
  score = ary %*% tnsr[["U"]][[2]]
  sd = apply(score, 2, sd)
  x = as.data.frame(sd)
  x$PC = rownames(x)
  x = x[order(x$sd, decreasing = T), ]
  num_PC = as.numeric(findPC::findPC(x$sd, r2))
  return(x$PC[1:num_PC])
}

keepPC <- function(ary, tnsr, num){
  ary = ary
  dim(ary) = c(7128*15, 209)
  score = ary %*% tnsr[["U"]][[2]]
  sd = apply(score, 2, sd)
  x = as.data.frame(sd)
  x$PC = rownames(x)
  x = x[order(x$sd, decreasing = T), ]
  return(as.numeric(x$PC[1:num]))
}

### TENSOR
set.seed(100)

# Create true label
true_label = as.data.frame(sub('.*:', '', colnames(dat[["Naive/Central Memory T cell_1:T/NK cell"]])))
colnames(true_label) = "true_label"
true_label$true_label[true_label$true_label == "Mod" | true_label$true_label == "Se"] <- "Mod&Se"
true_label$true_label[true_label$true_label == "Mi" | true_label$true_label == "HD"] <- "Mi&HD"


## SVM
pclist = keepPC(ary_scale1, dat1, 5)
d = cbind(dat1[["U"]][[2]][,pclist], true_label)
colnames(d) = c(paste('PC',1:5, sep = ''), "true_label")
d$true_label = as.factor(d$true_label)

# 10-folds cross validation
trctrl <- trainControl(method = "cv", number = 10, savePredictions = TRUE)
svm_tensor <- train(true_label ~., data = d, method = "svmPoly", trControl = trctrl, preProcess = c("center", "scale"))
print(svm_tensor)

## Random forest
rf_tensor <- train(true_label ~., data = d, method = "rf", trControl = trctrl, preProcess = c("center", "scale"))
print(rf_tensor)
#rf_tensor_pred <- predict(rf_tensor, tensor_testing)

### PCA
score = pca$x[,1:5]
pca_dat = cbind(score, true_label)
pca_dat$true_label = as.factor(pca_dat$true_label)

## SVM
svm_pca <- train(true_label ~., data = pca_dat, method = "svmPoly", trControl = trctrl, preProcess = c("center", "scale"))
print(svm_pca)

## Random forest
rf_pca <- train(true_label ~., data = pca_dat, method = "rf", trControl = trctrl, preProcess = c("center", "scale"))
print(rf_pca)
