# Pair-wise prediciton performance comparison between proposed methods and benchmark methods.

library(MASS)
library(dplyr)
library(glmnet)
library(Rfast)
library(sigmoid)
library(caret)
library(pROC)
library(sparsesvd)
library(xgboost)
library(tidyverse)

##### Data preparation
tsr = readRDS('/Users/sunxiaotan/Desktop/Tensor/Covid/tsr_covid.RData')
meta = readRDS('/Users/sunxiaotan/Desktop/Tensor/Covid/meta.rds')
rownames(meta) = meta$`Library Sample Code`
sub_meta = meta[,c('Library Sample Code', 'type')]

### Create true label
true_label = vector()
for (i in 1:ncol(tsr[[1]])) {
  true_label[i] = sub_meta[sub_meta$`Library Sample Code` == colnames(tsr[[1]])[i],]$type
}
true_label = as.data.frame(true_label)
colnames(true_label) = "true_label"
rownames(true_label) = sub(":.*", "", colnames(tsr[[1]]))
true_label = filter(true_label, true_label %in% c("Mod", "Se"))
for (i in 1:nrow(true_label)){
  if (true_label$true_label[i] == 'Mod'){
    true_label$encoder[i] = 1
  }
  else{
    true_label$encoder[i] = 0
  }
}

### Gene selection
# Combined matrix
combmat = do.call(cbind, tsr)

# Retain genes that are expressed in at least 1% of cells:subjects
nonzero = rowSums(combmat != 0)
combmat = cbind(combmat, nonzero)
summary(nonzero)
combmat = combmat[combmat[,ncol(combmat)] > 1814,]
combmat = combmat[,-ncol(combmat)]

# Retain half of the genes with highest coefficient of variance (CV)
cv = rowcvs(combmat)
summary(cv)
combmat = cbind(combmat, cv)
combmat = combmat[combmat[,ncol(combmat)] > 0.55516,]
combmat = combmat[,-ncol(combmat)]

# Standardization (gene-wise)
# Make sure genes are comparable --> control the effect from gene variation so we can isolate the effect from the sample variable.
combmat = t(combmat)
combmat = scale(combmat)
combmat = t(combmat)

# Convert back to tensor
ary = array(unlist(combmat), dim = c(4061, 243, 13))
filtered_tsr = list()
for (i in 1:dim(ary)[3]){
  filtered_tsr[[i]] = ary[,,i]
  colnames(filtered_tsr[[i]]) = colnames(combmat[,1:243])
  filtered_tsr[[i]] = as.data.frame(filtered_tsr[[i]])
  filtered_tsr[[i]] = select(filtered_tsr[[i]], rownames(true_label))
  filtered_tsr[[i]] = as.matrix(filtered_tsr[[i]])
  rownames(filtered_tsr[[i]]) = rownames(combmat)
}

### Make folds
set.seed(1)
ctrl = trainControl(method = 'cv', number = 10, index = createFolds(true_label$encoder, k = 10))
folds = ctrl$index

##### Baseline method -- XGBoost(tree)
matrix_dat = matrix(0, nrow = nrow(true_label), ncol = 4061*13)
for (i in 1:nrow(true_label)){
  subject_matrix = get_matrix(rownames(true_label[i,]), filtered_tsr, 4061, 13)
  subject_vector = c(subject_matrix)
  matrix_dat[i,] = subject_vector
}
auc_xgb = vector()

for (i in 1:10){
  test_dat = matrix_dat[folds[i][[1]],]
  test_label = true_label[folds[i][[1]],]$encoder
  train_dat = matrix_dat[-folds[i][[1]],]
  train_label = true_label[-folds[i][[1]],]$encoder
  dtrain = xgb.DMatrix(data = train_dat, label = train_label)
  dtest = xgb.DMatrix(data = test_dat, label = test_label)
  bst = xgboost(data = dtrain, objective = "binary:logistic", nrounds = 5)
  pred = predict(bst, dtest)
  auc_xgb = append(auc_xgb, auc(test_label, pred))
}
print(mean(auc_xgb))

##### Baseline method -- LASSO
auc_lasso = vector()
for (i in 1:10){
  test_dat = matrix_dat[folds[i][[1]],]
  test_label = true_label[folds[i][[1]],]$encoder
  train_dat = matrix_dat[-folds[i][[1]],]
  train_label = true_label[-folds[i][[1]],]$encoder
  lasso = cv.glmnet(train_dat, train_label, alpha = 1, family = "binomial", type.measure = "class", intercept = FALSE)
  lasso_pred = predict(lasso, s = "lambda.min", newx = test_dat)
  lasso_pred = as.vector(sigmoid(lasso_pred))
  auc_lasso = append(auc_lasso, auc(test_label, lasso_pred))
}
print(mean(auc_lasso))

##### Method 1: ISLET
A = list()
for (i in 1:10){
  seq = c(1,2,3,4,5,6,7,8,9,10)
  seq = seq[-i]
  train_idx = vector()
  for (j in 1:length(seq)){
    train_idx = append(train_idx, folds[seq[j]][[1]]) 
  }
  train = get_train_data(filtered_tsr, train_idx)
  train_label = get_label(train_idx)
  A[[i]] = Matrix_ISLET(train, train_label, 5)
}

### Testing
result_method1 = list()
for (i in 1:length(A)){
  test_idx = folds[i][[1]]
  test_label = get_label(test_idx)
  label = test_label$encoder
  label = as.matrix(label)
  test_dat = list()
  for (j in 1:length(test_idx)){
    test_dat[[j]] = get_matrix(rownames(test_label[j,]), filtered_tsr, 4061, 13)
  }
  probs = list()
  for (k in 1:length(test_dat)){
    probs[[k]] = sum(A[[i]]*test_dat[[k]])
  }
  probs = as.numeric(probs)
  probs = sigmoid(probs)
  result_method1[[i]] = cbind(probs, label)
}

## Compute AUC
auc_method1 = vector()
for (i in 1:length(result_method1)){
  auc_method1 = append(auc_method1, auc(result_method1[[i]][,2], result_method1[[i]][,1]))
}
print(mean(auc_method1))


##### Method 2: Cell proportion model -- XGBoost(tree)
cell_info = readRDS('/Users/sunxiaotan/Desktop/Tensor/Covid/cell_info.RData')

### Cell proportion model
cell_info = t(cell_info)
cp = merge(cell_info, sub_meta, by = 0, all = TRUE)
cp = cp[complete.cases(cp),]
rownames(cp) = cp$Row.names
cp = dplyr::select(cp, -c(1,15,16))
# Match row names
cp = cp[match(rownames(true_label), rownames(cp)),]

pred_method2 = list()
auc_method2 = vector()
for (i in 1:10){
  test_dat = as.matrix(cp[folds[i][[1]],])
  test_label = true_label[folds[i][[1]],]$encoder
  train_dat = as.matrix(cp[-folds[i][[1]],])
  train_label = true_label[-folds[i][[1]],]$encoder
  dtrain = xgb.DMatrix(data = train_dat, label = train_label)
  dtest = xgb.DMatrix(data = test_dat, label = test_label)
  bst = xgboost(data = dtrain, objective = "binary:logistic", nrounds = 5)
  pred = predict(bst, dtest)
  pred_method2[[i]] = pred
  auc_method2 = append(auc_method2, auc(test_label, pred))
}
print(mean(auc_method2))


##### Method 3: Combine gene and cell proportion info
auc_op = find_optim(result_method1, pred_method2, auc_method1, auc_method2)
print(mean(auc_op))

##### Method 4: Concatenate(ISLET)
A_concat = list()
for (i in 1:10){
  seq = c(1,2,3,4,5,6,7,8,9,10)
  seq = seq[-i]
  train_idx = vector()
  for (j in 1:length(seq)){
    train_idx = append(train_idx, folds[seq[j]][[1]]) 
  }
  train = get_train_data(filtered_tsr, train_idx)
  train_label = get_label(train_idx)
  A_concat[[i]] = Matrix_ISLET_Concat(train, train_label, 5)
}


### Testing
result_concat = list()
for (i in 1:length(A_concat)){
  test_idx = folds[i][[1]]
  test_label = get_label(test_idx)
  label = test_label$encoder
  label = as.matrix(label)
  test_dat = list()
  for (j in 1:length(test_idx)){
    m = get_matrix(rownames(test_label[j,]), filtered_tsr, 4061, 13)
    test_dat[[j]] = rbind(m, cell_info[rownames(test_label[j,]),])
  }
  probs = list()
  for (k in 1:length(test_dat)){
    probs[[k]] = sum(A_concat[[i]]*test_dat[[k]])
  }
  probs = as.numeric(probs)
  probs = sigmoid(probs)
  result_concat[[i]] = cbind(probs, label)
}

## Compute AUC
auc_concat = vector()
for (i in 1:length(result_concat)){
  auc_concat = append(auc_concat, auc(result_concat[[i]][,2], result_concat[[i]][,1]))
}
print(mean(auc_concat))
