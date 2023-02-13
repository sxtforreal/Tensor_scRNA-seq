# Implementation of ISLET algorithm in R

### Functions
get_train_data <- function(filtered_tsr, train_idx) {
  train = list()
  for (i in 1:length(filtered_tsr)){
    train[[i]] = filtered_tsr[[i]][,train_idx]
  }
  return(train)
}

get_label <- function(idx) {
  return(true_label[idx,])
}

get_matrix <- function(subject, tsr, gene_num, cell_num) {
  mat = matrix(nrow = gene_num, ncol = cell_num)
  for (i in 1:length(tsr)){
    mat[,i] = tsr[[i]][,subject]
  }
  return(mat)
}

Matrix_ISLET <- function(train, train_label, r) {
  # Step 1.1
  y = train_label$encoder
  p1 = nrow(train[[1]])
  p2 = length(train)
  n = ncol(train[[1]])
  A = matrix(0, nrow = p1, ncol = p2)
  X = list()
  for (i in 1:n){
    X[[i]] = get_matrix(rownames(train_label[i,]), train, p1, p2)
    A = A + train_label$encoder[i] * X[[i]]
  }
  A = A/n
  
  # Step 1.2
  U1 = sparsesvd(Matrix::Matrix(A, sparse = TRUE), rank = r)$u
  U2 = svd(t(A))$u[,1:r]
  
  # Step 1.3
  XB = matrix(nrow = n, ncol = r^2)
  XD_1 = matrix(nrow = n, ncol = (p1-r)*r)
  XD_2 = matrix(nrow = n, ncol = (p2-r)*r)
  U1_oc = MASS::Null(U1)
  U2_oc = MASS::Null(U2)
  t_U1 = t(U1)
  t_U1_oc = t(U1_oc)
  t_U2_oc = t(U2_oc)
  
  for (j in 1:length(X)){
    matrix_1 = X[[j]] %*% U2
    matrix_2 = t(X[[j]]) %*% U1
    XB[j,] = c(t_U1 %*% matrix_1)
    XD_1[j,] = c(t_U1_oc %*% matrix_1)
    XD_2[j,] = c(t_U2_oc %*% matrix_2)
  }
  
  # Step 1.4
  bh = cv.glmnet(XB, train_label$encoder, alpha = 0, family = "binomial", type.measure = "class", intercept = FALSE)
  B_h = coef(bh, s = "lambda.min")
  B_h = as.matrix(B_h)
  B_h = B_h[2:length(B_h)]
  B_hat = matrix(B_h, nrow = r, ncol = r)
  
  d1 = cv.glmnet(XD_1, train_label$encoder, alpha = 1, family = "binomial", type.measure = "class", intercept = FALSE)
  d1 = coef(d1, s = "lambda.min")
  D1_h = as.matrix(d1)
  D1_h = D1_h[2:length(D1_h)]
  D1_hat = matrix(D1_h, nrow = p1-r, ncol = r)
  
  d2 = cv.glmnet(XD_2, train_label$encoder, alpha = 0, family = "binomial", type.measure = "class", intercept = FALSE)
  d2 = coef(d2, s = "lambda.min")
  D2_h = as.matrix(d2)
  D2_h = D2_h[2:length(D2_h)]
  D2_hat = matrix(D2_h, nrow = p2-r, ncol = r)
  
  # Step 1.5
  L1_hat = (U1 %*% B_hat + U1_oc %*% D1_hat) %*% solve(B_hat)
  L2_hat = (U2 %*% t(B_hat) + U2_oc %*% D2_hat) %*% solve(t(B_hat))
  A_hat = L1_hat %*% B_hat %*% t(L2_hat)
  
  return(A_hat)
}

find_optim <- function(method1_pred, method2_pred, auc1, auc2) {
  auc_opt = vector()
  for(i in 1:10){
    test_idx = folds[i][[1]]
    test_label = get_label(test_idx)
    label = test_label$encoder
    auc_sum = auc1[i] + auc2[i]
    w1 = auc1[i]/auc_sum
    w2 = auc2[i]/auc_sum
    combined_prob = w1*method1_pred[[i]][,1] + w2*method2_pred[[i]]
    auc_opt = append(auc_opt, auc(label, combined_prob))
    }
  return(auc_opt)
}

Matrix_ISLET_Concat <- function(train, train_label, r) {
  # Step 1.1
  y = train_label$encoder
  p1 = nrow(train[[1]]) + 1
  p2 = length(train)
  n = ncol(train[[1]])
  A = matrix(0, nrow = p1, ncol = p2)
  X = list()
  for (i in 1:nrow(train_label)){
    m = get_matrix(rownames(train_label[i,]), train, p1-1, p2)
    X[[i]] = rbind(m, cell_info[rownames(train_label[i,]),])
    A = A + train_label$encoder[i] * X[[i]]
  }
  A = A/n
  
  # Step 1.2
  U1 = sparsesvd(Matrix::Matrix(A, sparse=TRUE), rank = r)$u
  U2 = svd(t(A))$u[,1:r]
  
  # Step 1.3
  XB = matrix(nrow = n, ncol = r^2)
  XD_1 = matrix(nrow = n, ncol = (p1-r)*r)
  XD_2 = matrix(nrow = n, ncol = (p2-r)*r)
  U1_oc = MASS::Null(U1)
  U2_oc = MASS::Null(U2)
  t_U1 = t(U1)
  t_U1_oc = t(U1_oc)
  t_U2_oc = t(U2_oc)
  
  for (j in 1:length(X)){
    matrix_1 = X[[j]] %*% U2
    matrix_2 = t(X[[j]]) %*% U1
    XB[j,] = c(t_U1 %*% matrix_1)
    XD_1[j,] = c(t_U1_oc %*% matrix_1)
    XD_2[j,] = c(t_U2_oc %*% matrix_2)
  }
  
  # Step 1.4
  bh = cv.glmnet(XB, train_label$encoder, alpha = 0, family = "binomial", type.measure = "class", intercept = FALSE)
  B_h = coef(bh, s = "lambda.min")
  B_h = as.matrix(B_h)
  B_h = B_h[2:length(B_h)]
  B_hat = matrix(B_h, nrow = r, ncol = r)
  
  d1 = cv.glmnet(XD_1, train_label$encoder, alpha = 1, family = "binomial", type.measure = "class", intercept = FALSE)
  d1 = coef(d1, s = "lambda.min")
  D1_h = as.matrix(d1)
  D1_h = D1_h[2:length(D1_h)]
  D1_hat = matrix(D1_h, nrow = p1-r, ncol = r)
  
  d2 = cv.glmnet(XD_2, train_label$encoder, alpha = 0, family = "binomial", type.measure = "class", intercept = FALSE)
  d2 = coef(d2, s = "lambda.min")
  D2_h = as.matrix(d2)
  D2_h = D2_h[2:length(D2_h)]
  D2_hat = matrix(D2_h, nrow = p2-r, ncol = r)
  
  # Step 1.5
  L1_hat = (U1 %*% B_hat + U1_oc %*% D1_hat) %*% solve(B_hat)
  L2_hat = (U2 %*% t(B_hat) + U2_oc %*% D2_hat) %*% solve(t(B_hat))
  A_hat = L1_hat %*% B_hat %*% t(L2_hat)
  
  return(A_hat)
}
