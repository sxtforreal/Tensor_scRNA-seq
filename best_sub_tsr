# Find the best sub-tensor in cases when raw tensor is not complete.

library(ggplot2)
library(ggExtra)
library(findPC)
library(tidyverse)
library(Matrix)
library(rTensor)
library(matrixcalc)
library(MASS)
library(parallel)
count = readRDS('/Users/sunxiaotan/Desktop/Tensor/Covid data/count.rds')
normcount = readRDS('/Users/sunxiaotan/Desktop/Tensor/Covid data/normcount.rds')
cluster = readRDS('/Users/sunxiaotan/Desktop/Tensor/Covid data/cluster.rds')
celltype = readRDS('/Users/sunxiaotan/Desktop/Tensor/Covid data/celltype.rds')
pb_count = readRDS('/Users/sunxiaotan/Desktop/Tensor/Covid data/pb_count.rds')
pb_norm = readRDS('/Users/sunxiaotan/Desktop/Tensor/Covid data/pb_norm.rds')
meta = readRDS('/Users/sunxiaotan/Desktop/Tensor/Covid data/meta.rds')

# Find unique subjects (429)
subjects <- unique(sub(":.*", "", count@Dimnames[[2]]))
cells <- names(pb_norm)

# Fill the gaps with zeros
fill_zeros <- function(x) {
  ary = array(numeric(),c(14881,429,31)) 
  for (i in 1:length(x)) {
    new_mat <- matrix(0, ncol = 429, nrow = 14881)
    colnames(new_mat) = subjects
    for (j in 1:ncol(new_mat)) {
      if (colnames(new_mat)[j] %in% colnames(x[[i]])) {
        new_mat[,j] = x[[i]][, colnames(new_mat)[j]]
      }
      else {
        next
      }
    }
    ary[,,i] = new_mat
  }
  return(ary)
}

# Use this after fill_zeros, provide the map of subject distribution in each cell type
subjects_distribution <- function(x) {
  df <- data.frame(matrix(ncol = dim(x)[2], nrow = dim(x)[3]))
  colnames(df) <- subjects
  rownames(df) <- cells
  for (i in 1:dim(x)[3]) {
    for (j in 1:dim(x)[2]) {
      if (sum(x[,j,i]) != 0) {
        df[i,j] = 1
      }
      else {
        df[i,j] = 0
      }
    }
  }
  o = df[,order(colSums(-df))]
  p = o[order(rowSums(-o)),]
  return(p)
}

# Find the column index of first zero value and suggest the best subset through computing the biggest rectangular area possible, use this after generating a dataset from subject_distribution
find_best_subset <- function(x) {
  first_idx = data.frame(matrix(ncol = 1, nrow = nrow(x)))
  rownames(first_idx) <- rownames(x)
  area = data.frame(matrix(ncol = 1, nrow = nrow(x)))
  colnames(area) <- c('area')
  rownames(area) <- rownames(x)
  for (i in 1:nrow(x)) {
    for (j in 1:ncol(x)) {
      if (x[i,j] == 1) {
        next
      }
      if (x[i,j] == 0) {
        first_idx[i,1] = j
        break
      }
    }
  }
  for (k in 1:nrow(first_idx)) {
    if (is.na(first_idx[k,1]) == T) {
      area[k,1] = k*ncol(x)
    }
    else {
      area[k,1] = k*(first_idx[k,1]-1)
    }
  }
  selected_cells = rownames(area)[1:which.max(area[,1])]
  selected_subjects = colnames(x)[1:(max(area[,1])/which.max(area[,1]))]
  newList <- list("cells" = selected_cells, "subjects" = selected_subjects)
  return(newList)
}

selected_dat <- function(x, o) {
  lst = list()
  for (i in 1:length(o$cells)) {
    lst[[i]] = x[[o$cells[i]]]
  }
  names(lst) = o$cells
  for (i in 1:length(lst)) {
    lst[[i]] = lst[[i]][,o$subjects]
  }
  return(lst)
}

rename_tsr <- function(x, meta) {
  v = vector()
  colname = colnames(x[[1]])
  for (i in 1:length(colname)) {
    v[[i]] = meta[meta$`Library Sample Code` == colname[i], ]$type
  }
  for (j in 1:length(x)) {
    colnames(x[[j]]) = v
  }
  return(x)
}

# Example
n = fill_zeros(pb_norm)
m = subjects_distribution(n)
o = find_best_subset(m)
p = m[1:14,1:267]
which(p==0,arr.ind = T)
# Remove c6, Subject 173,178,221,246,253
o$cells = o$cells[o$cells != "c6"]
o$subjects = o$subjects[! o$subjects %in% colnames(p[,c(173,178,221,246,253)])]
q = selected_dat(pb_norm, o)
s = rename_tsr(q, meta)
saveRDS(s, file = "/Users/sunxiaotan/Desktop/filtered_tsr.RData")
