# Tensor_scRNA-seq   
Apply tensor method on single cell RNA sequence data   

# Pipeline   
## Tensor:  
Convert the 14881*243*13 data into 14881*3159 matrix   
Retain genes that are expressed in at least 1%(31) of cells —> 14838*3159   
Retain half of the genes with highest CV —> 7419*3159   
Row-wise standardization   
Convert it back to 7419*243*13 array   
HOOI decomposition with r1,r2,r3 = 742,24,5   
Keep 10 columns from the factor matrix   
Train SVM and random forest   

## PCA(gene:cell type vs subject):   
Convert the  14881*243*13 data into 193453*243 matrix   
Retain genes that are expressed in at least 1%(2) of subjects —> 188903*243   
Retain half of the genes with highest CV —> 94451*243   
Transpose into 243*94451 matrix and do PCA with centering and scaling   
The score matrix is of size 243*243   
Keep first 10 PCs   
Train SVM and random forest         
