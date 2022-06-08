# Tensor_scRNA-seq   
Apply tensor method on single cell RNA sequence data   

# Pipeline   
### Tensor:   
Convert the 14881*209*15 original data into 14881*3135 matrix   
Retain genes that are expressed in at least 1%(31) of cells —> 14257*3135   
Retain half of the genes with highest CV —> 7128*3135   
Row-wise standardization   
Convert it back to 7128*209*15 array   
HOOI decomposition with r1,r2,r3 = 712,20,5//1068,20,5//1424,20,5//712,30,5//712,40,5   
Keep 10 columns from the factor matrix   
Train SVM and random forest   

### PCA(gene:cell type vs subject):   
Convert the 14881*209*15 original data into 223215*209 matrix   
Retain genes that are expressed in at least 1%(2) of subjects —> 197317*209   
Retain half of the genes with highest CV —> 98658*209   
Transpose into 209*98658 matrix and do PCA with centering and scaling   
The score matrix is of size 209*209   
Keep first 10 PCs   
Train SVM and random forest   
