# Ji Lab, Duke University, School of Medicine - Tensor analysis project

## Background
Precision medicine is a rapidly growing field that seeks to provide highly personalized healthcare treatments based on a patient's unique health data. The advancements in technology and the growing amount of health-related data sources such as patient demographics, health records, and single-cell RNA sequencing data have provided exciting potentials and challenges. The traditional data analysis methods cannot fully take advantage of the big data, and this is where tensor methods come into play. Tensor methods are designed to effectively analyze and extract information from large datasets and uncover latent structures that traditional data analysis methods might miss, especially the complex inter-relationships between multiple data sources.

More specifically, tensor-based statistical analysis is beneficial in many ways:
1. By integrating multiple data sources into a tensor representation, it captures complex relationships between different sources of data.
2. It allows for the simultaneous analysis of multiple data modalities, which provides a more complete picture of the underlying relationships between the data.
3. Improved interpretability and prediction performance as it considers the inter-relationships between data sources.
4. Better handling of missing data because it can leverage the relationships between different data sources to impute missing values.

In this project, we developed a pipeline to integrate multiple data sources into a high-dimensional tensor representation and perform downstream analyses including tensor decomposition and regression. The goal of this project is to predict the progression stage of Covid-19 from three patient-specific data sources: health records, gene expressions, and cell types. Our research may potentially contribute to the development of precision medicine and optimiized clinical resource allocation.

## Detail
The pipeline includes:
1. Data sources integration (matching, inputation, restriction).
2. Implemented the Importance Sketching Low-rank Estimation for Tensors (ISLET) algorithm to approximate the true low-rank structure of the tensor.
3. Applied Higher-Order Orthogonal Iteration (HOOI) tensor decomposition algorithm to identify major modes of variations in the latent structure across multiple data sources.
4. Trained predictive models using the factor matrices obtained from the tensor decomposition

## Software
R - 4.2.2

## References
1. Zhang, A.∗, Luo, Y.†, Raskutti, G., and Yuan, M. (2020). ISLET: fast and optimal low-rank tensor regression via importance sketchings. SIAM Journal on Mathematics of Data Science.
2. Zhang, A.∗ (2019). Cross: Efficient tensor completion. Annals of Statistics.
3. Zhang, A.∗ and Xia, D. (2018). Tensor SVD: Statistical and computational limits. IEEE Transactions on Information Theory.
4. Zhang, A.∗ and Han, R.† (2019). Optimal denoising and singular value decomposition for sparse high-dimensional high-order data. Journal of the American Statistical Association.
5. Han, R., Willett, R. and Zhang, A.∗ (2020). An optimal statistical and computational framework for generalized tensor estimation, Annals of Statistics.
