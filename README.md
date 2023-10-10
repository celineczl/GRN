# Gene Regulatory Network Inference

## Project Objectives
To develop and validate a computational method for inferring gene regulatory networks from genomic data using a non-local prior Bayesian model as a foundation for high-dimensional variable selection. Other statistical models and variable selection methods (lasso and mombf) are tested and compared.

## Set-up
### Github
`git clone `
Class 1: assume TF expression = TF activity
* So far, only focus on CDKN1A and its TFs
* testing method: 20 TFs for training, 10 TFs for testing
* calculate normalized root mean squared error (NRMSE) for evaluation

model 0 - univariate linear model
* uni_reg_cdkn1a.R

model 1 - multivariate linear model
* multi_reg_cdkn1a.R

model 2 - multivariate lasso model (glmnet)
* lasso_cdkn1a.R

model 3 - multivariate non-local prior model (mombf)
* mombf_cdkn1a.R

‌Class 2: infer TF activity from target genes' expression
* mean expression of target genes (DoRothEA) as activity of TF [no normalization performed]
* Summerized GTEX data [pick out all activities of the TF]
